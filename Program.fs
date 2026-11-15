module BasicNpt

open System
open System.Collections
open System.Threading.Tasks

[<Measure>]
type eV

[<Measure>]
type K

[<Measure>]
type J

[<Measure>]
type nm

[<Measure>]
type m

[<Measure>]
type Pa = J / m^3

module conversions =

    let jPerEV = 1.602_176_634e-19<J / eV>
    let nmPerM = 1e9<nm / m>

    let p1atm (pressure: float<Pa>) : float<eV / nm^3> =
        pressure / jPerEV / (nmPerM * nmPerM * nmPerM)

let kB: float<eV / K> = 8.617333262<eV / K> * (pown 10.0 -5)


[<Struct>]
type Position =
    { X: float<nm>
      Y: float<nm>
      Z: float<nm> }

type Position with
    member this.pbcDisplacement(other: Position, boxLength: float<nm>) : float<nm> =
        let deltaX = abs (this.X - other.X) % boxLength
        let deltaY = abs (this.Y - other.Y) % boxLength
        let deltaZ = abs (this.Z - other.Z) % boxLength

        (deltaX * deltaX) + (deltaY * deltaY) + (deltaZ * deltaZ) |> sqrt

type Positions = Position array

type EnergyPositions =
    { mutable Positions: Positions // needed for volume code
      mutable Energy: float<eV>
      VolumeMovePositionBuffer: Positions }

type Potentials =
    // length -> energy
    abstract member Pairwise: float<nm> -> float<eV>
    // number density -> energy
    abstract member Tail: float<nm^-3> -> float<eV>

type LJPotentials(ε: float<eV>, σ: float<nm>, rCutoff: float<nm>) =
    interface Potentials with
        member this.Pairwise(r: float<nm>) : float<eV> =
            (σ / r) |> fun m -> 4.0 * ε * ((pown m 12) - (pown m 6))

        member this.Tail(numberDensity: float<nm^-3>) : float<eV> =
            let a = (8.0 / 3.0) * Math.PI * numberDensity * ε * (σ * σ * σ)
            // I'd be curious as to if this would make any difference:
            //let b = (σ/rCutoff) |> fun ratio -> (1.0/3.0)*(pown ratio 9) - (pown ratio 3)
            let b = (1.0 / 3.0) * (pown (σ / rCutoff) 9) - (pown (σ / rCutoff) 3)
            a * b

[<Struct>]
type MultithreadingConstants =
    {
      // Couldn't figure out in closed form how to express '1/Nth of the unique i,j pairs, so
      // Throwing up my hands and precompputing this
      // This is a little inelegant, but shrug emoji
      ThreadIndexBins: (int * int) array array }

// Need to be careful how multithreaded options are created
type MultithreadingConstants with
    member this.getThreads() = this.ThreadIndexBins.Length


[<Struct>]
type SimulationConstants =
    { Molecules: int
      Cycle: int
      InitialBoxLength: float<nm>
      Pressure: float<eV / nm^3>
      Temperature: float<K>
      MolecularPotential: Potentials
      CutoffRadius: float<nm>
      MaximumTranslationMovement: float<nm>
      MaximumVolumeLogChange: float
      MultiThreading: MultithreadingConstants option }

type VolumeMoveOutcome =
    | Accepted of newBoxLength: float<nm>
    | Rejected

type TranslationMoveOutcome =
    | Accepted
    | Rejected

type SimulationStatistics() =
    let mutable volumeMoveAttempts = 0
    let mutable volumeMoveSuccesses = 0
    let mutable translationalMoveAttempts = 0
    let mutable translationalMoveSuccesses = 0
    
    member this.updateVolumeMoves(outcome: VolumeMoveOutcome) =
        volumeMoveAttempts <- volumeMoveAttempts + 1

        match outcome with
        | VolumeMoveOutcome.Accepted(_) -> volumeMoveSuccesses <- volumeMoveSuccesses + 1
        | VolumeMoveOutcome.Rejected -> ()

    // Not thread safe
    member this.updateTranslationalMoves(outcome: TranslationMoveOutcome) =
        translationalMoveAttempts <- translationalMoveAttempts + 1

        match outcome with
        | TranslationMoveOutcome.Accepted -> translationalMoveSuccesses <- translationalMoveSuccesses + 1
        | TranslationMoveOutcome.Rejected -> ()

    // Not thread safe
    member this.getMaybePVolumeSuccess() =
        if volumeMoveAttempts > 0 then
            (float32 volumeMoveSuccesses) / (float32 volumeMoveAttempts) |> Some
        else
            None

    member this.getMaybePTranslationalSuccess() =
        if translationalMoveAttempts > 0 then
            (float32 translationalMoveSuccesses) / (float32 translationalMoveAttempts)
            |> Some
        else
            None

let energy (boxLength: float<nm>) (sc: SimulationConstants) (positions: Positions) =
    
    let diagnosticList = System.Collections.Generic.List<float<eV>>()
    
    let pairwise =
        match sc.MultiThreading with
        | Some(mtc) ->
            let tasks =
                seq {
                    for bin in mtc.ThreadIndexBins do
                        Task.Run(fun () ->
                            let mutable energy = 0.0<eV>

                            for idxA, idxB in bin do
                                let displacement = positions[idxA].pbcDisplacement (positions[idxB], boxLength)
                                if (displacement <= sc.CutoffRadius) then
                                    let pairwisePotential = sc.MolecularPotential.Pairwise displacement
                                    // Here
                                    diagnosticList.Add(pairwisePotential)
                                    energy <- energy + pairwisePotential 
                            energy)
                }

            (Task.WhenAll tasks).GetAwaiter().GetResult() |> Array.sum
        | None ->
            let mutable energy = 0.0<eV>

            for idxA in 0 .. positions.Length - 1 do
                for idxB in idxA + 1 .. positions.Length - 1 do
                    let displacement = positions[idxA].pbcDisplacement (positions[idxB], boxLength)
                    let pairwisePotential = sc.MolecularPotential.Pairwise displacement
                    // Here
                    diagnosticList.Add(pairwisePotential)
                    energy <- energy + pairwisePotential

            energy

    let tailInteractions =
        (float positions.Length) / (boxLength * boxLength * boxLength)
        |> sc.MolecularPotential.Tail

    if (pairwise + tailInteractions).ToString().Equals "NaN" then
        failwith "checkpoint"
        
    pairwise + tailInteractions

let translationalMcMove
    selectedMolIdx
    (boxLength: float<nm>)
    (sc: SimulationConstants)
    (random: Random)
    (energyPositions: EnergyPositions)
    : TranslationMoveOutcome =

    let randX = random.NextDouble()
    let randY = random.NextDouble()
    let randZ = random.NextDouble()

    let translate (existingDimension: float<nm>) (rand: float) =
        let naiveTranslation =
            existingDimension + 2.0 * (rand - 0.5) * sc.MaximumTranslationMovement

        if (abs naiveTranslation) >= 2.0 * boxLength then
            let message =
                $"Illegal translation; boxLength:{boxLength}, randX: {randX}, randY: {randY}, randZ: {randZ}"

            failwith message
        else if naiveTranslation > 0.0<nm> then
            (naiveTranslation % boxLength)
        else
            boxLength + naiveTranslation % boxLength

    let originalMolecule = energyPositions.Positions[selectedMolIdx]
    let u0 = energyPositions.Energy

    let translatedMolecule =
        { X = translate originalMolecule.X randX
          Y = translate originalMolecule.Y randY
          Z = translate originalMolecule.Z randZ }

    energyPositions.Positions[selectedMolIdx] <- translatedMolecule
    let u' = energy boxLength sc energyPositions.Positions

    let randE = random.NextDouble()
    let power = (-1.0 * (u' - u0) / (kB * sc.Temperature))
    
    
    if randE < exp power then
        energyPositions.Energy <- u'
        TranslationMoveOutcome.Accepted
    else
        energyPositions.Positions[selectedMolIdx] <- originalMolecule
        TranslationMoveOutcome.Rejected


let volumeMcMove
    (sc: SimulationConstants)
    (random: Random)
    (boxLength: float<nm>)
    (energyPositions: EnergyPositions)
    : VolumeMoveOutcome =

    let v0 = boxLength * boxLength * boxLength
    let lnV' = log(v0/1.0<nm^3>) + (random.NextDouble() - 0.5) * sc.MaximumVolumeLogChange * 2.0
    let v' = 1.0<nm^3> * exp lnV'
    let deltaV = v' - v0
    let boxLength' = Math.Pow(v' / 1.0<nm^3>, 1.0 / 3.0) * 1.0<nm>
    // Would be interesting to see how expensive it is to reallocate entire new array each time
    let lengthRatio = boxLength' / boxLength
    // I suspect that multithreading is not actually performant but the binning should be handled in a pre-calculated
    // field on `sc`

    if boxLength' <= 2.0 * sc.CutoffRadius then
        failwith "Illegal state: box too small for sound periodic boundary condition"

    energyPositions.Positions
    |> Array.iteri (fun idx p ->
        let p' =
            { X = p.X * lengthRatio
              Y = p.Y * lengthRatio
              Z = p.Z * lengthRatio }

        energyPositions.VolumeMovePositionBuffer[idx] <- p')

    let beta = 1.0 / (kB * sc.Temperature)
    let u' = energy boxLength' sc energyPositions.VolumeMovePositionBuffer
    let deltaU = energyPositions.Energy - u'

    let power =
        (float (sc.Molecules + 1)) * log (v' / v0)
        - beta * deltaU
        - beta * sc.Pressure * deltaV

    let rand = random.NextDouble()

    if rand > exp power then
        VolumeMoveOutcome.Rejected
    else
        energyPositions.Positions <- energyPositions.VolumeMovePositionBuffer
        energyPositions.Energy <- u'
        VolumeMoveOutcome.Accepted boxLength'

let displayProgress cycle boxLength (stats: SimulationStatistics) =
    Console.WriteLine($"MC cycle: {cycle}")
    Console.WriteLine($"Box Length: {boxLength}")
    Console.WriteLine($"Translation Move success: {stats.getMaybePTranslationalSuccess ()}")
    Console.WriteLine($"Volume Move success: {stats.getMaybePVolumeSuccess ()}")


let monteCarloLoop (random: Random) sc (positions: Positions) =
    let mutable boxLength = sc.InitialBoxLength
    let stats = SimulationStatistics()

    let energyPositions =
        { Positions = positions
          Energy = energy boxLength sc positions
          VolumeMovePositionBuffer =
            Array.map
                (fun _ ->
                    { X = 0.0<nm>
                      Y = 0.0<nm>
                      Z = 0.0<nm> })
                positions }

    for cycle in 1 .. sc.Cycle do
        for step in 0 .. positions.Length do
            // Double check this
            let rand = random.Next(0, positions.Length + 1)

            if rand < positions.Length then
                translationalMcMove rand boxLength sc random energyPositions
                |> stats.updateTranslationalMoves
            else
                let outcome = volumeMcMove sc random boxLength energyPositions
                stats.updateVolumeMoves outcome

                match outcome with
                | VolumeMoveOutcome.Accepted newBoxLength -> boxLength <- newBoxLength
                | _ -> ()

        if cycle % 100 = 0 then
            displayProgress cycle boxLength stats

    (boxLength, stats)

let spacePositions (molecules: int) (boxLength: float<nm>) =
    let count = Math.Cbrt molecules |> ceil |> int
    let spacing = boxLength / (float count)

    seq {
        for xi in 0 .. count - 1 do
            for yi in 0 .. count - 1 do
                for zi in 0 .. count - 1 do
                    { X = (float xi) * spacing
                      Y = (float yi) * spacing
                      Z = (float zi) * spacing }
    }
    |> Seq.truncate molecules
    |> Seq.toArray

let getThreadIndexBins (positions: Positions) (cores: int) =
    seq {
        for i in 0 .. positions.Length - 1 do
            for j in i + 1 .. positions.Length - 1 do
                yield (i, j)
    }
    |> Seq.splitInto cores
    |> Seq.toArray



[<EntryPoint>]
let main _ =
    let cores = 12

    let temperature = 300.0<K>
    let cutoffRadius = 4.0<nm>
    let molecules = 125
    let initialBoxLength = 13.0<nm> //17.22<nm>

    let argonLjPotential =
        LJPotentials(120.0 * temperature * kB, 0.34<nm>, cutoffRadius)

    let sc: SimulationConstants =
        { Molecules = molecules
          Cycle = 1000
          InitialBoxLength = initialBoxLength
          Pressure = 0.000632420907<eV / nm^3>
          Temperature = 300.0<K>
          MolecularPotential = argonLjPotential
          CutoffRadius = cutoffRadius
          MaximumTranslationMovement = 2.0<nm>
          MaximumVolumeLogChange = 0.0298
          MultiThreading = None }

    let seed = 1001
    let random = Random(seed)


    let (boxLength, stats) =
        spacePositions molecules initialBoxLength |> monteCarloLoop random sc

    0
