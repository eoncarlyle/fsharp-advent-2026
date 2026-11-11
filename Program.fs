module BasicNpt

open System
open System.Collections
open System.Threading.Tasks

[<Measure>]
type eV

[<Measure>]
type K


[<Measure>]
type nm

let kB: float<eV / K> = 8.617333262<eV / K> * (10.0 ** -5)

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
    { Positions: Positions
      mutable Energy: float<eV> }

type Potentials =
    // length -> energy
    abstract member Pairwise: float<nm> -> float<eV>
    // number density -> energy
    abstract member Tail: float<nm^-3> -> float<eV>

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
    { Steps: int
      InitialBoxLength: float<nm>
      Pressure: float<eV / nm^3>
      Temperature: float<K>
      MolecularPotential: Potentials
      CutoffRadius: float<nm>
      MaximumVolumeLogChange: float
      MultiThreading: MultithreadingConstants option }

type SimulationStatistics() =
    let mutable volumeMoveAttempts = 0
    let mutable volumeMoveSuccesses = 0
    let mutable translationalMoveAttempts = 0
    let mutable translationalMoveSuccesses = 0

    member this.updateTranslationalMoves(successful: bool) =
        volumeMoveAttempts <- volumeMoveAttempts + 1

        if successful then
            (volumeMoveSuccesses <- volumeMoveAttempts + 1)

    // Not thread safe
    member this.updateVolumeMoves(successful: bool) =
        translationalMoveAttempts <- translationalMoveAttempts + 1

        if successful then
            (translationalMoveSuccesses <- translationalMoveSuccesses + 1)

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
    let pairwise =
        match sc.MultiThreading with
        | Some(mtc) ->
            let tasks =
                seq {
                    for bin in mtc.ThreadIndexBins do
                        Task.Run(fun () ->
                            let mutable energy = 0.0<eV>

                            for (idxA, idxB) in bin do
                                let displacement = positions[idxA].pbcDisplacement (positions[idxB], boxLength)
                                let pairwisePotential = sc.MolecularPotential.Pairwise displacement
                                energy <- energy + pairwisePotential

                            energy)
                }

            (Task.WhenAll tasks).GetAwaiter().GetResult() |> Array.sum
        | None ->
            let mutable energy = 0.0<eV>

            for idxA in 0 .. positions.Length do
                for idxB in idxA + 1 .. positions.Length do
                    let displacement = positions[idxA].pbcDisplacement (positions[idxB], boxLength)
                    let pairwisePotential = sc.MolecularPotential.Pairwise displacement
                    energy <- energy + pairwisePotential

            energy

    let tailInteractions =
        (float positions.Length) / (boxLength * boxLength * boxLength)
        |> sc.MolecularPotential.Tail

    pairwise + tailInteractions

let translationalMcMove
    selectedMolIdx
    (boxLength: float<nm>)
    (sc: SimulationConstants)
    (random: Random)
    (energyPositions: EnergyPositions)
    : bool =

    let randX = random.NextDouble()
    let randY = random.NextDouble()
    let randZ = random.NextDouble()

    let translate (existingDimension: float<nm>) (rand: float) =
        let naiveTranslation = existingDimension + (rand - 0.5) * sc.CutoffRadius

        if (abs naiveTranslation) >= 2.0 * boxLength then
            let message =
                $"Illegal translation; boxLength:{boxLength}, randX: {randX}, randY: {randY}, randZ: {randZ}"

            raise (Exception(message))
        else if naiveTranslation > 0.0<nm> then
            (naiveTranslation % boxLength)
        else
            boxLength + naiveTranslation % boxLength

    let originalMolecule = energyPositions.Positions[selectedMolIdx]
    let e0 = energyPositions.Energy

    let translatedMolecule =
        { X = translate originalMolecule.X randX
          Y = translate originalMolecule.Y randY
          Z = translate originalMolecule.Z randZ }

    energyPositions.Positions[selectedMolIdx] <- translatedMolecule
    let e' = energy boxLength sc energyPositions.Positions

    let randE = random.NextDouble()

    if randE < exp(-1.0 * (e' - e0) / (kB * sc.Temperature)) then
        energyPositions.Energy <- e'
        true
    else
        energyPositions.Positions[selectedMolIdx] <- originalMolecule
        false


let volumeMcMove (sc: SimulationConstants) (random: Random) (boxLength: float<nm>) (energyPositions: EnergyPositions) : float<nm> option =
    
    let v0 = boxLength * boxLength * boxLength
    let randlnV' = (random.NextDouble() - 0.5) * sc.MaximumVolumeLogChange * 2.0
    let deltaV = 1.0<nm^3> * exp randlnV'
    let v' = v0 + deltaV
    
    
    None

let monteCarloLoop (random: Random) sc (stats: SimulationStatistics) (positions: Positions) =
    let mutable boxLength = sc.InitialBoxLength

    let energyPositions =
        { Positions = positions
          Energy = energy boxLength sc positions }

    for step in 1 .. sc.Steps do
        for cycle in 0 .. positions.Length do
            // Need to do volume move
            let rand = random.Next(0, positions.Length + 1)

            if rand < positions.Length then
                translationalMcMove rand boxLength sc random energyPositions
                |> stats.updateTranslationalMoves
            else
                let maybeNewBoxLength = volumeMcMove sc random boxLength energyPositions

                match maybeNewBoxLength with
                | Some(newBoxLength) ->
                    stats.updateTranslationalMoves true
                    boxLength <- newBoxLength
                | None -> stats.updateTranslationalMoves false

            ()

[<EntryPoint>]
let main _ = 0
