open System
open System.Threading.Tasks

let kB = 8.617333262 * (10.0 ** -5)

[<Struct>]
type Position = { X: float; Y: float; Z: float }

type Position with
    member this.pbcDisplacement(other: Position, boxLength: float) : float =
        let ΔX = abs (this.X - other.X) % boxLength
        let ΔY = abs (this.Y - other.Y) % boxLength
        let ΔZ = abs (this.Z - other.Z) % boxLength
        ΔX ** 2 + ΔY ** 2 + ΔZ ** 2

type Positions = Position array

type Potentials =
    // length -> energy
    abstract member Pairwise: float -> float
    // number density -> energy
    abstract member Tail: float -> float

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
      InitialBoxLength: float // in nm
      Pressure: float // in eV/nm^3
      MolecularPotential: Potentials
      CutoffRadius: float //in nm
      MultiThreading: MultithreadingConstants option }

// Need to decide if doing npt or otherwise


let energy (boxLength: float) (sc: SimulationConstants) (positions: Positions) =
    let pairwise =
        match sc.MultiThreading with
        | Some(mtc) ->
            let a =
                seq {
                    for bin in mtc.ThreadIndexBins do
                        Task.Run(fun () ->
                            let mutable energy = 0.0

                            for (idxA, idxB) in bin do
                                let displacement = positions[idxA].pbcDisplacement (positions[idxB], boxLength)
                                let pairwisePotential = sc.MolecularPotential.Pairwise displacement
                                energy <- energy + pairwisePotential

                            energy)
                }

            (Task.WhenAll a).GetAwaiter().GetResult() |> Array.sum
        | None ->
            let mutable energy = 0.0

            for idxA in 0 .. positions.Length do
                for idxB in idxA + 1 .. positions.Length do
                    let displacement = positions[idxA].pbcDisplacement (positions[idxB], boxLength)
                    let pairwisePotential = sc.MolecularPotential.Pairwise displacement
                    energy <- energy + pairwisePotential

            energy

    let tailInteractions =
        (float positions.Length) / (boxLength ** 3) |> sc.MolecularPotential.Tail

    pairwise + tailInteractions

let translationalMcMove molIdx (boxLength: float) (sc: SimulationConstants) (positions: Positions) =

    ()

let volumeMcMove (sc: SimulationConstants) (positions: Positions) = ()

let monteCarloLoop (random: Random) sc (positions: Positions) =
    let mutable boxLength = sc.InitialBoxLength

    for step in 1 .. sc.Steps do
        for cycle in 0 .. positions.Length do
            // Need to do volume move
            let rand = random.Next(0, positions.Length + 1)

            if rand < positions.Length then
                translationalMcMove rand boxLength sc positions
            else
                volumeMcMove sc positions

            ()
