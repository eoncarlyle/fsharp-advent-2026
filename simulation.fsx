open System

[<Struct>]
type Position = { X: float; Y: float; Z: float }


type MolecularPotential =
    // length -> energy
    abstract member Potential: float -> float
    // number density -> energy
    abstract member TailPotential: float -> float

[<Struct>]
type SimulationParameters = {
    Cycles: int;
    InitialLength: float; // in nm
    Pressure: float; // in eV/nm^3
}


let rand = Random()

// Need to decide if doing npt or otherwise

let monteCarloMove molIdex (positions: Position array) =
    ()

let monteCarloLoop cycles (positions: Position array) =
    for _ in 1 .. cycles do
        // Need to do volume move
        let selectedMolIdx = rand.Next(0, positions.Length)
        monteCarloMove selectedMolIdx positions
