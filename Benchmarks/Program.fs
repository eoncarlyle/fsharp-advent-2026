module Benchmarks

open System
open BenchmarkDotNet.Attributes
open BenchmarkDotNet.Running
open BasicNpt
open System.IO

[<MemoryDiagnoser>]
type NptBenchmarks() =
    let temperature = 300.0<K>
    let cutoffRadius =0.34<nm> * 2.5
    let initialBoxLength = 13.0<nm>
    let seed = 1003

    let mutable simulation = Unchecked.defaultof<NpTSimulation>

    [<Params(125)>]
    member val molecules = 125 with get, set

    [<Params(1000)>]
    member val cycles = 1000 with get, set

    [<GlobalSetup>]
    member this.Setup() =
        let template = spacePositions this.molecules initialBoxLength

        let multiThreading = None

        let potential = LJPotentials(120.0 * temperature * kB, 0.34<nm>, cutoffRadius)

        let uncovertedPressure = 101325.0<J / m^3>
        let pressure = conversions.pressureConversion uncovertedPressure

        let argonLjPotential =
            LJPotentials(120.0 * temperature * kB, 0.34<nm>, cutoffRadius)

        simulation <-
            NpTSimulation (
              { Molecules = this.molecules
                Cycles = this.cycles
                InitialBoxLength = initialBoxLength
                Pressure = pressure
                Temperature = temperature
                MolecularPotential = argonLjPotential
                CutoffRadius = cutoffRadius
                MaximumTranslationMovement = 0.7<nm>
                MaximumVolumeLogChange = 0.0298 * 6.0
                MultiThreading = None
                StatisticsEvaluationInterval = 100
                FullStatisticsEvaluationCycles = this.cycles / 2 },
                None,
                seed,
                DirectoryInfo("/tmp"))


    [<Benchmark>]
    member this.MonteCarloLoop() =
        let positions = spacePositions this.molecules initialBoxLength
        simulation.monteCarloLoop positions

[<EntryPoint>]
let main argv =
    BenchmarkSwitcher.FromAssembly(typeof<NptBenchmarks>.Assembly).Run(argv) |> ignore
    0
