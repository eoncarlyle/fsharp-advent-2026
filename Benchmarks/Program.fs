module Benchmarks

open System
open BenchmarkDotNet.Attributes
open BenchmarkDotNet.Running
open BasicNpt

[<MemoryDiagnoser>]
type NptBenchmarks() =
    let temperature = 300.0<K>
    let cutoffRadius = 3.0<nm>
    let initialBoxLength = 13.0<nm>
    let seed = 1002

    let mutable simulation = Unchecked.defaultof<NpTSimulation>

    [<Params(125)>]
    member val Molecules = 125 with get, set

    [<Params(1000)>]
    member val Cycles = 1000 with get, set

    [<Params(1)>]
    member val Threads = 1 with get, set

    [<GlobalSetup>]
    member this.Setup() =
        let template = spacePositions this.Molecules initialBoxLength

        let multiThreading =
            if this.Threads = 0 then
                None
            else
                Some { ThreadIndexBins = getThreadIndexBins template this.Threads }

        let potential = LJPotentials(120.0 * temperature * kB, 0.34<nm>, cutoffRadius)

        simulation <-
            NpTSimulation
                { Molecules = this.Molecules
                  Cycle = this.Cycles
                  InitialBoxLength = initialBoxLength
                  Pressure = 0.000632420907<eV / nm^3>
                  Temperature = temperature
                  MolecularPotential = potential
                  CutoffRadius = cutoffRadius
                  MaximumTranslationMovement = 2.0<nm>
                  MaximumVolumeLogChange = 0.0298
                  MultiThreading = multiThreading }

    [<Benchmark>]
    member this.MonteCarloLoop() =
        let positions = spacePositions this.Molecules initialBoxLength
        simulation.monteCarloLoop (Random 1002) positions

[<EntryPoint>]
let main argv =
    BenchmarkSwitcher.FromAssembly(typeof<NptBenchmarks>.Assembly).Run(argv) |> ignore
    0
