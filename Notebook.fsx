#r "nuget: AngouriMath.FSharp"

[<Measure>]
type eV

[<Measure>]
type K

[<Measure>]
type nm

[<Measure>]
type m

[<Measure>]
type J

let kB: float<eV / K> = 8.617333262<eV / K> * (pown 10.0 -5)

let jPerEV = 1.602_176_634e-19<J / eV>
let nmPerM = 1e9<nm / m>

let p1atm: float<eV / nm^3> =
    101325.0<J / m^3> / jPerEV / (nmPerM * nmPerM * nmPerM)

let argonEpsilon = 120

let argonSigma = 0.34<nm>

let count = 150

let v0 = ((float count) * kB * 300.0<K>) / p1atm

let l0 = System.Math.Cbrt(v0 / 1.0<nm^3>) * 1.0<nm>

open AngouriMath.FSharp.Functions

// Normalised maximum volume difference from V0_n after change 
let v'_n: float = pown ((l0/1.0<nm>)*1.01) 3


let V0_n = v0 / 1.0<nm^3>

//printfn "%O" (solutions "a" "2a = 10")    
(solutions "lnv_max" $"{v'_n} = exp(ln({V0_n}) + lnv_max)").Evaled.ToString()
// ~ 0.0298