# 2026.11.01

Page 35, general `U_tot` for
$$
\mathcal{U}^{\mathrm{tot}} = \sum_{i<j} u_c(r_{ij}) + \frac{N\rho}{2} \int_{r_c}^{\infty} \mathrm{d}r \, u(r) \, 4\pi r^2
$$

where
$$
\frac{8}{3}\pi\rho\epsilon\sigma^3 \left[ \frac{1}{3}\left(\frac{\sigma}{r_c}\right)^9 - \left(\frac{\sigma}{r_c}\right)^3 \right]
$$

for number density `\rho`.

Other math that matters is that 10^6 Pa is 0.00624 eV nm^-3, which has much more convienent units.

Using the LJ paramaters of `\frac{\epsilon}/{k_B} = 120 K` and `\sigma = 0.34 \mathrm{nm}` for
Argon seems appropriate.

NPT algorithim shown on page 121

# 2026.11.03

Basic NPT-Ensemble Simulation

```
PROGRAM mc_npt                    basic NPT ensemble simulation

do icycl=1,ncycl                  perform ncycl MC cycles
  ran=ranf()*(npart+1)+1
  if (ran.le.npart) then
    call mcmove                   perform particle displacement
  else
    call mcvol                    perform volume change
  endif
  if (mod(icycl,nsamp).eq.0)
+   call sample                   sample averages
enddo
end
```

```
SUBROUTINE mcvol                          attempts to change
                                          the volume
call toterg(box,eno)                      total energy old conf.
vo=box**3                                 determine old volume
lnvn=log(vo)+(ranf()-0.5)*vmax            perform random walk in ln V
vn=exp(lnvn)
boxn=vn**(1/3)                            new box length
do i=1,npart
  x(i)=x(i)*boxn/box                      rescale center of mass
enddo
call toterg(boxn,enn)                     total energy new conf.
arg=-beta*((enn-eno)+p*(vn-vo)
+ -(npart+1)*log(vn/vo)/beta)             appropriate weight function!
if (ranf().gt.exp(arg)) then              acceptance rule (5.2.3)
  do i=1,npart                            REJECTED
    x(i)=x(i)*box/boxn                    restore the old positions
  enddo
endif
return
end
```

Basic metroplis displacement

```
SUBROUTINE mcmove                         attempts to displace a particle

o=int(ranf()*npart)+1                     select a particle at random
call ener(x(o),eno)                       energy old configuration
xn=x(o)+(ranf()-0.5)*delx                 give particle random displacement
call ener(xn,enn)                         energy new configuration
if (ranf().lt.exp(-beta                   acceptance rule (3.2.1)
+   *(enn-eno)) x(o)=xn                   accepted: replace x(o) by xn
return
end
```

Two places where I'm thinking 'I know there must be a better way to math this out' is a) equal distribution of unique
pairs to sum over and b) periodic boundary condition distances

I think there is probably a way to avoid the modulos, but I'm not positive.

Also, very annoying that the runtime does not provide a threadpool and I just have to trust it; I'd like to
deterministically prove that only a CPU number of threads are run.

┌───────────────┬──────────────┐
│ a │ a │
│ │ │
│ │ │
│ │ │
│ y │ y │
│ x │ x │
└───────────────└──────────────┘

Okay another note - it is frustrating that `**` is not unit aware. I also think that a written explanation of
`ThreadIndexBins` will make sense. If there were three molecules and two bins, the pairs would be
`[[(1, 2); (1, 3)]; [(2; 3)]`. For 4, `[[(1, 2); (1, 3); (1, 4)]; [(2, 3); (2, 4); (3, 4)]]]`.

# 2026.11.10

Measurements saved me:

```fsharp
let mistake = existingDimension * (rand - 0.5) * sc.CutoffRadius
let correct = existingDimension + (rand - 0.5) * sc.CutoffRadius 
```

I had an epiphany while calculating the `translationalMcMove` that there is no reason to recalculate the energy:
when the positions are first created they should be created with the starting energy.

# 2026.11.11

For volume moves it makes sense to preallocate the buffer where the newly translated particles will be placed into.
Otherwise lots of needless allocation

# 2026.11.12

I can only imagine that doing the volume normalization multithreaded isn't any more efficient, but that will be
something to actually test

# 2026.11.14

Did not konw that scientific notation was supported by F#!
Remark on `e037264`

That is a lot of `translate` allocations that can probably be avoided?

![img.png](img.png)

Attempting to optimise `this.Pairwise` because this is more interesting. Each run taken under profiler with sampling.

I saw two places in the CPU flamegraph that intrinsic `pown` was called. When I run these cold, off the bat they took
like

`d47e311`
25, 26, 25

`82c048d`
34, 34, 34

```
let a = pown (σ / r) 6
4.0 * ε * (a * a - a)
```

Mean: 10.024s, standard deviation: 0.081s

```
4.0 * ε * ((pown (σ / r) 12) - (pown (σ / r) 6))
```

At seed of 1003, these are the properly benchmarked results
New/fast: `Mean = 8.937 s, StdErr = 0.034 s `
Old/slow: `Mean = 9.711 s, StdErr = 0.016 s

The next thing to do is to a) more conviently launch multiple `simulation.monteCarloLoop` calls across different
threads? I'd think that the place to put this is `SimulationStatistics`, but I don't know what form this would take

# 2026.11.22

`c7855d6` and `f48cdcd` tested the simulation with and without tail correction respectively. There is no difference in
the (way too high) final box length.

So two problems - the existing box length issue and the tail correction is never called. But if it literally does not
impact the box length if the tail correction isn't called, then clearly there is some issue with the pairwise potential.

Am I not accounting correctly for something really basic like forgetting to count all pairs? It would make sense to do
pairwise energy the naive/unoptimised/obviously correct way first.

```
➜   grep -hoE 'Mean=[^,]+' * | sort -u
Mean=26.17925991752719
Mean=34.967065375678374
Mean=42.492491458571465
Mean=42.77739929873855
Mean=43.66653325913573
Mean=44.00921010552093
Mean=47.58311068642399
Mean=49.81561916482975
Mean=51.87762131231645
Mean=57.84137009566896
```

# 2026.11.27

Flipped sign + missing multilplication by _N_ for tail correction ended up
being my issue. Small and annoying issues and _hopefully_ the last problem
that we will have on the actual equation of state. Also before I forget,
benchmarking is done with

```
dotnet run -c Release --project Benchmarks
```

`c0da5a7` benchmarks
```
| Method         | molecules | cycles | Mean    | Error    | StdDev   | Median  | Allocated |
|--------------- |---------- |------- |--------:|---------:|---------:|--------:|----------:|
| MonteCarloLoop | 125       | 1000   | 6.279 s | 0.2968 s | 0.8611 s | 5.868 s | 149.73 KB |

```

Not any faster with changes - not likely to be literally slower though
`841468f`
```
| Method         | molecules | cycles | Mean    | Error    | StdDev   | Median  | Allocated |
|--------------- |---------- |------- |--------:|---------:|---------:|--------:|----------:|
| MonteCarloLoop | 125       | 1000   | 6.406 s | 0.3082 s | 0.8942 s | 5.977 s | 149.73 KB |
```