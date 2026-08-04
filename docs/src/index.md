```@meta
CurrentModule = AuCO2RR
```

# AuCO2RR.jl

Continuum simulation of CO₂ reduction on gold: a generalised Poisson–Nernst–Planck
(gPNP) electrolyte model coupled to a CatMAP microkinetic surface mechanism.

The package builds on

- [`LiquidElectrolytes.jl`](https://github.com/j-fu/LiquidElectrolytes.jl) — gPNP transport
  with finite ion size and solvation, IR compensation, CV and IV sweep drivers
- [`VoronoiFVM.jl`](https://github.com/j-fu/VoronoiFVM.jl) — finite volume discretisation
- [`CatmapInterface.jl`](https://github.com/j-fu/CatmapInterface.jl) — reads a CatMAP `.mkm`
  input and emits a symbolic microkinetic rate law

## What it computes

A one-dimensional cell from the electrode at $x = 0$ out to a bulk reservoir at $x = L$.
Seven transported species obey drift–diffusion with activity coefficients, coupled by a
five-reaction carbonate buffer network in the volume and by the CO₂-reduction
microkinetics at the electrode boundary.

Outputs are cyclic voltammograms, current–voltage curves, double-layer capacitance,
and concentration/pH profiles in space and time.

## Layout

| Path | Contents |
|:--|:--|
| `src/goldmodel.jl` | `GoldModel` — species layout, rate constants, buffer network, surface reaction, boundary conditions |
| `src/cv.jl` | CV sweep drivers (over $L$, over scan rate, over pH) |
| `src/iv.jl` | IV sweep drivers (over $L$, over CO₂ pressure, over pH) |
| `src/dlcap.jl` | double-layer capacitance |
| `src/sweeps_csv.jl` | CSV export |
| `src/AuCO2RR_plots/` | plotting layer (`cvplot.jl`, `ivplot.jl`, `capsplot.jl`) |
| `notebooks/` | Pluto notebooks driving the above |
| `data/catmap_CO2R_data/` | CatMAP `.mkm` input and formation energies |

## Building these docs

```julia
using Pkg
Pkg.activate("docs")
Pkg.add(["Documenter", "DocumenterMermaid"])
Pkg.develop(PackageSpec(path = "."))     # dev the package itself into the docs env
include("docs/make.jl")
```

The generated site lands in `docs/build/`. Open `docs/build/index.html`.

!!! note "No remote configured"
    `make.jl` sets `remotes = nothing`, so "edit on GitHub" and source links are
    disabled. Once the repository has a GitHub remote, drop that line, set
    `repo = Remotes.GitHub("USER", "Capacitance_Code")`, and re-enable the
    `deploydocs` block at the bottom of `make.jl`.

## Where to go next

The [Guide](@ref) documents the model itself — species indexing, the buffer network, how
the microkinetic turnover frequency becomes a boundary flux, and how "the current" is
defined. Read it before adding a plot or interpreting an absolute current.
