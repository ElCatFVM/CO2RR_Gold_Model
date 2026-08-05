```@meta
CurrentModule = AuCO2RR
```

# AuCO2RR.jl

Continuum simulation of CO₂ reduction on gold: a generalised Poisson–Nernst–Planck (gPNP)
electrolyte coupled to a CatMAP microkinetic surface mechanism.

The package builds on

- [`LiquidElectrolytes.jl`](https://github.com/j-fu/LiquidElectrolytes.jl) — gPNP transport
  with finite ion size and solvation, IR compensation, CV and IV sweep drivers
- [`VoronoiFVM.jl`](https://github.com/j-fu/VoronoiFVM.jl) — finite volume discretisation
- [`CatmapInterface.jl`](https://github.com/j-fu/CatmapInterface.jl) — reads a CatMAP `.mkm`
  input and emits a symbolic microkinetic rate law

The domain is one-dimensional: an electrode at $x = 0$, a bulk reservoir at $x = L$.

## Transport

Seven species are transported. With $c_i$ the concentration, $z_i$ the charge number and
$\gamma_i$ the activity coefficient, the generalised Nernst–Planck flux is

```math
\vec N_i = -D_i \left( \nabla c_i
   + c_i \nabla \left[ \frac{z_i F \phi}{RT} + \ln \gamma_i \right] \right)
```

so each species drifts along the electrostatic potential $\phi$ **and** along gradients of
its own activity coefficient. The $\ln\gamma_i$ term is where the electrolyte model
enters: finite ion size and solvation make $\gamma_i$ a function of the whole
composition, which couples the fluxes to one another.

Continuity carries the homogeneous reactions $R_i$,

```math
\partial_t c_i + \nabla \cdot \vec N_i = R_i ,
```

the potential follows Poisson's equation with the free charge density,

```math
-\nabla \cdot \left( \varepsilon \varepsilon_0 \nabla \phi \right) = F \sum_i z_i c_i ,
```

and a momentum balance fixes the pressure,

```math
\nabla p = -F \left( \sum_i z_i c_i \right) \nabla \phi .
```

## Homogeneous chemistry

Five reversible carbonate reactions act everywhere in the volume:

```math
\begin{aligned}
\ce{CO2 + OH^- &<=> HCO3^-} \\
\ce{HCO3^- + OH^- &<=> CO3^{2-} + H2O} \\
\ce{CO2 + H2O &<=> HCO3^- + H^+} \\
\ce{HCO3^- &<=> CO3^{2-} + H^+} \\
\ce{H2O &<=> H^+ + OH^-}
\end{aligned}
```

The alkaline and acidic routes are thermodynamically consistent
($K_\mathrm{b} = K_\mathrm{a}/K_\mathrm{w}$) but kinetically independent, so the network
can sit far from the water equilibrium locally — and at a driven electrode it does, by
several decades. Anything that infers one of H⁺/OH⁻ from the other through $K_\mathrm{w}$
is fine in the bulk and wrong at the surface.

## Surface chemistry

At $x = 0$ a CatMAP mechanism reduces CO₂ to CO in four elementary steps:

```math
\begin{aligned}
\ce{CO2 + * &<=> CO2*} \\
\ce{CO2* + H2O + e^- &<=> COOH* + OH^-} \\
\ce{COOH* + e^- &<=> CO* + OH^-} \\
\ce{CO* &<=> CO + *}
\end{aligned}
```

The rate law returns a turnover frequency per catalytic site; multiplying by the site
density converts it to a flux per unit electrode area. How that flux is booked onto the
transported species is not a formality — see
[Booking the proton stoichiometry](@ref).

## Repository layout

| Path | Contents |
|:--|:--|
| `src/goldmodel.jl` | `GoldModel` — species layout, rate constants, buffer network, surface reaction, boundary conditions |
| `src/cv.jl` | CV sweep drivers (over $L$, scan rate, pH) |
| `src/iv.jl` | IV sweep drivers (over $L$, CO₂ pressure, pH) |
| `src/dlcap.jl` | double-layer capacitance |
| `src/sweeps_csv.jl` | CSV export |
| `src/AuCO2RR_plots/` | figure layer (`struct.jl`, `cvplot.jl`, `ivplot.jl`, `capsplot.jl`) |
| `notebooks/` | Pluto notebooks driving the above |
| `data/catmap_CO2R_data/` | CatMAP `.mkm` input and formation energies |
| `data/Langmuir_CV_data/` | digitised experimental CVs used as plot backgrounds |

## Validity

- **Potential window.** The CatMAP input declares
  `descriptor_ranges = [[-1.5, 0.0], [298, 298]]`. Sweeping to $+0.8$ V vs. SHE
  extrapolates the surface energetics about 0.8 V past where they were parameterised —
  read the anodic branch as extrapolation.
- **Positivity.** Check it rather than assume it. Scan rate, boundary-layer thickness $L$
  and the anodic vertex all move the boundary between a well-posed problem and one whose
  discrete solution carries negative concentrations, which have no physical reading.
- **Absolute currents.** Linear in the site density $S$, and dependent on which species
  and electron count the current is read from. Curve *shapes* are far more robust than
  magnitudes.

## Building these docs

```julia
using Pkg
Pkg.activate("docs")
Pkg.add(["Documenter", "DocumenterMermaid"])
Pkg.develop(PackageSpec(path = "."))     # dev the package itself into the docs env
include("docs/make.jl")
```

The site lands in `docs/build/`. For a live preview:

```julia
using LiveServer
@async serve(dir = "docs/build")     # re-run include("docs/make.jl") after each edit
```

!!! note "No remote configured"
    `make.jl` sets `remotes = nothing` and disables `edit_link` / `repolink`, so there are
    no "edit on GitHub" links. Once the repository has a GitHub remote, drop those, set
    `repo = Remotes.GitHub("USER", "Capacitance_Code")` and re-enable the `deploydocs`
    block at the bottom of `make.jl`.

```@docs
AuCO2RR
```
