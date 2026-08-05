```@meta
CurrentModule = AuCO2RR
```

# Internal API

`src/goldmodel.jl`, submodule `GoldModel`: everything that defines the physical problem.
Also the parts of the repository that are easy to trip over.

## Types

| Name | What it is |
|:--|:--|
| `SpeciesLayout` | Fixes the index of every unknown plus the `ibuffer*` / `isurface*` slice bounds. See [Notations](@ref). |
| `ReactionData` | All physical constants — buffer rate constants, `pH`, `T`, Henry coefficients, site density `S`, gap capacitance `C_gap`, PZC, reference-electrode position. |
| `BulkSpecies` | One transported species: charge, diffusivity, bulk concentration, effective radius, solvation number, plot colour. |

## Species bookkeeping

| Name | What it does |
|:--|:--|
| `make_species_dict` | Maps display names (`"OH⁻"`, `"CO₂"`, …) to indices. |
| `make_species_dict_catmap` | Maps CatMAP names (`"OH_g"`, `"CO2_aq"`, …) to indices; orders the microkinetic state vector. |
| `make_eneutral` | Builds K⁺ at whatever concentration makes the bulk electroneutral. |

## Electrolyte assembly

| Name | What it does |
|:--|:--|
| `elydata_Au` | Assembles the `ElectrolyteData` for the gold cell from the species list, activity model, constants and IR-compensation mode. |
| `elydata_NaClO₄`, `elydata_NaF` | Two-species reference electrolytes for capacitance benchmarks. |

## Activity coefficient models

Selected by the `γ_select` keyword of `create_model`.

| Name | Model |
|:--|:--|
| `DGML_γ!` | Dreyer–Guhlke–Müller–Landstorfer: solvation and volume exclusion. |
| `Stefan_γ!` | Pure volume exclusion — every species gets the same factor. |
| `Potassium_γ!` | Volume exclusion for K⁺ only; everything else ideal. |

## Reaction networks

| Name | What it does |
|:--|:--|
| `buffer_system` | Builds the five-reaction carbonate network and emits its rate law. |
| `surface_reaction` | Parses the CatMAP `.mkm`, builds the microkinetic network, emits its rate law and the parameter index map. |
| `calc_QBL_local` | Surface charge from the local pressure, when the Robin boundary condition is not selected. |

## Entry point

`create_model` picks the activity model and ion radii, builds the species list, closes
over the volume and boundary reactions, and returns
`(bcondition, reaction, elydata, bulk, bulknames, bulkcolors, species_dict)`.

Four closures live inside it and are not separately callable:

| Closure | Role |
|:--|:--|
| `reaction` | The volume buffer term, evaluated at every node. |
| `we_breactions` | The electrode boundary term — surface charge, local pH, the microkinetic call, and the turnover-frequency-to-flux conversion. |
| `pnp_bcondition` | Boundary conditions handed to `PNPSystem`; also decides whether the potential is imposed directly or left to the IR-compensation machinery. |
| `pb_bcondition` | Poisson–Boltzmann variant, used by the capacitance calculations. |

## From turnover frequency to boundary flux

The microkinetic rate law returns a **turnover frequency**, per catalytic site per second.
Coverages want exactly that; aqueous species want a flux per unit electrode area, so they
are multiplied by the site density:

```math
\underbrace{r}_{\mathrm{s^{-1}}} \times \underbrace{S}_{\mathrm{mol\,sites\,m^{-2}}}
   = \underbrace{j}_{\mathrm{mol\,m^{-2}\,s^{-1}}}
```

With `ReactionData().S ≈ 1.6e-8` mol m⁻² — about $10^{16}$ sites per m², roughly 0.07 % of
the gold surface atom density. $S$ therefore encodes an *active site fraction*, not a
geometric count.

!!! important "`S` scales the whole faradaic current"
    Steady-state coverages do not depend on `S`, so $I_\mathrm{F} = n_e F r S$ is linear
    in it. Matching an absolute measured current is, to a large extent, fitting `S`.

## Booking the proton stoichiometry

The `.mkm` writes both electron transfers with water as the proton donor, so read
literally the reverse (anodic) direction consumes free OH⁻. The electrolyte cannot supply
that:

| Reservoir | Concentration | Max diffusive supply | vs. demand $\sim 10^{-4}$ mol m⁻² s⁻¹ |
|:--|--:|--:|--:|
| OH⁻ | $6\times10^{-5}$ mol m⁻³ | $3\times10^{-9}$ | 1/37 000 |
| H⁺ | $1.6\times10^{-4}$ mol m⁻³ | $1.5\times10^{-8}$ | 1/8 200 |
| H₂O | $5.6\times10^{4}$ mol m⁻³ | $5.5\times10^{-1}$ | 4 600× surplus |

A boundary sink that does not vanish as its own species is depleted admits **no
non-negative solution** once it exceeds the maximum local supply; the Newton solve then
converges to a negative concentration. That is a property of the discrete system, not of
the solver settings — refining the time step, tightening tolerances, refining the grid and
slowing the buffer all leave it unchanged.

`we_breactions` therefore books the stoichiometry on whichever species is being
*produced*, branch-free so the sparsity tracer can walk both paths:

```julia
r_oh  = f[iohminus] * S
f[iohminus] = min(r_oh, zero(r_oh))    # cathodic: water reduced   → OH⁻ produced
f[ihplus]  -= max(r_oh, zero(r_oh))    # anodic:   proton released → H⁺ produced
```

Charge balance is unaffected — one ion per electron either way — and the split between H⁺
and OH⁻ is left to the autoprotolysis reaction already present in the buffer network.
Where water *is* locally equilibrated the two formulations are identical, so the cathodic
branch is unchanged bit for bit.

!!! warning "This is a mechanism choice"
    Microscopic reversibility ties the reverse of a water-donor step to OH⁻ consumption.
    Routing the anodic direction to H⁺ production is the *proton route*, which is not what
    the `.mkm` text says. Either rewrite steps 2–3 there as proton transfers, or feed the
    transported OH⁻ activity into the rate law and accept a smaller anodic current.

## Repository traps

!!! warning "Three files are never loaded"
    `AuCO2RR.jl` includes `goldmodel.jl`, `sweeps_csv.jl`, `AuCO2RR_plots/`, `cv.jl`,
    `dlcap.jl` and `iv.jl`; `AuCO2RR_plots.jl` includes `struct.jl`, `cvplot.jl`,
    `ivplot.jl` and `capsplot.jl`. Everything else is **not** part of the module.

    | File | Contents |
    |:--|:--|
    | `src/overwritten_functions.jl` | older copies of package functions |
    | `src/AuCO2RR_plots/overwritten_plots.jl` | older `plot_cv_current`, `plot_conc_time_electrode`, `plot_pressure_varied_sweep`, `capsplot_κ` |
    | `src/AuCO2RR_plots/koper_figure_plots.jl` | `plot_koper_fig1/3/5` — direct plots of the Langmuir CV data |

!!! warning "Two names are defined twice"
    `cv.jl` and `iv.jl` both define `sweep_over_L_cv` and `scanrate_varied_sweep`, and
    `iv.jl` is included last, so **its** definitions win. The `cv.jl` versions never run.
