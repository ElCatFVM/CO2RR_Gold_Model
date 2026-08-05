```@meta
CurrentModule = AuCO2RR
```

# Standard calculations

The sweep drivers and the CSV export. Each builds a `PNPSystem` from an electrolyte and
hands it to `LiquidElectrolytes`' `cvsweep` / `ivsweep`, varying one parameter at a time.

Every calculation starts from a model and a grid:

```julia
using AuCO2RR, LiquidElectrolytes, ExtendableGrids, LessUnitful
@unitfactors μm

L    = 1000 * μm
X    = ExtendableGrids.geomspace(0, L, 1.0e-7 * μm, 0.01 * L)
grid = ExtendableGrids.simplexgrid(X)

m = GoldModel.create_model(; use_md_hydrated = false, γ_select = "Stefan",
                             ircompensation = NoIRCompensation())

sawtooth = SawTooth(scanrate = 0.05, vmin = -1.2, vmax = 0.8,
                    scanup = false, vstart = 0.0; tstart = 0.0)
```

A single voltammogram is then one call:

```julia
cv = sweep(m, grid, sawtooth; nperiods = 1)
```

## `src/cv.jl` — cyclic voltammetry

| Name | What it does |
|:--|:--|
| `sweep` | Thin wrapper: build a `PNPSystem` from a copied electrolyte, run one `cvsweep`. The one to reach for when sweeping nothing. |
| `sweep_over_L_cv` | CV over a list of boundary-layer thicknesses, building a fresh grid per `L`. **Shadowed by the `iv.jl` definition.** |
| `scanrate_varied_sweep` | CV over a list of scan rates. **Shadowed by the `iv.jl` definition.** |
| `run_pH_sweep` | CV over a list of bulk pH values, optionally re-balancing a counter-ion to keep the bulk electroneutral. |
| `cvsweep_compensated_over_L` | CV over `L` with the ohmic drop compensated by a resistance computed from conductivity and electrode area. |
| `cvsweep_odr_over_L` | CV over `L` on pre-built grids, rescaling `Ru` per `L` when `OhmicDropEstimation` is active. |
| `blthickness` | Outermost position at which a species still deviates from bulk, maximised over all stored times — the diffusion-layer thickness. |
| `blthickness_t` | Same criterion at every stored time; returns `(times, δ)`. |

!!! note "Two caveats in the shadowed scan-rate driver"
    `cv.jl`'s `scanrate_varied_sweep` rebuilds the `SawTooth` without `vstart` / `tstart`,
    so the reconstructed protocol starts from a different potential than the original, and
    it accepts `sweep_kwargs` without forwarding them to `cvsweep`, so solver settings
    passed to it are silently dropped. Both are moot in practice because `iv.jl`'s
    definition wins, but the same two mistakes are easy to reintroduce.

## `src/iv.jl` — steady-state current–voltage

| Name | What it does |
|:--|:--|
| `pressure_varied_sweep` | Scales one species' bulk concentration over a list of pressures and runs the supplied sweep function for each. Returns `(p, result)` pairs — the input shape the pressure plots expect. |
| `sweep_over_L_cv` | Sweep over `L` with a caller-supplied sweep function. Overrides the `cv.jl` definition. |
| `scanrate_varied_sweep` | Sweep over scan rate. Overrides the `cv.jl` definition. |
| `ivsweep_over_L` | IV curves over a list of boundary-layer thicknesses. |
| `simulate_CO2R` | One IV sweep for a given grid and electrolyte. |
| `simulate_CO2R_dir` | As above, writing results into a directory. |
| `_pressure_colname` | Formats a pressure into a CSV column name. |

Because `pressure_varied_sweep` takes the sweep function as an argument, it drives either
a CV or an IV run — pass a closure that builds a `PNPSystem` and calls `cvsweep` to get
pressure-resolved voltammograms.

## `src/dlcap.jl` — double-layer capacitance

| Name | What it does |
|:--|:--|
| `capscalc` | Capacitance over a voltage range for several bulk molarities, with a flag selecting the gold model or a plain reference electrolyte. |

## `src/sweeps_csv.jl` — CSV export

| Name | What it does |
|:--|:--|
| `sweepcomparedir` | Path helper into `data/sweepcompare`. |
| `filename` | Builds a result filename from the sweep parameters. |
| `ivsweep_csv` | Runs an IV sweep and writes voltage and current to CSV. |
| `export_scanrate_varied_species_csv_long` | Long-format CSV of concentrations against scan rate. |
| `export_cv_profile_csv` | CSV of one CV's spatial profile. |
| `export_pressure_varied_species_csv_long` | Long-format CSV of concentrations against CO₂ pressure. |

!!! warning "The exported current is the OH⁻ flux"
    `ivsweep_csv` hard-codes the current column to `currents(result, iohminus)`. OH⁻ takes
    part in the buffer network, so its boundary flux is not purely faradaic, and it is
    also the species whose stoichiometry the electrode reaction re-routes. The plotting
    layer moved to CO for exactly this reason — see
    [Current and voltage definitions](@ref) — but the CSV export has not followed.
