```@meta
CurrentModule = AuCO2RR
```

# Cyclic voltammetry walkthrough

A guided tour of `notebooks/CyclicVoltammetry_notebook.jl`, in the order the notebook runs.
Each section states what the figure is for, the call that produces it, and what to check
before believing it.

Figures are dropped into `docs/src/assets/` and referenced below. A placeholder that shows a
broken image means that file has not been exported yet — the `![](...)` line names the file
each section expects.

!!! info "Run the notebook, not this page"
    Nothing here is executed at build time. The code blocks are transcriptions of the
    notebook cells, kept so the page reads on its own; the notebook is the source of truth.

## Setup

Every figure on this page starts from one model, one grid and one potential protocol.

```julia
using AuCO2RR, LiquidElectrolytes, ExtendableGrids, LessUnitful, CairoMakie
using AuCO2RR.AuCO2RR_plots
@unitfactors μm

L    = 1000 * μm
X    = ExtendableGrids.geomspace(0, L, 1.0e-7 * μm, 0.01 * L)
grid = ExtendableGrids.simplexgrid(X)

sawtooth = SawTooth(scanrate = 0.05, vmin = -1.2, vmax = 0.8,
                    scanup = false, vstart = 0.0; tstart = 0.0)
```

Two electrolytes are built, differing only in how the ohmic drop is treated. Every
compensated-vs-uncompensated figure below compares these two:

```julia
elystruct_unc = GoldModel.create_model(; ircompensation = NoIRCompensation())
elystruct_odr = GoldModel.create_model(;
    ircompensation = OhmicDropEstimation(Ru, ico, 2, factor))
```

!!! warning "Style changes need an explicit call"
    The Makie theme is applied once, from `AuCO2RR_plots.__init__`. Revise does not re-run
    `__init__`, so after editing `src/AuCO2RR_plots/struct.jl` the session still holds the
    old theme. Call `AuCO2RR_plots.apply_electrochemistry_style!()` and re-plot.

## Choosing the time grid

The solver's adaptive stepping is the right default for a single run and the wrong one for
any figure that compares two runs: two sweeps of the *same* protocol land on different time
points, so a difference between them is dominated by the grid rather than by the physics.

Pin the step whenever runs are compared:

```julia
Δt_cv = 0.05
cv = sweep(elystruct_odr, grid, sawtooth; nperiods,
           Δt_min = Δt_cv, Δt_max = Δt_cv,
           abstol = 1.0e-12, reltol = 1.0e-12)
```

For a family swept over scan rate this has a cost that is easy to miss. A cycle covers a
fixed voltage distance, so its duration scales as `1/scanrate`, and a fixed `Δt` makes the
slowest run the most expensive *and* the most over-resolved:

| scan rate (V s⁻¹) | cycle (s) | steps at `Δt = 0.05` | V per step |
|--:|--:|--:|--:|
| 0.005 | 800 | 16 000 | 0.00025 |
| 0.01 | 400 | 8 000 | 0.0005 |
| 0.02 | 200 | 4 000 | 0.001 |
| 0.05 | 80 | 1 600 | 0.0025 |
| 0.1 | 40 | 800 | 0.005 |

To give every member of a scan-rate family the same *voltage* resolution instead, scale the
step with the rate — `Δt = Δu / scanrate`. Use that for scan-rate figures and the fixed `Δt`
for compensation comparisons, where a shared time grid is the whole point.

## The five-panel summary

The standard single-run figure: everything measured at the electrode against time, stacked
on one time axis.

```julia
plot_cv_summary(cv, elystruct_odr)
```

![Five-panel CV summary](assets/cv_summary.png)

*Caption placeholder — describe the run: scan rate, `L`, compensation mode, `Δt`.*

Reading order, top to bottom:

| Panel | Shows | Watch for |
|:--|:--|:--|
| (a) | Surface concentrations, log scale | Any species reaching the floor of the axis. Negative concentrations are clipped by the log scale, not fixed by it. |
| (b) | CO partial current | The cathodic peak and, an order of magnitude smaller, the anodic one. |
| (c) | Electrode potential | The applied protocol is dotted behind it; the gap between the two is the compensation term. |
| (d) | Surface pH | Follows (b) with a lag set by buffer transport. |
| (e) | Reaction quotients over their equilibria | `Q/K = 1` is the dashed line. Excursions mark where the buffer is driven out of equilibrium. |

!!! info "Current is read from CO"
    CO is the only species in this model whose boundary flux is purely faradaic — every
    other transported species also takes part in the carbonate buffer, so its flux carries
    buffer-driven transport that is not current. See
    [Current and voltage definitions](plots.md#Current-and-voltage-definitions).

## IR compensation

### Compensated against uncompensated

```julia
plot_ircomp_compare(cv_odr, elystruct_odr;
                    reference = cv_unc, reference_model = elystruct_unc)
```

![Driving force, electrode potential and ohmic drop](assets/ircomp_compare.png)

*Caption placeholder — state `Ru`, the compensation factor and both electrolytes.*

Panels are selected with `panels`; the default is `(:driving, :metal_time, :ircomp)`.

!!! warning "`reference_model` is not optional"
    Read a reference run through the *compensated* model and `result.voltages` means
    something different than it does for that run — the reference then plots as a straight
    line and looks like a physical difference that is not there. Pass the reference's own
    model.

### The residual

The two currents overlap almost everywhere, so the informative quantity is their difference:

```julia
panel_time_current_diff!(f, f[1, 1], cv_odr, cv_unc, elystruct_odr;
                         include_capacitive = false)
```

![Faradaic current difference](assets/ircomp_diff.png)

*Caption placeholder — note that this is the faradaic difference only.*

`include_capacitive = false` is deliberate: the compensated run carries a capacitive term
the uncompensated one effectively does not, so including it makes the difference mostly that
term rather than a difference in the reaction.

### Sweeping the compensation factor

```julia
results = cvsweep_over_ircompfactor(elystruct_odr, grid, sawtooth;
                                    factors = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9],
                                    Δt_fixed = 0.05,
                                    abstol = 1.0e-12, reltol = 1.0e-12)
```

![Compensation factor family](assets/ircomp_factors.png)

*Caption placeholder — list the factors and the fixed `Δt`.*

!!! info "`factor = 0` reproduces the uncompensated run exactly"
    At `factor = 0` the PNP subsystem is identical to the uncompensated one — `ϕ_DL = 0`
    gives the same boundary condition and the same surface reactions, and the extra unknowns
    are slaved with no feedback. Any residual difference is the time grid. With
    `Δt_min = Δt_max = 0.05` and `abstol = reltol = 1e-12` on both, the relative difference
    measured on this model is `8.7e-13`. If you see more than that, the two runs are not on
    the same grid.

    This is the cheapest available check that a compensation change did what you think it
    did — run it before trusting any other factor.

## Parameter families

Both family figures use the broken voltage axis: the reductive window on the left, the
oxidative one on the right, with the featureless region between them cut out. The anodic
feature is one to two orders of magnitude below the reduction peak and is unreadable on a
single axis.

### Scan rate

```julia
plot_scanrate_sweeps_split(SR_vec, scanrates)
```

![Scan-rate family](assets/scanrate_split.png)

*Caption placeholder — quote both windows, both current ranges and `anodic_gain`.*

### CO₂ partial pressure

```julia
pressure_varied_cvsweep_split(P_recs)
```

![Pressure family](assets/pressure_split.png)

*Caption placeholder — as above, plus the pressures.*

!!! warning "Nothing is comparable across the seam"
    The two panels share neither a current scale nor a continuous voltage axis. Heights and
    slopes may not be compared across the break, and `anodic_gain ≠ 1` magnifies the right
    panel further — which is why any gain other than 1 is written into its title.

## Comparison with experiment

```julia
plot_exp_sim_cvsweep_split(P_recs)
```

![Measured against simulated CVs](assets/exp_sim_split.png)

*Caption placeholder — cite the measurement source and state `sim_gain`.*

Experiment occupies the top row, theory the bottom, both on the same broken axis and with
linked voltage axes. Each row keeps its own current scale: the two do not have to agree in
magnitude for their shapes to be compared, which is what the figure is for.

Colour means pressure and only pressure — the rows are told apart by depth of ink
(`CMAP_PRESSURE_EXP` against `CMAP_PRESSURE`, both listed in [Plots](plots.md)) so a given
pressure keeps its hue in both.

!!! warning "`sim_gain` defaults to 1"
    The older `plot_combined_exp_sim_ivc` multiplies its simulated current by a hard-coded
    factor of 2 with no comment. That factor is not carried over. Pass `sim_gain = 2` to
    reproduce the old figure, and say so in the caption if you do.

## Known open points

Carry these into any caption written from this notebook.

- **Mechanism mismatch.** The boundary reaction in `goldmodel.jl` books its proton
  stoichiometry on the species being produced, which makes water the donor. The `.mkm` input
  still describes the OH⁻-acceptor route. The two are not equivalent under microscopic
  reversibility — see [Internal API](internal.md#Internal-API).
- **Site density.** The turnover-to-flux conversion uses `S ≈ 1.6e-8 mol m⁻²`, about 0.07 %
  of the gold surface atoms. `I_F` scales linearly with `S`, so any absolute current quoted
  from these figures inherits that choice.
- **CSV export lags the plots.** `ivsweep_csv` still exports the OH⁻ flux as its current
  column while the plotting layer has moved to CO.
