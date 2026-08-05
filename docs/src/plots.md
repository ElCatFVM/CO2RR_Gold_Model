```@meta
CurrentModule = AuCO2RR
```

# Plots

`src/AuCO2RR_plots/`. The module loads Makie, re-declares the unit factors for the figure
layer, and includes `struct.jl`, `cvplot.jl`, `ivplot.jl` and `capsplot.jl`.

## `struct.jl` — shared style

Everything here exists so that no figure hard-codes a style value. Reach for one of these
constants instead of a literal when adding a plot.

| Name | What it is |
|:--|:--|
| `electrochemistry_theme` | The shared Makie `Theme` — fonts, spine and tick weights, label sizes, grid off. |
| `apply_electrochemistry_style!` / `__init__` | Applies the theme globally at load, so figures that skip `with_theme` still match. |
| `PLOT_FONT_FAMILY` | Single sans-serif family for all text. |
| `LW_LINE`, `LW_HIGHLIGHT`, `LW_EXP`, `LW_GUIDE` | Line widths by role. These replaced seven numeric values scattered under six different keyword names. |
| `CMAP_SCANRATE`, `CMAP_PRESSURE` | Colour sequences for a swept family of curves. Index as `CMAP[t]` with `t in range(0, 1, length = n)`. |
| `lab_time`, `lab_voltage`, `lab_current` | Pre-built rich-text axis labels: italic quantity, roman unit. |
| `powlab` | `10^n` as rich text, for log tick labels. |

## `cvplot.jl` — CV figures

### Current and voltage definitions

Every CV figure goes through these four. See [Current and voltage definitions](@ref) for why
`species` defaults to CO.

| Name | What it does |
|:--|:--|
| `faradaic_current` | `sgn · n_e · currents(result, species)`. |
| `capacitive_current` | Reads `result.j_cap`, which exists in every IR-compensation mode. |
| `cv_current` | Faradaic plus, optionally, capacitive. |
| `cv_abscissa` | Returns the voltage values **and** their axis label, for `:applied`, `:reaction_plane` or `:dl`. |

Alongside them: the species index constants `ikplus` … `ico`, and the two extra voltage
labels `lab_voltage_rp` and `lab_voltage_dl`.

### Single CV curves

| Name | What it does |
|:--|:--|
| `plot_cv_current` | One CV curve, colour-graded along the sweep. |
| `plot_cv_current_variedL` | CV curves for a dict of boundary-layer thicknesses. |
| `plot_cv_current_dict` | As above for an arbitrary keyed dict. |
| `plot_cv_over_L` | IV curves over `L` with a cutoff voltage. |
| `CV_overlay_currents` | Faradaic, capacitive and total current for several results on one axis. |
| `CV_total_current` | Capacitive current alone. |
| `CV_dsp_cap_result` | Displacement versus capacitive current, optionally their difference. |
| `plot_cv_total_current` | Faradaic plus capacitive for one result. |
| `plot_cv_total_current_tot` | As above with the three components drawn separately and a blended total colour. |

### Scan-rate families

| Name | What it does |
|:--|:--|
| `plot_scanrate_sweeps` | All scan rates on one axis, with an optional highlighted member. |
| `plot_scanrate_sweeps_cv` | Same, pastel colour ramp. |
| `plot_scanrate_sweeps_cv_2` | Same with the capacitive term included and per-curve annotations. |
| `plot_scanrate_sweeps_cap` | Capacitive current only. |
| `plot_scanrate_sweeps_split` | Broken-axis version: reduction and oxidation windows side by side, each with its own current scale. |
| `plot_cv_scanrate_grid` | 3×3 grid — faradaic, capacitive and total current for three scan rates. |
| `plot_cv_scanrate_grid_unc` | 1×3 grid of total current only. |

### CO₂ pressure families

| Name | What it does |
|:--|:--|
| `plot_pressure_varied_sweep` | Simulation curves over an experimental CSV background. |
| `pressure_varied_cvsweep` | Simulation curves only. |
| `pressure_varied_cvsweep_split` | Broken-axis version; also the shared implementation behind `plot_scanrate_sweeps_split`, which differs only in its legend labelling. |
| `plot_pressure_varied_sweep_ivc` | Pastel overlay with an inset legend. |
| `plot_combined_exp_sim_ivc` | Experiment and simulation on one axis. |

### Concentration and pH

| Name | What it does |
|:--|:--|
| `plot_conc_time_electrode` | All species at the electrode against time, log scale. |
| `plot_activity_time_electrode` | Same for activities rather than concentrations. |
| `plot_conc_profile_logx` | Spatial profile on a log distance axis. |
| `plot_conc_profile_with_delta` | Profile with the diffusion-layer thickness marked. |
| `plot_co2_profiles` | CO₂ profiles at several times. |
| `co2_log_contour` | CO₂ contour over distance and time. |
| `plot_7species_contours` | One contour panel per species. |
| `cv_conc_gif` | Animated concentration profile over a CV. |
| `plottsol` | Generic contour of any unknown from a transient solution. |

### Multi-panel building blocks

Used to assemble the five-panel CV summary.

| Name | Panel |
|:--|:--|
| `panel_conc_time!` | (a) concentrations at the electrode |
| `panel_time_current!` | (b) current |
| `panel_time_voltage!` | (c) applied potential |
| `panel_time_ph!` | (d) surface pH |
| `QoverK` | (e) reaction quotient over equilibrium constant for the two buffer reactions |
| `panel_co2_log_contour!`, `panel_log_contour!` | contour panels |

### Comparison with experiment

| Name | What it does |
|:--|:--|
| `plot_cv_model_vs_koper_facets` | Model against the Koper facet data. |
| `plot_iv_with_experiment` | IV curve with an experimental CSV overlay. |
| `plot_pH_varied_sweep` | Curves for a list of bulk pH values. |
| `plot_time_voltage_and_dt` | Applied protocol and solver time step against time — a solver diagnostic rather than a result plot. |

## `ivplot.jl` — IV figures

| Name | What it does |
|:--|:--|
| `iv_curve_axis` | IV curve with an experimental overlay, log current axis. |
| `conc_vs_voltage_axis` | Electrode concentrations against voltage. |
| `conc_vs_voltage_axis_compare` | As above for two results side by side. |
| `activity_vs_voltage_axis` | Activities against voltage. |
| `electrode_activity_vs_voltage` | Extracts the activities, then calls the above. |
| `plot1d`, `plot1d_makie` | Spatial profile at one voltage. |
| `plot1d_movie` | Animated spatial profile over the sweep. |
| `addplot_ax!`, `addplot_solution!`, `addplot_df!` | Add a solution or a dataframe to an existing axis. |
| `plotcurr_over_L` | Current against voltage for several boundary-layer thicknesses. |
| `plot_iv_with_ringe_refs` | IV curve against the Ringe reference data. |

!!! note "Still on the old current convention"
    The IV figures take `abs.(currents(...))` for their log axes and carry no electron
    count. That is deliberate for a log plot, but it means their magnitudes are not
    directly comparable with the CV figures, which use [`cv_current`](@ref).

## `capsplot.jl` — capacitance figures

| Name | What it does |
|:--|:--|
| `capsplot` | Capacitance against voltage. |
| `capsplot_v` | Capacitance for a named result set. |
| `capsplot_κ` | Capacitance for a range of solvation numbers. |
| `capsplot_fixed` | Fixed-axis variant for figure panels. |
| `capsplot_with_csv` | Capacitance with an experimental CSV overlay. |
| `plot_caps_comparison` | Several capacitance curves on one axis. |
| `plot_κ_sweep_with_refs` | Solvation-number sweep against reference data. |
| `overlay_csv_on_axis!` | Adds a CSV dataset to an existing axis. |
