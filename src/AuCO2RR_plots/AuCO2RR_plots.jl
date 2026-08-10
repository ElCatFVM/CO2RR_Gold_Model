module AuCO2RR_plots
using LiquidElectrolytes
using VoronoiFVM
using ExtendableGrids, GridVisualize
using LessUnitful
using CairoMakie
using CSV
using DataFrames
using Printf
using LessUnitful
using DelimitedFiles
using CairoMakie, Colors
using LinearAlgebra
# `panel_time_current_diff!` resamples one sweep onto the other's time grid: two adaptive
# runs of the same protocol do not share time points, so a difference needs interpolation.
using Interpolations

@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;


include("cvplot.jl")
export plot_cv_current, plot_time_voltage_and_dt, plot_conc_time_electrode, plot_cv_current, plot_conc_time_electrode, plot_conc_profile_logx, cv_conc_gif, plot_cv_model_vs_koper_facets, plot_pressure_varied_sweep, plot_pH_varied_sweep, plot_iv_with_experiment, plot_scanrate_sweeps, plot_conc_profile_with_delta, pressure_varied_cvsweep
# added from row_interaction_script.jl (batch 1: CV current / scanrate family)
export CV_dsp_cap_result, CV_total_current, plot_cv_total_current_tot, plot_scanrate_sweeps_cv, plot_scanrate_sweeps_cv_2, plot_scanrate_sweeps_cap, plot_cv_scanrate_grid, plot_cv_scanrate_grid_unc
# added from row_interaction_script.jl (batch 2: contour / multi-panel / overlay family)
export plot_co2_profiles, co2_log_contour, plot_activity_time_electrode, CV_overlay_currents, plot_cv_current_variedL, plot_pressure_varied_sweep_ivc, plot_cv_current_dict, panel_conc_time!, panel_time_current!, panel_time_voltage!, panel_co2_log_contour!, panel_time_ph!, plot_7species_contours, panel_log_contour!, QoverK, plot_cv_total_current, plot_combined_exp_sim_ivc, powlab
# single source of truth for the CV current / voltage definitions
export faradaic_current, capacitive_current, cv_current, cv_abscissa
export pressure_varied_cvsweep_split, plot_scanrate_sweeps_split
export plot_exp_sim_cvsweep_split
# shared broken-axis building blocks, used by all three split plotters
export split_axis_pair!, split_series!, annotate_split!, split_widths!
export reaction_plane_potential, electrode_potential, compensation_potential, driving_force
export double_layer_capacitance, plot_ircomp_compare, plot_cv_summary
export panel_time_current_diff!
export lab_time, lab_voltage, lab_current, lab_current_co, lab_pressure, lab_scanrate
export lab_current_F, lab_current_C, lab_current_tot
export lab_current_abs, lab_distance, lab_conc_surface, lab_activity_surface
export lab_capacitance, SPECIES_COLORS
# CSV-driven polarization panels
export species_columns, panel_polarization!, panel_species_vs_voltage!, panel_dlcap!
export plot_iv_summary, plot_iv_summary_from_result, species_rich
export species_table, polarization_table, surface_conc_table, surface_activity_table
export lab_time_short, lab_voltage_short, lab_current_short


include("capsplot.jl")
export plot_caps_comparison, capsplot_v, capsplot, capsplot_κ, plot_κ_sweep_with_refs

include("ivplot.jl")
export iv_curve_axis, conc_vs_voltage_axis, addplot_ax!, plot1d_makie, conc_vs_voltage_axis_compare, ivsweep_over_L, addplot, plot1d, plotcurr_over_L, plot_iv_with_ringe_refs
# added from row_interaction_script.jl (batch 3: electrode activity vs voltage)
export electrode_activity_vs_voltage, activity_vs_voltage_axis

include("struct.jl")
export electrochemistry_theme
# Revise does not re-run `__init__`, so this has to be callable by hand after a style edit.
export apply_electrochemistry_style!

end
