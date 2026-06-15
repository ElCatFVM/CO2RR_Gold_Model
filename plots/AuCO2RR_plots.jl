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

    @unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;


    include("cvplot.jl")
    export plot_cv_current, plot_time_voltage_and_dt, plot_conc_time_electrode, plot_cv_current, plot_conc_time_electrode, plot_conc_profile_logx, cv_conc_gif, plot_cv_model_vs_koper_facets, plot_pressure_varied_sweep, plot_pH_varied_sweep, plot_iv_with_experiment, plot_scanrate_sweeps, plot_conc_profile_with_delta, pressure_varied_cvsweep 


    include("capsplot.jl")
    export plot_caps_comparison, capsplot_v, capsplot, capsplot_κ, plot_κ_sweep_with_refs

    include("ivplot.jl")
    export iv_curve_axis, conc_vs_voltage_axis, addplot_ax!, plot1d_makie, conc_vs_voltage_axis_compare, ivsweep_over_L, addplot, plot1d, plotcurr_over_L, plot_iv_with_ringe_refs

    include("struct.jl")
    export electrochemistry_theme

end # module