"""
Placeholder for a short summary about {PKGNAME}.
"""
module AuCO2RR
    using LiquidElectrolytes
    using VoronoiFVM
    using ExtendableGrids, ExtendableGrids
    using LessUnitful
    using DelimitedFiles
    using CairoMakie, Colors

    @unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;

    include(joinpath(@__DIR__, "..", "plots", "AuCO2RR_plots.jl"))
    export AuCO2RR_plots

    include("cv.jl")
    export sweep_over_L_cv


    include("dlcap.jl")
    export capscalc

    include("iv.jl")
    export pressure_varied_sweep, sweep_over_L_c, scanrate_varied_sweep, run_pH_sweep
 

end # module