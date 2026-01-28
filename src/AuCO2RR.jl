"""
Placeholder for a short summary about {PKGNAME}.
"""
module AuCO2RR
    using LiquidElectrolytes
    using VoronoiFVM
    using ExtendableGrids
    using LessUnitful


    include(joinpath(@__DIR__, "..", "plots", "AuCO2RR_plots.jl"))
    export AuCO2RR_plots

    include("cv.jl")

    include("dlcap.jl")
    export capscalc

    include("iv.jl")
    export pressure_varied_sweep, sweep_over_L, scanrate_varied_sweep, run_pH_sweep
 

end # module