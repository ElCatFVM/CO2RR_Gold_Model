"""
Continuum simulation of CO₂ reduction on gold.

A one-dimensional generalised Poisson–Nernst–Planck electrolyte — finite ion size,
solvation, activity coefficients — coupled to a CatMAP microkinetic surface mechanism.
Seven transported species run from the electrode at `x = 0` out to a bulk reservoir at
`x = L`, linked by a five-reaction carbonate buffer in the volume and by the CO₂-reduction
microkinetics at the boundary.

Produces cyclic voltammograms, current–voltage curves, double-layer capacitance and
space/time concentration and pH profiles.

Submodules: [`GoldModel`](@ref) builds the electrolyte and boundary conditions,
`AuCO2RR_plots` holds the figure layer.
"""
module AuCO2RR
using LiquidElectrolytes
using VoronoiFVM
using ExtendableGrids, ExtendableGrids
using LessUnitful
using DelimitedFiles, DataFrames, CSV
using DrWatson
using CairoMakie, Colors
using Printf   # added for _pressure_colname / filename (moved from row_interaction_script.jl)


include("goldmodel.jl")
export GoldModel

include("sweeps_csv.jl")
export ivsweep_csv, sweepcomparedir
# added from row_interaction_script.jl (batch 6: CSV export utilities)
export export_scanrate_varied_species_csv_long, export_cv_profile_csv, export_pressure_varied_species_csv_long

@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;

include("AuCO2RR_plots/AuCO2RR_plots.jl")
export AuCO2RR_plots

include("cv.jl")
export sweep_over_L_cv
# added from row_interaction_script.jl (batch 4: extra CV sweep drivers)
export sweep, cvsweep_compensated_over_L, cvsweep_odr_over_L, blthickness, blthickness_t


include("dlcap.jl")
export capscalc

include("iv.jl")
export pressure_varied_sweep, sweep_over_L_c, scanrate_varied_sweep, run_pH_sweep
export ivsweep_over_L
# added from row_interaction_script.jl (batch 5: pressure colname + simulate drivers)
export _pressure_colname, simulate_CO2R, simulate_CO2R_dir

end # module
