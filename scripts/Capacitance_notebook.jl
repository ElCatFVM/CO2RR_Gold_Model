### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 71af95dc-6924-11f1-8bc7-0dcb68d28002
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	using Revise
    using LiquidElectrolytes
	using CatmapInterface
	using Catalyst#: unknowns
	using VoronoiFVM
	using LessUnitful
	using ExtendableGrids, GridVisualize
	using DelimitedFiles
	using Interpolations
	using PlutoUI, HypertextLiteral
	using PreallocationTools
	using Latexify
	using Catalyst
	using Printf
	using Test
	using LinearAlgebra
	using Colors
	using FileIO
	using CSV, DataFrames
	if isdefined(Main,:PlutoRunner)
        using CairoMakie	
   		default_plotter!(CairoMakie)
 		CairoMakie.activate!(type="svg")
    end
end;

# ╔═╡ cb404d54-7a38-4886-b8da-97e20a6a3ce2
begin
	using AuCO2RR
	include(joinpath(pkgdir(AuCO2RR), "plots", "AuCO2RR_plots.jl"))
	using .AuCO2RR_plots
end

# ╔═╡ c2fc186a-d0b6-42a4-8c84-97d03f0cea97
begin
	const Γ_we 		= 1
	const Γ_bulk 	= 2
	const L = 80 * ufac"μm"
end

# ╔═╡ 3d068002-b866-4957-91a3-262df3e2eb37
begin
    hmin = 1.0e-6 	* ufac"μm"
    hmax = 1.0 		* ufac"μm"
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ 34573291-da44-4424-ba27-90d02bd4a151
elystruct = AuCO2RR.GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = :none)

# ╔═╡ d34b83d6-1a18-4696-9600-9a42ef3d2f71
begin	
	sys_pnp = PNPSystem(grid; bcondition = elystruct.bcondition, celldata = deepcopy(elystruct.elydata), elystruct.reaction)
	result_pnp = capscalc(sys_pnp, false; vrange = range(-1.0, 1.0, length = 201))
end

# ╔═╡ 10f73411-78fe-4dab-b572-06dad2031c40
begin
		fig_pnp, ax_pnp = AuCO2RR_plots.capsplot_fixed(result_pnp, "Poisson-Nernst-Planck", xlimits_L=(-1.0, 1.0), ylimits_L=(0, 200))
		fig_pnp
end

# ╔═╡ Cell order:
# ╠═71af95dc-6924-11f1-8bc7-0dcb68d28002
# ╠═cb404d54-7a38-4886-b8da-97e20a6a3ce2
# ╠═c2fc186a-d0b6-42a4-8c84-97d03f0cea97
# ╠═3d068002-b866-4957-91a3-262df3e2eb37
# ╠═34573291-da44-4424-ba27-90d02bd4a151
# ╠═d34b83d6-1a18-4696-9600-9a42ef3d2f71
# ╠═10f73411-78fe-4dab-b572-06dad2031c40
