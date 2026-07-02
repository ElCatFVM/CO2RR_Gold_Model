### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 68e52c7a-6922-11f1-8634-f109a9649529
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

# ╔═╡ 0a53e71e-5f80-46d0-bb44-2a89f1f914fa
begin
	using AuCO2RR
	include(joinpath(pkgdir(AuCO2RR), "plots", "AuCO2RR_plots.jl"))
	using .AuCO2RR_plots
end

# ╔═╡ 0a5f93b7-4f7c-4049-8a27-dceddc40da97
begin
	const Γ_we 		= 1
	const Γ_bulk 	= 2
	const L = 80 * ufac"μm"
end

# ╔═╡ 7f013b2a-9ab4-4ab4-8129-28c365a7a13d
begin
    hmin = 1.0e-6 	* ufac"μm"
    hmax = 1.0 		* ufac"μm"
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ 0cc0e806-7a04-4254-9ecb-831df5ae7328
elystruct = GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = :none)

# ╔═╡ 4ea455cd-5cd9-40f0-b75f-e7f7f0e7ad96
solver_control = (; max_round 	= 4,
					maxiters 	= 20,
              		tol_round 	= 1.0e-9,
              		verbose 	= "a",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)

# ╔═╡ 124f8548-0619-4911-b6e0-58d35ac768bd
function simulate_CO2R(grid, celldata, bcondition, reaction; voltages = (-1.5:0.1:0.0) * ufac"V", kwargs...)
	kwargs 	 	= merge(solver_control, kwargs) 
    cell        = PNPSystem(grid; bcondition=bcondition, reaction=reaction, celldata)
	ivresult    = ivsweep(cell; voltages, store_solutions=true, kwargs...)

	cell, ivresult
end;

# ╔═╡ 4c313b9a-fc0d-4585-9c4b-dcc48b8664d6
cell, ivresult = simulate_CO2R(grid, elystruct.elydata, elystruct.bcondition, elystruct.reaction)

# ╔═╡ bf704270-29bd-4782-9682-b03a70c74af7
findfirst(==("CO₂"), elystruct.bulknames)

# ╔═╡ 7195c5eb-e956-4b43-ad71-201f3e6e7467
function iv_curve_axis(ivresult;
    cutoff = -0.4,
    showlegend = true,
    species = iohminus,
    data_dir = "../data/catmap_CO2R_data",
)
    v_all = ivresult.voltages
    mask  = v_all .< cutoff
    volts = v_all[mask]
    I_sim = abs.(currents(ivresult, species))[mask] .* ufac"cm^2/mA"

    table2 = readdlm(joinpath(data_dir, "Ringe-theorical.csv"),    ',', Float64, '\n')
    table3 = readdlm(joinpath(data_dir, "Ringe-experimental.csv"), ',', Float64, '\n')
    df2_v = table2[:, 1]; df2_I = max.(abs.(table2[:, 2]), eps(Float64))
    df3_v = table3[:, 1]; df3_I = max.(abs.(table3[:, 2]), eps(Float64))

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (960, 540))
        a = Axis(f[1, 1];
            xlabel = L"\textbf{Voltage}\ \;\; U\ \text{(V vs. SHE)}",
            ylabel = L"\textbf{Partial Current}\ \;\; I_{CO}\ \text{(mA/cm}^2\text{)}",
            yscale = log10,
            limits = ((-1.5, -0.50), (1e-11, 1e2)),
        )
        return f, a
    end

    xt = [-1.4, -1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])
    ax.yticks = (10.0 .^ (0:-3:-9),
                 [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"])

    ax.xlabelfont     = :regular
    ax.ylabelfont     = :regular
    ax.xticklabelfont = :regular
    ax.yticklabelfont = :regular
    ax.xlabelsize     = 26
    ax.ylabelsize     = 26
    ax.xticklabelsize = 24
    ax.yticklabelsize = 24
    ax.spinewidth     = 5.5
    ax.xtickwidth     = 2.0
    ax.ytickwidth     = 2.0
    ax.xticksize      = 8
    ax.yticksize      = 8
    ax.xlabelpadding  = 10
    ax.ylabelpadding  = 10
    ax.xgridvisible   = false
    ax.ygridvisible   = false

    lines!(ax, volts, I_sim;
        color     = "#008b3f",
        linewidth = 8,
        label     = "LiquidElectrolyte.jl"
    )
	lines!(ax, df2_v, df2_I;
	    color     = ("#0083fe", 0.6),
	    linewidth = 10,
	    label     = "CatINT"
	)
    scatter!(ax, df3_v, df3_I;
        marker     = :circle,
        markersize = 15,
        color      = "#e52c40",
        label      = "Experiment"
    )

    if showlegend
        axislegend(ax;
            position     = :rt,
            labelsize    = 24,
            titlesize    = 24,
            labelfont    = :regular,
            titlefont    = :regular,
            framevisible = false,
            patchsize    = (40, 20)
        )
    end

    return fig
end

# ╔═╡ 9a7ec635-851f-4d76-9a6f-e9731aa35626
iv_curve_axis(ivresult; cutoff=-0.4, showlegend=true, species = findfirst(==("OH⁻"), elystruct.bulknames))

# ╔═╡ Cell order:
# ╠═68e52c7a-6922-11f1-8634-f109a9649529
# ╠═0a53e71e-5f80-46d0-bb44-2a89f1f914fa
# ╠═0a5f93b7-4f7c-4049-8a27-dceddc40da97
# ╠═7f013b2a-9ab4-4ab4-8129-28c365a7a13d
# ╠═0cc0e806-7a04-4254-9ecb-831df5ae7328
# ╠═4ea455cd-5cd9-40f0-b75f-e7f7f0e7ad96
# ╠═124f8548-0619-4911-b6e0-58d35ac768bd
# ╠═4c313b9a-fc0d-4585-9c4b-dcc48b8664d6
# ╠═bf704270-29bd-4782-9682-b03a70c74af7
# ╠═9a7ec635-851f-4d76-9a6f-e9731aa35626
# ╟─7195c5eb-e956-4b43-ad71-201f3e6e7467
