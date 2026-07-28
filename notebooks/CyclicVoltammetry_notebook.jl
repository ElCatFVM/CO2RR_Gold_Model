### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 8d20515c-54c6-11f1-aeac-bdc0c7e51b18
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	using Revise
    using LiquidElectrolytes
	using AuCO2RR, AuCO2RR.AuCO2RR_plots
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

# ╔═╡ 8e524886-b0fb-49b0-b6fa-f3bfcc2a1bd7
TableOfContents()

# ╔═╡ ea90b4a5-6709-4a85-87de-b71fe87e71f7
md"""
## Setup
"""

# ╔═╡ 4b663702-b895-4a84-b065-0673e287b0ca
md"""
### Units
"""

# ╔═╡ 0854d2c7-1b61-46b0-aa0a-496cdc8439f2
begin
	@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;
	@phconstants N_A c_0 k_B e h ε_0 R
	
end

# ╔═╡ b24b7697-2c64-4e7a-8f7f-87cab6add489
md"""
### Data
"""

# ╔═╡ 96af1c36-0fab-4c1d-bad1-96b98a927d66
begin
	const Γ_we 		= 1
	const Γ_bulk 	= 2
	const L_blinf = 2500 * μm
end

# ╔═╡ d011050d-c66e-4e8a-a173-f55679403a90
begin
	sawtooth = SawTooth(
	        scanrate = 1,
#	        scanrate = 0.05,
	        vmin     = -1.2, 
	        vmax     = 0.8,
	        scanup   = false,
			vstart = 0.0; tstart = 0.0
	    )
	const nperiods = 1
end

# ╔═╡ 035cf151-b62a-42ee-8b03-38e68fc4e4b3
begin
    #Vmax = 2 * V
    L = 7500 * μm
    hmin = 1.0e-7 	* μm
    hmax = 0.1*L
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ bcfc1e78-c1e9-4538-890a-35e74bfc4060
elystruct_unc = GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = :none)

# ╔═╡ 2fe95a9f-202e-418a-b422-cef1ef013c9d
md"""
### Set IR Compensation
"""

# ╔═╡ 8626e977-a77f-4646-be99-33e713dbf593
ircompensation=:pseudopotentiostat

# ╔═╡ 0963720a-e310-45b2-92a5-a9e5bc3e6888
begin
    elystruct_odr = GoldModel.create_model(; use_md_hydrated = false, γ_select = "Stefan", ircompensation)
    ely_odr = elystruct_odr.elydata
    ely_odr.Ru           = L / LiquidElectrolytes.conductivity(ely_odr, ely_odr.c_bulk)
    
    elystruct_odr
end

# ╔═╡ f4f59329-e817-495a-9e83-1ab53e7738a9
function sweep(model, grid, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
    celldata = deepcopy(model.elydata)   # copy the electrolyte → each run is independent
    pnpcell = PNPSystem(grid; bcondition = model.bcondition, celldata = celldata, reaction = model.reaction)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
		kwargs...
    )

end

# ╔═╡ 590a17bd-c98d-4b23-9139-04b550c95efe
pnpcell_odr = PNPSystem(grid; bcondition = elystruct_odr.bcondition, celldata = elystruct_odr.elydata, reaction = elystruct_odr.reaction)

# ╔═╡ 1e59ac64-17d4-42a0-befc-32961332c4c7
pnpcell_unc = PNPSystem(grid; bcondition = elystruct_unc.bcondition, celldata = elystruct_unc.elydata, reaction = elystruct_unc.reaction)

# ╔═╡ e2eb4160-550a-4a07-b4c8-66863aaab46f
cv_odr = sweep(elystruct_odr, grid, sawtooth; nperiods = nperiods)

# ╔═╡ 47515ef3-b6aa-49d0-b4fc-b471cb947aa1
cv_unc = sweep(elystruct_unc, grid, sawtooth; nperiods = nperiods)

# ╔═╡ f77ec140-e091-46c1-8bc6-1c00f3880e83
AuCO2RR_plots.plot_conc_time_electrode(cv_odr, elystruct_odr)

# ╔═╡ bd9b5c58-0375-4ce2-aa55-c74b921aa050
let
    result = cv_odr
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 1260))       
        ax1, leg1 = AuCO2RR_plots.panel_conc_time!(f, f[1, 2], result, elystruct_odr;     xlabel = "")
        ax2 = AuCO2RR_plots.panel_time_current!(f, f[2, 2], result, elystruct_odr;      xlabel = "")
        ax3 = AuCO2RR_plots.panel_time_voltage!(f, f[3, 2], result;             xlabel = "")
        ax4 = AuCO2RR_plots.panel_time_ph!(f, f[4, 2], result;                  xlabel = "")   # (d) pH
        #ax5, hm = panel_co2_log_contour!(f, f[5, 2], f[5, 3], result, X, bulk; L_val = L)  # (e)
        ax5, leg5 = AuCO2RR_plots.QoverK(f, f[5, 2], result,elystruct_odr; legend_pos = f[5, 3])       # (f) Q/K

        labels = ["(a)", "(b)", "(c)", "(d)", "(e)"]
        for i in 1:5
            Label(f[i, 1], labels[i],
                fontsize = 25, font = :bold, valign = :top,
                padding = (0, -20, -10, 0))
        end

        hidexdecorations!(ax1, grid = false)
        hidexdecorations!(ax2, grid = false)
        hidexdecorations!(ax3, grid = false)
        hidexdecorations!(ax4, grid = false)
        #hidexdecorations!(ax5, grid = false)     # hide x-labels on the contour too (only ax6 shows the x-axis)

        ax3.yticks = LinearTicks(3)
        ax4.yticks = LinearTicks(4)
        ylims!(ax1, 1e-14, 1e2)

        linkxaxes!(ax1, ax2, ax3, ax4, ax5)
        rowgap!(f.layout, 15)
        for r in 1:5
			if r == 1
				rowsize!(f.layout, r, 350)
			else
            	rowsize!(f.layout, r, Relative(1/6))   # split into 6
			end
        end
        f
    end
    fig
end

# ╔═╡ 01f688a7-3265-4bc5-9b74-6eedfc36476d
begin
	grid_dict = Dict{Float64, Any}()
	for Lv in [1000, 2500, 5000, 10000, 15000, 20000] .* μm
	    Xg = ExtendableGrids.geomspace(0, Lv, 1.0e-7*μm, Lv*0.1)
	    grid_dict[Lv] = ExtendableGrids.simplexgrid(Xg)
	end
	grid_dict
end

# ╔═╡ cf4713e6-706c-484f-b398-8ba6cc33561b
begin
	results = cvsweep_odr_over_L(
	    elystruct_odr.elydata,       
	    grid_dict,
	    elystruct_odr.bcondition,
	    elystruct_odr.reaction,
	    sawtooth;
		ircompensation,
	    nperiods,
		Δu_opt=0.025,
		Δt_min=1.0e-7
	)
end

# ╔═╡ affdc880-3a70-429a-bcfb-a53cf7ed1f31
let
	fig=AuCO2RR_plots.plot_cv_current_variedL(results, elystruct_odr; species = 6)
	CairoMakie.save("cv-$(ircompensation).png",fig)
	fig
end

# ╔═╡ 6301f323-16d8-4805-bbcf-4a055bca2d59
blthickness(grid, elystruct_unc.elydata, cv_unc.tsol; species = 5) / μm

# ╔═╡ 6f56453e-384b-4be6-9298-8101d12d8b79
blthickness(grid, elystruct_odr.elydata, cv_odr.tsol; species = 5) / μm

# ╔═╡ 9bd12311-fbe8-436b-8c14-18edb5744929
function panel_time_current_diff!(fig, panel_pos, result_1, result_2, m;
                             redox_species = nothing,
                             include_capacitive = true,
                             scale = cm^2/mA,
                             sign = -1,
                             color_gradient = false,
                             lw = 4,
                             label_1 = "odr",
                             label_2 = "unc",
                             xlims = nothing,          # ★ (xmin, xmax) or nothing
                             ylims = nothing)          # ★ (ymin, ymax) or nothing
    model = m.elydata
    xlabel = "Curr"
    ax = Axis(panel_pos; xlabel = xlabel,
              ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"))
    redox = Dict(model.cspecies[5] => 2)
    function calc_total_current(res)
        n_t = length(res.times)
        I_F = zeros(n_t)
        for (idx, n_e) in redox
            I_F .+= n_e .* currents(res, idx)
        end
        I_C = zeros(n_t)
        if include_capacitive && model.ircompensation == :ohmicdrop
            I_C = [u[model.icc, 1] for u in res.tsol[1:end-1]]
        end
        return res.times, sign .* (I_F .+ I_C) .* scale
    end
    t1, I_1 = calc_total_current(result_1)
    t2, I_2 = calc_total_current(result_2)
    # the two sweeps have different time grids, so interpolate I_2 onto t1
    itp2   = Interpolations.linear_interpolation(t2, I_2; extrapolation_bc = Interpolations.Line())
    I_diff = I_1 .- itp2.(t1)
    lines!(ax, t1, I_diff; color = :red, linewidth = lw)

    # ★ set axis limits (only when given)
    xlims === nothing || xlims!(ax, xlims...)
    ylims === nothing || ylims!(ax, ylims...)

    return ax
end

# ╔═╡ 96eb220e-d88c-4f3c-872e-174938b19a8c
let
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 400))
        panel_time_current_diff!(
            f, f[1, 1],
            cv_odr,            # result_1  (label_1 = "odr")
            cv_unc,            # result_2  (label_2 = "unc")
            elystruct_odr;     # m  → m.elydata (cspecies, capacitive-term mode)
            include_capacitive = false,   # ★ recommend false when the modes differ (see note below)
			#xlims = (0, 10),
         #   ylims = (-1e-8, 1e-8)
        )
        f
    end
    fig
end

# ╔═╡ Cell order:
# ╠═8d20515c-54c6-11f1-aeac-bdc0c7e51b18
# ╠═8e524886-b0fb-49b0-b6fa-f3bfcc2a1bd7
# ╟─ea90b4a5-6709-4a85-87de-b71fe87e71f7
# ╟─4b663702-b895-4a84-b065-0673e287b0ca
# ╠═0854d2c7-1b61-46b0-aa0a-496cdc8439f2
# ╟─b24b7697-2c64-4e7a-8f7f-87cab6add489
# ╠═96af1c36-0fab-4c1d-bad1-96b98a927d66
# ╠═d011050d-c66e-4e8a-a173-f55679403a90
# ╠═035cf151-b62a-42ee-8b03-38e68fc4e4b3
# ╠═bcfc1e78-c1e9-4538-890a-35e74bfc4060
# ╟─2fe95a9f-202e-418a-b422-cef1ef013c9d
# ╠═8626e977-a77f-4646-be99-33e713dbf593
# ╠═0963720a-e310-45b2-92a5-a9e5bc3e6888
# ╠═f4f59329-e817-495a-9e83-1ab53e7738a9
# ╠═590a17bd-c98d-4b23-9139-04b550c95efe
# ╠═1e59ac64-17d4-42a0-befc-32961332c4c7
# ╠═e2eb4160-550a-4a07-b4c8-66863aaab46f
# ╠═47515ef3-b6aa-49d0-b4fc-b471cb947aa1
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═bd9b5c58-0375-4ce2-aa55-c74b921aa050
# ╠═96eb220e-d88c-4f3c-872e-174938b19a8c
# ╠═01f688a7-3265-4bc5-9b74-6eedfc36476d
# ╠═cf4713e6-706c-484f-b398-8ba6cc33561b
# ╠═affdc880-3a70-429a-bcfb-a53cf7ed1f31
# ╠═6301f323-16d8-4805-bbcf-4a055bca2d59
# ╠═6f56453e-384b-4be6-9298-8101d12d8b79
# ╠═9bd12311-fbe8-436b-8c14-18edb5744929
