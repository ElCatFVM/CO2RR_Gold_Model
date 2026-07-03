### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 8d20515c-54c6-11f1-aeac-bdc0c7e51b18
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

# ╔═╡ cd3522a0-904f-4b89-b454-fcd62ed8de52
begin
	using AuCO2RR
	include(joinpath(pkgdir(AuCO2RR), "plots", "AuCO2RR_plots.jl"))
	using .AuCO2RR_plots
end

# ╔═╡ 5471cc8c-f4e9-41bb-a0ec-a8c5125a2b4d
Pkg.status()

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
	        scanrate = 0.05,
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
    L = 2500 * μm
    hmin = 1.0e-6 	* μm
    hmax = 50 * μm
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ 0963720a-e310-45b2-92a5-a9e5bc3e6888


# ╔═╡ e7e23a3a-2383-42e0-b5cd-0b315ad7d3e9
begin
	elystruct = GoldModel.create_model(; use_md_hydrated=false, γ_select="Stefan", ircompensation=:none)
	
	ely = elystruct.elydata
	ely.x_ref        = [X[end], 0.0, 0.0]    
	ely.ircompfactor = 0                  
	R_u = L / conductivity(ely, ely.c_bulk)   
	
	ely_pts = deepcopy(ely); ely_pts.ircompensation = :pseudopotentiostat     
	ely_odr = deepcopy(ely); ely_odr.ircompensation = :ohmicdrop;
	ely_odr.Ru = R_u 
	
end

# ╔═╡ f4f59329-e817-495a-9e83-1ab53e7738a9
function sweep(model, grid, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
	#elystruct.elydata.ircompfactor = 1
	
    celldata = deepcopy(model)
    pnpcell = PNPSystem(grid; bcondition = model.bcondition, celldata = model.elydata, reaction = model.reaction)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
		kwargs...
    )

end

# ╔═╡ 590a17bd-c98d-4b23-9139-04b550c95efe
pnpcell = PNPSystem(grid; bcondition = elystruct.bcondition, celldata = elystruct.elydata, reaction = elystruct.reaction)

# ╔═╡ fd933939-647a-46a4-9322-53cebb07b8d9
unc = sweep(merge(elystruct, (elydata = ely,     )), grid, sawtooth; nperiods)

# ╔═╡ d70e5a81-9b9e-43ce-8809-a5aaf00b1bd9
ely_pts

# ╔═╡ d00584d6-f876-4eb5-860d-c32bd21d51f6
pts = sweep(merge(elystruct, (elydata = ely_pts, )), grid, sawtooth; nperiods)

# ╔═╡ e535b3d1-90d4-4b4c-979e-5e51a8f35e7d
ely_odr

# ╔═╡ ab5edc9a-c650-49c9-a71b-f15287b67bd3
odr = sweep(merge(elystruct, (elydata = ely_odr, )), grid, sawtooth; nperiods)

# ╔═╡ f77ec140-e091-46c1-8bc6-1c00f3880e83
AuCO2RR_plots.plot_conc_time_electrode(unc, elystruct)

# ╔═╡ 144c4dec-4cdf-440c-94da-38f6fb25742d


# ╔═╡ bd9b5c58-0375-4ce2-aa55-c74b921aa050
let
    result = unc
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 1260))       
        ax1, leg1 = AuCO2RR_plots.panel_conc_time!(f, f[1, 2], result, elystruct;     xlabel = "")
        ax2 = AuCO2RR_plots.panel_time_current!(f, f[2, 2], result, elystruct;      xlabel = "")
        ax3 = AuCO2RR_plots.panel_time_voltage!(f, f[3, 2], result;             xlabel = "")
        ax4 = AuCO2RR_plots.panel_time_ph!(f, f[4, 2], result;                  xlabel = "")   # (d) pH
        #ax5, hm = panel_co2_log_contour!(f, f[5, 2], f[5, 3], result, X, bulk; L_val = L)  # (e)
        ax5, leg5 = AuCO2RR_plots.QoverK(f, f[5, 2], result,elystruct; legend_pos = f[5, 3])       # (f) Q/K

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
        #hidexdecorations!(ax5, grid = false)     # contour도 x라벨 숨김 (ax6만 x축 표시)

        ax3.yticks = LinearTicks(3)
        ax4.yticks = LinearTicks(4)
        ylims!(ax1, 1e-14, 1e2)

        linkxaxes!(ax1, ax2, ax3, ax4, ax5)
        rowgap!(f.layout, 15)
        for r in 1:5
			if r == 1
				rowsize!(f.layout, r, 350)
			else
            	rowsize!(f.layout, r, Relative(1/6))   # 6등분
			end
        end
        f
    end
    fig
end

# ╔═╡ 57fa330b-0ab4-4910-8e8e-5a3a305d298e
begin
	@show maximum(abs.(currents(odr, 7) .- currents(unc, 7)))
	@show maximum(abs.(currents(pts, 7) .- currents(unc, 7)))
	@show maximum(abs.(odr.voltages .- unc.voltages))
	
	@show R_u
	@show maximum(abs.(currents(unc, 7)))         
	
end

# ╔═╡ e0373a0d-2577-4c18-a6a7-8db361cf506a
@show elystruct.elydata.Ru elystruct.elydata.ircompfactor elystruct.elydata.x_ref

# ╔═╡ ca1c75c0-e005-40a3-acf3-0c96675b22ef
function panel_time_current_comp!(fig, panel_pos, result_1, result_2, m;
                             redox_species = nothing,
                             include_capacitive = true,
                             scale = cm^2/mA,
                             sign = -1,
                             color_gradient = false, 
                             lw = 4,
                             xlabel = lab_time,
                             label_1 = "odr", 
                             label_2 = "unc")  
    model = m.elydata

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
        if include_capacitive
            ely = model
            if ely.ircompensation == :ohmicdrop
                icc = ely.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in res.tsol[1:end-1]]
            end
        end
        return sign .* (I_F .+ I_C) .* scale
    end

    I_1 = calc_total_current(result_1)
    I_2 = calc_total_current(result_2)

    color_1 = color_gradient ? :skyblue : :black
    color_2 = color_gradient ? :orange : :red      

    line1 = lines!(ax, result_1.times, I_1; color = color_1, linewidth = lw, label = label_1)
    line2 = lines!(ax, result_2.times, I_2; color = color_2, linewidth = lw, linestyle = :dash, label = label_2) 

    axislegend(ax; position = :rt) 

    return ax
end

# ╔═╡ 59ce784f-15b8-4880-a0a5-c5d78f5ce60a
function panel_time_current_diff!(fig, panel_pos, result_1, result_2, m;
                             redox_species = nothing,
                             include_capacitive = true,
                             scale = cm^2/mA,
                             sign = -1,
                             color_gradient = false, 
                             lw = 4,
                             xlabel = lab_time,
                             label_1 = "odr", 
                             label_2 = "unc")  
    model = m.elydata

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
        if include_capacitive
            ely = model
            if ely.ircompensation == :ohmicdrop
                icc = ely.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in res.tsol[1:end-1]]
            end
        end
        return sign .* (I_F .+ I_C) .* scale
    end

    I_1 = calc_total_current(result_1)
    I_2 = calc_total_current(result_2)

	I_diff = I_1 .- I_2
	
    color_1 = color_gradient ? :skyblue : :black
    color_2 = color_gradient ? :orange : :red       

    #line1 = lines!(ax, result_1.times, I_1; color = color_1, linewidth = lw, label = label_1)
    #line2 = lines!(ax, result_2.times, I_2; color = color_2, linewidth = lw, linestyle = :dash, label = label_2) 
	lines!(ax, result_1.times, I_diff; color = :red, linewidth = lw)

    #axislegend(ax; position = :rt) 

    return ax
end

# ╔═╡ 2879a1f3-b1e7-47f0-99f6-ccb44b62572d
let
    result = unc
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 800))      
		ax1 = panel_time_current_comp!(f, f[1,1], odr, unc, elystruct; redox_species = 5, xlabel = "")
		ax2 = panel_time_current_diff!(f, f[2,1], odr, unc, elystruct; redox_species = 5, xlabel = "time")
        hidexdecorations!(ax1, grid = false)

		f
	end
	fig
end

# ╔═╡ Cell order:
# ╠═8d20515c-54c6-11f1-aeac-bdc0c7e51b18
# ╠═cd3522a0-904f-4b89-b454-fcd62ed8de52
# ╟─5471cc8c-f4e9-41bb-a0ec-a8c5125a2b4d
# ╟─ea90b4a5-6709-4a85-87de-b71fe87e71f7
# ╟─4b663702-b895-4a84-b065-0673e287b0ca
# ╠═0854d2c7-1b61-46b0-aa0a-496cdc8439f2
# ╟─b24b7697-2c64-4e7a-8f7f-87cab6add489
# ╠═96af1c36-0fab-4c1d-bad1-96b98a927d66
# ╠═d011050d-c66e-4e8a-a173-f55679403a90
# ╠═035cf151-b62a-42ee-8b03-38e68fc4e4b3
# ╠═0963720a-e310-45b2-92a5-a9e5bc3e6888
# ╠═e7e23a3a-2383-42e0-b5cd-0b315ad7d3e9
# ╠═f4f59329-e817-495a-9e83-1ab53e7738a9
# ╠═590a17bd-c98d-4b23-9139-04b550c95efe
# ╠═fd933939-647a-46a4-9322-53cebb07b8d9
# ╠═d70e5a81-9b9e-43ce-8809-a5aaf00b1bd9
# ╠═d00584d6-f876-4eb5-860d-c32bd21d51f6
# ╠═e535b3d1-90d4-4b4c-979e-5e51a8f35e7d
# ╠═ab5edc9a-c650-49c9-a71b-f15287b67bd3
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═144c4dec-4cdf-440c-94da-38f6fb25742d
# ╠═bd9b5c58-0375-4ce2-aa55-c74b921aa050
# ╠═57fa330b-0ab4-4910-8e8e-5a3a305d298e
# ╠═e0373a0d-2577-4c18-a6a7-8db361cf506a
# ╠═2879a1f3-b1e7-47f0-99f6-ccb44b62572d
# ╠═ca1c75c0-e005-40a3-acf3-0c96675b22ef
# ╠═59ce784f-15b8-4880-a0a5-c5d78f5ce60a
