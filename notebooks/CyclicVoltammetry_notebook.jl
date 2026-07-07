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
elystruct = GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = :ohmicdrop)

# ╔═╡ f4f59329-e817-495a-9e83-1ab53e7738a9
function sweep(model, grid, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
	elystruct.elydata.ircompfactor = 0.00
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

# ╔═╡ e2eb4160-550a-4a07-b4c8-66863aaab46f
odr_cv_pnp = sweep(elystruct, grid, sawtooth; nperiods = nperiods)

# ╔═╡ f77ec140-e091-46c1-8bc6-1c00f3880e83
AuCO2RR_plots.plot_conc_time_electrode(odr_cv_pnp, elystruct)

# ╔═╡ 144c4dec-4cdf-440c-94da-38f6fb25742d


# ╔═╡ bd9b5c58-0375-4ce2-aa55-c74b921aa050
let
    result = odr_cv_pnp
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

# ╔═╡ Cell order:
# ╠═8d20515c-54c6-11f1-aeac-bdc0c7e51b18
# ╟─ea90b4a5-6709-4a85-87de-b71fe87e71f7
# ╟─4b663702-b895-4a84-b065-0673e287b0ca
# ╠═0854d2c7-1b61-46b0-aa0a-496cdc8439f2
# ╟─b24b7697-2c64-4e7a-8f7f-87cab6add489
# ╠═96af1c36-0fab-4c1d-bad1-96b98a927d66
# ╠═d011050d-c66e-4e8a-a173-f55679403a90
# ╠═035cf151-b62a-42ee-8b03-38e68fc4e4b3
# ╠═0963720a-e310-45b2-92a5-a9e5bc3e6888
# ╠═f4f59329-e817-495a-9e83-1ab53e7738a9
# ╠═590a17bd-c98d-4b23-9139-04b550c95efe
# ╠═e2eb4160-550a-4a07-b4c8-66863aaab46f
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═144c4dec-4cdf-440c-94da-38f6fb25742d
# ╠═bd9b5c58-0375-4ce2-aa55-c74b921aa050
