### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 8d20515c-54c6-11f1-aeac-bdc0c7e51b18
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	using Revise
    using LiquidElectrolytes
	using AuCO2RR, AuCO2RR.AuCO2RR_plots
	using VoronoiFVM
	using LessUnitful
	using ExtendableGrids, GridVisualize
	using DelimitedFiles
	using Interpolations
	using PlutoUI, HypertextLiteral
	using Latexify
	using Printf
	using Test
	using LinearAlgebra
	using Colors
	using FileIO
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
	        scanrate = 0.05,
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
    L = 1000 * μm
    hmin = 1.0e-6 	* μm
    hmax = 0.05*L
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ bcfc1e78-c1e9-4538-890a-35e74bfc4060
elystruct_unc = GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = NoIRCompensation())

# ╔═╡ 530ed232-5e1a-4c0a-85c2-6bea12824c71
σ=conductivity(elystruct_unc.elydata, elystruct_unc.elydata.c_bulk)

# ╔═╡ 2fe95a9f-202e-418a-b422-cef1ef013c9d
md"""
### Set IR compensation mode
Choose between
- `:none`: No IR compensation
- `:pseudopotentiostat`: "measure" voltage at point `x_ref` and and add this to applied voltage at 0
- `:ohmicdrop`: Estimate voltage to add to applied voltage R_u from current of species `ircompspecies` and uncompensated resistance `Ru`.

"""

# ╔═╡ 8626e977-a77f-4646-be99-33e713dbf593
begin
    ircomp=:ohmicdrop
    
if  ircomp==:none
        ircompensation=NoIRCompensation()
elseif ircomp==:pseudopotentiostat
        ircompensation=PseudoPotentiostat()
elseif ircomp==:ohmicdrop
    ircompensation=OhmicDropEstimation(
                          Ru = L / σ,
                          species=5,
                          ne=2,
        factor=0.0
    )
else
    error("wrong value of ircompensation: $ircomp")
end
end

# ╔═╡ 0963720a-e310-45b2-92a5-a9e5bc3e6888
begin
    elystruct_odr = GoldModel.create_model(; use_md_hydrated = false, γ_select = "Stefan", ircompensation)
    elystruct_odr
end

# ╔═╡ a6561efa-29d5-498d-a66f-b8fffe1bc293
begin
    elystruct_irc = GoldModel.create_model(; use_md_hydrated = false, γ_select = "Stefan", ircompensation = PseudoPotentiostat())
    elystruct_irc
end

# ╔═╡ f4f59329-e817-495a-9e83-1ab53e7738a9
function sweep(model, grid, sawtooth; nperiods = 1,kwargs...)
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
# ╠═╡ show_logs = false
cv_odr = sweep(elystruct_odr, grid, sawtooth; nperiods = nperiods)

# ╔═╡ 47515ef3-b6aa-49d0-b4fc-b471cb947aa1
# ╠═╡ show_logs = false
cv_unc = sweep(elystruct_unc, grid, sawtooth; nperiods = nperiods)

# ╔═╡ ba48cb67-a7af-4d65-b65f-b55d9c29f54c
# ╠═╡ show_logs = false
cv_irc = sweep(elystruct_irc, grid, sawtooth; nperiods = nperiods)

# ╔═╡ f77ec140-e091-46c1-8bc6-1c00f3880e83
AuCO2RR_plots.plot_conc_time_electrode(cv_odr, elystruct_odr)

# ╔═╡ bd9b5c58-0375-4ce2-aa55-c74b921aa050
let
    result = cv_irc
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 1260))       
        ax1, leg1 = AuCO2RR_plots.panel_conc_time!(f, f[1, 2], result, elystruct_odr;     xlabel = "")
        ax2 = AuCO2RR_plots.panel_time_current!(f, f[2, 2], result, elystruct_odr;      xlabel = "")
        ax3 = AuCO2RR_plots.panel_time_voltage!(f, f[3, 2], result; sawtooth = sawtooth, xlabel = "")
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

# ╔═╡ 3ad8349a-f95f-4ba8-98c4-b404f98137aa
plot_7species_contours(cv_odr, X, elystruct_unc)

# ╔═╡ 34ab3828-53ef-41f9-bcd7-95755b0a8c8b
plot_potential_split(cv_unc, elystruct_unc; scanrate = 0.05)

# ╔═╡ 36238c29-eea8-47b4-b5a3-17905b3483a8
extrema(currents(cv_odr, 7) ./ currents(cv_odr, 5))   # CO / CO₂

# ╔═╡ 59ec47b9-2a47-4350-b708-ca7a980bbdb2
md"""
## CV Solution
"""

# ╔═╡ 01f688a7-3265-4bc5-9b74-6eedfc36476d
begin
	grid_dict = Dict{Float64, Any}()
	for Lv in [100, 300, 500, 1000, 3000, 5000] .* μm
	    Xg = ExtendableGrids.geomspace(0, Lv, 1.0e-7*μm, Lv*0.1)
	    grid_dict[Lv] = ExtendableGrids.simplexgrid(Xg)
	end
	grid_dict
end

# ╔═╡ cf4713e6-706c-484f-b398-8ba6cc33561b
begin
	results = cvsweep_odr_over_L(
	    elystruct_irc.elydata,       
	    grid_dict,
	    elystruct_irc.bcondition,
	    elystruct_irc.reaction,
	    sawtooth;
		unknown_storage=:dense,
	    nperiods,
		Δu_opt=0.025,
	#	Δt_min=1.0e-8,
	#	damp_initial=0.5
	)
end

# ╔═╡ 5cc3a435-8fae-4eb9-8e2f-2de72e2e807a
md"""
## CV PLot
"""

# ╔═╡ affdc880-3a70-429a-bcfb-a53cf7ed1f31
let
	fig=plot_cv_current_variedL(results, elystruct_odr)
	CairoMakie.save("cv-$(ircomp).png",fig)
	fig
end

# ╔═╡ 2f128f25-d629-44b5-bcb7-d3b9d59362ea
Lmax=sort(keys(grid_dict))[end]

# ╔═╡ 0fb39632-67a5-4780-a3d4-71c51bd9a3c3
function fixed(fig)
    @htl("""
    <div style="
        display: inline-block;
        align-self: flex-start;
        flex: 0 0 auto;
    ">
        $(fig)
    </div>
    """)
end

# ╔═╡ b5abb13d-b6e4-4a47-a4b0-d099588b1b60
plots=[AuCO2RR_plots.plottsol(grid_dict[L],
					   elystruct_odr.elydata, 
					   results[L].tsol;
							  figscale=L/Lmax,
							  stride=3, # xscale=log10 is slow
					   species=5, levels=15)|>fixed for L in sort(keys(grid_dict))];

# ╔═╡ c2e9573b-2105-4930-9f17-5a3bf3fe7058
PlutoUI.ExperimentalLayout.vbox(plots)

# ╔═╡ 1ea52fc4-0921-43ea-88e4-04fadee047f9
[
	L=>blthickness(grid_dict[L], elystruct_odr.elydata, results[L].tsol; species=5, atol=5.0e-2)/μm
	for L in sort(keys(grid_dict))]

# ╔═╡ 6301f323-16d8-4805-bbcf-4a055bca2d59
blthickness(grid, elystruct_unc.elydata, cv_unc.tsol; species = 5) / μm

# ╔═╡ 6f56453e-384b-4be6-9298-8101d12d8b79
blthickness(grid, elystruct_odr.elydata, cv_odr.tsol; species = 5) / μm

# ╔═╡ 69bcb95c-69a5-4daf-9646-2313c8452018
begin
    factors = [0.0, 0.3, 0.5, 0.7, 0.9]
    ely0 = elystruct_odr.elydata
    
    splitruns = Dict{Any, Any}()
    for f in factors
        # copy() so `redoxreaction` survives — create_model bound we_breactions into it,
        # and a freshly constructed OhmicDropEstimation would not have it
        ircompensation=OhmicDropEstimation(
                          Ru = L / σ,
                          species=5,
                          ne=2,
        factor=f
        )
        
        cd = copy(ely0; ircompensation)
        pnp = PNPSystem(grid; bcondition = elystruct_odr.bcondition, celldata = cd,
                        reaction = elystruct_odr.reaction, unknown_storage = :dense)
        splitruns[f] = LiquidElectrolytes.cvsweep(pnp; voltages = sawtooth, nperiods,
                                                  store_solutions = true, Δu_opt = 0.025)
    end
    splitruns["unc"] = cv_unc          # 무보정 기준선
    splitruns
end


# ╔═╡ f32f2d90-a774-4485-a850-b15a7d1268da
let
    fig = Figure(size = (700, 520))
    ax = Axis(fig[1, 1]; xlabel = AuCO2RR_plots.lab_voltage,
              ylabel = rich(rich("ϕ", font = :italic), "(0)  (V)"),
              title = "potential divider vs. compensation")
    ablines!(ax, 0, 1; color = (:black, 0.45), linestyle = :dot, linewidth = 2)

    ks = [0.0, 0.3, 0.5, 0.7, 0.9]
    for (i, f) in enumerate(ks)
        r = splitruns[f]
        lines!(ax, r.sawtooth, r.voltages;
               color = AuCO2RR_plots.CMAP_SCANRATE[(i - 1) / (length(ks) - 1)],
               linewidth = 4, label = "f = $f")
    end
    axislegend(ax; position = :lt, framevisible = false)
    fig
end


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
    lines!(ax, t1, I_diff; color = colorant"#C4844C", linewidth = lw)

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

# ╔═╡ 743e985e-1fa3-426f-b2c0-cdb52517443d
md"""
### scanrate
"""

# ╔═╡ 535e8412-7e45-4be2-9533-b1ee1245d3e1
begin
    scanrates = [0.005, 0.01, 0.02, 0.05, 0.1]
    SR_vec = [sweep(elystruct_odr, grid,
                    SawTooth(scanrate = sr, vmin = -1.2, vmax = 0.8,
                             scanup = false, vstart = 0.0; tstart = 0.0);
                    nperiods, Δu_opt = 0.025) for sr in scanrates]
end

# ╔═╡ d7b50e64-d0c9-4551-a732-a7a1dea445a4
plot_scanrate_sweeps(SR_vec, scanrates)

# ╔═╡ 4c520a35-1109-472d-a800-d42133a7ef94
plot_scanrate_sweeps_split(SR_vec, scanrates)

# ╔═╡ d291fbb9-0cc4-4f6b-a7ed-4e3f103f4d11
md"""
### Pressure
"""

# ╔═╡ 2333ae9c-0f59-4a87-90ed-bc5e8fba9f76
begin
    cvfun = ely -> LiquidElectrolytes.cvsweep(
        PNPSystem(grid; bcondition = elystruct_odr.bcondition,
                  celldata = ely, reaction = elystruct_odr.reaction,
                  unknown_storage = :dense);
        voltages = sawtooth, nperiods, store_solutions = true, Δu_opt = 0.025)

    P_recs = pressure_varied_sweep(elystruct_odr.elydata, cvfun;
                                   Pvec = [0.2, 0.4, 0.6, 0.8, 1.0], ispec = 5)
end

# ╔═╡ bde0975f-957f-4e8c-8c15-59e7f1a09400
pressure_varied_cvsweep(P_recs)

# ╔═╡ ad3c4119-cd97-4933-820f-84eee96cbc07
pressure_varied_cvsweep_split(P_recs)

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
# ╠═530ed232-5e1a-4c0a-85c2-6bea12824c71
# ╟─2fe95a9f-202e-418a-b422-cef1ef013c9d
# ╠═8626e977-a77f-4646-be99-33e713dbf593
# ╠═0963720a-e310-45b2-92a5-a9e5bc3e6888
# ╠═a6561efa-29d5-498d-a66f-b8fffe1bc293
# ╠═f4f59329-e817-495a-9e83-1ab53e7738a9
# ╠═590a17bd-c98d-4b23-9139-04b550c95efe
# ╠═1e59ac64-17d4-42a0-befc-32961332c4c7
# ╠═e2eb4160-550a-4a07-b4c8-66863aaab46f
# ╠═47515ef3-b6aa-49d0-b4fc-b471cb947aa1
# ╠═ba48cb67-a7af-4d65-b65f-b55d9c29f54c
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═bd9b5c58-0375-4ce2-aa55-c74b921aa050
# ╠═3ad8349a-f95f-4ba8-98c4-b404f98137aa
# ╠═34ab3828-53ef-41f9-bcd7-95755b0a8c8b
# ╠═36238c29-eea8-47b4-b5a3-17905b3483a8
# ╠═96eb220e-d88c-4f3c-872e-174938b19a8c
# ╟─59ec47b9-2a47-4350-b708-ca7a980bbdb2
# ╠═01f688a7-3265-4bc5-9b74-6eedfc36476d
# ╠═cf4713e6-706c-484f-b398-8ba6cc33561b
# ╟─5cc3a435-8fae-4eb9-8e2f-2de72e2e807a
# ╠═affdc880-3a70-429a-bcfb-a53cf7ed1f31
# ╠═2f128f25-d629-44b5-bcb7-d3b9d59362ea
# ╠═0fb39632-67a5-4780-a3d4-71c51bd9a3c3
# ╠═b5abb13d-b6e4-4a47-a4b0-d099588b1b60
# ╠═c2e9573b-2105-4930-9f17-5a3bf3fe7058
# ╠═1ea52fc4-0921-43ea-88e4-04fadee047f9
# ╠═6301f323-16d8-4805-bbcf-4a055bca2d59
# ╠═6f56453e-384b-4be6-9298-8101d12d8b79
# ╠═69bcb95c-69a5-4daf-9646-2313c8452018
# ╠═f32f2d90-a774-4485-a850-b15a7d1268da
# ╠═9bd12311-fbe8-436b-8c14-18edb5744929
# ╠═743e985e-1fa3-426f-b2c0-cdb52517443d
# ╠═535e8412-7e45-4be2-9533-b1ee1245d3e1
# ╠═d7b50e64-d0c9-4551-a732-a7a1dea445a4
# ╠═4c520a35-1109-472d-a800-d42133a7ef94
# ╠═d291fbb9-0cc4-4f6b-a7ed-4e3f103f4d11
# ╠═2333ae9c-0f59-4a87-90ed-bc5e8fba9f76
# ╠═bde0975f-957f-4e8c-8c15-59e7f1a09400
# ╠═ad3c4119-cd97-4933-820f-84eee96cbc07
