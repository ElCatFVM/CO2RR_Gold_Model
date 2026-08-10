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
    L = 2500 * μm
    hmin = 1.0e-6 	* μm
    #hmax = 0.05*L
    hmax = 0.1*L
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
        factor=0.85
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

# ╔═╡ 57b5d6b5-cf65-4997-903d-f8b08e983b8b
#Δt_min = Δt_cv, Δt_max = Δt_cv, abstol = 1.0e-12, reltol = 1.0e-12

# ╔═╡ e2eb4160-550a-4a07-b4c8-66863aaab46f
cv_odr = sweep(elystruct_odr, grid, sawtooth; nperiods)

# ╔═╡ 7bf09e84-8ad2-4b3b-92bd-8627d2c0dedc
Δt_cv = 0.05

# ╔═╡ 47515ef3-b6aa-49d0-b4fc-b471cb947aa1
cv_unc = sweep(elystruct_unc, grid, sawtooth; nperiods)

# ╔═╡ ba48cb67-a7af-4d65-b65f-b55d9c29f54c
cv_irc = sweep(elystruct_irc, grid, sawtooth; nperiods)

# ╔═╡ f77ec140-e091-46c1-8bc6-1c00f3880e83
AuCO2RR_plots.plot_conc_time_electrode(cv_odr, elystruct_odr)

# ╔═╡ bd9b5c58-0375-4ce2-aa55-c74b921aa050
plot_cv_summary(cv_odr, elystruct_odr)	

# ╔═╡ 3ad8349a-f95f-4ba8-98c4-b404f98137aa
plot#_7species_contours(cv_odr, X, elystruct_odr)

# ╔═╡ 36238c29-eea8-47b4-b5a3-17905b3483a8
extrema(currents(cv_odr, 7) ./ currents(cv_odr, 5))   # CO / CO₂

# ╔═╡ 96eb220e-d88c-4f3c-872e-174938b19a8c
let
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 400))
        panel_time_current_diff!(
            f, f[1, 1],
            cv_odr,            # result_1  (label_1 = "odr")
            cv_unc,            # result_2  (label_2 = "unc")
            elystruct_odr;     # m  → m.elydata (cspecies, capacitive-term mode)
            include_capacitive = false,   # ★ recommend false when the modes differ 
            title = "IR Compensation Effect"
        )
        f
    end
    fig
end

# ╔═╡ 5412784d-68b9-4ac1-bc87-0051e1da0ce8
tsol = (; Δt_min = 0.05, Δt_max = 0.05,   # 적응 제어 무력화 → 같은 격자
          abstol  = 1.0e-12,
          reltol  = 1.0e-12)


# ╔═╡ 59ec47b9-2a47-4350-b708-ca7a980bbdb2
md"""
## CV Solution
"""

# ╔═╡ 01f688a7-3265-4bc5-9b74-6eedfc36476d
begin
	grid_dict = Dict{Float64, Any}()
	for Lv in [100, 500, 1000, 2500] .* μm
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

# ╔═╡ de959fd3-0a8f-4f18-924c-6c9785408a04
cv_unc_fixed = sweep(elystruct_unc, grid, sawtooth; nperiods, tsol...)

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
    factors = [0.0, 0.1, 0.5, 0.9]

    ely0    = elystruct_odr.elydata

    #same = (; Δt_min = 0.05, Δt_max = 0.05, abstol = 1.0e-12, reltol = 1.0e-12)

    splitruns = Dict{Float64, Any}()
    for f in factors
        cd  = copy(ely0; ircompensation = copy(ely0.ircompensation; factor = f))
        pnp = PNPSystem(grid; bcondition = elystruct_odr.bcondition, celldata = cd,
                        reaction = elystruct_odr.reaction, unknown_storage = :dense)
        splitruns[f] = LiquidElectrolytes.cvsweep(pnp; voltages = sawtooth, nperiods,
                                                  store_solutions = true, same...)
    end
    splitruns
end

# ╔═╡ 0fb39632-67a5-4780-a3d4-71c51bd9a3c3
begin
    cd  = copy(ely0; ircompensation = copy(ely0.ircompensation; factor = 0.0))
    pnp = PNPSystem(grid; bcondition = elystruct_odr.bcondition, celldata = cd,
                    reaction = elystruct_odr.reaction, unknown_storage = :dense)
    cv_f0 = LiquidElectrolytes.cvsweep(pnp; voltages = sawtooth, nperiods,
                                       store_solutions = true, tsol...)
end

# ╔═╡ 56baa4d6-6838-40d3-9649-c2b151fd9755
begin
	I0 = AuCO2RR_plots.cv_current(cv_f0)       .* (cm^2/mA)
	Iu = AuCO2RR_plots.cv_current(cv_unc_fixed) .* (cm^2/mA)
	length(I0) == length(Iu), maximum(abs, I0 .- Iu) / maximum(abs, Iu)
end

# ╔═╡ bb8cd305-1bb7-4d4e-8158-014a290e2325
plot_ircomp_compare(splitruns, elystruct_odr;
                                  reference = cv_unc_fixed,
                                  reference_model = elystruct_unc,
                                  panels = (:driving, :metal_time, :ircomp, :cv),
                                  layout = (2, 2))

# ╔═╡ 68527236-79a2-4064-ace0-ddd046401c13


# ╔═╡ 743e985e-1fa3-426f-b2c0-cdb52517443d
md"""
### scanrate
"""

# ╔═╡ 535e8412-7e45-4be2-9533-b1ee1245d3e1
begin
    scanrates = [0.005, 0.05, 0.5, 5.0]
    SR_vec = [sweep(elystruct_odr, grid,
                    SawTooth(scanrate = sr, vmin = -1.2, vmax = 1.2,
                             scanup = false, vstart = 0.0; tstart = 0.0);
                    nperiods,
                    Δu_opt = 0.01,
                    abstol = 1.0e-12, reltol = 1.0e-12) for sr in scanrates]
end


# ╔═╡ d7b50e64-d0c9-4551-a732-a7a1dea445a4
plot_scanrate_sweeps(SR_vec, scanrates)

# ╔═╡ 4c520a35-1109-472d-a800-d42133a7ef94
plot_scanrate_sweeps_split(SR_vec, scanrates; ured = (-1.3, -0.6),
        uox = (-0.2, 1.45),)

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
        voltages = sawtooth, nperiods, store_solutions = true,
        Δt_min = Δt_cv, Δt_max = Δt_cv,
        abstol = 1.0e-12, reltol = 1.0e-12)

    P_recs = pressure_varied_sweep(elystruct_odr.elydata, cvfun;
                                   Pvec = [0.1, 0.5, 1.0], ispec = 5)
end

# ╔═╡ bde0975f-957f-4e8c-8c15-59e7f1a09400
pressure_varied_cvsweep(P_recs)

# ╔═╡ ad3c4119-cd97-4933-820f-84eee96cbc07
pressure_varied_cvsweep_split(P_recs)

# ╔═╡ 46f1c6f3-d36e-442b-9936-3dd4dec397b3
rec = first(r for (p, r) in P_recs if p == 0.1)

# ╔═╡ 9530079d-d4e1-4712-9c8d-e01ec26a1d79
plot_cv_summary(rec, elystruct_odr)

# ╔═╡ 7620826e-c562-41e2-aa2b-e717843f9187
plot_cv_scanrate_grid(SR_vec, m; scanrates = scanrates)

# ╔═╡ 9849e70d-47f9-4b04-aafb-54033c674b54
plot_exp_sim_cvsweep_split(P_recs; annotate_x = 0.93)

# ╔═╡ 9a5e3eed-0f4e-43b6-b7b5-e28c62a30a4b
begin
    ispec = 5                                    # CO₂
    ely_01 = at_pressure(elystruct_odr.elydata, 0.1; ispec)

    SR_01 = [ LiquidElectrolytes.cvsweep(
                  PNPSystem(grid;
                            bcondition = elystruct_odr.bcondition,
                            celldata = deepcopy(ely_01),
                            reaction = elystruct_odr.reaction,
                            unknown_storage = :dense);
                  voltages = SawTooth(scanrate = sr, vmin = -1.2, vmax = 0.8,
                                      scanup = false, vstart = 0.0; tstart = 0.0),
                  nperiods, store_solutions = true,
                  Δu_opt = 0.002, abstol = 1.0e-12, reltol = 1.0e-12)
              for sr in [0.01, 0.05, 0.2] ]
end

# ╔═╡ f91ff7ae-334b-4263-858a-3fa0892911e8
plot_cv_scanrate_grid(SR_01, elystruct_odr; scanrates = [0.01, 0.05, 0.2])


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
# ╠═57b5d6b5-cf65-4997-903d-f8b08e983b8b
# ╠═e2eb4160-550a-4a07-b4c8-66863aaab46f
# ╠═7bf09e84-8ad2-4b3b-92bd-8627d2c0dedc
# ╠═47515ef3-b6aa-49d0-b4fc-b471cb947aa1
# ╠═ba48cb67-a7af-4d65-b65f-b55d9c29f54c
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═bd9b5c58-0375-4ce2-aa55-c74b921aa050
# ╠═3ad8349a-f95f-4ba8-98c4-b404f98137aa
# ╠═36238c29-eea8-47b4-b5a3-17905b3483a8
# ╠═96eb220e-d88c-4f3c-872e-174938b19a8c
# ╠═5412784d-68b9-4ac1-bc87-0051e1da0ce8
# ╟─59ec47b9-2a47-4350-b708-ca7a980bbdb2
# ╠═01f688a7-3265-4bc5-9b74-6eedfc36476d
# ╠═cf4713e6-706c-484f-b398-8ba6cc33561b
# ╟─5cc3a435-8fae-4eb9-8e2f-2de72e2e807a
# ╠═affdc880-3a70-429a-bcfb-a53cf7ed1f31
# ╠═2f128f25-d629-44b5-bcb7-d3b9d59362ea
# ╠═de959fd3-0a8f-4f18-924c-6c9785408a04
# ╠═0fb39632-67a5-4780-a3d4-71c51bd9a3c3
# ╠═56baa4d6-6838-40d3-9649-c2b151fd9755
# ╠═b5abb13d-b6e4-4a47-a4b0-d099588b1b60
# ╠═c2e9573b-2105-4930-9f17-5a3bf3fe7058
# ╠═1ea52fc4-0921-43ea-88e4-04fadee047f9
# ╠═6301f323-16d8-4805-bbcf-4a055bca2d59
# ╠═6f56453e-384b-4be6-9298-8101d12d8b79
# ╠═69bcb95c-69a5-4daf-9646-2313c8452018
# ╠═bb8cd305-1bb7-4d4e-8158-014a290e2325
# ╠═68527236-79a2-4064-ace0-ddd046401c13
# ╠═743e985e-1fa3-426f-b2c0-cdb52517443d
# ╠═535e8412-7e45-4be2-9533-b1ee1245d3e1
# ╠═d7b50e64-d0c9-4551-a732-a7a1dea445a4
# ╠═4c520a35-1109-472d-a800-d42133a7ef94
# ╠═d291fbb9-0cc4-4f6b-a7ed-4e3f103f4d11
# ╠═2333ae9c-0f59-4a87-90ed-bc5e8fba9f76
# ╠═bde0975f-957f-4e8c-8c15-59e7f1a09400
# ╠═ad3c4119-cd97-4933-820f-84eee96cbc07
# ╠═46f1c6f3-d36e-442b-9936-3dd4dec397b3
# ╠═9530079d-d4e1-4712-9c8d-e01ec26a1d79
# ╠═7620826e-c562-41e2-aa2b-e717843f9187
# ╠═9849e70d-47f9-4b04-aafb-54033c674b54
# ╠═9a5e3eed-0f4e-43b6-b7b5-e28c62a30a4b
# ╠═f91ff7ae-334b-4263-858a-3fa0892911e8
