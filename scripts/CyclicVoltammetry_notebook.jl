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
    #L = user_input_model.L * μm
    hmin_blinf = 1.0e-6 	* μm
    hmax_blinf = L_blinf * 0.02 
    X_blinf = ExtendableGrids.geomspace(0, L_blinf, hmin_blinf, hmax_blinf)
    grid_blinf = ExtendableGrids.simplexgrid(X_blinf)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ 0963720a-e310-45b2-92a5-a9e5bc3e6888
elystruct = GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = :none)

# ╔═╡ f4f59329-e817-495a-9e83-1ab53e7738a9
function sweep(model, grid, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
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

# ╔═╡ 144c4dec-4cdf-440c-94da-38f6fb25742d
function plot_conc_time_electrode(result, elystruct;
                                  nspecies=7,
                                  scale=(mol/dm^3))

    names  = elystruct.bulknames
    colors = elystruct.bulkcolors

    times = result.tsol.t
    nt    = length(times)

    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    iOH = findfirst(isequal("OH⁻"), names)
    iH  = (iOH === nothing) ? 1 : iOH # fallback
    I   = currents(result, iH) .* (cm^2/mA)

    cols = RGBf(0.3, 0.5, 1.0) 

    fig, ax_conc, ax_current = with_theme(AuCO2RR_plots.electrochemistry_theme()) do
        f = Figure(size = (960, 540))
        
        a_conc = Axis(f[1, 1],
            xlabel = L"\text{time / s}",
            ylabel = L"\mathbf{c_{i}^\ddagger}\;(\mathrm{M})",
            yscale = log10,
            limits = ((times[1] - (times[end] / 200), times[end] + (times[end] / 100)), (1e-12, 1e4)),
            rightspinevisible = false
        )
        
        yt_vals = 10.0 .^ (4:-4:-12)
        yt_lbls = [L"10^{4}", L"10^{0}", L"10^{-4}", L"10^{-8}", L"10^{-12}"]
        a_conc.yticks = (yt_vals, yt_lbls)

        # 2. Right axis: current density (linear scale)
        a_current = Axis(f[1, 1],
            ylabel = L"\text{Current density}\;(\mathrm{mA/cm^{2}})",
            yaxisposition = :right,
            ygridvisible = false,
            leftspinevisible = false,
            rightspinecolor = cols,
            ylabelcolor = cols,
            yticklabelcolor = cols,
            ytickcolor = cols,
            yticks = LinearTicks(4) 
        )

        return f, a_conc, a_current
    end

    # sync x-axis between both axes
    linkxaxes!(ax_conc, ax_current)

    # 1. Concentration profile (left axis)
    for i in 1:nspecies
        y = conc[i, :]
        y_fixed = max.(y, eps(Float64)) # floor for log-scale safety
        lines!(ax_conc, times, y_fixed; color=colors[i])
    end

    # 2. Current density profile (right axis)
    lines!(ax_current, result.times, I ./ 2; color=cols, linewidth = 5.5)

    return fig
end

# ╔═╡ f77ec140-e091-46c1-8bc6-1c00f3880e83
plot_conc_time_electrode(odr_cv_pnp, elystruct)

# ╔═╡ 98223a9d-e7a4-4e9a-b781-61486c1a61ad
elystruct.species_dict

# ╔═╡ 5077d283-8067-4bee-9e94-14cb3f759e9d
begin
	powlab(n) = rich("10", superscript(string(n)))
	lab_time  = rich(rich("t", font=:italic), "  (s)")

	function panel_conc_time!(fig, panel_pos, result, elystruct; nspecies=7, scale=mol/dm^3, lw=4,
	                          xlabel = lab_time, legend_pos = nothing)
	    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
	                 "#C0392B","#27AE60","#2980B9"]
	    colors = sp_colors[1:nspecies]
    	names  = elystruct.bulknames
	    times  = result.tsol.t
	    nt     = length(times)
	    conc   = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

	    ax_c = Axis(panel_pos;
			        xlabel = xlabel,
			        ylabel = rich(rich("c", font=:italic),
			                      subscript(rich("α", font=:italic)),
			                      superscript("‡"), "  (M)"),
			        yscale = log10,
			        yminorticksvisible = true,
			        yminorticks = IntervalsBetween(9)
				   )

	    ax_c.yticks = (10.0 .^ (4:-4:-12),
	                   [powlab(4), powlab(0), powlab(-4), powlab(-8), powlab(-12)])

	    for i in 1:nspecies
	        lines!(ax_c, times, max.(conc[i, :], eps(Float64));
	               color = colors[i], linewidth = lw)
	    end

	    legend_elements = [ [LineElement(color = colors[i], linewidth = lw)] for i in 1:nspecies ]
	    legend_labels   = [ rich(string(names[i]), color = colors[i]) for i in 1:nspecies ]

	    leg = Legend(legend_pos === nothing ? fig[1, 3] : legend_pos,
	                 legend_elements, legend_labels;
	                 framevisible = false,
	                 labelsize    = 20)

	    return ax_c, leg
	end

	function panel_time_current!(fig, panel_pos, result, elystruct;
                             redox_species = nothing,
                             include_capacitive = true,
                             scale = cm^2/mA,
                             sign = -1,
                             color_gradient = true,
                             lw = 4,
                             xlabel = lab_time)

    ax = Axis(panel_pos; xlabel = xlabel,
              ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"))

    redox = Dict(elystruct.species_dict["CO₂"] => 2)

    n_t = length(result.times)
    I_F = zeros(n_t)
    for (idx, n_e) in redox
        I_F .+= n_e .* currents(result, idx)
    end

    I_C = zeros(n_t)
    if include_capacitive
        ely = elystruct
        if ely.elydata.ircompensation == :ohmicdrop
            icc = ely.icc
            node_we = 1
            I_C = [u[icc, node_we] for u in result.tsol[1:end-1]]
        end
    end

    I = sign .* (I_F .+ I_C) .* scale

    cols = color_gradient ? :skyblue : :black
    lines!(ax, result.times, I; color = cols, linewidth = lw)

    return ax
end

	function panel_time_voltage!(fig, panel_pos, result; lw=5, xlabel = lab_time)
	    ax = Axis(panel_pos; xlabel = xlabel,
	              ylabel = rich(rich("U", font=:italic), "  (V vs. SHE)"))
	    lines!(ax, result.times, result.voltages;
	           color = parse(Colorant, "#D7C2F0"), linewidth = lw)
	    return ax
	end

	function panel_co2_log_contour!(fig, panel_pos, cbar_pos, result, X, bulk;
	                                 scale=mol/dm^3, num_levels=24, L_val=nothing)
	    co2_idx = findfirst(s -> s.name == "CO₂", bulk)
	    times = result.tsol.t

	    log_c_bulk = log10(result.tsol[co2_idx, end, 1] / scale)

	    M = [log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / scale, 1e-12))
	         for ix in 1:length(X), it in 1:length(times)]

	    c_min = 0.0
	    c_max = log10(result.tsol[co2_idx, end, 1] / scale) - log10(1e-5)

	    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

	    ax = Axis(panel_pos;
	        xlabel = lab_time,
	        ylabel = rich(rich("x", font=:italic), "  (m)"),
	        yscale = log10,
	        yminorticksvisible = true,
	        yminorticks = IntervalsBetween(9),
	        yticks = (10.0 .^ (-12:3:-6),
	          [powlab(-12), powlab(-9), powlab(-6)])
				 )

	    hm = heatmap!(ax, times, X .+ 1e-12, M';
	        colorrange = (c_min, c_max),
	        colormap = discrete_cmap,
	        interpolate = false)

	    Colorbar(cbar_pos, hm;
	        label = rich("log", subscript("10"), "(", rich("c", font=:italic),
	                     subscript("bulk"), ") − log", subscript("10"), "(",
	                     rich("c", font=:italic),
	                     subscript(rich("CO", subscript("2"))), ")"),
	        ticklabelsize = 20, labelsize = 20)
	    return ax, hm
	end

	function panel_time_ph!(fig, panel_pos, result; ihplus=2, scale=mol/dm^3, lw=5,
	                        xlabel = lab_time, color = parse(Colorant, "#7BB661"))
	    times = result.tsol.t
	    nt    = length(times)
	    cH    = [result.tsol[ihplus, 1, t] / scale for t in 1:nt]   # 표면(node 1) H⁺ 농도 (M)
	    pH    = -log10.(max.(cH, eps(Float64)))
	    ax = Axis(panel_pos; xlabel = xlabel, ylabel = rich("pH"))
	    lines!(ax, times, pH; color = color, linewidth = lw)
	    return ax
	end
	function QoverK(fig, panel_pos, result; scale=mol/dm^3, lw=4,
	                xlabel = lab_time, legend_pos = nothing)
		
	    times = result.tsol.t
	    nt    = length(times)

		ihco3 		= elystruct.species_dict["HCO₃⁻"]
		iohminus 	= elystruct.species_dict["OH⁻"]
		ico3 		= elystruct.species_dict["CO₃²⁻"]
		ico2 		= elystruct.species_dict["CO₂"]

		
	    # Equilibrium Constant
	    EqK_hco3 = (elystruct.elydata.c_bulk[ihco3] * elystruct.elydata.c_bulk[iohminus]) /
	           		elystruct.elydata.c_bulk[ico3]
		EqK_co3 = (elystruct.elydata.c_bulk[ico2] * elystruct.elydata.c_bulk[iohminus]) /
	           elystruct.elydata.c_bulk[ihco3]
		
	    # reaction quotient Q(t)
	    cco3  = [result.tsol[ico3,     1, t] for t in 1:nt]
	    cohm  = [result.tsol[iohminus, 1, t] for t in 1:nt]
	    chco3 = [result.tsol[ihco3,    1, t] for t in 1:nt]
		cco2  = [result.tsol[ico2,     1, t] for t in 1:nt]
		
	    Qt_hco3    = (chco3 .* cohm) ./ cco3
		Qt_co3     = (cco2 .* cohm) ./ chco3
		
	    ax = Axis(panel_pos;
	        xlabel = xlabel,
	        ylabel = rich(rich("Q", font=:italic), " / ", rich("K", font=:italic)),
	        yscale = log10,
	        yticks = (10.0 .^ (-6:3:6),
		          [powlab(-6), powlab(-3), powlab(0), powlab(3), powlab(6)]),
						        yminorticksvisible = true,
				        yminorticks = IntervalsBetween(27)	 
				 )
	    hlines!(ax, [1e0]; color = :black, linestyle = :dash, linewidth = 2)
		
	    l1 = lines!(ax, times, max.(Qt_hco3 ./ EqK_hco3, eps(Float64));
	                color = "#2980B9", linewidth = lw)
	    l2 = lines!(ax, times, max.(Qt_co3  ./ EqK_co3,  eps(Float64));
	                color = "#E67E22", linewidth = lw)
		
		
		ylims!(ax, low = 1e-7, high = 1e7)
		
		leg = Legend(legend_pos,
	                 [l1, l2],
	                 [rich("HCO", subscript("3"), superscript("-"), " ⇌ CO", subscript("3"), superscript("2-")),
	                  rich("CO", subscript("2"), " ⇌ HCO", subscript("3"), superscript("-"))],
	                 framevisible = false)
	
	    return ax, leg
	end

	
end

# ╔═╡ bd9b5c58-0375-4ce2-aa55-c74b921aa050
let
    result = odr_cv_pnp
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 1260))       
        ax1, leg1 = panel_conc_time!(f, f[1, 2], result, elystruct;     xlabel = "")
        ax2 = panel_time_current!(f, f[2, 2], result, elystruct;      xlabel = "")
        ax3 = panel_time_voltage!(f, f[3, 2], result;             xlabel = "")
        ax4 = panel_time_ph!(f, f[4, 2], result;                  xlabel = "")   # (d) pH
        #ax5, hm = panel_co2_log_contour!(f, f[5, 2], f[5, 3], result, X, bulk; L_val = L)  # (e)
        ax5, leg5 = QoverK(f, f[5, 2], result; legend_pos = f[5, 3])       # (f) Q/K

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
# ╠═f4f59329-e817-495a-9e83-1ab53e7738a9
# ╠═590a17bd-c98d-4b23-9139-04b550c95efe
# ╠═e2eb4160-550a-4a07-b4c8-66863aaab46f
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═144c4dec-4cdf-440c-94da-38f6fb25742d
# ╠═bd9b5c58-0375-4ce2-aa55-c74b921aa050
# ╠═98223a9d-e7a4-4e9a-b781-61486c1a61ad
# ╠═5077d283-8067-4bee-9e94-14cb3f759e9d
