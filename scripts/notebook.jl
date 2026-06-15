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
	const L = 2500 * μm
end

# ╔═╡ d011050d-c66e-4e8a-a173-f55679403a90
begin
	sawtooth = SawTooth(
	        scanrate = 0.05,
	        vmin     = -1.2, 
	        vmax     = 1.2,
	        scanup   = false,
			vstart = 0.0; tstart = 0.0
	    )
	const nperiods = 2
end

# ╔═╡ 035cf151-b62a-42ee-8b03-38e68fc4e4b3
begin
    #Vmax = 2 * V
    #L = user_input_model.L * μm
    hmin = 1.0e-6 	* μm
    hmax = L * 0.02 
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ 0963720a-e310-45b2-92a5-a9e5bc3e6888
elystruct = GoldModel.create_model(;use_md_hydrated = false, γ_select = "Stefan", ircompensation = :none)

# ╔═╡ f4f59329-e817-495a-9e83-1ab53e7738a9
function sweep(model, grid, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
    celldata = deepcopy(model)
    pnpcell = PNPSystem(grid; bcondition = model.bcondition, celldata = model.elydata, reaction = model.reaction)
    return result = LiquidElectrolytes.cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
		kwargs...
    )

end

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

# ╔═╡ 5077d283-8067-4bee-9e94-14cb3f759e9d
# I'll reconstruct after the publishing...

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
# ╟─035cf151-b62a-42ee-8b03-38e68fc4e4b3
# ╠═0963720a-e310-45b2-92a5-a9e5bc3e6888
# ╠═f4f59329-e817-495a-9e83-1ab53e7738a9
# ╠═e2eb4160-550a-4a07-b4c8-66863aaab46f
# ╠═f77ec140-e091-46c1-8bc6-1c00f3880e83
# ╠═144c4dec-4cdf-440c-94da-38f6fb25742d
# ╠═5077d283-8067-4bee-9e94-14cb3f759e9d
