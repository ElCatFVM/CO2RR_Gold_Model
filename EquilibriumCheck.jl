### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 91ac9e35-71eb-4570-bef7-f63c67ce3881
begin
    using Pkg
    Pkg.activate(@__DIR__)
	using Revise
    using LiquidElectrolytes
	using CatmapInterface
	using Catalyst: unknowns
	using VoronoiFVM
	using LessUnitful
	using ExtendableGrids, GridVisualize
	using DelimitedFiles
	using PlutoUI, HypertextLiteral
	using PreallocationTools
	using Latexify
	using Catalyst
	using Printf
	using Test
	using Colors
	using FileIO
	using CSV, DataFrames
	if isdefined(Main,:PlutoRunner)
        using CairoMakie	
   		default_plotter!(CairoMakie)
 		CairoMakie.activate!(type="svg")
    end
end;

# ╔═╡ 3ac837b8-559b-41c2-8f83-1331839dcf7e
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
end

# ╔═╡ a94bc4e1-506f-4e40-bfe8-1ce7e6093974
pkgdir(LiquidElectrolytes)

# ╔═╡ beae1479-1c0f-4a55-86e1-ad2b50174c83
md"""
## Setup
"""

# ╔═╡ ab2184fc-0279-46d9-9ee4-88fe3e732789
md"""
### Units
"""

# ╔═╡ 7316901c-d85d-48e9-87dc-3614ab3d81a5
@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm;

# ╔═╡ 6b7cfe87-8190-40a5-8d25-e39ef8d55db5
md"""
### Data
"""

# ╔═╡ 5a146a44-03dc-45f3-ae15-993d11c2edac
begin
	@phconstants N_A c_0 k_B e h
	const F = N_A * e

	const voltages = (-1.5:0.1:-0.0) * V
	
	# geometrical constants
	const Γ_we 		= 1
	const Γ_bulk 	= 2
	
	# bulk constants
	const pH 		= 6.8
	const T 		= 298.0 * K
	const Hcp_CO  	= 9.7e-6 * mol/(m^3 * Pa)
    const Hcp_CO2 	= 3.3e-4 * mol/(m^3 * Pa)
	# species not involved in reactions
	const ikplus 	= 1
	# species involved in buffer reactions but not in surface reactions
	const ibufferstart = 2
	const ihplus 	= 2
	const ihco3 	= 3
    const ico3 		= 4
	# species involved in buffer reactions and surface reactions
    const isurfacestart = 5
	const ico2 		= 5
    const iohminus 	= 6
	const ibufferend = 6
	# species involved in surface reactions but not in buffer reactions
	const ico  		= 7
	const nc 		= 7
	## reaction rate constants for bulk reactions
	### CO2 + OH- <=> HCO3-
	const kbe1 = 4.44e7 / (mol/dm^3)
	const kbf1 = 5.93e3 / (mol/dm^3) / s
	const kbr1 = kbf1 / kbe1
	### HCO3- + OH- <=> CO3-- + H2O
	const kbe2 = 4.66e3 / (mol/dm^3)
	const kbf2 = 1.0e8 / (mol/dm^3) / s
	const kbr2 = kbf2 / kbe2
	### CO2 + H20 <=> HCO3- + H+
    const kae1 = 4.44e-7 * (mol/dm^3)
    const kaf1 = 3.7e-2 / s
    const kar1 = kaf1 / kae1
    ### HCO3- <=> CO3-- + H+ 
    const kae2 = 4.66e-5 / (mol/dm^3)
    const kaf2 = 59.44e3 / (mol/dm^3) / s
    const kar2 = kaf2 / kae2
	### autoprotolyse
    const kwe  = 1.0e-14 * (mol/dm^3)^2
    const kwf  = 2.4e-5 * (mol/dm^3) / s
    const kwr  = kwf / kwe
	const aH₂O = 1.0 #* mol/dm^3

	# surface constants
	const S 		= 9.61e-5 / N_A * (1.0e10)^2 * mol/m^2
	const C_gap 	= 20 * μF/cm^2
    const ϕ_pzc 	= 0.16 * V
	const ico_t 	= 8
	const icooh_t 	= 9
	const ico2_t 	= 10
	const isurfaceend = 10
	const na 		= 3 # CO_t, CO2_t, COOH_t
	const M0 		= 18.0153 * ufac"g/mol"
	const v0        = N_A * (8.2 * Å)^3 # 1 / (55.4 * ufac"M") #
	
	const species_dict = Dict(
		"K⁺" => ikplus,
		"HCO₃⁻" => ihco3,
		"CO₃²⁻" => ico3,
		"CO₂" => ico2,
		"OH⁻" => iohminus, 
		"H⁺" => ihplus,
		"CO" => ico,
		"CO_t" => ico_t,
		"COOH_t" => icooh_t,
		"CO2_t" => ico2_t,
	)

	const species_dict_catmap = Dict(
		"OH_g" => iohminus,
		"CO2_aq" => ico2,
		"CO_aq" => ico,
		"CO_t" => ico_t,
		"COOH_t" => icooh_t,
		"CO2_t" => ico2_t,
	)
end;

# ╔═╡ 00947475-c96e-4ecc-a1ef-5be5e3e3c864
begin
	@kwdef struct BulkSpecies
		name::String 
		z::Int
		D::Float64
		c_bulk::Union{Nothing, Float64}
		κ::Float64
		a::Float64
		v::Float64
		M::Float64
		color::Symbol
	end
	function BulkSpecies(;name, z, c_bulk=nothing, D, κ=0.0, a=0.0, color)
		D *= m^2/s
		c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol/dm^3
		a *= Å
		v = N_A * a^3
		M = M0 * v
		BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
	end
	function make_eneutral(bulk_species::Vector{BulkSpecies}
						  ;name, z, D, κ=0.0, a=0.0, v=N_A*(a*Å)^3, M=M0*v, color)
		a *= Å
		c_bulk = -mapreduce(x -> x.c_bulk * x.z, +, bulk_species)/z
		BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
	end
end;

# ╔═╡ 4b64e168-5fe9-4202-9657-0d4afc237ddc
md"""
### Reaction Description
"""

# ╔═╡ de2c826d-6c05-47cf-b5f5-44a00ea9889c
md"""
#### Buffer System
"""

# ╔═╡ d8f00649-e2ed-4bdd-853f-05268f0d5353
md"""
Consider the bicarbonate buffer system in base and acid as well as autoprotolysis of water:
"""

# ╔═╡ 47b36c81-b57e-4dd0-a22f-999e4fd3ac9f
begin
	@variables t
	@species HCO₃⁻(t), CO₃²⁻(t), CO₂(t), OH⁻(t),H⁺(t)
	@parameters γHCO₃⁻ γCO₃²⁻ γCO₂ γOH⁻ γH⁺
	buffer_rn = @reaction_network buffer begin
		($kbf1 * γCO₂ * γOH⁻, $kbr1 * γHCO₃⁻), CO₂ + OH⁻ <--> HCO₃⁻
		($kbf2 * γHCO₃⁻ * γOH⁻, $kbr2 * $aH₂O * γCO₃²⁻), HCO₃⁻ + OH⁻ <--> CO₃²⁻
		($kaf1 * $aH₂O * γCO₂, $kar1 * γHCO₃⁻ * γH⁺), CO₂ <--> HCO₃⁻ + H⁺
		($kaf2 * γHCO₃⁻, $kar2 * γCO₃²⁻ * γH⁺), HCO₃⁻ <--> CO₃²⁻ + H⁺
		($kwf * $aH₂O, $kwr * γH⁺ * γOH⁻), ∅ <--> H⁺ + OH⁻
	end
	odesys_buffer = convert(ODESystem, buffer_rn; combinatoric_ratelaws=false)
	const f_buffer! = CatmapInterface.generate_function(
		buffer_rn; 
		dvs = [H⁺, HCO₃⁻, CO₃²⁻, CO₂, OH⁻], 
		ps = [γH⁺, γHCO₃⁻, γCO₃²⁻, γCO₂, γOH⁻]
	)
	latexify(buffer_rn; env=:chemical)
end

# ╔═╡ 1e877f17-0219-45f1-b640-3a25ae085dbd
latexify(odesys_buffer) # hide

# ╔═╡ 8912f990-6b02-467a-bd11-92f94818b1c7
md"""
#### Surface Reactions
"""

# ╔═╡ a8157cc1-1761-4b11-a37c-9e12a9ca695e
md"""
A microkinetic modeling approach is taken:

The reaction mechanism for the $CO_2$ reduction is divided into four elementary reactions at the electrode surface:

1. Adsorption of $CO_2$ molecules at the oxygen atoms
${CO_2}_{(aq)} + * \rightleftharpoons {CO_2 *}_{(ad)}$

2. First proton-coupled electron transfer
${CO_2*}_{(ad)} + H_2O_{(l)} + e^- \rightleftharpoons COOH*_{(ad)} + OH^-_{(aq)}$

3. Second proton-coupled electron transfer
$COOH*_{(ad)} + e^- \rightleftharpoons CO*_{(ad)} + OH^{-}_{(aq)}$
with the transition state: $*CO-OH^{TS}$

4. Desorption of $CO$
$*CO_{(ad)} \rightleftharpoons CO_{(aq)} + *$
"""

# ╔═╡ 6b5cf93c-0df3-4a18-8786-502361736838
begin
	catmap_params 		= CatmapInterface.parse_catmap_input("catmap_CO2R_data/catmap_CO2R_template.mkm")
	rn 					= create_reaction_network(catmap_params)
	odesys 				= convert(ODESystem, rn; combinatoric_ratelaws=false)
	odesys 				= CatmapInterface.liquidize(odesys, catmap_params)
	vars 				= Catalyst.unknowns(odesys)
	const f_microkinetics! 	= CatmapInterface.generate_function(
		odesys;
		dvs = sort(vars, by=x->species_dict_catmap[string(operation(x))])
	)
	const paramsidx = paramsmap(odesys)
	latexify(odesys)
end

# ╔═╡ d2c0642d-dfa5-4a76-bd36-ac4a735a3299
md"""
##### Reaction Rates
"""

# ╔═╡ 06d45088-ab8b-4e5d-931d-b58701bf8464
md"""
For the calculation of the reaction rates a __mean field approach__ based is applied.

The acitivities of H⁺, OH⁻ and e⁻ are set to be zero. The dependence on the rates on the pH-value, applied voltage and surface charges are contained in the reaction rate constants.

The surface charging relation $σ = σ(U)$ is given by the Robin boundary condition

$σ(U) = C_{gap} (ϕ_{we} - ϕ_{pzc} - ϕ^\ddagger)$

where the gap capacitance between the working electrode and the reaction plane ($\ddagger$) is given by $C_{gap} = 20~μF/cm^2$. The potential of zero current is measured to be $ϕ_{pzc} = 0.16~V$.

__Question is the pH-dependence only in the reaction rate constants (i.e. activity of OH⁻ must be set to 0)?__
"""

# ╔═╡ d0093605-0e35-4888-a93c-8456c698e6f0
md"""
##### Reaction Rate Constants
"""

# ╔═╡ f0b5d356-6b97-4878-98de-bee5f380d41a
md"""
The details of the rate constant calculations can be found in the documentation of CatMAP and the CatMAP specification for this model: "./models/Au/catmap_CO2R_template.mkm". Here only a few import points are highlighted.

According to the transition state theory the kinetic model assumes that the dynamics of a reaction from an educt complex to a product complex can be projected onto a path parametrized by a reaction coordinate. On the reaction path an activated complex can be identified. The rate constant is interpreted according to the Eyring equation
```math
k = κ \frac{k_B T}{h} exp\left(-ΔG^\ddagger\right)
```
where the $κ$ is a transmission coefficient and $ΔG^\ddagger$ is the change in the Gibbs free energy between the educt and the transition state for canonical ensembles in the respective states. For the transition state only degrees of freedom perpendicular to the reaction coordinate are considered and treated as harmonic degrees of freedom.

Only the second concerted proton-electron transfer (CPET) is assumed to pass through a distinct activated complex. For the other reactions no kinetic barrier is assumed and $ΔG^\ddagger$ is interpreted as the change in the relative Gibbs free energy of formation between the educt and product state $Δᵣμ^{⦵}$ in case that the reaction does not occur spontaneously $Δᵣμ^{⦵} > 0$.
The (entropy-free) pre-exponential factors $κ k_B T/ h$ are approximated as $10^{13}$ for all but the adsorption of CO where it set to $10^8$.
"""

# ╔═╡ e3eda42f-e2f3-4c10-81c4-610246ca528d
md"""
The Gibbs free energies are based on potential energies calculated by applying DFT-methods. In order to obtain thermodynamic values ab-initio statistical models are used.
"""

# ╔═╡ 7b87aa2a-dbaf-441c-9ad7-444abf15f664
md"""
The free energies $ΔG_f$ of the surface species are corrected according to the (excess) surface charge density $σ$ by a fitted quadratic model:

$ΔG_f(σ) = a_σ~σ + b_σ~σ^2$
"""

# ╔═╡ 2b9d9bfd-d660-4b4d-8f0b-b6b5bcc0dbfa
begin
    Vmax = 2 * V

    L = 80 * μm

    hmin = 1.0e-6 	* μm

    hmax = 1.0 		* μm 

    X = ExtendableGrids.geomspace(0, L, hmin, hmax)

    grid = ExtendableGrids.simplexgrid(X)
end;

# ╔═╡ 6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
md"""
### Electrolyte Data
"""

# ╔═╡ 848b7aeb-968f-4116-8038-b61276f02b6c
elydata_NaClO₄ = ElectrolyteData(
 		z = [-1, 1],
		κ = [8.0, 8.0],
		#v = [25.0, 25.0],
		c_bulk = [0.5, 0.5],
		ε = 26.0,
		#vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
)

# ╔═╡ f18dc873-1c9d-46d3-9596-92d28705e894
elydata_NaF = ElectrolyteData(
 		z = [-1, 1],
		κ = [25.0, 25.0],
		#v = [25.0, 25.0],
		c_bulk = [0.5, 0.5],
		ε = 26.0,
		#vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
		#vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
		v0 = 18.048 * ufac"cm^3" / ufac"mol",		
)

# ╔═╡ 53f12821-7d8d-4971-87fd-ad4689ec62a5
md"""
### γ(activity coefficients) function
"""

# ╔═╡ e1e0ca0f-7f88-40f0-850e-590b25da0331
md"""
    DGML_gamma!(γ, c, p, electrolyte)

Activity coefficients according to Dreyer, Guhlke, Müller, Landstorfer.

```math
γ_i = \exp\left(\frac{\tilde v_i p}{RT}\right) \left(\frac{\bar c}{c_0}\right)^{m_i} \frac{1}{v_0\bar c}
```

**Input**:
- c: vector of concentrations
- p: pressure
- electrolyte: instance of `ElectrolyteData`

**Output**: 
- γ (mutated) activity coefficients
"""

# ╔═╡ 0db74a70-af86-492c-affb-9de62ffe4455
function DGML_γ!(γ, c, p, electrolyte)
	
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
    for ic in cspecies
        γ[ic] = rexp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc))
    end
    return γ
end

# ╔═╡ 4f388fe0-6bc8-4a29-bccc-fa725e62c6a7
md"""
    Stefan_gamma!(γ, c, p, electrolyte)

Activity coefficients according to Nat Commun 11, 33 __(2020)__.

```math
γ_i = \frac{1}{1-\Sigma_1^i c_i v_i}
```

**Input**:
- c_i: concentration of species i 
- v_i: volume of species i
- electrolyte: instance of `ElectrolyteData`

**Output**: 
- γ (mutated) activity coefficients
"""

# ╔═╡ 5d179c52-43d7-4bcb-a2df-93c5806876fa
function Stefan_γ!(γ, c, p, electrolyte)
	
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp, c_bulk, v, nc) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
    for ic in cspecies
        γ[ic] = 1.0 / (1 - sum(c_bulk[i] * v[i] for i in 1:nc) / (mol/dm^3))
    end
    return γ
end

# ╔═╡ 5f17b4f7-54d6-4ad0-9886-252854840a80
function activity_coefficient!(
    γ::AbstractVector,
	u,
    data,
    mode::Function;
    idxs::AbstractVector{<:Integer} = 1:data.nc
)
    (; v, ip, pscale, p_bulk, M, M0, v0, κ, RT, nc, Mrel, tildev, v0, cspecies, rexp) = data

    if mode == Stefan_γ!
        # Ringe et al. approach (volume fraction-based)
        total_volfrac = sum(data.c_bulk[i] * v[i] for i in idxs)
        γ[idxs] .= 1.0 / (1 - total_volfrac / (mol / dm^3))

    elseif mode == DGML_γ!
        # Dreyer et al. approach
        p = u[ip] * pscale - p_bulk
        c0, barc = c0_barc(u, data)

        for ic in idxs
            barv = v[ic] + κ[ic] * v0
            γ[ic] = rexp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc))

        end
    else
        error("Unsupported mode: must be Stefan_γ! or DGML_γ!")
    end

    return γ
end


# ╔═╡ 2a20d9be-6c1e-4c1f-8bb6-a7693800732d
md"""
## Double Layer Capacitance
"""

# ╔═╡ 4f7ec19d-cd60-4c2b-a766-7557caa471c0
md""" 
### System Setup
"""

# ╔═╡ 952a26ce-2610-48cc-9158-eda816da3a1c
molarities = [0.005, 0.01, 0.05, 0.1, 0.5] 

# ╔═╡ 9a4e01d9-f469-4427-bf4c-883adb67ae24
function pb_bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data
	
    ## Dirichlet ϕ=ϕ_we at Γ_we
    #boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
    #boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_bulk, value = data.ϕ_bulk)
    #boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)

	## Robin ϕ=dϕ₀/dx
	boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))


    return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ d76d8413-c019-4728-b182-7f7cb78dede4
md"""
### Result Plots
"""

# ╔═╡ 783e2058-c720-4f31-8e51-7c313813924c
md"""
##### κ(Solvation Number) Plots
"""

# ╔═╡ 9598e2c6-521e-4f8d-82d8-a836809736f3
md"""
### Plotting Functions
"""

# ╔═╡ f2043f2c-f3c8-4b0f-944c-7b55624dac08
function capsplot_v(vis, result_named)
    color = [:magenta, :blue]
    for (i, (name, res)) in enumerate(result_named)
		voltages = res[1].voltage_range[1:201]
        caps = vec(res[1].dlcaps)[1:201]
        scalarplot!(
            vis, voltages, caps / (μF / cm^2);
            limits = (-1, 100), xlimits = (-1.1, 1.1),
            color = color[i], clear = false, label = name,
            markershape = :none, yscale = 10,
            xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)"
        )
    end
    return vis
end;

# ╔═╡ 8bfdf2f5-c80a-4ce0-a8e1-b315affffb5f
begin
	
	#computed capacitance plot(Fig 13 & Fig 14)
	Landstorfer_NaF_5mM = CSV.read("Landstorfer_data/Landstorfer_NaF_0.005M.csv", DataFrame);
	Landstorfer_NaF_100mM = CSV.read("Landstorfer_data/Landstorfer_NaF_0.1M.csv", DataFrame);
	Landstorfer_NaClO₄_100mM = CSV.read("Landstorfer_data/Landstorfer_NaClO4_0.1M.csv", DataFrame);
	Landstorfer_NaClO₄_5mM = CSV.read("Landstorfer_data/Landstorfer_NaClO4_0.005M.csv", DataFrame);
	
	#Solvation number plot(Fig 7)
	Landstorfer_Low_κ = CSV.read("Landstorfer_data/Landstorfer_kappa0.csv", DataFrame);
	Landstorfer_High_κ = CSV.read("Landstorfer_data/Landstorfer_kappa40.csv", DataFrame);
end;

# ╔═╡ 5c808c71-6094-49d7-8215-e88262f34e1f
md"""
## Cyclic Voltammetry
"""

# ╔═╡ da8390d1-47e8-451f-b12b-45b8aca7b6ec
md"""
### System Setup
"""

# ╔═╡ 39c8ef0d-aac2-4c7f-8004-4166c460ebc5
# ╠═╡ disabled = true
#=╠═╡
nnpresult = sweep(model; eneutral = true, tunnel = false)
  ╠═╡ =#

# ╔═╡ eb920b6e-86a6-4dd6-8e66-6b7e27d81257
md"""
### Result Plots
"""

# ╔═╡ 3f30fae0-18d4-4e5f-9618-cbd9852d7857
function zstr(z::Int)
    if z == -1
        return "-"
    elseif z < -1
        return "$(-z)-"
    elseif z == 0
        return ""
    elseif z == 1
        return "+"
    elseif z > 1
        return "$(z)+"
    end
end


# ╔═╡ bb00b5bb-326e-47f9-a4f4-e7b4f29dd1f2
md"""
### Plotting Functions
"""

# ╔═╡ 82baec54-048a-4851-8808-dd8c31eae4b9
function plot_concentration_profile(result, X, bulk; filename = "concentration_profile.gif")
    species = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    XX = X[2:end] / nm
    ru_all = result.tsol[:, :, :] ./ (mol / dm^3)
    voltages = result.voltages

    vis = GridVisualizer(;
        size = (650, 400),
        clear = true,
        legend = :rt,
        limits = (-25, 2),
        xlimits = (XX[1], maximum(XX)),
        xlabel = "x / nm",
        ylabel = "log c / (mol/dm³)",
        xscale = :log,
    )

    movie(vis; file=filename, framerate=10) do vis
        for it in 1:length(voltages)
            title = "V = $(round(voltages[it], digits=2)) V"
            for i in 1: nc
                yvals = log10.([c > 0 ? c : NaN for c in ru_all[i, 2:end, it]])
                scalarplot!(
                    vis,
                    XX,
                    yvals,
                    color = colors[i],
                    label = species[i],
                    clear = (i == 1),
                    title = title,
                )
            end
            reveal(vis)
        end
    end

    return isdefined(Main, :PlutoRunner) ? LocalResource(filename) : nothing
end

# ╔═╡ b7cb5183-65e8-4ee8-af86-2bedd11daecc
function CVPlot!(result, model)
    ic = model.cspecies
    fig = Figure(size = (650, 400))
    ax = Axis(fig[1, 1], 
              limits = ((-1.0, 1.5),(-0.50, 7)),
			  ylabel = "Current Density(mA/cm²)",
			  xlabel = "Voltage (ϕ-ϕₚ)"
			 )
	
    total_current = zero(currents(result, ic[1]))
    #for s in ic
    #    total_current .+= (currents(result, s) * mA / cm^2)
    #end
	total_current = currents(result, ico)

	
    lines!(ax, result.voltages, total_current,
           color = RGBf.(range(0, 1, length(result.voltages)), 0.0, 0.0))

    fig
end

# ╔═╡ 842b074b-f808-48d8-8dc5-110ddd907f90
md"""
## IV Results
"""

# ╔═╡ 31298257-d35a-4f6f-8a76-ff00d5361ced
md"""
### System Setup
"""

# ╔═╡ 72269ec4-a56e-46d9-85c8-0dd8ccaf43e1
solver_control = (; max_round 	= 4,
					maxiters 	= 20,
              		tol_round 	= 1.0e-9,
              		verbose 	= "a",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)

# ╔═╡ 7a02463d-cfd9-4648-af53-f1e65d46733f
md"""
### Result Plots
"""

# ╔═╡ 114d2324-5289-4e44-8d77-736a9bdec365
md"""
Show only pH: $(@bind useonly_pH PlutoUI.CheckBox(default=false))
"""

# ╔═╡ c1d2305e-fb8b-4845-a414-08fff84aa9b0
md"""
### Plotting Functions
"""

# ╔═╡ d5ab1a28-3a60-49d9-bb3e-ca589b1c79fd
begin
	curr(J, ix) = [F * abs(j[ix]) for j in J]
	
	function plotcurr(result; df = nothing)
	    scale = 1 / (mol / dm^3)
	    volts = result.voltages[result.voltages .< -0.4]
	    vis = GridVisualizer(;
	                         size = (600, 400),
	                         tilte = "IV Curve",
	                         xlabel = "Φ_WE/(V vs. SHE)",
	                         ylabel = "I/(mA/cm²)",
	                         legend = :lb,
							 yscale = :log,
		)
							 
	    scalarplot!(vis,
	                volts,
	                curr(result.j_we, iohminus)[result.voltages .< -0.4] .* cm^2/mA;
	                color = :green,
	                clear = false,
	                linestyle = :solid,
	                label = "e⁻, we")
		if !isnothing(df)
			scalarplot!(vis,
						df[:voltage],
						df[:current],
						clear = false,
						linewidth = 0,
						markershape = :cross,
						markersize = 8,
						markevery = 1,
						color = :red,
						label = "Ringe et. al")
		end
		
	    reveal(vis)
	end
end

# ╔═╡ 686ac3dc-c191-4575-ba0c-d4c2551474b5
md"""
### Regression Test
"""

# ╔═╡ d1ab199f-1a40-4377-bca3-7f72f3cde3a9
md"""
Run regression tests $(@bind runregtest PlutoUI.CheckBox())
"""

# ╔═╡ de144adb-a467-4077-8cb1-d86462f56110
html"""<hr>"""

# ╔═╡ d0985ca6-fef5-4b67-9ad6-f51d84b595b4
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 8ae53b8a-0fb3-4c1c-8e5f-a3782a85141c
begin
    hrule() = html"""<hr>"""
    function highlight(mdstring, color)
        htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""
    end

    macro important_str(s)
        :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        :(highlight(Markdown.parse($s), "#ccffcc"))
    end

    html"""
        <style>
	     pluto-log-dot-sizer  { max-width: 655px;}
         pluto-log-dot.Stdout { background: #002000;
	                            color: #10f080;
                                border: 6px solid #b7b7b7;
                                min-width: 18em;
                                max-height: 300px;
                                width: 675px;
                                    overflow: auto;
 	                           }
	
    </style>
"""
end

# ╔═╡ 9d7d4d68-c9cc-4a42-a99d-ae25a1ab554c
html"""<hr>"""

# ╔═╡ e5fc814f-a8e1-41ef-b81a-c3b0839a2f87
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	left: 1rem;
            	top: $(top)px;
            	width: 350px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;

# ╔═╡ ae50877a-20da-4b1a-acd5-cc0a73e428da
floataside(
    @bind user_input confirm(
        PlutoUI.combine() do Child
	md"""
	##### __User Data__  
	``CO_2 + H_2O + 2e^- \leftrightharpoons CO + 2OH^-``
			
	- ``V_\mathrm{min}``: $(Child("vmin", NumberField(-2.0:0.1:0.0; default = -1.5)))  ``V_\mathrm{max}``: $(Child("vmax", NumberField(0.0:0.1:2.0; default = 1.5)))
	- ``z_R``: $(Child("zR", NumberField(-2:2; default = -1)))   
	  ``n``: $(Child("n", NumberField(0:2; default = 2)))
	
	- Scan rate ``(V/s)``: $(Child("scanrate", TextField(6; default = "0.3")))  
	- Periods: $(Child("nperiods", NumberField(1:10; default = 3)))
	
	- ``L``: 80 μm  
	- `Double64`: $(Child("double64", CheckBox()))  
	- `tunnel`: $(Child("tunnel", CheckBox()))
	
	---
	
	##### __Parameter Set__
	
	- __Other ions__  
	  ``a``: $(Child("at", NumberField(0.0:0.1:20.0; default = 8.2)))  
	  ``κ``: $(Child("κt", NumberField(0.0:0.1:20.0; default = 8.0)))
	
	- __Cation__  
	  ``a``: $(Child("ak", NumberField(0.0:0.1:20.0; default = 8.2)))  
	  ``κ``: $(Child("κk", NumberField(0.0:0.1:20.0; default = 8.0)))
	
	---
	
	##### __Model Selection__  
	- Model: $(Child("model_choice", Select(["Gold_Model", "Landstorfer_NaClO₄ model", "Landstorfer_NaF model"])))
	
	---
	
	##### __Activity Coefficient__  
	- Mode: $(Child("mode", Select(["Stefan_γ", "DGML_γ"])))
	"""
	        end;
	        label = "Submit"
	    );
	    top = 50
	)


# ╔═╡ ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
begin
	at = user_input[:at]
	κt = user_input[:κt]
	ak = user_input[:ak]
	κk = user_input[:κk]
	const bulk = let 
		bulk = [
				BulkSpecies(;name = "HCO₃⁻", 
							z = -1, 
							D = 1.185e-9, 
							c_bulk = 0.091, 
							a = at, 
							κ = κt, 
							color = :brown
				),
				BulkSpecies(;name = "CO₃²⁻",
							z = -2, 
							D = 0.923e-9, 
							c_bulk = 2.68e-5,
							a = at,  
							κ = κt, 
							color = :violet
				),
				BulkSpecies(;name = "CO₂",
							z = 0, 
							D = 1.91e-9, 
							c_bulk = 0.033, 
							a = at, 
							κ = 0, 
							color=:red
				),
				BulkSpecies(;name = "OH⁻",
							z = -1, 
							D = 5.273e-9, 
							c_bulk = 10^(pH-14), 
							a = at, 
							κ = κt, 
							color = :green
				),
				BulkSpecies(;name = "H⁺", 
							z = 1, 
							D = 9.310e-9, 
							c_bulk = 10^(-pH), 
							a = at, 
							κ = κt, 
							color = :gray
				),
				BulkSpecies(;name="CO",
							z = 0,
							D = 2.23e-9,
							c_bulk = 0.0,
							a = at, 
							κ = 0,  
							color=:blue
				)
		]
		push!(bulk, make_eneutral(bulk;name="K⁺", 
									   z = 1, 
									   D = 1.957e-9,
		
									   a = ak, 
									   κ = κk, 
									   color = :orange
								  )
		)
		sort(bulk, by=x->species_dict[x.name])
	end
end;

# ╔═╡ 06f52599-7006-4a5c-ba86-0b668b6952c9
begin
	function create_markdown(bulk::Vector{BulkSpecies})
		table = """
| Name | z | D | c_bulk | a | v | M | κ | color |
|------|---|---|--------|---|---|---|---|-------|
"""
		for sp in bulk
			(; name, z, D, c_bulk, a, v, M, κ, color) = sp
			table *= @sprintf("| %s | %i | %1.3e | %1.3e | %1.3e | %1.3e | %1.3e | %1.2f | %s |\n", name, z, D, c_bulk, a, v, M, κ, color) 
		end
		Markdown.parse(table)
	end
	create_markdown(bulk)
end

# ╔═╡ f8e6c01b-e64e-4fa2-a84a-e20f9b60287c
function capsplot_κ(vis, sys; n::Int=7)
    color = [RGB(0, 0, (i/n)) for i in 1:n]
    dls = LiquidElectrolytes.DLCapSweepResult[]
    κ_values = [0, 1, 5, 10, 20, 30, 40]

    sys = deepcopy(sys)
    κ_original = deepcopy(electrolytedata(sys).κ)

    for j in 1:n
        electrolytedata(sys).κ .= κ_values[j]
        electrolytedata(sys).c_bulk .= [0.05, 0.05] * ufac"mol/dm^3"

        try
            result = caps(sys)
            push!(dls, result)
            
            scalarplot!(
                vis,
                result.voltages,
                result.dlcaps / (μF / cm^2),
			    linestyle = :solid,
                color = color[j],
                clear = false,
                label = "κ = $(κ_values[j])"
            )
        catch e
            @warn "caps failed at κ=$(κ_values[j])" exception=e
        end
    end
	electrolytedata(sys)
	
    electrolytedata(sys).κ .= κt # 복구
end

# ╔═╡ 8c367e8f-df43-4f21-bef0-55060f36f44e
function conc_time_func(result, scan)
    species = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    times = result.tsol.t
    nt = length(times)
    nspecies = 7

    conc_at_electrode = [result.tsol[i, 1, t] / (mol / dm^3) for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = (650, 400))
    ax = Axis(
        fig[1, 1];
        xlabel = L"time / s",
        ylabel = L"c_{i,\,\text{electrode}} / (mol/dm^3)",
        yscale = log10
    )

    for i in 1:nspecies
    	yvals = conc_at_electrode[i, :]
    	yvals_safe = [c > 0 ? c : NaN for c in yvals]
    	lines!(ax, times, yvals_safe, color = colors[i], label = species[i])
	end


    Legend(fig[1, 2], ax; labelsize = 10, backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5))
    fig
end


# ╔═╡ 2ce5aa45-4aa5-4c2a-a608-f581266e55f0
begin
	function addplot(vis, sol, vshow)
		species = getproperty.(bulk, :name)
		colors = getproperty.(bulk, :color)
		
		scale = 1.0 / (mol / dm^3)
	    title = @sprintf("Φ_we=%+1.2f [V vs. SHE]", vshow)
	
		if useonly_pH
			i = findfirst(isequal("H⁺"), species)
			scalarplot!(vis, 
					    grid.components[XCoordinates] .+ 1.0e-14, 
					    log10.(sol[ihplus, :] * scale), 
					    color = colors[i],
					    label = species[i],
					    clear = true,
						title = title)
		else
			scalarplot!(vis, 
						grid.components[XCoordinates] .+ 1.0e-14, 
						log10.(sol[1, :] * scale), 
						color = colors[1],
						label = species[1],
						clear = true,
						title = title)
			for ia = 2:nc			
				scalarplot!(vis, 
						    grid.components[XCoordinates] .+ 1.0e-14, 
						    log10.(sol[ia, :] * scale), 
						    color = colors[ia],
						    label = species[ia],
						    clear = false,)
			end
		end
	end
	function plot1d(result, celldata, vshow; df_compare = nothing)
		tsol 	= LiquidElectrolytes.voltages_solutions(result)
		vis 	= GridVisualizer(;
								 size 	= (600, 300),
								 clear 	= true,
								 legend 	= :rt,
								 limits 	= (-14, 2),
								 xlimits    = (10e-12, 80 * μm),
								 xlabel 	= "Distance from electrode [m]",
	 							 ylabel 	= "log c(aᵢ)", 
								 xscale 	= :log,)
	    addplot(vis, tsol(vshow), vshow)
		if !isnothing(df_compare)
			addplot(vis, df_compare)
		end
		reveal(vis)
	end

	function plot1d(result, celldata)
    	tsol  	= LiquidElectrolytes.voltages_solutions(result)
		vis  	= GridVisualizer(; 
								 size 	= (600, 300),
								 clear 	= true,
							 	 legend = :rt,
								 limits = (-14, 2),
								 xlimits= (10e-12, 80 * μm),
								 xlabel = "Distance from electrode [m]",
 								 ylabel = "log c(aᵢ)", 
								 xscale = :log,)
	
		vrange = result.voltages[end:-5:1]
		movie(vis, file="concentrations.gif", framerate=3) do vis
		for vshow_it in vrange
			addplot(vis, tsol(vshow_it), vshow_it)
			reveal(vis)
		end
		end
		isdefined(Main, :PlutoRunner) && LocalResource("concentrations.gif")
	end
end

# ╔═╡ 9d814b85-a5b6-42e5-abf4-15500bbdb717
begin
	γ_key = user_input.mode 
	γ_mode = γ_key == "Stefan_γ" ? Stefan_γ! : DGML_γ!
end;

# ╔═╡ 8a1047fa-e483-40d9-8904-7576f30acfb4
begin
	const γ_cache = DiffCache(zeros(nc), 12)
	
	function reaction(
		f, 
		u::VoronoiFVM.NodeUnknowns{Tv, Tc, Tp, Ti}, 
		node, 
		data
	) where {Tv, Tc, Tp, Ti}  
		
		(; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, pscale, p_bulk) = data


		γ = get_tmp(γ_cache, u[ico2])
		# compute activity coefficients according to the approach in Ringe et al.
		activity_coefficient!(γ, u, data, γ_mode)
		
		# compute activity coefficients according to the approach in Dreyer et al.
		# p = u[ip] * pscale-p_bulk
    	# c0, bar_c = c0_barc(u, data)
		# for ic in 1:nc
		#	Mrel = M[ic] / M0
		#	barv = v[ic] + κ[ic] * v0
		#	tildev = barv - Mrel * v0
		# 	γ[ic] = exp(tildev * p / (RT)) * (bar_c / c0)^Mrel*(1/bar_c) /v0
		# end


		
		
		@views f_buffer!(
			f[ibufferstart:ibufferend], 
			u[ibufferstart:ibufferend],
			γ[ibufferstart:ibufferend],
			nothing
		)
		nothing
	end
end;

# ╔═╡ 91113083-d80e-4528-be41-82d10f6860fc
begin
	const ps_cache = DiffCache(zeros(8), 13)
	const us_cache = DiffCache(zeros(isurfaceend-isurfacestart+1), 13)
	
	function we_breactions(f, 
			u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
			bnode, 
			data
		) where {Tval, Tv, Tc, Tp, Ti}
		(; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we) = data
				
		γ = get_tmp(γ_cache, u[ico2])
		γ_co2 	= activity_coefficient!(γ, u, data, γ_mode; idxs=[ico2])[ico2]
		γ_co 	= activity_coefficient!(γ, u, data, γ_mode; idxs=[ico])[ico]
		#γ_co 	= 1.0
		σ 			= C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)
		local_pH 	= -log10(u[ihplus] / (mol/dm^3))

	
		ps = get_tmp(ps_cache, u[iϕ])
		#for (p, default_value) in odesys.defaults
		#	ps[paramsidx[p]] = default_value
		#	println(default_value)
		#end

		ps = get_tmp(ps_cache, u[iϕ])
		ps[paramsidx[Symbolics.rename(odesys.σ, :σ)]] = σ
		ps[paramsidx[Symbolics.rename(odesys.γCO2_aq, :γCO2_aq)]] = γ_co2 
		ps[paramsidx[Symbolics.rename(odesys.aH2O_g, :aH2O_g)]] = aH₂O 
		ps[paramsidx[Symbolics.rename(odesys.ϕ, :ϕ)]] = u[iϕ] 
		ps[paramsidx[Symbolics.rename(odesys.ϕ_we, :ϕ_we)]] = ϕ_we 
		ps[paramsidx[Symbolics.rename(odesys.local_pH, :local_pH)]] = local_pH 
		ps[paramsidx[Symbolics.rename(odesys.γCO_aq, :γCO_aq)]] = γ_co 
		ps[paramsidx[Symbolics.rename(odesys.βCOOHΔH2OΔele_t, :βCOOHΔH2OΔele_t)]] = 0.59 


	    #println[1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))]

		
		@views f_microkinetics!(
			f[isurfacestart:isurfaceend], 
			u[isurfacestart:isurfaceend],
			ps,
			nothing
		)
		#if ϕ_we == -0.9
		#	println(u[isurfacestart:isurfaceend])
		#end

		
		global coeff = u[isurfacestart:isurfaceend]
		# conversion from turnover frequency (appropriate for change in coverage) to 
		# production rate (per unit area) (approprite for change in concentration)
		# by S = number of free catalyst sites in mole per unit area
		f[ico2] *= S
		f[iohminus] *= S
		f[ico] *= S
		f[ikplus] *= S
	end
end

# ╔═╡ e510bce3-d33f-47bb-98d6-121eee8f2252
elydata_Gold = ElectrolyteData(;
                               	nc = size(bulk)[1],
								na    = na,
								z     = getproperty.(bulk, :z),
							  	D     = getproperty.(bulk, :D),
							  	T     = T,
							  	eneutral=false,
							  	κ     = getproperty.(bulk, :κ),
	                            c_bulk= getproperty.(bulk, :c_bulk),
							    v0 	  = v0,
								v     = getproperty.(bulk, :v),
								M0 	  = M0,
								M     = getproperty.(bulk, :M),
							  	Γ_we  = Γ_we,
							  	Γ_bulk= Γ_bulk,
							   	actcoeff! = γ_mode
							   )

# ╔═╡ 924f8f5d-2cb0-4381-a522-509ff4c002b6
begin
	model_key = user_input[:model_choice]

	model = model_key == "Gold_Model" ? elydata_Gold :
	        model_key == "Landstorfer_NaClO₄ model" ? elydata_NaClO₄ :
	        model_key == "Landstorfer_NaF model" ? elydata_NaF :
	        error("Unknown model choice: $model_key")
end

# ╔═╡ dc203e95-7763-4b13-8408-038b933c5c9c
function pnp_bcondition(
	f,
	u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
	bnode,
	data
) where {Tval, Tv, Tc, Tp, Ti}
	
	(; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

	bulkbcondition(f, u, bnode, data; region = Γ_bulk)

    ## Dirichlet ϕ=ϕ_we at Γ_we
    #boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)	
	
	# Robin b.c. for the Poisson equation
	boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))

	if bnode.region == Γ_we && model == elydata_Gold
			we_breactions(f, u, bnode, data)
	end
	nothing
end;

# ╔═╡ 84d1270b-8df5-4d5d-a153-da4ffdb1d283
function simulate_CO2R(grid, celldata; voltages = (-1.5:0.1:0.0) * V, kwargs...)
    kwargs 	 	= merge(solver_control, kwargs) 
    cell        = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)
	ivresult    = ivsweep(cell; voltages, store_solutions=true, kwargs...)

	cell, ivresult
end;

# ╔═╡ 45451a14-17b9-4754-b56b-c9b8e8cce1b4
begin
	reaction_arg = model == elydata_Gold ? (; reaction=reaction) : NamedTuple()
	sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = model, reaction_arg...)

end

# ╔═╡ 675ab1d0-a4e8-44a1-9d16-c2293e802868
sys_pb = PBSystem(grid; celldata = deepcopy(model), bcondition = pb_bcondition)

# ╔═╡ f0aecbbc-3c8c-4984-8704-fd79f986beb2
begin
	try
		vis_κ = GridVisualizer(Plotter = CairoMakie, legend = :lt, size = (650, 650))
		Low_κ = scalarplot!(
		    vis_κ,
		    Landstorfer_Low_κ.voltages,
		    Landstorfer_Low_κ.dlcaps,
		    color = :black,
		    linestyle = :dot,
		    label = "κ = 0",
			xlimits = (-0.5, 0.5),
	        xlabel = "φ / (V vs φ_pzc)",
	        ylabel = "dlcaps / (μF / cm²)",
		    clear = true,  
		)
		High_κ = scalarplot!(
		    vis_κ,
		    Landstorfer_High_κ.voltages,
		    Landstorfer_High_κ.dlcaps,
		    color = :blue,
		    linestyle = :dot,
		    label = "κ = 40",
		    clear = false,
		)
		
		capsplot_κ(vis_κ, sys_pb)
		reveal(vis_κ)
	catch
	end
end

# ╔═╡ 5df5ee46-b0d3-47a8-835b-b3f21a3cab34
is_Landstorfer = model != elydata_Gold

# ╔═╡ 36e756a9-4d9b-40ef-9e37-d86f1194cc51
function capscalc(sys, molarities)
    result = []
	vrange = range(-1, 1, length = 201)
	if is_Landstorfer
		for imol in 1:length(molarities)
		    if !isa(sys, AbstractElectrochemicalSystem)
		     	data = sys.physics.data
		        set_molarity!(data, molarities[imol])
		      	t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 		101)
		    else
		        data = sys.vfvmsys.physics.data
		        data.c_bulk .= molarities[imol] * ufac"mol/dm^3"
		    	t = @elapsed r = dlcapsweep(
		                sys,
		     	        voltages = range(-1, 1, length = 201)
		        )
		        volts = vrange
		        caps = r.dlcaps
		    end
		    cdl0 = dlcap0(data)
		    @info "elapsed=$(t)"
		    push!(result, (voltage_range = volts, dlcaps = caps, cdl0 = cdl0, molarity = 	molarities[imol]))
		end
	else
	   	if !isa(sys, AbstractElectrochemicalSystem)
	   		data = sys.physics.data
	        t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
		else
	        data = sys.vfvmsys.physics.data
	        t = @elapsed r = dlcapsweep(
	    	        sys,
	                voltages = range(-1, 1, length = 201)
	        )
			volts = vrange
			caps = r.dlcaps
	    end
	    cdl0 = dlcap0(data)
	    @info "elapsed=$(t)"
	    push!(result, (voltage_range = volts, dlcaps = caps, cdl0 = cdl0))
	end
    return result
end


# ╔═╡ 6c0608c4-a785-471a-8e03-16a5b16508b8
result_pb = capscalc(sys_pb, molarities)

# ╔═╡ 1ff59725-ba1e-4309-9b8c-19a30adb36db
result_pnp = capscalc(sys_pnp, molarities)

# ╔═╡ 50ccc291-3625-4639-afd1-5209899d904e
let
	if is_Landstorfer
		f = Figure()
	    result = result_pb
	    l = 1 / length(result)
		ϕ0_pzc = 0.972
		
		ax = Axis(f[1, 1], xlabel="φ / (V vs φ_pzc)", ylabel="dlcaps / (μF / cm²)", title="CSV Plot")
	
		if model_key == "Landstorfer_NaClO₄ model"
			Low_c0 = lines!(ax, Landstorfer_NaClO₄_5mM.voltages .+ ϕ0_pzc, Landstorfer_NaClO₄_5mM.dlcaps, color = :darkblue, linestyle = :dash)
			High_c0 = lines!(ax, Landstorfer_NaClO₄_100mM.voltages .+ ϕ0_pzc, Landstorfer_NaClO₄_100mM.dlcaps, color = :red, linestyle = :dash)
		else	
			Low_c0 = lines!(ax, Landstorfer_NaF_5mM.voltages .+ ϕ0_pzc, Landstorfer_NaF_5mM.dlcaps, color = :darkblue, linestyle = :dash)
			High_c0 = lines!(ax, Landstorfer_NaF_100mM.voltages .+ ϕ0_pzc, Landstorfer_NaF_100mM.dlcaps, color = :red, linestyle = :dash)
		end
		
	    k = []  
	    legend_labels = [model_key*"\t 5mM", model_key*"\t 100mM"]  
		
		for i in 1:length(result)
		    c = RGB(i * l, 0.0, 1 - i * l)
		    v = result_pb[i].voltage_range
		    cdl = result_pb[i].dlcaps / (μF / cm^2)
		    
		    minlength = min(length(v), length(cdl))
		    push!(k, lines!(ax, v[1:minlength], cdl[1:minlength], color=c, label="Result $i"))
		
		    m = molarities[i]
		    push!(legend_labels, "LiquidElectrolyte $m M")
		end
	    
		Legend(f[1, 1], [Low_c0, High_c0, k...], legend_labels, halign = :left, valign =:top, tellheight = false, tellwidth = false, framevisible = false)  
	    
		f 
	else
		results = [
		    ("Poisson-Boltzmann", result_pb),
		    ("Poisson-Nernst-Planck", result_pnp)
		]
		plots = []
		
		for (name, result) in results
		    try
		        if !isempty(result)
		            push!(plots, (name, result))
		        end
		    catch
		        continue
		    end
		end
		
	    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, title="Compare PB & PNP")
	    capsplot_v(vis, plots)
	    reveal(vis)
	end
end

# ╔═╡ b21c8394-f847-477e-ac8f-713398b81166
function capsplot(vis, result, title)
    if is_Landstorfer
        hmol = 1 / length(result)
        for imol in 1:length(result)
            c = RGB(imol * hmol, 0, 1 - imol * hmol)
            voltages = result[imol].voltage_range[1:201]
            caps = vec(result[imol].dlcaps)[1:201]

            scalarplot!(
                vis, voltages, caps / (μF / cm^2);
                color = c, clear = false, label = "$(result[imol].molarity)M",
                markershape = :none, title = title,
                xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)"
            )

            scalarplot!(
                vis, [0], [result[imol].cdl0] / (μF / cm^2);
                clear = false, markershape = :circle, markersize = 8, label = ""
            )
        end
    else
        voltages = result[1].voltage_range[1:201]
        caps = vec(result[1].dlcaps)[1:201]

        scalarplot!(
            vis, voltages, caps / (μF / cm^2);
            limits = (-1, 100), xlimits = (-1.1, 1.1),
            color = :green, clear = false, label = "$title",
            title = title, markershape = :none, yscale = 10,
            xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)"
        )
        # scalarplot!(
        #     vis, [0], [result[1].cdl0] / (μF / cm^2);
        #     clear = false, markershape = :circle, markersize = 8, label = ""
        # )
    end
    return vis
end;

# ╔═╡ 4f991d6d-3a3f-45d8-b2e0-662c5292251c
let
    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, layout = (1, 2), size = 	(650, 350))
    capsplot(vis[1, 1], result_pb, "Poisson-Boltzmann")
    capsplot(vis[1, 2], result_pnp, "Poisson-Nernst-Planck")

    reveal(vis)
end

# ╔═╡ 11b12556-5b61-42c2-a911-4ea98a0a1e85
cell, result = simulate_CO2R(grid, model; voltages)

# ╔═╡ 659091d3-60b2-4158-80e2-cd28a492e870
(~, default_index) = findmin(abs, result.voltages .+ 0.9 * ufac"V");

# ╔═╡ 15fadfc2-3cf8-4fda-9aed-a79c602b1d51
plot1d(result, celldata)

# ╔═╡ 1cd669ac-05eb-48b2-b457-8c395cd5807d
let
	table = readdlm("./catmap_CO2R_data/IV-Ringe-digitized.csv", ',', Float64, '\n')
	df = Dict(:voltage => table[:,1], :current => table[:,2])
	plotcurr(result; df=df)
end

# ╔═╡ 32eb1122-5013-4a8e-be54-18a30c151515
if runregtest
	sresult = load("./data/regressionresults.jld2")["regressionresults"]

	vidxs_result = [findfirst(isequal(v), result.voltages) for v in voltages[1:end-1]]
	vidxs_sresult = [findfirst(isequal(v), sresult.voltages) for v in voltages[1:end-1]]
	
	if any(isnothing.(vidxs_sresult))
		throw(ArgumentError("For the full regression test use the applied voltages  -1.5:0.1:0.0"))
	end

	@testset begin
	@testset "Concentrations" begin
		@testset "$(bulk[ia].name)" for ia in 1:nc
			@testset "U=$(result.voltages[vidx_result])" for (vidx_result, vidx_sresult) in zip(vidxs_result, vidxs_sresult)	
				@test all(isapprox(
					result.solutions[vidx_result][ia,:], sresult.solutions[vidx_sresult][ia,:], 
					rtol = 1.0e-5
				))
			end
		end
	end

	@testset "Currents" begin
		for (vidx_result, vidx_sresult) in zip(vidxs_result, vidxs_sresult)
			for (j_result, j_sresult) in zip(result.j_we[vidx_result][iohminus], 
											sresult.j_we[vidx_result][iohminus])
				@test isapprox(j_result, j_sresult, rtol=1.0e-5)
			end
		end
	end
	end
end;

# ╔═╡ b4aaf070-d4ab-409a-b1e8-f5469b9f398b
begin 
	sawtooth = SawTooth(
        scanrate = parse(Float64, user_input.scanrate),
        vmin = user_input.vmin , vmax = user_input.vmax
    )
	    const nperiods = user_input.nperiods

end

# ╔═╡ e50fe651-11d4-45ee-89dd-371a7fbc097e
function sweep(pnpdata; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    celldata.eneutral = eneutral
	reaction_arg = model == elydata_Gold ? (; reaction) : NamedTuple()
    pnpcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = model, reaction_arg...)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
    )

end

# ╔═╡ 9fb47b83-a853-4316-bb8d-30e65b16ef78
pnpresult = sweep(model; eneutral = false, tunnel = false)

# ╔═╡ cb9b0158-f17d-4136-8994-360f3078c7df
let
    fig = Figure(size = (600, 200))
    ax = Axis(fig[1, 1], yscale = log10)
    T = pnpresult.times
    #lines!(ax,T, voltages.(T))
    lines!(ax, T[2:end], T[2:end] - T[1:(end - 1)])
    fig
end

# ╔═╡ 52a5bbd2-0278-4d92-95f1-367797f636e6
CVPlot!(pnpresult, model)

# ╔═╡ 4e894347-2ce6-4c5f-a06e-7f1af1983bbc
conc_time_func(pnpresult, sawtooth)

# ╔═╡ 2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
plot_concentration_profile(pnpresult, X, bulk)

# ╔═╡ 315dd351-9d68-48f1-aa7a-8f43f3dec6ac
floataside(
    md"""
    __Input Voltage Index:__ $(@bind vindex PlutoUI.Slider(1:5:length(result.voltages), default=default_index))

    __Input Time Index:__ $(@bind it PlutoUI.Slider(1:length(pnpresult.tsol.t)-1, show_value=false))
    """,
    top = 700
)


# ╔═╡ b64c0d67-016d-4bca-9fae-150cf50efc77
let
    species = getproperty.(bulk, :name)
	colors = getproperty.(bulk, :color)
	
	
	ZR = zstr(user_input.zR)
    ZO = zstr(user_input.zR + user_input.n)

    XX = (X[2:end] / nm)
    uu = pnpresult.tsol[:, :, it + 1]
    ru = uu / (mol / dm^3)

    fig = Figure(size = (650, 400))
    ax1 = Axis(
        fig[1, 1];
        #       xlabel=L"x/nm",
        ylabel = L"c/(mol/dm^3)",
        xscale = log10,
        yscale = log10,
        title = "V=$(pnpresult.voltages[it])"
    )
    xlims!(ax1, XX[1], L / nm)
    ylims!(ax1, 1.0e-25, 1.0e2)

	for i in 1:nc
	    yvals = [c > 0 ? log10(c) : NaN for c in ru[i, 2:end]]
	    lines!(ax1, XX, 10 .^ yvals, color = colors[i], linestyle = :solid, label = species[i])
	end


    Legend(
        fig[1, 2], ax1; labelsize = 10,
        backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)

    ) 
    fig
	#println(species, colors)
end

# ╔═╡ c4876d26-e841-4e28-8303-131d4635fc23
md"""
Potential at the working electrode 
$(vshow = result.voltages[vindex]; @sprintf("%+1.4f", vshow))
"""

# ╔═╡ 5dd1a1e6-7db1-479e-a684-accec53ce06a
plot1d(result, celldata, vshow)

# ╔═╡ 7454f68a-64dc-4676-b2b2-ed8fcb35d81e
floataside(
	md"""
	**Voltage:** $(round(result.voltages[vindex], digits=3)) V    
	**Time:** $(round(pnpresult.tsol.t[it+1], digits=3)) s
	""",
	top = 649
)

# ╔═╡ Cell order:
# ╠═91ac9e35-71eb-4570-bef7-f63c67ce3881
# ╠═a94bc4e1-506f-4e40-bfe8-1ce7e6093974
# ╟─beae1479-1c0f-4a55-86e1-ad2b50174c83
# ╟─ab2184fc-0279-46d9-9ee4-88fe3e732789
# ╠═7316901c-d85d-48e9-87dc-3614ab3d81a5
# ╟─6b7cfe87-8190-40a5-8d25-e39ef8d55db5
# ╠═5a146a44-03dc-45f3-ae15-993d11c2edac
# ╠═00947475-c96e-4ecc-a1ef-5be5e3e3c864
# ╠═ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
# ╟─06f52599-7006-4a5c-ba86-0b668b6952c9
# ╟─4b64e168-5fe9-4202-9657-0d4afc237ddc
# ╟─de2c826d-6c05-47cf-b5f5-44a00ea9889c
# ╟─d8f00649-e2ed-4bdd-853f-05268f0d5353
# ╠═47b36c81-b57e-4dd0-a22f-999e4fd3ac9f
# ╠═1e877f17-0219-45f1-b640-3a25ae085dbd
# ╠═8a1047fa-e483-40d9-8904-7576f30acfb4
# ╟─8912f990-6b02-467a-bd11-92f94818b1c7
# ╟─a8157cc1-1761-4b11-a37c-9e12a9ca695e
# ╠═6b5cf93c-0df3-4a18-8786-502361736838
# ╟─d2c0642d-dfa5-4a76-bd36-ac4a735a3299
# ╟─06d45088-ab8b-4e5d-931d-b58701bf8464
# ╠═91113083-d80e-4528-be41-82d10f6860fc
# ╟─d0093605-0e35-4888-a93c-8456c698e6f0
# ╟─f0b5d356-6b97-4878-98de-bee5f380d41a
# ╟─e3eda42f-e2f3-4c10-81c4-610246ca528d
# ╟─7b87aa2a-dbaf-441c-9ad7-444abf15f664
# ╠═2b9d9bfd-d660-4b4d-8f0b-b6b5bcc0dbfa
# ╟─6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
# ╠═e510bce3-d33f-47bb-98d6-121eee8f2252
# ╠═848b7aeb-968f-4116-8038-b61276f02b6c
# ╠═f18dc873-1c9d-46d3-9596-92d28705e894
# ╟─53f12821-7d8d-4971-87fd-ad4689ec62a5
# ╟─e1e0ca0f-7f88-40f0-850e-590b25da0331
# ╠═0db74a70-af86-492c-affb-9de62ffe4455
# ╟─4f388fe0-6bc8-4a29-bccc-fa725e62c6a7
# ╟─5d179c52-43d7-4bcb-a2df-93c5806876fa
# ╠═9d814b85-a5b6-42e5-abf4-15500bbdb717
# ╠═5f17b4f7-54d6-4ad0-9886-252854840a80
# ╟─2a20d9be-6c1e-4c1f-8bb6-a7693800732d
# ╟─4f7ec19d-cd60-4c2b-a766-7557caa471c0
# ╠═952a26ce-2610-48cc-9158-eda816da3a1c
# ╠═924f8f5d-2cb0-4381-a522-509ff4c002b6
# ╠═dc203e95-7763-4b13-8408-038b933c5c9c
# ╠═45451a14-17b9-4754-b56b-c9b8e8cce1b4
# ╠═9a4e01d9-f469-4427-bf4c-883adb67ae24
# ╠═675ab1d0-a4e8-44a1-9d16-c2293e802868
# ╠═6c0608c4-a785-471a-8e03-16a5b16508b8
# ╠═1ff59725-ba1e-4309-9b8c-19a30adb36db
# ╟─d76d8413-c019-4728-b182-7f7cb78dede4
# ╠═4f991d6d-3a3f-45d8-b2e0-662c5292251c
# ╠═50ccc291-3625-4639-afd1-5209899d904e
# ╠═5df5ee46-b0d3-47a8-835b-b3f21a3cab34
# ╟─783e2058-c720-4f31-8e51-7c313813924c
# ╠═f8e6c01b-e64e-4fa2-a84a-e20f9b60287c
# ╠═f0aecbbc-3c8c-4984-8704-fd79f986beb2
# ╟─9598e2c6-521e-4f8d-82d8-a836809736f3
# ╠═f2043f2c-f3c8-4b0f-944c-7b55624dac08
# ╠═36e756a9-4d9b-40ef-9e37-d86f1194cc51
# ╠═8bfdf2f5-c80a-4ce0-a8e1-b315affffb5f
# ╠═b21c8394-f847-477e-ac8f-713398b81166
# ╟─5c808c71-6094-49d7-8215-e88262f34e1f
# ╟─da8390d1-47e8-451f-b12b-45b8aca7b6ec
# ╠═9fb47b83-a853-4316-bb8d-30e65b16ef78
# ╠═39c8ef0d-aac2-4c7f-8004-4166c460ebc5
# ╟─eb920b6e-86a6-4dd6-8e66-6b7e27d81257
# ╠═cb9b0158-f17d-4136-8994-360f3078c7df
# ╠═52a5bbd2-0278-4d92-95f1-367797f636e6
# ╠═4e894347-2ce6-4c5f-a06e-7f1af1983bbc
# ╠═b64c0d67-016d-4bca-9fae-150cf50efc77
# ╠═2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
# ╟─3f30fae0-18d4-4e5f-9618-cbd9852d7857
# ╟─bb00b5bb-326e-47f9-a4f4-e7b4f29dd1f2
# ╠═b4aaf070-d4ab-409a-b1e8-f5469b9f398b
# ╠═e50fe651-11d4-45ee-89dd-371a7fbc097e
# ╠═82baec54-048a-4851-8808-dd8c31eae4b9
# ╠═8c367e8f-df43-4f21-bef0-55060f36f44e
# ╠═b7cb5183-65e8-4ee8-af86-2bedd11daecc
# ╟─842b074b-f808-48d8-8dc5-110ddd907f90
# ╟─31298257-d35a-4f6f-8a76-ff00d5361ced
# ╠═72269ec4-a56e-46d9-85c8-0dd8ccaf43e1
# ╠═84d1270b-8df5-4d5d-a153-da4ffdb1d283
# ╠═11b12556-5b61-42c2-a911-4ea98a0a1e85
# ╟─7a02463d-cfd9-4648-af53-f1e65d46733f
# ╟─114d2324-5289-4e44-8d77-736a9bdec365
# ╟─659091d3-60b2-4158-80e2-cd28a492e870
# ╟─c4876d26-e841-4e28-8303-131d4635fc23
# ╠═5dd1a1e6-7db1-479e-a684-accec53ce06a
# ╠═15fadfc2-3cf8-4fda-9aed-a79c602b1d51
# ╠═1cd669ac-05eb-48b2-b457-8c395cd5807d
# ╟─c1d2305e-fb8b-4845-a414-08fff84aa9b0
# ╠═2ce5aa45-4aa5-4c2a-a608-f581266e55f0
# ╠═d5ab1a28-3a60-49d9-bb3e-ca589b1c79fd
# ╟─686ac3dc-c191-4575-ba0c-d4c2551474b5
# ╟─d1ab199f-1a40-4377-bca3-7f72f3cde3a9
# ╠═32eb1122-5013-4a8e-be54-18a30c151515
# ╟─de144adb-a467-4077-8cb1-d86462f56110
# ╟─d0985ca6-fef5-4b67-9ad6-f51d84b595b4
# ╟─8ae53b8a-0fb3-4c1c-8e5f-a3782a85141c
# ╟─9d7d4d68-c9cc-4a42-a99d-ae25a1ab554c
# ╠═e5fc814f-a8e1-41ef-b81a-c3b0839a2f87
# ╠═315dd351-9d68-48f1-aa7a-8f43f3dec6ac
# ╠═7454f68a-64dc-4676-b2b2-ed8fcb35d81e
# ╠═ae50877a-20da-4b1a-acd5-cc0a73e428da
# ╠═3ac837b8-559b-41c2-8f83-1331839dcf7e
