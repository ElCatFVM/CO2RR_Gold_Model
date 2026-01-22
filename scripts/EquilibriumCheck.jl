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
    Pkg.activate(joinpath(@__DIR__, ".."))
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

# ╔═╡ 68fe205c-9dd0-441b-9f12-3ddc12ec0a0d
begin
	include(joinpath(@__DIR__, "..", "src", "Cyclic_Voltammetry.jl"))
	include(joinpath(@__DIR__, "..", "plots", "cvplot.jl"))
end;

# ╔═╡ a94bc4e1-506f-4e40-bfe8-1ce7e6093974
pkgdir(CatmapInterface)

# ╔═╡ 22f2574c-abbe-4a06-89e9-14635fb30932
pkgdir(LiquidElectrolytes)

# ╔═╡ bd8134d8-5a69-486e-8429-7cf810b3ccbe
Pkg.status()

# ╔═╡ beae1479-1c0f-4a55-86e1-ad2b50174c83
md"""
## Setup
"""

# ╔═╡ ab2184fc-0279-46d9-9ee4-88fe3e732789
md"""
### Units
"""

# ╔═╡ 7316901c-d85d-48e9-87dc-3614ab3d81a5
@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;

# ╔═╡ 6b7cfe87-8190-40a5-8d25-e39ef8d55db5
md"""
### Data
"""

# ╔═╡ 5a146a44-03dc-45f3-ae15-993d11c2edac
begin
	@phconstants N_A c_0 k_B e h
	const F = N_A * e

	const voltages = (-1.15:0.1:-0.0) * V
	
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
	const v0        = N_A * (8.2 * Å)^3#18.048 * ufac"cm^3/mol"# # 1 / (55.4 * ufac"M") #
	
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

# ╔═╡ a1b895c2-d3d1-4afd-bd85-da575effef1c
@which boundary_neumann!

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
	function BulkSpecies(;name, z, c_bulk=nothing, a, D, κ=nothing, color)
		D *= m^2/s
		c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol/dm^3
		a *= Å #8.2 * Å#8.2 * Å            # (v0/N_A)^(1/3)
		v = N_A * (a * Å)^3 #v0 * (κ + 1)
		M = M0 * v
		BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
	end
	function make_eneutral(bulk_species::Vector{BulkSpecies}
						  ;name, z, D, κ=0.0, a=8.2, v = v0 * (κ * abs(z) + 1)
, M=M0*v, color)
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

# ╔═╡ a72aa512-941c-4815-b69f-eaa9d94e6b06
# ╠═╡ disabled = true
#=╠═╡
bulk
  ╠═╡ =#

# ╔═╡ 68491839-1398-4807-9a50-07c87bb349cf


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
	catmap_params 		= CatmapInterface.parse_catmap_input("../data/catmap_CO2R_data/catmap_CO2R_template.mkm")
	rn 					= create_reaction_network(catmap_params; symbolic_formation_energies = false)
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

# ╔═╡ 42efe97d-bf43-42b3-8ff0-e3d98a7597f1
ex_grid = ExtendableGrids.simplexgrid(ExtendableGrids.geomspace(0, 4000* μm, 1.0e-6 	* μm, 1.0	* μm))

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
	
    (; Mrel, tildev, v0, v, RT, v0, cspecies, rexp, c_bulk) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
	
    for ic in cspecies
        γ[ic] = rexp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc))
    end

	#for ic in cspecies
    #    γ[ic] = rexp(tildev[ic] * p / RT) * (barc)^(Mrel[ic]-1) 
    #end
    #for ic in cspecies
    #    γ[ic] = v0 * ( 1.0 / (1 - sum(c_bulk[i] * v[i] for i in 1:nc) / (mol/dm^3)))
    #end
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
        γ[ic] = (1.0 / (1 - sum(c[i] * v[i] for i in 1:nc))) # (mol/dm^3)))
    end
    return γ
end

# ╔═╡ 161a810d-c05e-42ad-97ab-131059d6784a
function Potassium_γ!(γ, c, p, electrolyte)
	
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp, c_bulk, v, nc) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
	γ .= 0 
    γ[ikplus] = 1.0 / (1 - v[ikplus] * c[ikplus] / (mol/dm^3))
    
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
molarities = [0.005, 0.1]#, 0.05, 0.1, 0.5] 

# ╔═╡ 4656ee04-ae86-442f-b37c-c5563170f992
md"""
Run double layer curve $(@bind double_layer_curve PlutoUI.CheckBox())
"""

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
		voltages = res[1].voltage_range[1:end]
        caps = vec(res[1].dlcaps)[1:end]
        scalarplot!(
            vis, voltages, caps / (μF / cm^2);
            limits = (-1, 400), xlimits = (-1.1, 1.1),
            color = color[i], clear = false, label = name,
            markershape = :none, yscale = 10,
            xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)"
        )
    end
    return vis
end;

# ╔═╡ 0e61a0f7-1611-4eb3-8fda-3f807a4ffca2
typeof(odesys)

# ╔═╡ 8bfdf2f5-c80a-4ce0-a8e1-b315affffb5f
begin
	
	#computed capacitance plot(Fig 13 & Fig 14)
	Landstorfer_NaF_5mM = CSV.read("../data/E.Acta_Cap_data/Landstorfer_NaF_0.005M.csv", DataFrame);
	Landstorfer_NaF_100mM = CSV.read("../data/E.Acta_Cap_data/Landstorfer_NaF_0.1M.csv", DataFrame);
	Landstorfer_NaClO₄_100mM = CSV.read("../data/E.Acta_Cap_data/Landstorfer_NaClO4_0.1M.csv", DataFrame);
	Landstorfer_NaClO₄_5mM = CSV.read("../data/E.Acta_Cap_data/Landstorfer_NaClO4_0.005M.csv", DataFrame);
	
	#Solvation number plot(Fig 7)
	Landstorfer_Low_κ = CSV.read("../data/E.Acta_Cap_data/Landstorfer_kappa0.csv", DataFrame);
	Landstorfer_High_κ = CSV.read("../data/E.Acta_Cap_data/Landstorfer_kappa40.csv", DataFrame);
end;

# ╔═╡ 44258eea-f114-4dfe-aa61-1e2cac31baa4
function capsplot(vis, ::Nothing, title)
    scalarplot!(
        vis, [NaN], [NaN];
        clear = false,
        markershape = :none,
        label = "(no data)",
        title = title,
        xlabel = L"φ / (V vs φ_{pzc})",
        ylabel = L"dlcaps / (μF / cm²)"
        # xlimits = (-1.1, 1.1), ylimits = (-1, 250)
    )
    return vis
end

# ╔═╡ 5c808c71-6094-49d7-8215-e88262f34e1f
md"""
## Cyclic Voltammetry
"""

# ╔═╡ da8390d1-47e8-451f-b12b-45b8aca7b6ec
md"""
### System Setup
"""

# ╔═╡ ef7212fc-a3d0-4784-b901-219204b79dc0
md"""
#### General CV function
"""

# ╔═╡ b95160b5-18f7-49d9-80be-9159abd2dcd1
md"""
Run general cyclic voltammetry curve $(@bind CV PlutoUI.CheckBox())
"""

# ╔═╡ 39c8ef0d-aac2-4c7f-8004-4166c460ebc5
# ╠═╡ disabled = true
#=╠═╡
nnpresult = sweep(model; eneutral = true, tunnel = false)
  ╠═╡ =#

# ╔═╡ fba53a72-5d27-4db7-8453-41200db29481
md"""
#### RDE-rpm-dependance function
"""

# ╔═╡ 4af42420-4630-42a8-8941-6702c16beb99
md"""
Run RDE-CV calculation $(@bind CV_RPM_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 83defc30-2532-490f-ac19-f883f06e6e2d
md"""
##### RDE-rpm-export function
"""

# ╔═╡ 57a4947d-8de3-42de-b17a-9b769ca0589d
function export_iv_csv_each(Pressure_vec, P, currents, iohminus; dir="iv_8000_csv", xlim=nothing)
    mkpath(dir)
    sanitize(s) = replace(string(s), r"[^\w\.\-]+" => "_")

    for (j, rec) in enumerate(Pressure_vec)
        x = rec.voltages
        y = currents(rec, iohminus) .* cm^2/mA  # convert to mA/cm²

        # optional x-range filter
        if xlim !== nothing
            xmin, xmax = xlim
            mask = (xmin .<= x) .& (x .<= xmax)
            x, y = x[mask], y[mask]
        end

        # label & filename
        label = j == 1 ? @sprintf("%s_sat", string(P[1])) : @sprintf("%s_atm", string(P[j]))
        fname = joinpath(dir, @sprintf("IV_%02d_%s.csv", j, sanitize(label)))

        # DataFrame and write
        df = DataFrame(φ = float.(x), I = float.(y))  # ensure plain Float64
        CSV.write(fname, df)
    end
    nothing
end

# ╔═╡ b2546d23-825e-4356-a640-3fd53852cdcf
# ╠═╡ disabled = true
#=╠═╡
export_iv_csv_each(Pressure_vec, P, currents, iohminus)
  ╠═╡ =#

# ╔═╡ 2d950a96-9404-4b78-9db2-42ae2a4c44bf
md"""
#### pressure varied CV function
"""

# ╔═╡ dfd42e6e-a98e-4759-b988-52dc5a793f15
md"""
Run CO2 pressure varied CV calculation $(@bind CO2_pressure_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 12a4df23-3c63-41d6-bd50-d7209e423cb8
# ╠═╡ disabled = true
#=╠═╡
begin
	if CO2_pressure_varied_checkbox
	    P = [0.0, 0.1, 0.2, 0.3, 0.5, 0.6, 1.0]
	    base_CO2 = elydata_Gold.c_bulk[5]
		Pressure_vec = []
	    for p in P
	        ely_pressure = deepcopy(elydata_Gold)   
			ely_pressure.c_bulk[5] = base_CO2 .* p
	
	        pnp_rec = sweep(ely_pressure; eneutral=true, tunnel=false)
	
	        push!(Pressure_vec, pnp_rec)              
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 05c8e2fa-8650-49d0-a190-b865c0fb3261
md"""
Run bicarbonate pressure varied CV calculation $(@bind CO3_pressure_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 3ba47c4f-b8f8-41e1-bada-8146062b099e
ihco3

# ╔═╡ 3ec78a69-a7b3-4c32-82cb-5b863c88798a
# ╠═╡ disabled = true
#=╠═╡
begin
	if CO3_pressure_varied_checkbox
	    PHCO3 = [0.01, 0.05, 0.1, 0.5, 1]
	    base_CO3 = elydata_Gold.c_bulk[ihco3]
		Pr_HCO3 = []
	    for p in PHCO3
	        ely_pressure = deepcopy(elydata_Gold)   
			ely_pressure.c_bulk[ihco3] = base_CO3 .* p
	
	        pnp_rec = sweep(ely_pressure; eneutral=true, tunnel=false)
	
	        push!(Pr_HCO3, pnp_rec)              
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 25eb8aa3-697e-4538-9472-ceea45fbfbd9
md"""
#### pH varied CV function
"""

# ╔═╡ 11892724-1851-46f2-802d-4da45127b0af
md"""
Run pH varied CV calculation $(@bind pH_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ b3649b03-25cd-4f3e-99c4-85e24ddd3d11
md"""
#### Koper paper Figure 5 CV function
"""

# ╔═╡ fee347ff-5401-4540-a1ce-fc2e8ff0ce63
# ╠═╡ disabled = true
#=╠═╡
begin
    c_bulk2 = [0.2, 1, 5]
    F5 = [8.9, 9.3, 8.8]
    
    F5_vec = Any[]

    for p in 1:length(c_bulk2)
        ely_pH = deepcopy(elydata_Gold)
        ely_pH.c_bulk[5] = elydata_Gold.c_bulk[5] * 0
        ely_pH.c_bulk .*= c_bulk2[p]
        ely_pH.c_bulk[2] = 10^(-F5[p])
        ely_pH.c_bulk[6] = 10^(F5[p]-14)
        pnp_pH = sweep(ely_pH; eneutral=false, tunnel=false)

        push!(F5_vec, pnp_pH)
    end
end
  ╠═╡ =#

# ╔═╡ c048e472-3983-4279-bf60-82784baa145e
md"""
#### Scan Rate CV function
"""

# ╔═╡ 3bdaab98-c0f7-46af-86b7-d68374e8a5d0
md"""
Run scan rate varied CV calculation $(@bind scan_rate_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ d38c2b43-4d8b-4be7-8d77-5a30da384541
# ╠═╡ disabled = true
#=╠═╡
function sweep2(pnpdata, sawtooth; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    celldata.eneutral = eneutral
	reaction_arg = model == elydata_Gold ? (; reaction) : NamedTuple()
    pnpcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = pnpdata, reaction_arg)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
    )

end
  ╠═╡ =#

# ╔═╡ 1f085f56-e0ee-4cb5-a37e-eb82ef3d7589
#=╠═╡
begin
	if scan_rate_varied_checkbox
	    scanrates = [0.002, 0.003, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 10.0] 
	
	    sweep_vec = Vector{Any}(undef, length(scanrates))
	
	    for (i, sr) in pairs(scanrates)
	        sawtooth = SawTooth(
	            scanrate = sr,
	            vmin     = user_input_cv.vmin,
	            vmax     = user_input_cv.vmax,
	            scanup   = user_input_cv.scanup
	        )
	        sweep_vec[i] = sweep2(elydata_Gold, sawtooth; eneutral = false, tunnel = 
								  false)
	    end
	end
end

  ╠═╡ =#

# ╔═╡ eb920b6e-86a6-4dd6-8e66-6b7e27d81257
md"""
### Result Plots
"""

# ╔═╡ 79018ef0-6ab1-4420-9a52-8f8e2812fd40
md"""
#### Time-Voltage plots
"""

# ╔═╡ 278dd577-2d1f-4608-aee4-f7de466cf736
md"""
#### Time-Concentration plots
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


# ╔═╡ 7dd05779-3ffc-471c-9ae9-4bb00b45b7e8
md"""
#### Facet plots
"""

# ╔═╡ de2baeae-eaf6-4565-9ed1-f2eb8c666839
md"""
#### Pressure plots
"""

# ╔═╡ 04790584-5822-460a-be5c-c9efb3bc26b5
#=╠═╡
let
    try
    fig = Figure(size = (1600, 900))
    ax = Axis(fig[1, 1];
        xlabel = L"φ (V vs SHE)",
        ylabel = L"I (mA/cm²)",
        #yscale = log10
    )

    n = length(Pressure_vec)

	cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]

    plots = Makie.AbstractPlot[] 
    labels = String[]
    for (j, rec) in enumerate(Pressure_vec)
        label2 = j == 1 ? "$(P[1])\t\t sat" : "$(P[j])\t pCO2(atm)"
        push!(labels, label2)
        line = lines!(ax, rec.voltages, ((currents(rec, iohminus) .* cm^2/mA)); color = cols[j])
        push!(plots, line)
    end
    Legend(fig[1, 2], plots, labels, "Theoretical"; framevisible = true)
    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 780e8faa-e346-45cb-81f1-34df2a99bc17
#=╠═╡
let
    try
	    fig = Figure(size = (1600, 900))
	    ax = Axis(fig[1, 1];
	        xlabel = L"φ (V vs SHE)",
	        ylabel = L"I (mA/cm²)",
			limits = ((-0.5, 1.5),(-2e-14, 2e-7)),	
	        #yscale = log10
	    )
	
	    n = length(Pr_HCO3)
	
		cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]
	
	    plots = Makie.AbstractPlot[] 
	    labels = String[]
	    for (j, rec) in enumerate(Pr_HCO3)
	        label2 = "$(j)\t pCO2(atm)"
	        push!(labels, label2)
	        line = lines!(ax, rec.voltages, ((currents(rec, ico) .* cm^2/mA)); color = cols[j])
	        push!(plots, line)
	    end
	    Legend(fig[1, 2], plots, labels, "Theoretical"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 8cbc4ced-7ac2-4def-97cf-ff056c1dcb4a
let
	try
	    fig = Figure(size = (1600, 900))
	    ax = Axis(fig[1, 1], 
			      ylabel = L"I (mA/cm²)", 
			      xlabel = L"φ (V vs SHE)",
				  yscale = log10,
			     )

	    keys_sorted = sort(collect(keys(Lresult)))
	    idx = 11

	    L = keys_sorted[idx]
	    rec = Lresult[L]

	    line = lines!(ax, rec.voltages, abs.(currents(rec, ico) .* cm^2/mA))
	    Legend(fig[1, 2], [line], ["L = $(L) μm"], "Theoretical";
	        framevisible = true,
	        labelsize = 12,
	        titlesize = 13,
	        patchsize = (15, 5),
	        rowgap = 1, colgap = 1,
	        patchlineattrs = (linewidth = 30,)
	    )

	    fig
	catch e
	    if e isa UndefVarError
	        # normal → skip
	    else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end


# ╔═╡ 0106756b-594d-4fdb-81b5-cf0739898521
cv_fig3, ax3 = plot_koper_fig3()

# ╔═╡ 4be2d50a-7289-44a5-a4a1-2f736c466f4f
cv_fig3

# ╔═╡ 8fc7877e-c4e4-40d1-a720-7806f7dbde0a
#=╠═╡
let
	try
	    fig = Figure(size = (1600, 900))
	       ax = Axis(fig[1, 1],
	        xlabel = L"φ (V vs SHE)",
	        ylabel = L"I (mA/cm²)",
	        yscale = log10,
	        yminorticksvisible = true,  
	        yminorticks = IntervalsBetween(5),
			#limits = ((-1.3, -0.7),(1e-4, 1e2))
	    )
	
	
	    # Experimental Data Plotting based on M.T.M Koper
	    raw = CSV.read("Langmuir 2021, 37, 5707−5716/Figure_5.csv", DataFrame; header=false)
	    pres = vec(Matrix(raw[1:1, :]))
	    sub = Matrix(raw[4:end, :])
	    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
	    num_df = DataFrame(num, :auto)
	    npairs = size(num_df, 2) ÷ 2
	    pink, pblue = RGB(0.0, 0.7, 0.8), RGB(0.2, 0.5, 0.0)
	    cols1 = [RGB(pink.r + t*(pblue.r-pink.r),
	                 pink.g + t*(pblue.g-pink.g),
	                 pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]
	
	    plot_objs1 = []
	    labels1 = String[]
	    for j in 1:npairs
	        xcol, ycol = 2j - 1, 2j
	        label = j == 1 ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
	        push!(labels1, label)
	        line = lines!(ax, num_df[!, xcol], (abs.(num_df[!, ycol])); color = cols1[j])
	        #line = lines!(ax, num_df[!, xcol], ((num_df[!, ycol])); color = cols1[j])
	        push!(plot_objs1, line)
	    end
	    Legend(fig[1, 2], plot_objs1, labels1, "Experimental"; framevisible = true)
	
	    # Theoretical Data Plotting based on `LiquidElectrolytes.jl`
	    plot_objs2 = []
	    labels2 = String[]
	    for (j, rec) in enumerate(F5_vec)
	        label2 = j == 1 ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
	        push!(labels2, label2)
	        line = lines!(ax, rec.voltages, (abs.(currents(rec, iohminus) .* cm^2/mA)); color = cols1[j])
			#line = lines!(ax, rec.voltages, ((currents(rec, iohminus) .* cm^2/mA)); color = cols2[j])
	        push!(plot_objs2, line)
	    end
	    Legend(fig[1, 3], plot_objs2, labels2, "Theoretical"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 58ac8edc-2432-4054-88d8-52dafe0a2a61
#=╠═╡
let
	try
		table = readdlm("../data/catmap_CO2R_data/IV-Ringe-digitized.csv", ',', Float64, '\n')
		raw = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
	
		x_exp_all = table[:, 1]
		y_exp_all = table[:, 2]
		
		mask_exp = (y_exp_all .> 0) .& isfinite.(y_exp_all) .& isfinite.(x_exp_all)
		x_exp = x_exp_all[mask_exp]
		y_exp = y_exp_all[mask_exp]
		
		volts_mask = ivresult.voltages .< -0.4
		x_iv  = ivresult.voltages[volts_mask]
		y_iv0 = abs.(currents(ivresult, iohminus))[volts_mask] .* cm^2/mA
		mask_iv = (y_iv0 .> 0) .& isfinite.(y_iv0) .& isfinite.(x_iv)
		x_iv = x_iv[mask_iv]; y_iv = y_iv0[mask_iv]
		
		fig = Figure(size = (900, 550))
		ax  = Axis(fig[1, 1];
				   xlabel = L"\phi_{we} \, (\mathrm{V \; vs \; SHE})",
				   ylabel = L"I \; (\mathrm{mA/cm^2})",
				   #yscale = log10,
				   #yminorticksvisible = true,
				   #yminorticks = IntervalsBetween(10),
				   #limits = ((-1.3, -0.4),(1e-12, 1e-7))
				  )
		
		#scatter!(ax, x_exp, y_exp; markersize=8, marker=:cross, color=:red, label="Ringe et al.")
		
		#lines!(ax, x_iv, y_iv; color=:green, label="e⁻, we")
	
		
		u = 7
		cols = [RGB(0.2 + 0.6*(i/length(sweep_vec)), 
					0.3 + 0.5*(1-i/length(sweep_vec)), 
					0.8 - 0.7*(i/length(sweep_vec)))
				for i in 1:length(sweep_vec)]
	    plot_objs = []
	    labels = String[]
	    for (j, rec) in enumerate(sweep_vec)
	        label = "$(scanrates[j])\t\t "
	        push!(labels, label)
	        line = lines!(ax, rec.voltages, ((currents(rec, iohminus) .* 
				cm^2/mA));linewidth = 1, color = cols[j], label = "$(scanrates[j])\t\t V/s")
	        push!(plot_objs, line)
	    end
		pres = vec(Matrix(raw[1:1, :])) 
		sub = Matrix(raw[4:end, :]) 
		num = map(x -> x === missing ? NaN : parse(Float64, x), sub) 
		num_df = DataFrame(num, :auto)
		npairs = size(num_df, 2) ÷ 2
		
		js    = u:npairs
		xcols = 2 .* js .- 1
		ycols = 2 .* js
		
		xs = [@view num_df[!, i] for i in xcols]
		ys = [@view num_df[!, i] for i in ycols]
		
		lines!.(Ref(ax), xs, [(y) for y in ys], color = RGB(0.0, 0.0, 0.7), linestyle = :dot, linewidth = 1, label = "CV Experimental")
	
		Legend(fig[1, 2], ax, "Legend"; framevisible=true)
		fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
  ╠═╡ =#

# ╔═╡ fbe4aca2-6a47-4457-98bb-588a5cde0ed5
md"""
#### Position at 0.99 Cbulk at a Given Time
"""

# ╔═╡ f918dc11-80e4-4223-8f02-3d5f0a10f8e5
function x99_at_time(tsol, ico2, t_index, x, c_bulk)
    c = tsol[ico2, :, t_index] ./ (mol/dm^3)

    idx = findfirst(ci -> ci ≥ 0.99c_bulk, c)

    return isnothing(idx) ? NaN : x[idx]
end

# ╔═╡ 1753c20f-9b53-4120-a8c8-e2b086f46f44
md"""
#### RDE_CV plots
"""

# ╔═╡ b1e64332-95a4-46a5-a45d-457c26e3fc67
const target_time = 20

# ╔═╡ 48029647-f162-459b-8824-fbf652d127f7
let
    try
        L_values = sort(collect(keys(RDE_result)))

        conc_vals = Float64[]

        for L in L_values
            rec = RDE_result[L]              # NamedTuple(δ, result)
            result_L = rec.result            # ✅ CVSweepResult

            times = result_L.tsol.t
            idx   = argmin(abs.(times .- target_time))

            c_ico2 = result_L.tsol.u[idx][ico2, 1] / (mol / dm^3)   # ✅ 안전하게 u로 접근
            push!(conc_vals, c_ico2)
        end

        L_labels = [ @sprintf("%0.4f", L) for L in L_values ]

        fig = Figure(size = (650, 400))
        ax = Axis(fig[1, 1];
            xlabel = L"L / \mu\mathrm{m}",
            ylabel = L"c_{\mathrm{CO_2}}(x=0)\ /(\mathrm{mol}/\mathrm{dm}^3)",
            title  = L"\mathrm{CO_2}\ \text{at electrode vs } L \text{ at } t \approx %$target_time",
            xticklabelrotation = π/4,
            xticks = (L_values, L_labels),
			xscale = log10
        )

        scatter!(ax, L_values, conc_vals; markersize = 8)
        fig

    catch e
        if e isa UndefVarError
            # normal → skip
        else
            println("⚠️ Error occurred: ", e)
            println(stacktrace(catch_backtrace()))
        end
    end
end


# ╔═╡ 4116166d-5f82-4d9b-80fb-c8035b9b6ade
let
    try 
		rpms = sort(collect(keys(RDE_result)))
	
	    δ_values        = Float64[] 
	    I_pos_peaks     = Float64[]  
	    I_neg_peaks     = Float64[]  
	    I_at_voltage    = Float64[]  
	    target_time = 40.0 #0.05 * 24 = -1.2V, mass transport limits currents     
	
	    for rpm in rpms
	        rec = RDE_result[rpm].result
	        δ   = RDE_result[rpm].δ        
	        δ_um = δ * 1e6                
	
	        push!(δ_values, δ_um)
	
	        # 전체 CV current
	        I = currents(rec, ico) .* (cm^2/mA)
	
	        # 1) positive / negative peak
	        push!(I_pos_peaks, maximum(I))
	        push!(I_neg_peaks, minimum(I))
	        times = rec.times            
	        idx   = argmin(abs.(times .- target_time))
	
	        push!(I_at_voltage, I[idx])      
	    end
	
	    δ_labels = [ @sprintf("%0.1f", δ) for δ in δ_values ]
	
	    fig = Figure(size = (1200, 800))
	    ax  = Axis(fig[1, 1],
	               xlabel = L"\text{L}\;(\mu m)",
	               ylabel = L"I_{peak} \; (mA/cm^2)",
	               title  = @sprintf("Peak current vs boundary layer thickness (t = %.1f s 포함)", target_time),
				   xticklabelrotation = π/4,
	               xticks = (δ_values, δ_labels),
	               xscale = log10
	              )
	
	    # peak currents
	    plot1 = scatterlines!(ax, δ_values, I_pos_peaks; marker = :circle)
	    plot2 = scatterlines!(ax, δ_values, I_neg_peaks; marker = :utriangle)
	
		plot3 = scatterlines!(ax, δ_values, I_at_voltage; marker = :diamond)
	
	   	Legend(fig[1, 2], [plot1, plot2, plot3], ["Positive", "Negative", "I(t)"])
	
	
	   	fig
	catch e
	   	if e isa UndefVarError
			# normal case → skip
	   	else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 9a92d4a9-f489-4bba-9361-03cad3ea12e1
md"""
#### Concentration plots
"""

# ╔═╡ 6a11b8e7-ed7f-4972-a4d5-d713e045ee1c
let
    fig = Figure()
    ax  = Axis(fig[1, 1],
        xlabel = L"\phi \, (\mathrm{V \; vs \; SHE})",
        ylabel = L"I \; (\mathrm{mA/cm^2})",
       # limits = ((-1.25, 0.8), (-0.5, 0.07))
    )

	csv_path = "../data/Langmuir_CV_data/Figure_5.csv"  
	unit_scale = cm^2/mA 

    raw = CSV.read(csv_path, DataFrame; header = false)
    pres_labels = vec(Matrix(raw[1:1, :])) 
    sub  = Matrix(raw[4:end, :])
    num  = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)

    npairs = size(num_df, 2) ÷ 2

	colors = [RGB(0.9 - (0.05 * i/npairs), 0.15 * (1 - i/npairs), 0.3 * i/npairs) for i in 1:npairs]


    plots = Plot[] 
    labels = String[]

    for j in 1:npairs
        xcol, ycol = 2j - 1, 2j
		label = "$(pres_labels[min(2j, length(pres_labels))]) pH"
        push!(labels, label)

		h = lines!(ax, num_df[!, xcol], log.(abs.(num_df[!, ycol])); color = colors[j])
        push!(plots, h)
    end
    Legend(fig[1, 2], plots,  labels,  "Experimental"; framevisible = true)

	fig
end

# ╔═╡ fe1e2a72-4772-4482-88da-f9e5f90e928a
md"""
#### Scan Rate CV plots
"""

# ╔═╡ d94ec33c-3d9d-4d70-b0e1-e3d861a62821
#=╠═╡
let
	try
	    fig = Figure(size = (1600, 900))
	    ax = Axis(fig[1, 1], ylabel = L"I (mA/cm²)", xlabel = L"φ (V vs SHE)")
		
		cols = [RGB(0.2 + 0.6*(i/length(sweep_vec)), 
					0.3 + 0.5*(1-i/length(sweep_vec)), 
					0.8 - 0.7*(i/length(sweep_vec)))
				for i in 1:length(sweep_vec)]
	    plot_objs = []
	    labels = String[]
	    for (j, rec) in enumerate(sweep_vec)
   		label = "$(scanrates[j])\t\t "
    	push!(labels, label)
	    color_j = cols[j]
	    lw_j = 1
	
	    if j == length(sweep_vec) - 4
	        color_j = RGB(1, 0.2, 0.2)   
	        lw_j = 4                       
	    end
	
	    line = lines!(ax,
	                  rec.voltages,
	                  (currents(rec, iohminus) .* cm^2/mA);
	                  linewidth = lw_j,
	                  color = color_j)
	
	    push!(plot_objs, line)
		end
	    Legend(fig[1, 2], plot_objs, labels, "Scan Rates (V/s)"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 333492ec-9016-44c5-9059-e3cb42c05a89
#=╠═╡
let
	try
	    fig = Figure(size = (1600, 900))
	    ax = Axis(fig[1, 1], ylabel = L"I (mA/cm²)", xlabel = L"φ (V vs SHE)", 
				  #limits = ((-0.8, 1.5),(-2e-8, 10))), 
				  yscale = log10
				 )
		
		cols = [RGB(0.2 + 0.6*(i/length(sweep_vec)), 
					0.3 + 0.5*(1-i/length(sweep_vec)), 
					0.8 - 0.7*(i/length(sweep_vec)))
				for i in 1:length(sweep_vec)]
	    plot_objs = []
	    labels = String[]
	    for (j, rec) in enumerate(sweep_vec)
   		label = "$(scanrates[j])\t\t "
    	push!(labels, label)
	    color_j = cols[j]
	    lw_j = 1
	
	    if j == length(sweep_vec) - 4
	        color_j = RGB(1, 0.2, 0.2)   
	        lw_j = 4                       
	    end
	
	    line = lines!(ax,
	                  rec.voltages,
	                  abs.(currents(rec, iohminus) .* cm^2/mA);
	                  linewidth = lw_j,
	                  color = color_j)
	
	    push!(plot_objs, line)
		end
	    Legend(fig[1, 2], plot_objs, labels, "Scan Rates (V/s)"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
  ╠═╡ =#

# ╔═╡ 9e2b6a47-a113-4107-acdf-3901e0578898
#=╠═╡
let
    try
        target_voltage = -1.2

        fig = Figure(size = (1600, 900))
        ax = Axis(fig[1, 1],
                  ylabel = L"I \; (mA/cm^2)",
                  xlabel = L"v \; (V/s)",
                  title  = L"I(\varphi = %$(target_voltage) \text{ V vs SHE}) \text{ vs. scan rate}"
        )

        I_at_V = Float64[]
        cols = [RGB(0.2 + 0.6*(i/length(sweep_vec)),
                    0.3 + 0.5*(1-i/length(sweep_vec)),
                    0.8 - 0.7*(i/length(sweep_vec)))
                for i in 1:length(sweep_vec)]

        labels    = String[]
        plot_objs = Any[]

        for (j, rec) in enumerate(sweep_vec)
            idx = argmin(((rec.voltages .- target_voltage)))

            Ival = currents(rec, ico)[idx] .* (cm^2/mA)
            push!(I_at_V, Ival)

            color_j = cols[j]
            lw_j    = 2
            marker_j = :circle

            if j == length(sweep_vec) - 4
                color_j  = RGB(1, 0.2, 0.2)  
                lw_j     = 4
                marker_j = :utriangle
            end

            line = scatter!(ax,
                            ([scanrates[j]]),
                            ([Ival]);
                            markersize = 10,
                            marker     = marker_j,
                            color      = color_j)

            push!(plot_objs, line)
            push!(labels, "$(scanrates[j])")
        end

        Legend(fig[1, 2], plot_objs, labels,
               "Scan Rates (V/s)";
               framevisible = true)

        fig
    catch e
        if e isa UndefVarError
            # normal case → skip
        else
            println("⚠️ Error occurred: ", e)
            println(stacktrace(catch_backtrace()))
        end
    end
end

  ╠═╡ =#

# ╔═╡ be26b92a-14e2-45bc-bb6f-a2664e2e3cd9
md"""
#### Electrode concentration_Boundary Thickness Layer Function
"""

# ╔═╡ 61be3485-960b-42f8-82e6-71e213a5c9a1
let
	try
	    δ_keys = sort(collect(keys(Lresult)))
	    δ_um   = Float64.(δ_keys)
	
	    I_pos_peaks  = Float64[]
	    I_neg_peaks  = Float64[]
	    I_at_voltage = Float64[]
	
	    target_time = 85.0
	
	    rec_ref = Lresult[δ_keys[1]]
	    idx_ref = argmin(abs.(rec_ref.times .- target_time))
	    V_target = rec_ref.voltages[idx_ref]   
	
	    for δ in δ_keys
	        rec = Lresult[δ]
	
	        I = currents(rec, ico) .* (cm^2/mA)
	
	        push!(I_pos_peaks, maximum(I))
	        push!(I_neg_peaks, minimum(I))
	
	        times = rec.times
	        idx   = argmin(abs.(times .- target_time))
	        push!(I_at_voltage, I[idx])
	    end
	
	    δ_labels = [ @sprintf("%0.1f", δ) for δ in δ_um ]
	
	    fig = Figure(size = (1200, 800))
	    ax  = Axis(fig[1, 1],
	        xlabel = L"\text{L}\;(\mu\mathrm{m})",
	        ylabel = L"I_{\text{peak}} \; (\mathrm{mA}/\mathrm{cm}^2)",
	        title  = @sprintf("Peak current vs boundary layer thickness (V = %.2f V, t = %.1f s) [Scan Rate = 50 mV/s]",
	                          round(V_target, digits=2), target_time),
	        xticklabelrotation = π/4,
	        xticks = (δ_um, δ_labels),
	    )
	
	    plot1 = scatterlines!(ax, δ_um, I_pos_peaks;  marker = :circle)
	    plot2 = scatterlines!(ax, δ_um, I_neg_peaks;  marker = :utriangle)
	    plot3 = scatterlines!(ax, δ_um, I_at_voltage; marker = :diamond)
	
	    Legend(fig[1, 2],
	           [plot1, plot2, plot3],
	           ["Positive peak", "Negative peak",
	            @sprintf("I at V = %.2f V", round(V_target, digits=2))])
	
	    fig
	catch e
	    if e isa UndefVarError
	        # normal → skip
	    else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
		
end

# ╔═╡ def960de-f74a-4ca8-9d95-8af4e0240b60
let
	try
    fig = Figure(size = (1600, 900))
    ax = Axis(fig[1, 1], 
			  ylabel = L"\log|I| (mA/cm²)", 
			  xlabel = L"φ (V vs SHE)",
			  #limits = ((-2, -1.4),(10e-1, 1.0e1)),
			  #yscale = log10
			 )	

    keys_sorted = sort(collect(keys(Lresult)))
	cols2 = [RGB(0.3 + 0.7*(i/length(Lresult)), 
				0.1 + 0.7*(1-i/length(Lresult)), 
				0.9 - 0.6*(i/length(Lresult)))
				for i in 1:length(Lresult)]   
	plot_objs2 = []
    labels2 = String[]

    for (j, L) in enumerate(keys_sorted)
        rec = Lresult[L] 
        line = lines!(ax, rec.voltages, (currents(rec, ico) .* cm^2/mA);
                      color = cols2[j])
        push!(plot_objs2, line)
        push!(labels2, "L = $(L) μm")
    end
	#vspan!(ax,  -1.2,  -1.19, color=(colorant"#000000"))
    #Legend(fig[1, 2], plot_objs2, labels2, "Theoretical"; framevisible = true)
		Legend(fig[1, 2], plot_objs2, labels2, "Theoretical";
		    framevisible = true,
		    labelsize = 12,      
		    titlesize = 13,      
		    patchsize = (15, 5), 
		    rowgap = 1, colgap = 1,
		    patchlineattrs = (linewidth = 30,) 
		)

    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 89520d6a-7a44-41f6-92ba-3d9416ac2047
md"""
#### Constant time_Boundary Thickness Layer Function
"""

# ╔═╡ d3493ce8-85d1-4132-b3e0-4ec35ac9d36d
md"""
#### Target-time CV function
"""

# ╔═╡ f90190cc-555d-47e1-a2cb-99e5d78d4ff5
let
    try
        fig = Figure(size = (1600, 900))
        ax = Axis(fig[1, 1], 
                  ylabel = L"\log|I| (mA/cm²)", 
                  xlabel = L"φ (V vs SHE)",
                  #limits = ((-2, -1.4),(10e-1, 1.0e1)),
                  #yscale = log10
        )	

        keys_sorted = sort(collect(keys(Lresult)))
        cols2 = [RGB(0.3 + 0.7*(i/length(Lresult)), 
                     0.1 + 0.7*(1-i/length(Lresult)), 
                     0.9 - 0.6*(i/length(Lresult)))
                 for i in 1:length(Lresult)]   

        plot_objs2 = Vector{Any}()
        labels2    = String[]

        for (j, L) in enumerate(keys_sorted)
            rec = Lresult[L]

            # 전체 CV 곡선
            I = currents(rec, ico) .* (cm^2/mA)
            φ = rec.voltages

            line = lines!(ax, φ, I; color = cols2[j])
            push!(plot_objs2, line)
            push!(labels2, "L = $(L) μm")

            times = rec.times
            idx   = argmin(abs.(times .- target_time))

            φ_t = φ[idx]
            I_t = I[idx]

            scatter!(ax, [φ_t], [I_t];
                     markersize = 15,
                     marker = :circle,
                     color = cols2[j])
        end

        Legend(fig[1, 2], plot_objs2, labels2, "Theoretical";
            framevisible   = true,
            labelsize      = 12,      
            titlesize      = 13,      
            patchsize      = (15, 5), 
            rowgap         = 1, 
            colgap         = 1,
            patchlineattrs = (linewidth = 30,) 
        )

        fig
    catch e
        if e isa UndefVarError
            # normal case → skip
        else
            println("⚠️ Error occurred: ", e)
            println(stacktrace(catch_backtrace()))
        end
    end
end


# ╔═╡ dadf76f0-cbea-4c34-a142-41e120679674
function voltage_at_time(result::CVSweepResult, t_input)
    idx = argmin(abs.(result.times .- t_input))  
    return result.voltages[idx]
end

# ╔═╡ 91242a8c-c09b-402c-a0ea-40b8e3e26ae7
if @isdefined Lresult
    v = voltage_at_time(Lresult[2108], 76.0)
else
    @info "Lresult undefined — skipping"
end

# ╔═╡ bb00b5bb-326e-47f9-a4f4-e7b4f29dd1f2
md"""
### Plotting Functions
"""

# ╔═╡ b7cb5183-65e8-4ee8-af86-2bedd11daecc
function CVPlot!(result, model)
    ic = model.cspecies
    fig = Figure(size = (650, 400))
    ax = Axis(fig[1, 1], 
 			 # limits = ((-1.25, 0.8),(-0.0000002, 0.0000002)),
              ylabel = L"I (mA/cm²)",
              xlabel = L"φ (V vs SHE)",
			  #yscale = log10
			 )
	
    total_current = zero(currents(result, ic[1]))
    #for s in ic
    #    total_current .+= (currents(result, s) * mA / cm^2)
    #end
	total_current = (currents(result, ico)) .* cm^2/mA

	
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

# ╔═╡ b976ab43-69f1-47a0-b2c6-c63e1c15cdb4
md"""
Boundary Layer thickness varied polarization curve: $(@bind Ldependancy PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 7a02463d-cfd9-4648-af53-f1e65d46733f
md"""
### Result Plots
"""

# ╔═╡ 114d2324-5289-4e44-8d77-736a9bdec365
md"""
Show only pH: $(@bind useonly_pH PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 180c12b0-d410-4a5f-97bb-226e6624a39b
catmap_params

# ╔═╡ f8b5dc8f-1f41-4600-825e-2f9653f2d925
md"""
### **Polarization Curve**
"""

# ╔═╡ f672a256-641a-478e-b0aa-2df6e68b4d86
md"
#### Boundary Layer Thickness varied Polarization Curve
"

# ╔═╡ 904ac4c2-50a8-4f70-8050-a0a1d4a448fa
md"
### Voltage Concentration Plots
"

# ╔═╡ c1d2305e-fb8b-4845-a414-08fff84aa9b0
md"""
### Plotting Functions
"""

# ╔═╡ a81dd9a4-7938-4a72-b3d2-1780e8ecd536
function plotcurr_over_L(results::Dict{Int,Any};
    species=ico, cutoff=-0.4, title="IV vs L"
)
    vis = GridVisualizer(;
        size   = (800, 500),
        title  = title,
        xlabel = L"\phi_{we}\;(\mathrm{V\;vs\;SHE})",
        ylabel = L"I\;(\mathrm{mA/cm^2})",
        legend = :rt,
        yscale = :log,
    )

    items = sort(collect(results); by=first)
    n = length(items)

    cols = Makie.resample_cmap(:cool, n)   
    for (k, (L, ivres)) in enumerate(items)
        volts = ivres.voltages
        mask  = volts .< cutoff

        scalarplot!(vis,
            volts[mask],
            abs.(currents(ivres, species))[mask] .* cm^2/mA;
            clear = false,
            label = "L = $(L) μm",
            color = cols[k],              
        )
    end

    reveal(vis)
end

# ╔═╡ d5ab1a28-3a60-49d9-bb3e-ca589b1c79fd
begin
	curr(J, ix) = [F * abs(j[ix]) for j in J]
	
	function plotcurr(result; df = nothing)
	    scale = 1 / (mol / dm^3)
	    volts = result.voltages[result.voltages .< -0.4]
	    vis = GridVisualizer(;
	                         size = (600, 400),
	                         tilte = "IV Curve",
	                         xlabel = L"\phi_{we} \, (\mathrm{V \; vs \; SHE})",
	                         ylabel = L"I / (\mathrm{mA/cm^2})",       
	                         legend = :lb,
							 yscale = :log,
		)
							 
	    scalarplot!(vis,
	                volts,
	                #abs.(currents(ivresult, iohminus))[result.voltages .< -0.4] .* cm^2/mA;
					curr(result.j_we, ico)[result.voltages .< -0.4] .* cm^2/mA;
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

# ╔═╡ d377af90-9b7d-4fd9-8bc8-d93e58617542
function iv_curve_axis(ivresult;
    cutoff = -0.4,
    showlegend = true,
    #title = "IV Curve",
    iohminus = iohminus,
    data_dir = "./catmap_CO2R_data",
)

    # ---- data from simulation ----
    v_all = ivresult.voltages
    mask  = v_all .< cutoff
    volts = v_all[mask]

    # currents(ivresult, iohminus) 
    I_sim = abs.(currents(ivresult, iohminus))[mask] .* (cm^2/mA)
    #I_sim = max.(I_sim, eps(Float64))

    # ---- load csvs ----
    table  = readdlm(joinpath(data_dir, "IV-Ringe-digitized.csv"), ',', Float64, '\n')
    df_v   = table[:, 1]
    df_I   = abs.(table[:, 2])
    df_I   = max.(df_I, eps(Float64))

    table2 = readdlm(joinpath(data_dir, "Ringe-theorical.csv"), ',', Float64, '\n')
    df2_v  = table2[:, 1]
    df2_I  = abs.(table2[:, 2])
    df2_I  = max.(df2_I, eps(Float64))

    table3 = readdlm(joinpath(data_dir, "Ringe-experimental.csv"), ',', Float64, '\n')
    df3_v  = table3[:, 1]
    df3_I  = abs.(table3[:, 2])
    df3_I  = max.(df3_I, eps(Float64))

    # ---- FIGURE STYLE (match conc_vs_voltage_axis) ----
    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ \mathrm{(V)}}",
        ylabel = L"\mathbf{I}\;(\mathrm{mA/cm^2})",
        yscale = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e2)),
    )

	xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    ax.spinewidth = 5.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25

    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xlabelpadding = 10
    ax.ylabelpadding = 10
    ax.xlabelfont = :bold

    # ---- PLOTS ----
    # simulation line
    lines!(ax, volts, I_sim; color=:green, linewidth=5, label="e⁻, we")


	"""	
    # digitized (cross)
    scatter!(ax, df_v, df_I;
        marker = :xcross,
        markersize = 12,
        color = :red,
        label = "Ringe et. al"
    )
	"""

    # theoretical (circle)
    scatter!(ax, df2_v, df2_I;
        marker = :circle,
        markersize = 8,
        color = :blue,
        label = "Ringe et. al : Theoretical"
    )

    # experimental (triangle up)
    scatter!(ax, df3_v, df3_I;
        marker = :utriangle,
        markersize = 12,
        color = :magenta,
        label = "Ringe et. al : Experimental"
    )

    showlegend && axislegend(ax, position=:rt)

    return fig
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
TableOfContents(title="📚 Table of Contents", indent=true, depth=5, aside=true)

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

# ╔═╡ 6a9fad5b-4964-4e12-b131-8cb2628d1ab3
floataside(
    @bind user_input_cv confirm(
        PlutoUI.combine() do Child
			md"""
			###### __User Data__  
			``CO_2 + H_2O + 2e^- \leftrightharpoons CO + 2OH^-``
					
			- ``V_\mathrm{min}``: $(Child("vmin", NumberField(-2.0:0.1:2.0; default = -1.2)))  ``V_\mathrm{max}``: $(Child("vmax", NumberField(-2.0:0.1:2.0; default = 0.8)))
			- ``z_R``: $(Child("zR", NumberField(-2:2; default = -1)))   
			  ``n``: $(Child("n", NumberField(0:2; default = 2)))
			
			- Scan rate ``(V/s)``: $(Child("scanrate", TextField(6; default = "0.05")))  
			- Periods: $(Child("nperiods", NumberField(1:10; default = 1)))
			- `Double64`: $(Child("double64", CheckBox()))  
			- `tunnel`: $(Child("tunnel", CheckBox()))
			- `scanup`: $(Child("scanup", CheckBox()))
			"""
		end;
	        label = "Submit"
	);
	top = 50
)

# ╔═╡ b4aaf070-d4ab-409a-b1e8-f5469b9f398b
begin 
	sawtooth = SawTooth(
        scanrate = parse(Float64, user_input_cv.scanrate),
        vmin = user_input_cv.vmin , vmax = user_input_cv.vmax,
		scanup = user_input_cv.scanup
    )
	    const nperiods = user_input_cv.nperiods
end

# ╔═╡ e1ef1e83-c472-4267-8450-38c65f48d3dc
let
    try
        fig = Figure(size = (1920, 1080))
        ax = Axis(fig[1, 1],
            ylabel = L"I (mA/cm^2)",
            xlabel = L"\phi\ (V\ vs\ SHE)",
            title  = @sprintf("pH = 9 | v=%.1f mV/s", sawtooth.scanrate * 1e3),
            titlesize = 30,
        )

        rpms = sort(collect(keys(RDE_result)))

        cols = [RGB(0.3 + 0.7*(i/length(rpms)),
                    0.1 + 0.7*(1 - i/length(rpms)),
                    0.9 - 0.6*(i/length(rpms)))
                for i in 1:length(rpms)]

        plot_objs = Any[]
        labels    = String[]

        for (j, rpm) in enumerate(rpms)
            recNT = RDE_result[rpm]
            rec   = recNT.result
            δ     = recNT.δ

            ϕ = rec.voltages
            I = currents(rec, ico) .* (cm^2/mA)

            line = lines!(ax, ϕ, I; color = cols[j], linewidth=2)
            push!(plot_objs, line)

            # --- target_time에 해당하는 점 찍기 ---
            t = rec.tsol.t
            n = min(length(t), length(ϕ), length(I))   # 길이 mismatch 방어
            idx = argmin(abs.(t[1:n] .- target_time))

            scatter!(ax, [ϕ[idx]], [I[idx]];
                     markersize = 12,
                     marker = :circle,
                     color = cols[j])

            # --- legend label ---
            rpm_str = @sprintf("%7.4f", rpm)
            δ_um    = δ * 1e6
            δ_str   = @sprintf("%9.6f", δ_um)
            push!(labels, "rpm=$(rpm_str), δ=$(δ_str) μm")
        end

        Legend(fig[1, 2], plot_objs, labels, "Theoretical";
            framevisible = true,
            labelsize = 12,
            titlesize = 13,
            patchsize = (15, 5),
            rowgap = 1, colgap = 1,
            patchlineattrs = (linewidth = 6,)
        )

        fig
    catch e
        if e isa UndefVarError
            # normal → skip
        else
            println("⚠️ Error occurred: ", e)
            println(stacktrace(catch_backtrace()))
        end
    end
end


# ╔═╡ d75725cd-0ef6-421f-be56-f559312e73b6
floataside(
    @bind user_input_ion confirm(
        PlutoUI.combine() do Child
            md"""
            ###### __Parameter Set__			
            - __Other ions__  
              ``κ``: $(Child("κt", NumberField(0.0:0.1:20.0; default = 0.0)))  
              ``a``: $(Child("at", NumberField(0.0:0.1:20.0; default = 0.0)))
            - __Cation__  
              ``κ``: $(Child("κk", NumberField(0.0:0.1:20.0; default = 0.0)))  
              ``a``: $(Child("ak", NumberField(0.0:0.1:20.0; default = 8.2)))
            """
        end;
        label = "Submit"
    );
    top = 350
)

# ╔═╡ ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
begin
	at = user_input_ion[:at]
	κt = user_input_ion[:κt]
	ak = user_input_ion[:ak]
	κk = user_input_ion[:κk]
	const bulk = let 
		bulk = [
				BulkSpecies(;name = "HCO₃⁻", 
							z = -1, 
							D = 1.185e-9, 
							c_bulk = 0.091, 
							#v = v0*(κt+1),
							a = at, #15.6,#at,
							κ = κt, #4, #κt, 
							color = :brown
				),
				BulkSpecies(;name = "CO₃²⁻",
							z = -2, 
							D = 0.923e-9, 
							c_bulk = 2.68e-6,
							#v = v0*(κt+1), 
							a = at, #17.8,
							κ = κt, #7, #κt, 
							color = :violet
				),
				BulkSpecies(;name = "CO₂",
							z = 0, 
							D = 1.91e-9, 
							# = 0.033
							c_bulk = 0.033, 
							#v =v0, 
							a = at, #17,#at,
							κ = κt, #0, #κt, 
							color=:red
				),
				BulkSpecies(;name = "OH⁻",
							z = -1, 
							D = 5.273e-9, 
							c_bulk = 10^(pH-14), 
							#v = v0*(κt+1), 
							a = at, #13.3,#at,
							κ = κt, #3, #κt, 
							color = :green
				),
				BulkSpecies(;name = "H⁺", 
							z = 1, 
							D = 9.310e-9, 
							c_bulk = 10^(-pH), 
							#v = v0*(κt+1), 
							a = at, #10, #at,
							κ = κt, #4, #κt, 
							color = :gray
				),
				BulkSpecies(;name="CO",
							z = 0,
							D = 2.23e-9,
							c_bulk = 0.0,
							#v = v0, 
							a = at, #140,
							κ = κt, #0, # κt,  
							color=:blue
				)
		]
		push!(bulk, make_eneutral(bulk;name="K⁺", 
									   z = 1, 
									   D = 1.957e-9,
									   a = ak, # 13.3, #ak,
									   #v = v0*(κt+1), 
									   κ = κk, #4, 
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
	
    electrolytedata(sys).κ .= κt 
end

# ╔═╡ 59f05654-a8de-4e17-b2ad-60a7ac64e122
let
	try
	    L_values = sort(collect(keys(RDE_result)))
	    species = getproperty.(bulk, :name)
	
	    conc_vals = Float64[]
	
	    for L in L_values
	        result_L = RDE_result[L]
	
	        times = result_L.tsol.t
	        idx = argmin(abs.(times .- target_time))
	
	        c_ico2 = result_L.tsol[ico2, 1, idx] / (mol / dm^3)
	        push!(conc_vals, c_ico2)
	    end
	
	    L_labels = [ @sprintf("%0.1f", L) for L in L_values ]
	
	    fig = Figure(size = (650, 400))
	    ax = Axis(fig[1, 1];
	        xlabel = L"L / \mu\mathrm{m}",
	        ylabel = L"c_{\mathrm{CO_2}}(x=0) / (\mathrm{mol}/\mathrm{dm}^3)",
	        title  = L"\mathrm{CO_2}\ \text{concentration at electrode vs. } L \text{ at } t = %$target_time",
	        xticklabelrotation = π/4,
	        xticks = (L_values, L_labels),
	    )
	
	    scatter!(ax, L_values, conc_vals; markersize = 8, color = :green)
	
	    fig
	catch e
		if e isa UndefVarError
	        # normal → skip
	    else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end


# ╔═╡ 13645622-ed9d-4843-9b83-e67313473c35
function plot_pnp_summary(res;
    ispec::Int,
    icurr_t::Int = ispec,     
    icurr_cv::Int = icurr_t,   
    ix_e::Int = 1,
    bulk = bulk,
    fig_size = (950, 1500),
)
	_get_name(bulk, i)  = string(getproperty(bulk[i], :name))
	_get_color(bulk, i) = getproperty(bulk[i], :color)

	
    # ---- shared pulls (single source) ----
    pnp = res.pnpresult
    t   = pnp.tsol.t
    u   = pnp.tsol.u
    ϕ   = pnp.voltages

    name  = _get_name(bulk, ispec)
    col   = _get_color(bulk, ispec)

    # ---- data extraction ----
    It = currents(pnp, icurr_t) .* (cm^2/mA)
    n1 = min(length(t), length(It))
    t1 = t[1:n1]
    It = It[1:n1]

    Icv = currents(pnp, icurr_t) .* (cm^2/mA)
    n2  = min(length(ϕ), length(Icv))
    ϕ2  = ϕ[1:n2]
    Icv = Icv[1:n2]

    ce = [u[k][ispec, ix_e] for k in eachindex(t)] ./ (mol/dm^3)

    # ---- plotting (3 axes in one figure) ----
    fig = Figure(size = fig_size)

    ax1 = Axis(fig[1, 1],
        xlabel = "t / s",
        ylabel = L"I\ (mA/cm^2)",
        title  = "I(t)  |  icurr=$(icurr_t)  |  $(name)"
    )
    lines!(ax1, t1, It; linewidth=2, color=col)

    ax2 = Axis(fig[2, 1],
        xlabel = "t / s",
        ylabel = L"c(x=0)\ (mol/dm^3)",
        title  = "c at electrode  |  ispec=$(ispec), ix_e=$(ix_e)  |  $(name)"
    )
    lines!(ax2, t, ce; linewidth=2, color=col)

	ax3 = Axis(fig[3, 1],
        xlabel = L"\phi\ (V\ vs\ SHE)",
        ylabel = L"I\ (mA/cm^2)",
        title  = "CV  |  icurr=$(icurr_cv)  |  $(name)"
    )
    lines!(ax3, ϕ2, Icv; linewidth=2, color=col)

    fig
end

# ╔═╡ 323abbb8-6f34-4b8b-839c-e0682fed1971
plot_pnp_summary(res; ispec = iohminus)

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
		#limits = ((0, 24),(1e-30, 1e2)),
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


# ╔═╡ 4e894347-2ce6-4c5f-a06e-7f1af1983bbc
begin
	try
		conc_time_func(res.pnpresult, sawtooth)
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ b14d67ca-5f24-4d8a-9334-6e072e2b39eb
function conc_time_vs_L(Lresult, ico; nspecies = 7)
    L_values = sort(collect(keys(Lresult)))

    first_result = Lresult[L_values[1]]
    times = first_result.tsol.t
    nt = length(times)

    species = getproperty.(bulk, :name)
    colors =  resample_cmap(:winter, length(Lresult))

    fig = Figure(size = (1290, 960))
    ax = Axis(
        fig[1, 1];
        xlabel = L"time / s",
        ylabel = L"c_{i,\,\text{electrode}} / (mol/dm^3)",
        yscale = log10,
        title  = L"c_{\mathrm{%$(species[ico])}} \text{ at electrode for different } L"
    )

    plot_objs = AbstractPlot[]
    labels    = String[]

    for (j, L) in enumerate(L_values)
        result_L = Lresult[L]

        conc_at_electrode_ico = [
            result_L.tsol[ico, 1, t] / (mol / dm^3) for t in 1:nt
        ]

        yvals_safe = [c > 0 ? c : NaN for c in conc_at_electrode_ico]

        line = lines!(ax, times, yvals_safe;
                      color = colors[j],
					  linewidth = 1
        )

        push!(plot_objs, line)
        push!(labels, @sprintf("L = %.1f μm", L ))
    end

    Legend(fig[1, 2], plot_objs, labels;
           labelsize = 10,
           backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5),
           title = "Boundary-layer thickness L")

    fig
end


# ╔═╡ 3b41341e-d174-4b2f-8c19-068cb84ba571
if @isdefined Lresult
    conc_time_vs_L(Lresult, ico; nspecies = 7)
else
    @info "Lresult undefined — skipping"
end

# ╔═╡ deb15672-0e35-4855-b1c0-b2c0b7e78d41
function conc_time_vs_RPM(RDE_result, ico; nspecies = 7)
    rpms = sort(collect(keys(RDE_result)))

    species = getproperty.(bulk, :name)
    colors  = resample_cmap(:winter, length(rpms))

    fig = Figure(size = (1290, 960))
    ax  = Axis(
        fig[1, 1];
        xlabel = L"time / s",
        ylabel = L"c_{i,\,\text{electrode}} / (mol/dm^3)",
        yscale = log10,
        title  = L"c_{\mathrm{%$(species[ico])}} \text{ at electrode for different RDE rotation rates}"
    )

    plot_objs = AbstractPlot[]
    labels    = String[]

    for (j, rpm) in enumerate(rpms)
        entry = RDE_result[rpm]     
        rec   = entry.result        
        δ     = entry.δ              
        δ_um  = δ * 1e6

        times = rec.tsol.t
        nt    = length(times)

        conc_at_electrode_ico = [
            rec.tsol[ico, 1, t] / (mol / dm^3) for t in 1:nt
        ]

        yvals_safe = [c > 0 ? c : NaN for c in conc_at_electrode_ico]

        line = lines!(
            ax, times, yvals_safe;
            color     = colors[j],
            linewidth = 1,
        )

        push!(plot_objs, line)
        push!(labels, @sprintf("ω = %4d rpm  (L = %.1f μm)", rpm, δ_um))
    end

    Legend(
        fig[1, 2], plot_objs, labels;
        labelsize       = 10,
        backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5),
        title           = "RDE rotation & boundary-layer thickness"
    )

    fig
end


# ╔═╡ 3f50a881-337f-4578-b0ff-a143439b0a6d
if @isdefined RDE_result
	conc_time_vs_RPM(RDE_result, ico; nspecies = 7)
else
    @info "RDE_result undefined — skipping"
end

# ╔═╡ b48b4acb-ed25-4d9b-bee6-2e316ecb44c3
function conc_x_vs_L_at_time(Lresult, ico; target_time = 20.0)

	L_values = sort(collect(keys(Lresult)))

    first_result = Lresult[L_values[1]]
    times = first_result.tsol.t
    nt    = length(times)

    t_idx = argmin(abs.(times .- target_time))

    species = getproperty.(bulk, :name)
    colors  = resample_cmap(:autumn1, length(L_values))

    fig = Figure(size = (1290, 960))
    ax  = Axis(
        fig[1, 1];
        xlabel = L"x \; (\mu m)",
        ylabel = L"c_{i}(x, t^*) \; (mol/dm^3)",
        title  = L"c_{\mathrm{%$(species[ico])}}(x, t^*) \text{ for different } L \; \text{(}t^* = %$(round(target_time, digits=2)) \text{ s)}",
		
    )
	xlims!(ax, -1, 540)
    plot_objs = AbstractPlot[]
    labels    = String[]
	
    for (j, L) in enumerate(L_values)
        result_L = Lresult[L]

        nx = size(result_L.tsol, 2)

        xvals_um  = range(0, L; length = nx)

        conc_profile = [
            result_L.tsol[ico, ix, t_idx] / (mol / dm^3) for ix in 1:nx
        ]

        line = lines!(
            ax, xvals_um, conc_profile;
            color     = colors[j],
            linewidth = 1,
        )

        push!(plot_objs, line)
        push!(labels, @sprintf("L = %.1f μm", L))
    end

    Legend(
        fig[1, 2], plot_objs, labels;
        labelsize       = 10,
        backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5),
        title           = @sprintf("Profiles at t = %.1f s", target_time)
    )

    fig
end


# ╔═╡ 666c55e5-f7f5-4f83-b3a3-ea6632ca5a86
if @isdefined Lresult
    conc_x_vs_L_at_time(Lresult, ico; target_time = 40.0)
else
    @info "Lresult undefined — skipping"
end

# ╔═╡ 0a665fc0-1230-4978-8ebd-e551c595e857
function conc_x_vs_RPM_at_time(RDE_result, ico; target_time = 24.0)
    rpms = sort(collect(keys(RDE_result)))

    first_rec = RDE_result[rpms[1]].result
    times = first_rec.tsol.t
    t_idx = argmin(abs.(times .- target_time))

    species = getproperty.(bulk, :name)
    colors  = resample_cmap(:autumn1, length(rpms))

    fig = Figure(size = (1290, 960))
    ax  = Axis(
        fig[1, 1];
        xlabel = L"x \; (\mu m)",
        ylabel = L"c_{i}(x, t^*) \; (mol/dm^3)",
        title  = L"c_{\mathrm{%$(species[ico])}}(x, t^*) \text{ for different rpm}"
    )
	xlims!(ax, -1, 100)

    plot_objs = AbstractPlot[]
    labels    = String[]

    for (j, rpm) in enumerate(rpms)
        rec = RDE_result[rpm].result
        δ   = RDE_result[rpm].δ        
        nx  = size(rec.tsol, 2)

        xvals_m  = range(0, δ; length = nx)
        xvals_um = collect(xvals_m .* 1e6)

        conc_profile = [
            rec.tsol[ico, ix, t_idx] / (mol / dm^3) for ix in 1:nx
        ]

        line = lines!(
            ax, xvals_um, conc_profile;
            color = colors[j],
            linewidth = 1,
        )

        push!(plot_objs, line)
        push!(labels, @sprintf("ω = %4d rpm (L = %.1f μm)", rpm, δ * 1e6))
    end

    Legend(
        fig[1, 2], plot_objs, labels;
        labelsize       = 10,
        backgroundcolor = RGBA(1,1,1,0.5),
        title           = @sprintf("Profiles at t = %.1f s", target_time)
    )

    fig
end


# ╔═╡ 960ee06e-b15f-4890-912c-851691e5c1a7
begin
    species = getproperty.(bulk, :name)
    colors  = getproperty.(bulk, :color)
    nspecies = length(species)

    # --- catmap CSV ---
    df = CSV.read("catmap_CO2R_data/voltage-conc.csv", DataFrame)
    rename!(df, Dict(names(df)[1]=>:Index, names(df)[2]=>:Voltage, names(df)[3]=>:Concentration))

    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = L"\text{Voltage}\ U\ \mathrm{vs.}\ \text{SHE}\ (V)",
        ylabel = L"c_i^{+}\;(\mathrm{M})",
        yscale  = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    ax.spinewidth = 2.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 30
    ax.ylabelsize = 30
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20

    for g in groupby(df, :Index)
        idx = Int(first(g.Index))
        p = sortperm(g.Voltage)

        V = g.Voltage[p]
        y = max.(g.Concentration[p], eps(Float64))  

        if 1 <= idx <= nspecies
            lines!(ax, V, y; color=colors[idx], linewidth=3, label=species[idx])
        else
            lines!(ax, V, y; linewidth=3, label="Index = $idx")
        end
    end

    axislegend(ax; position=:rt) 
    fig
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
			for (j_result, j_sresult) in zip(result.j_we[vidx_result][ico], 
											sresult.j_we[vidx_result][ico])
				@test isapprox(j_result, j_sresult, rtol=1.0e-5)
			end
		end
	end
	end
end;

# ╔═╡ ab0e28f4-4310-4dcc-817e-81b9e45fd501
floataside(
    @bind user_input_model confirm(
        PlutoUI.combine() do Child
			md"""
			###### __Model Selection__  
			- Model: $(Child("model_choice", Select(["Gold_Model", "Landstorfer_NaClO₄ model", "Landstorfer_NaF model", "Toy model"])))
			---
			###### __Boundary layer thickness__   
			- ``δ``: $(Child("L", NumberField(1:80000; default = 80))) μm  			
			---
			
			###### __Activity Coefficient__  
			- LiquidElectrolyte.Mode: $(Child("Lmode", Select(["DMGL_γ", "Stefan_γ"])))
			- NoteBook.Mode: $(Child("Nmode", Select(["DMGL_γ", "Stefan_γ", "Potassium_γ"])))
			- Boundary Condition : $(Child("BC_Select", Select(["Robin", "Neumann", "Dirichlet"])))
			"""
	    end;
	    label = "Submit"
	);
	top = 505
)

# ╔═╡ 2b9d9bfd-d660-4b4d-8f0b-b6b5bcc0dbfa
begin
    Vmax = 2 * V

    L = user_input_model.L * μm

    hmin = 1.0e-6 	* μm

    hmax = 1.0	* μm 

    X = ExtendableGrids.geomspace(0, L, hmin, hmax)

    grid = ExtendableGrids.simplexgrid(X)
end;

# ╔═╡ 7b032dac-97fb-4fd7-adc7-0dbd0e34d1a5
let
	try
	    # --- helpers ---
	    lastU(sol) = hasproperty(sol, :u) ? (sol.u isa AbstractVector ? sol.u[end] : sol.u) : sol
		
	    function to_inival(sys, u_raw)
	        Ut = VoronoiFVM.unknowns(sys)                 # template
	        u  = u_raw
	        if size(u) != size(Ut)
	            if size(u,1) == size(Ut,1) && size(u,2) ≥ size(Ut,2)
	                u = u[:, 1:size(Ut,2)]
	            elseif reverse(size(u)) == size(Ut)
	                u = permutedims(u)
	            else
	                error("Shape mismatch: u=$(size(u)) vs Ut=$(size(Ut))")
	            end
	        end
	        U0 = copy(Ut);  U0 .= u
	        return U0
	    end
	
	    # --- settings ---
	    tindex = 88
	 
		Vhold  = round(res.pnpresult.voltages[tindex], digits = 2) #-0.87 * ufac"V"
		thold  = round(res.pnpresult.tsol.t[tindex], digits = 2)
	    ispec  = ico2
	
	    # --- x axis ---
	    x  = grid.components[XCoordinates]
	    xx = x ./ μm
	
	    # --- CV snapshot ---
	    u_cv = res.pnpresult.tsol.u[tindex]               # Matrix
	    c_cv = @view u_cv[ispec, :]
	
	    # --- HOLD end (transient) ---
	    u_hold = lastU(res.hold_tsol)                     # Matrix
	    c_hold = @view u_hold[ispec, :]
	
	    # --- IV / steady at the same Vhold ---
	    pnpcell = res.pnpcell
	    sys     = pnpcell.vfvmsys
	    cdata   = celldata(pnpcell)
	
	    function pre_iv(sol, t)
	        LiquidElectrolytes.working_electrode_voltage!(cdata, Vhold)
	    end
	
	    control_iv = SolverControl(; verbose="", handle_exceptions=true, damp_initial=0.1)
	
	    Uinit = to_inival(sys, u_hold)
	
	    iv_sol = solve(sys; inival=Uinit, control=control_iv, pre=pre_iv)  # steady solve
	    u_iv   = lastU(iv_sol)
	    c_iv   = @view u_iv[ispec, :]
	
	    # --- plotting (log10 mol/dm^3) ---
	    to_moldm3(v) = v .* 1e-3
	    y(v) = [vv > 0 ? log10(vv) : NaN for vv in to_moldm3(v)]
	
	    fig = Figure(size=(780, 450))
	    ax  = Axis(fig[1,1],
	        xlabel="x / μm",
	        ylabel=L"\log_{10} c \; (mol/dm^3)",
	        title="CO2 @ $Vhold V, $thold s: CV snapshot vs HOLD end vs IV(steady)",
			xscale = log10,
			#limits = ((1e-20, 1000), (1e-50, 1e2)), 
	    )
	    lines!(ax, xx, y(c_cv),   label="CV (t=$thold)", linewidth=2)
	    lines!(ax, xx, y(c_hold), label="HOLD end", linewidth=2)
	    lines!(ax, xx, y(c_iv),   label="IV (steady)", linewidth=2, linestyle = :dash)
	    axislegend(ax; position=:rb)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 3a7f5ebe-aed6-4edb-84de-a99de0755453
function cv_conc_gif(pnpresult; file="concentrations_cv.gif", framerate=10, step=5)
    tsol = pnpresult.tsol ./ (mol / dm^3)
    nvar, nx, nt = size(tsol)

    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = min(nvar, length(species), length(colors))

    nplot = min(nx, length(X))
    xx = (X[1:nplot] .+ 1e-14)               
    t_indices = 1:step:nt

    fig = Figure(size=(650, 400))
    ax  = Axis(fig[1, 1];
        xlabel = "Distance from electrode [m]",
        ylabel = "log c(aᵢ)",
        xscale = log10,
        limits = ((1e-12, L), (-14, 2)),
    )

    ys = [Observable(fill(NaN, nplot)) for _ in 1:nspecies]
    for i in 1:nspecies
        lines!(ax, xx, ys[i]; color=colors[i], label=species[i])
    end
    axislegend(ax)

    record(fig, file, t_indices; framerate=framerate) do ti
        ax.title = "t = $(round(pnpresult.tsol.t[ti], digits=4)) s | ϕ = $(round(pnpresult.voltages[ti], digits=2)) V"
        for i in 1:nspecies
            conc   = tsol[i, 1:nplot, ti]
            ys[i][] = [c > 0 ? log10(c) : NaN for c in conc]
        end
    end

    return LocalResource(abspath(file))
end

# ╔═╡ a7537912-16d8-4312-a1a7-513695ad86de
function conc_vs_voltage_axis_compare(result; useonly_pH=false, showlegend=true,
                              catmap_csv::Union{Nothing,String}=nothing)

    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    tsol  = LiquidElectrolytes.voltages_solutions(result)
    vgrid = result.voltages

    xcoords    = grid.components[XCoordinates]
    ielectrode = argmin(xcoords)

    scale = 1.0 / (mol / dm^3)
    nv = length(vgrid)
    conc_electrode = fill(NaN, nspecies, nv)

    for (j, v) in enumerate(vgrid)
        sol = tsol(v)
        if sol === nothing
            @warn "No solution available at voltage $v; skipping."
            continue
        end
        @inbounds for ia in 1:nspecies
            conc_electrode[ia, j] = sol[ia, ielectrode] * scale
        end
    end

    # --- FIGURE STYLE (keep exactly) ---
    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = L"\text{Voltage}\ U\ \mathrm{vs.}\ \text{SHE}\ (V)",
        ylabel = L"c_i^{+}\;(\mathrm{M})",
        yscale  = log10,
        limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    ax.spinewidth = 2.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 30
    ax.ylabelsize = 30
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20
	ax.xgridvisible = false
	ax.ygridvisible = false
	

    # --- SIMULATION: solid thick ---
    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")

        y = max.(conc_electrode[iH, :], eps(Float64))
        lines!(ax, vgrid, y; color=colors[iH], linewidth=3, label=species[iH])
    else
        for ia in 1:nspecies
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=3, label=species[ia])
        end
    end

     for g in groupby(df, :Index)
        idx = Int(first(g.Index))
        p = sortperm(g.Voltage)

        V = g.Voltage[p]
        y = max.(g.Concentration[p], eps(Float64)) 

        if 1 <= idx <= nspecies
            lines!(ax, V, y; color=colors[idx], linewidth=2, label=species[idx], linestyle = :dashdot)
        else
            lines!(ax, V, y; linewidth=3, label="Index = $idx")
        end
    end

    axislegend(ax; position=:rt)  
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
						    clear = false)
			end
		end
	end

	function addplot(vis, df::DataFrame, vshow)
		species  = getproperty.(bulk, :name)
		colors   = getproperty.(bulk, :color)
		nspecies = length(species)

		nms = names(df)
		@assert length(nms) ≥ 3
		rename!(df, Dict(nms[1]=>:Index, nms[2]=>:Distance, nms[3]=>:Concentration))

		idxs = unique(skipmissing(df.Index))

		for idx_raw in idxs
			idx = try
				Int(idx_raw)
			catch
				try
					parse(Int, String(idx_raw))
				catch
					continue
				end
			end

			(1 <= idx <= nspecies) || continue

			mask = (df.Index .== idx_raw)
			x = df.Distance[mask]
			y = df.Concentration[mask]

			x = Float64.(x)
			if maximum(x) > 1e-3
				x .*= 1e-6
			end
			x .+= 1e-14

			y = max.(Float64.(y), eps(Float64))

			p = sortperm(x)

			scalarplot!(vis, x[p], log10.(y[p]);
					    color = colors[idx],
					    clear = false,
					    label = "")
		end
	end

	function plot1d(result, celldata, vshow; df_compare = nothing)
		tsol = LiquidElectrolytes.voltages_solutions(result)
		vis  = GridVisualizer(;
							  size    = (600, 300),
							  clear   = true,
							  legend  = :rt,
							  limits  = (-11, 1),
							  xlimits = (10e-12, L*1.2),
							  xlabel  = "Distance from electrode [m]",
							  ylabel  = "log c(aᵢ)",
							  xscale  = :log)

		addplot(vis, tsol(vshow), vshow)
		if !isnothing(df_compare)
			addplot(vis, df_compare, vshow)
		end
		reveal(vis)
	end

	function plot1d(result, celldata)
		tsol = LiquidElectrolytes.voltages_solutions(result)
		vis  = GridVisualizer(;
							  size    = (650, 400),
							  clear   = true,
							  legend  = :rt,
							  limits  = (-11, 1),
							  xlimits = (10e-12, L*1.2),
							  xlabel  = "Distance from electrode [m]",
							  ylabel  = "log c(aᵢ)",
							  xscale  = :log)

		vrange = result.voltages[end:-5:1]
		movie(vis, file="concentrations.gif", framerate=2) do vis
			for vshow_it in vrange
				addplot(vis, tsol(vshow_it), vshow_it)
				reveal(vis)
			end
		end
		isdefined(Main, :PlutoRunner) && LocalResource("concentrations.gif")
	end
end

# ╔═╡ c10697b7-6e67-4a5b-937e-09d97ca5b7f8
function conc_vs_voltage_axis(result; useonly_pH=false, showlegend=false)
    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    tsol  = LiquidElectrolytes.voltages_solutions(result)
    vgrid = result.voltages

    xcoords    = grid.components[XCoordinates]
    ielectrode = argmin(xcoords)

    scale = 1.0 / (mol / dm^3)
    nv = length(vgrid)
    conc_electrode = fill(NaN, nspecies, nv)

    for (j, v) in enumerate(vgrid)
        sol = tsol(v)
        if sol === nothing
            @warn "No solution available at voltage $v; skipping."
            continue
        end
        @inbounds for ia in 1:nspecies
            conc_electrode[ia, j] = sol[ia, ielectrode] * scale
        end
    end

    # --- FIGURE STYLE (match screenshot) ---
    fig = Figure(size=(960, 540))
    ax = Axis(fig[1, 1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{c_i^{+}}\;(\mathrm{M})",
        yscale  = log10,
       # limits = ((-1.25, -0.50), (1e-11, 1e1)),
    )

    # x ticks: -1.2, -1.0, -0.8, -0.6
    xt = [-1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])

    # y ticks: 10^0, 10^-3, 10^-6, 10^-9
    yt_vals = 10.0 .^ (0:-3:-9)
    yt_lbls = [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"]
    ax.yticks = (yt_vals, yt_lbls)

    # thick spines like screenshot
    ax.spinewidth = 5.5
    ax.xtickwidth = 2.0
    ax.ytickwidth = 2.0
    ax.xticksize  = 8
    ax.yticksize  = 8
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
	ax.xgridvisible = false
	ax.ygridvisible = false
	ax.xlabelpadding = 10
	ax.ylabelpadding = 10
	ax.xlabelfont = :bold

 #Species text 	
	# CO	
	#text!(ax, -1.15, 0.3, text=L"\mathrm{K^+}", color=colors[ikplus], fontsize=24, font = "sans-bold")
	#text!(ax, -0.90, 0.000005, text=L"\mathrm{H^+}", color=colors[ihplus], fontsize=24, font = "sans-bold")
	#text!(ax, -1.05, 5e-11, text=L"\mathrm{CO_3^{2-}}", color=colors[ico3], fontsize=24, font = "sans-bold")	
	#text!(ax, -1.0, 2.5e-7, text=L"\mathrm{HCO_3^-}", color=colors[ihco3], fontsize=24, font = "sans-bold") 
	#text!(ax, -1.2, 5.2e-6, text=L"\mathrm{CO_2}", color=colors[ico2], fontsize=24, font = "sans-bold")
	#text!(ax, -0.90, 7e-9, text=L"\mathrm{OH^-}", color=colors[iohminus], fontsize=24, font = "sans-bold")
	#text!(ax, -1.15, 0.00024, text=L"\mathrm{CO}", color=colors[ico], fontsize=24, font = "sans-bold")

	
    # --- PLOT ---
    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")

        y = max.(conc_electrode[iH, :], eps(Float64))  # log축 보호
        lines!(ax, vgrid, y; color=colors[iH], linewidth=5, label=species[iH])
    else
        for ia in 1:nspecies
            y = max.(conc_electrode[ia, :], eps(Float64))
            lines!(ax, vgrid, y; color=colors[ia], linewidth=5, label=species[ia])
        end
    end

    showlegend && axislegend(ax, position=:rt)
    #display(fig)
	@show minimum(xcoords) maximum(xcoords) xcoords[ielectrode]
    return fig
end

# ╔═╡ 754ab149-7a65-4227-ba8b-d7d48e0092b8
begin
	function addplot_ax!(ax, sol, vshow; clear=true)
	    species = getproperty.(bulk, :name)
	    colors  = getproperty.(bulk, :color)
	
	    scale = 1.0 / (mol / dm^3)
	    title = @sprintf("Φ_we=%+1.2f [V vs. SHE]", vshow)
	
	    x = grid.components[XCoordinates] .+ 1.0e-14
	
	    clear && empty!(ax)        
	    #ax.title = title          
	
	    if useonly_pH
	        i = findfirst(isequal("H⁺"), species)
	        i === nothing && error("H⁺ not found in species list.")
	
	        y = log10.(sol[ihplus, :] .* scale)
	        lines!(ax, x, y; color = colors[i], linewidth=5, label = species[i])
	    else
	        # first species (clear=true 느낌)
	        y1 = log10.(sol[1, :] .* scale)
	        lines!(ax, x, y1; color = colors[1], linewidth=5, label = species[1])
	
	        # rest
	        for ia in 2:nc
	            y = log10.(sol[ia, :] .* scale)
	            lines!(ax, x, y; color = colors[ia], linewidth=5, label = species[ia])
	        end
	    end
	
	    return ax
	end
	function plot1d_makie(result, celldata, vshow; df_compare=nothing)
	    tsol = LiquidElectrolytes.voltages_solutions(result)
	
	    fig = Figure(size=(960, 540))
	    ax  = Axis(fig[1, 1];
	        xlabel = "Distance from electrode [m]",
	        ylabel = L"\log_{10} c(a_i)",
	        xscale = log10,
	        limits = ((10e-12, L*1.2), (-11, 1)),
	    )

	    ax.spinewidth = 5.5
	    ax.xtickwidth = 2.0
	    ax.ytickwidth = 2.0
	    ax.xticksize  = 8
	    ax.yticksize  = 8
	    ax.xlabelsize = 25
	    ax.ylabelsize = 25
	    ax.xticklabelsize = 25
	    ax.yticklabelsize = 25
		ax.xgridvisible = false
		ax.ygridvisible = false
		ax.xlabelpadding = 10
		ax.ylabelpadding = 10
		#ax.xlabelfont = :bold
	
	    addplot_ax!(ax, tsol(vshow), vshow; clear=true)

		#text!(ax, 1e-10, -0.5, text=L"\mathrm{K^+}", color=colors[ikplus], fontsize=24, font = "sans-bold")
		#text!(ax, 4e-11, -6.6, text=L"\mathrm{H}^+", color=colors[ihplus], fontsize=24, font = "sans-bold")
		#text!(ax, 1.5e-9, -5.2, text=L"\mathrm{CO_3^{2-}}", color=colors[ico3], fontsize=24, font = "sans-bold")	
		#text!(ax, 7e-10, -3, text=L"\mathrm{HCO_3^-}", color=colors[ihco3], fontsize=24, font = "sans-bold") 
		#text!(ax, 4e-11, -2.4, text=L"\mathrm{CO_2}", color=colors[ico2], fontsize=24, font = "sans-bold")
		#text!(ax, 4e-11, -9, text=L"\mathrm{OH^-}", color=colors[iohminus], fontsize=24, font = "sans-bold")
		#text!(ax, 4e-11, -4.2, text=L"\mathrm{CO}", color=colors[ico], fontsize=24, font = "sans-bold")	


		
	    if df_compare !== nothing
	        addplot_ax!(ax, df_compare, vshow; clear=false)
	    end
	
	    #axislegend(ax, position=:rt)
	    return fig
	end
end

# ╔═╡ 9d814b85-a5b6-42e5-abf4-15500bbdb717
begin
	Lγ_key = user_input_model.Lmode 
	Lγ_mode = Lγ_key == "Stefan_γ" ? Stefan_γ! : DGML_γ!

	Nγ_key = user_input_model.Nmode
	Nγ_mode = Nγ_key == "Stefan_γ" ? Stefan_γ! : Nγ_key == "DMGL_γ" ? DGML_γ! : Potassium_γ! 
end;

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
							   	actcoeff! = Lγ_mode
							   )

# ╔═╡ 12235c3c-18f2-4fc7-95ef-800f71783036
elydata_toy = elydata_Gold

# ╔═╡ 60698d90-e67d-4d85-bc00-6935c95b5a69
println(nc, elydata_Gold.cspecies)

# ╔═╡ 5f17b4f7-54d6-4ad0-9886-252854840a80
function activity_coefficient!(
    γ::AbstractVector,
	u,
    data,
    mode::Function;
)
    (; v, ip, pscale, p_bulk, M, M0, v0, κ, RT, nc, Mrel, tildev, v0, cspecies, rexp) = data

    if Nγ_mode == Stefan_γ!
        # Ringe et al. approach (volume fraction-based)
    	for ic in cspecies
        	γ[ic] = (1.0 / (1 - sum(u[i] * v[i] for i in 1:nc)))# / (mol/dm^3)))
    	end
		 #.= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))
    elseif Nγ_mode == DGML_γ!
		
        # Dreyer et al. approach
        p = u[ip] * pscale - p_bulk
        c0, barc = c0_barc(u, data)
        for ic in cspecies
            γ[ic] = rexp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc)) #* (1 /barc)
		end

    else
		γ .= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))
	end

    return γ
end


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
		activity_coefficient!(γ, u, data, Nγ_mode)
		
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
	const ps_cache = DiffCache(zeros(18), 13)
	const us_cache = DiffCache(zeros(isurfaceend-isurfacestart+1), 13)
	
	function we_breactions(f, 
			u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
			bnode, 
			data
		) where {Tval, Tv, Tc, Tp, Ti}
		(; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we) = data
				
		γ = get_tmp(γ_cache, u[ico2])
		γ_co2 	 = activity_coefficient!(γ, u, data, Nγ_mode)[ico2]
		γ_co 	 = activity_coefficient!(γ, u, data, Nγ_mode)[ico]
		σ 			= C_gap * (ϕ_we - ϕ_pzc) #- u[iϕ]
		local_pH 	= -log10(u[ihplus] * γ[ihplus] / (mol/dm^3))

	
		#for (p, default_value) in odesys.defaults
		#	ps[paramsidx[p]] = default_value
		#	println(default_value)
		#end

		ps = get_tmp(ps_cache, u[iϕ])
		ps[paramsidx[Symbolics.rename(odesys.σ, :σ)]] = σ
		ps[paramsidx[Symbolics.rename(odesys.γCO2_aq, :γCO2_aq)]] = γ_co2 
		ps[paramsidx[Symbolics.rename(odesys.aH2O_g, :aH2O_g)]] = aH₂O 
		ps[paramsidx[Symbolics.rename(odesys.ϕ, :ϕ)]] = 0 #u[iϕ] 
		ps[paramsidx[Symbolics.rename(odesys.ϕ_we, :ϕ_we)]] = ϕ_we 
		ps[paramsidx[Symbolics.rename(odesys.local_pH, :local_pH)]] = local_pH 
		ps[paramsidx[Symbolics.rename(odesys.γCO_aq, :γCO_aq)]] = γ_co 
		ps[paramsidx[Symbolics.rename(odesys.βCOOHΔH2OΔele_t, :βCOOHΔH2OΔele_t)]] = 0.59 
		#ps[paramsidx[Symbolics.rename(odesys.ECO2_g, :ECO2_g)]] = 0.0 
		#ps[paramsidx[Symbolics.rename(odesys.ECO2_t, :ECO2_t)]] = 0.65*e
		#ps[paramsidx[Symbolics.rename(odesys.ECOOHΔH2OΔele_t, :ECOOHΔH2OΔele_t)]] = 0.95*e
		#ps[paramsidx[Symbolics.rename(odesys.E_t, :E_t)]] = 0.0 
		#ps[paramsidx[Symbolics.rename(odesys.Eele_g, :Eele_g)]] = 0.0 
		#ps[paramsidx[Symbolics.rename(odesys.ECO_g, :ECO_g)]] = 0.270185*e 
		#ps[paramsidx[Symbolics.rename(odesys.ECO_t, :ECO_t)]] = -0.02145*e 
		#ps[paramsidx[Symbolics.rename(odesys.ECOOH_t, :ECOOH_t)]] = 0.1282*e
		#ps[paramsidx[Symbolics.rename(odesys.EH2O_g, :EH2O_g)]] = 0.0 

	    #println[1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))]
		#ps[paramsidx[Symbolics.rename(odesys.ECOOH_g, :ECOOH_g)]] = 0.0 

		
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

# ╔═╡ 924f8f5d-2cb0-4381-a522-509ff4c002b6
begin
	model_key = user_input_model[:model_choice]

	model = model_key == "Gold_Model" ? elydata_Gold :
	        model_key == "Landstorfer_NaClO₄ model" ? elydata_NaClO₄ :
	        model_key == "Landstorfer_NaF model" ? elydata_NaF :
			model_key == "Toy model" ? elydata_toy :
	        error("Unknown model choice: $model_key")
end

# ╔═╡ 5df5ee46-b0d3-47a8-835b-b3f21a3cab34
is_Landstorfer = model != elydata_Gold

# ╔═╡ 36e756a9-4d9b-40ef-9e37-d86f1194cc51
function capscalc(sys; molarities =molarities)
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
	                voltages = vrange
	        )
			volts = r.voltages
			caps = r.dlcaps
	    end
	    cdl0 = dlcap0(data)
	    @info "elapsed=$(t)"
	    push!(result, (voltage_range = volts, dlcaps = caps, cdl0 = cdl0))
	end
    return result
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
        voltages = result[1].voltage_range
        caps = vec(result[1].dlcaps)

        scalarplot!(
            vis, voltages, caps / (μF / cm^2) ;
            limits = (-1, 250), xlimits = (-1.1, 1.1),
            color = :green, clear = false, label = "$title",
            title = title, markershape = :none, xlabel = L"φ / (V vs φ_{pzc})", ylabel = L"dlcaps / (μF / cm²)"
        )
        # scalarplot!(
        #     vis, [0], [result[1].cdl0] / (μF / cm^2);
        #     clear = false, markershape = :circle, markersize = 8, label = ""
        # )
    end
    return vis
end;

# ╔═╡ 52a5bbd2-0278-4d92-95f1-367797f636e6
begin
	try
		CVPlot!(res.pnpresult, model)
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 27075c18-4da9-42f4-b5a8-d36bc7b4930d
let
    ic = model.cspecies
    fig = Figure(size = (650, 400))
    ax = Axis(fig[1, 1], 
              ylabel = "Current Density (mA/cm²)",
              xlabel = "Voltage (ϕ-ϕₚ)"
    )
    
    colors = (:pink, :skyblue, :lightgreen) 
    
    raw_df = CSV.read("../data/Langmuir_CV_data/Figure_1.csv", DataFrame; header=false)

    facet_row = collect(raw_df[1, :])
    datatype_row = collect(raw_df[2, :])

	numeric_data = [
	    parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
	    for i in 3:nrow(raw_df)
	]
	num_df = DataFrame(hcat(numeric_data...)', names(raw_df))

    # Au(110)
    lines!(ax, num_df[!, 1], num_df[!, 2], color = colors[1], label = facet_row[1])
    # Au(111)
    lines!(ax, num_df[!, 3], num_df[!, 4], color = colors[2], label = facet_row[3])
    # Au(100)
    lines!(ax, num_df[!, 5], num_df[!, 6], color = colors[3], label = facet_row[5])

    axislegend(ax, position = :rt)
    fig
end

# ╔═╡ 138948ff-5e55-4d9a-86f8-542e21024964
let
    ic = model.cspecies
    fig = Figure(size = (650, 400))
    ax = Axis(fig[1, 1], 
              ylabel = "Current Density (mA/cm²)",
              xlabel = "Voltage (ϕ-ϕₚ)"
    )
    
    colors = (:pink, :skyblue, :lightgreen) 
    
    raw_df = CSV.read("../data/Langmuir_CV_data/Figure_1.csv", DataFrame; header=false)

    facet_row = collect(raw_df[1, :])
    datatype_row = collect(raw_df[2, :])

	numeric_data = [
	    parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
	    for i in 3:nrow(raw_df)
	]
	num_df = DataFrame(hcat(numeric_data...)', names(raw_df))

    # Au(110)
    lines!(ax, num_df[!, 1], num_df[!, 2], color = colors[1], label = facet_row[1])
    # Au(111)
    #lines!(ax, num_df[!, 3], num_df[!, 4], color = colors[2], label = facet_row[3])
    # Au(100)
    #lines!(ax, num_df[!, 5], num_df[!, 6], color = colors[3], label = facet_row[5])

    axislegend(ax, position = :rt)
    fig
end

# ╔═╡ dc203e95-7763-4b13-8408-038b933c5c9c
function pnp_bcondition(
	f,
	u, #::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
	bnode,
	data
)# where
	#{Tval, Tv, Tc, Tp, Ti}
	
	(; Γ_we, Γ_bulk, ϕ_we, iϕ, ϕ_bulk, ip, p_bulk, c_bulk, cspecies) = data
	
	if user_input_model.BC_Select == "Dirichlet"
		boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)

	elseif user_input_model.BC_Select == "Robin"
		boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))	


	else
		#continue
	#boundary_neumann!(f, u, bnode; species = ic, region, value = 0)
 	#boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))	
 	#for ic in cspecies 
		#if ic == ico2 || ic == ico
       	#	boundary_neumann!(f, u, bnode; species = ic, region = Γ_we, value = 0)
		#else
		#boundary_dirichlet!(f, u, bnode; species = ic, region = Γ_we, value = c_bulk[ic])
		#end
   #end
	end
		
		
	if bnode.region == Γ_we && model == elydata_Gold
		we_breactions(f, u, bnode, data)
	end


	
	return bulkbcondition(f, u, bnode, data; region = Γ_bulk)

end;

# ╔═╡ 3ef57b7d-ec19-46bc-a881-0506cf5167f3
begin
	if double_layer_curve
		#pnp
		reaction_arg = model == elydata_Gold ? (reaction) : NamedTuple()
		sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = deepcopy(model), reaction_arg)

		result_pnp = capscalc(sys_pnp)
	else
		result_pnp = nothing
	end
end

# ╔═╡ e50fe651-11d4-45ee-89dd-371a7fbc097e
function sweep(pnpdata; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    celldata.eneutral = eneutral
	#reaction_arg = model == elydata_Gold ? (; reaction) : NamedTuple()
    pnpcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = pnpdata)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
    )

end

# ╔═╡ 9fb47b83-a853-4316-bb8d-30e65b16ef78
if CV
	pnpresult = sweep(model; eneutral = false, tunnel = false)
end

# ╔═╡ f9dade9f-8431-48a6-a2ee-2c88f178e76e
let
	try
	    fig = Figure()
	    ax = Axis(fig[1, 1])
	    T = 0:1.0:pnpresult.tsol.t[end]
	    lines!(ax, T, sawtooth.(T))
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ cb9b0158-f17d-4136-8994-360f3078c7df
let
	try
	    fig = Figure(size = (600, 200))
	    ax = Axis(fig[1, 1], yscale = log10)
	    T = pnpresult.times
	    #lines!(ax,T, voltages.(T))
	    lines!(ax, T[2:end], T[2:end] - T[1:(end - 1)])
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ c62ab378-0988-4fa5-b21d-5e1622c63c87
let
	try
		ic = model.cspecies
	    fig = Figure(size = (650, 400))
	    ax = Axis(fig[1, 1], 
	 			  #limits = ((-0.5, 0.9),(-2e-20, 2e-20)),
	              ylabel = L"I (mA/cm²)",
	              xlabel = L"φ (V vs SHE)",
				 )
		
	    total_current = zero(currents(pnpresult, ic[1]))
	    #for s in ic
	    #    total_current .+= (currents(result, s) * mA / cm^2)
	    #end
		total_current = ((currents(pnpresult, ico) .* cm^2/mA))
	
		
	    lines!(ax, pnpresult.voltages, total_current,
	           color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))
	
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ c23741d2-3ece-41b8-8a2f-16743425c5cd
begin
	try
		conc_time_func(pnpresult, sawtooth)
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
let
    try
        tsol = pnpresult.tsol ./ (mol / dm^3)
        nvar, nx, nt = size(tsol)

        species = getproperty.(bulk, :name)
        colors  = getproperty.(bulk, :color)
        nspecies = min(nvar, length(species), length(colors))

        nplot = min(nx, length(X))
        t_index = 1

        xx = X[1:nplot] .+ 1e-14   # [m]

        fig = Figure(size = (650, 400))
        ax = Axis(fig[1, 1];
            xlabel = "Distance from electrode [m]",
            ylabel = "log c(aᵢ)",
            xscale = log10,                
            limits = ((1e-12, 1e-3), (-14, 1)), 
            title  = "t = $(round(pnpresult.tsol.t[t_index], digits=4)) s | ϕ = $(round(pnpresult.voltages[t_index], digits=2)) V",
			
        )

        for i in 1:nspecies
            conc  = tsol[i, 1:nplot, t_index]
            yvals = [c > 0 ? log10(c) : NaN for c in conc]
            lines!(ax, xx, yvals; color = colors[i], label = species[i])
        end

        axislegend(ax;position=:rt)
        fig

    catch e
        println("⚠️ Error occurred: ", e)
        println(stacktrace(catch_backtrace()))
    end
end


# ╔═╡ 2420382d-227a-4063-9450-1f1726df018e
cv_conc_gif(pnpresult; file="concentrations_cv.gif", framerate=8, step=3)

# ╔═╡ 3d661549-a8d2-40b0-add8-b186193f90fe
let
	try
	    ic = model.cspecies
	    fig = Figure(size = (1050, 650))
	    ax = Axis(fig[1, 1],
				  #limits = ((-1.3, 1.2),(-6, 2)),
	              ylabel = L"I (mA/cm²)",
	              xlabel = L"φ (V vs SHE)"
	    )
	    
	    colors = (:pink, :skyblue, :lightgreen) 
	    
	    # === Gold Model ===
	    total_current = currents(pnpresult, iohminus).* cm^2/mA
	    gold_line = lines!(ax, pnpresult.voltages, total_current,                        color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))
	    labels1 = ["CO2RR Gold Model"]
	
	    # === Koper Data ===
	    raw_df = CSV.read("Langmuir 2021, 37, 5707−5716/Figure_1.csv", DataFrame; header=false)
	    facet_row = collect(raw_df[1, :])
	    
	    numeric_data = [
	        parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
	        for i in 3:nrow(raw_df)
	    ]
	    num_df = DataFrame(hcat(numeric_data...)', names(raw_df))
	
	    # === Facet data ===
	    facet_lines = [
	        lines!(ax, num_df[!, 1].-0.4, num_df[!, 2], color = colors[1]),
	        lines!(ax, num_df[!, 3].-0.4, num_df[!, 4], color = colors[2]),
	        lines!(ax, num_df[!, 5].-0.4, num_df[!, 6], color = colors[3])
	    ]
	    labels2 = [facet_row[1], facet_row[3], facet_row[5]]
	
	    Legend(fig[1, 2],
	        [[gold_line], facet_lines],   
	        [labels1, labels2],           
	        ["Model", "Koper\nFacets"];          
	    )
	
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 81c4e515-89b5-4ecf-8437-070e5a51cb4c
let
	try
	    fig = Figure(size = (1600, 900))
	    ax = Axis(fig[1, 1], ylabel = L"I (mA/cm²)", xlabel = L"φ (V vs SHE)")
	
		total_current = currents(pnpresult, ico) .* (cm^2/mA)
	    gold_line = lines!(ax, pnpresult.voltages, total_current,
	                       color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))
		scatter!(ax, pnpresult.voltages, total_current, markersize = 12,
	                       color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))
	    # Experimental Data Plotting based on M.T.M Koper
	    raw = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
	    pres = vec(Matrix(raw[1:1, :]))
	    sub = Matrix(raw[4:end, :])
	    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
	    num_df = DataFrame(num, :auto)
	    npairs = size(num_df, 2) ÷ 2
	    pink, pblue = RGB(1.0, 0.7, 0.8), RGB(0.2, 0.5, 1.0)
	    cols1 = [RGB(pink.r + t*(pblue.r-pink.r),
	                 pink.g + t*(pblue.g-pink.g),
	                 pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]
	
	    plot_objs1 = []
	    labels1 = String[]
	    for j in 1:npairs
	        xcol, ycol = 2j - 1, 2j
	        label = j == 1 ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
	        push!(labels1, label)
	        line = lines!(ax, num_df[!, xcol], ((num_df[!, ycol])); color = cols1[j])
	        push!(plot_objs1, line)
	    end
	    Legend(fig[1, 2], plot_objs1, labels1, "Experimental"; framevisible = true)
		fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 69e6f136-e31c-4d71-a9a7-ea9ca6530669
let
    try
        tsol = pnpresult.tsol ./ (mol / dm^3)
        nvar, nx, nt = size(tsol)

        species = getproperty.(bulk, :name)
        colors  = getproperty.(bulk, :color)
        nspecies = min(nvar, length(species), length(colors))

        nplot = min(nx, length(X))
        xx = X[1:nplot] ./ μm

        ispec = 5       
        frac  = 0.99          
        t_indices = unique(clamp.([1, 10, 50, 100, 140, nt], 1, nt))  

        c_bulk = tsol[ispec, nplot, 1]
        if !(c_bulk > 0)
            error("c_bulk is not positive (c_bulk=$(c_bulk)). Try different ispec or far-field index.")
        end

        function x_at_frac(tsol, ispec, t_index, xx, c_bulk, frac)
            c = tsol[ispec, 1:length(xx), t_index]
            target = frac * c_bulk

            idx = findfirst(ci -> (ci ≥ target), c)
            return isnothing(idx) ? NaN : xx[idx]
        end

        δs = Float64[]
        ts = Float64[]
        for ti in t_indices
            δ = x_at_frac(tsol, ispec, ti, xx, c_bulk, frac)
            push!(δs, δ)
            push!(ts, pnpresult.tsol.t[ti])
            @info "t=$(round(ts[end], digits=6)) s → δ$(Int(round(frac*100)))=$(round(δ, digits=6)) μm (c_bulk≈$(round(c_bulk, digits=6)) mol/dm^3)"
        end

        t_index = min(140, nt)

        fig = Figure(size = (800, 420))
        ax  = Axis(fig[1, 1],
                   xlabel = "x / μm",
                   ylabel = L"\log_{10} c_i\; (mol/dm^3)",
                   title  = "t = $(round(pnpresult.tsol.t[t_index], digits=4)) s  |  δ$(Int(round(frac*100))) for $(species[ispec])")

        conc = tsol[ispec, 1:nplot, t_index]
        yvals = [c > 0 ? log10(c) : NaN for c in conc]
        lines!(ax, xx, yvals; color = colors[ispec], linewidth = 2, label = species[ispec])

        δ_here = x_at_frac(tsol, ispec, t_index, xx, c_bulk, frac)
        if isfinite(δ_here)
            vlines!(ax, [δ_here]; linestyle = :dash, linewidth = 2)
            text!(ax, δ_here, maximum(skipmissing(yvals));
                  text = "  δ$(Int(round(frac*100)))≈$(round(δ_here, digits=4)) μm",
                  align = (:left, :top))
        end

        axislegend(ax; position = :rb)
        fig

    catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
    end
end


# ╔═╡ e5956bb0-a33a-488d-906e-fb5a7e2473a9
let
	try
	    ic = model.cspecies
	    fig = Figure(size = (1050, 650))
	    ax = Axis(fig[1, 1], 
	  			  #limits = ((-1.25, 0.8),(-0.02, 0.1
										 
										 #)),
	              ylabel = L"I (mA/cm²)",
	              xlabel = L"φ (V vs SHE)"
	    )
	    colors = [RGB(i/3, 0, 1-(i/3)) for i in 1:3]
	
	    total_current = currents(pnpresult, iohminus) .* (cm^2/mA)
	    gold_line = lines!(ax, pnpresult.voltages, total_current,
	                       color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))
	    labels1 = ["CO2RR Gold Model"]
	
	    raw_df = CSV.read("Langmuir 2021, 37, 5707−5716/Figure_5.csv", DataFrame; header=false)
	    pH_row = collect(raw_df[1, :])
	    electrolyte_row = collect(raw_df[2, :])
	
	    numeric_data = [
	        parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
	        for i in 4:nrow(raw_df)
	    ]
	    num_df = DataFrame(hcat(numeric_data...)', names(raw_df))
	
	    conc_lines = [
	        lines!(ax, num_df[!, 1], num_df[!, 2], color = colors[1]),
	        lines!(ax, num_df[!, 3], num_df[!, 4], color = colors[2]),
	        lines!(ax, num_df[!, 5], num_df[!, 6], color = colors[3])
	    ]
	    labels2 = [
	        electrolyte_row[1]*"\t"*pH_row[2]*"pH",
	        electrolyte_row[3]*"\t\t"*pH_row[4]*"pH",
	        electrolyte_row[5]*"\t\t"*pH_row[6]*"pH"
	    ]
	
	    Legend(fig[1, 2],
	        [[gold_line], conc_lines], 
	        [labels1, labels2],         
	        ["Model", "Koper\nElectrolyte"];   
	    )
	
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ a42b1afb-86f2-4a11-8328-c726b614aaba
begin
	if pH_varied_checkbox
		pH_var = [3.0, 4.0, 5.0, 6.0, 6.8, 7.0, 8.0, 9.0]	
	    pH_vec = Any[]  
	    for p in pH_var
	        ely_pressure = deepcopy(elydata_Gold)   
	        ely_pressure.c_bulk[2] = 10.0.^(-p)
			ely_pressure.c_bulk[6] = 10.0.^(-14+p)
		#	ely_pressure.c_bulk[5] = base_CO2 .* p
	
	        pnp_rec = sweep(ely_pressure; eneutral=true, tunnel=false)
	
	        push!(pH_vec, pnp_rec)              
	    end
	end
end

# ╔═╡ 0607672c-9177-4717-8ddf-e07a5dd82ec4
let
	try
	    fig = Figure(size = (1600, 900))
	    ax = Axis(fig[1, 1], ylabel = L"I (mA/cm²)", xlabel = L"φ (V vs SHE)")
	
		ntheo = min(length(pH_var), length(pH_vec))
	    cols2 = [RGB(0.5 - 0.1*(i/ntheo), 0.5 - 0.3*(i/ntheo), 0.4 + 0.7*(i/ntheo)) for i in 1:ntheo]
	
	    plot_objs2 = Makie.AbstractPlot[]
	    labels2 = String[]
	
	    for (j, (p, rec)) in enumerate(zip(pH_var[1:ntheo], pH_vec[1:ntheo]))
	        label2 = (j == 1) ? "pH=$(p)" : "⋅ pH=$(p)"
	        push!(labels2, label2)
	        I = currents(rec, ico) .* cm^2/mA
	        push!(plot_objs2, lines!(ax, rec.voltages, I; color = cols2[j]))
	    end
	
	    Legend(fig[1, 2], plot_objs2, labels2, "Theoretical"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ d4fb4803-1c1b-4fd7-a782-112777f55be0
let
	try
	    fig = Figure(size = (1600, 900))
	       ax = Axis(fig[1, 1],
	        xlabel = L"φ (V vs SHE)",
	        ylabel = L"I (mA/cm²)",
	        #yscale = log10,
	       # yminorticksvisible = true,  
	       # yminorticks = IntervalsBetween(5),
			#limits = ((-0.3, 1.0),(1e-4, 1.5e-1))
	    )
	
	
	    # Experimental Data Plotting based on M.T.M Koper
	    raw = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
	    pres = vec(Matrix(raw[1:1, :]))
	    sub = Matrix(raw[4:end, :])
	    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
	    num_df = DataFrame(num, :auto)
	    npairs = size(num_df, 2) ÷ 2
	    pink, pblue = RGB(1.0, 0.7, 0.8), RGB(0.2, 0.5, 1.0)
	    cols1 = [RGB(pink.r + t*(pblue.r-pink.r),
	                 pink.g + t*(pblue.g-pink.g),
	                 pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]
	
	    plot_objs1 = []
	    labels1 = String[]
	    for j in 1:npairs
	        xcol, ycol = 2j - 1, 2j
	       # label = j == 1 ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
	       # push!(labels1, label)
	        line = lines!(ax, num_df[!, xcol], ((num_df[!, ycol])); color = cols1[j])
	        #line = lines!(ax, num_df[!, xcol], ((num_df[!, ycol])); color = cols1[j])
	        push!(plot_objs1, line)
	    end
	   # Legend(fig[1, 2], plot_objs1, "Experimental"; framevisible = true)
	
	    # Theoretical Data Plotting based on `LiquidElectrolytes.jl`
	    cols2 = [RGB(1 - i/length(pH_vec), 0, i/length(pH_vec)) for i in 0:length(pH_vec)]
	    plot_objs2 = []
	    labels2 = String[]
	    for (j, rec) in enumerate(pH_vec)
	       # label2 = j == 1 ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
	      #  push!(labels2, label2)
	        line = scatterlines!(ax, rec.voltages, ((currents(rec, iohminus) .* cm^2/mA)); color = cols2[j])
			#line = lines!(ax, rec.voltages, ((currents(rec, iohminus) .* cm^2/mA)); color = cols2[j])
	        push!(plot_objs2, line)
	    end
	   # Legend(fig[1, 3], plot_objs2, labels2, "Theoretical"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ df5b1bfb-ce96-4d32-abdb-a6fcecb195a1
let
	try
	    fig = Figure(size = (1600, 900))
	    ax  = Axis(fig[1, 1],
	        xlabel = L"\phi \, (\mathrm{V \; vs \; SHE})",
	        ylabel = L"I \; (\mathrm{mA/cm^2})",
	       # limits = ((-1.25, 0.8), (-0.5, 0.07))
	    )
	
		csv_path = "../data/Langmuir_CV_data/Figure_5.csv"  
		unit_scale = cm^2/mA 
	
	    #Experimental (Koper, Langmuir 2021, Fig.3)
	    raw = CSV.read(csv_path, DataFrame; header = false)
	    pres_labels = vec(Matrix(raw[1:1, :])) 
	    sub  = Matrix(raw[4:end, :])
	    num  = map(x -> x === missing ? NaN : parse(Float64, x), sub)
	    num_df = DataFrame(num, :auto)
	
	    npairs = size(num_df, 2) ÷ 2
	
	    exp_colors = [RGB(0.9 - (0.05 * i/npairs), 0.8 * (1 - i/npairs), 0.7 + 0.3 * i/npairs) for i in 1:npairs]
	
	
	    exp_plots = Plot[] 
	    exp_labels = String[]
	
	    for j in 1:npairs
	        xcol, ycol = 2j - 1, 2j
			label = "$(pres_labels[min(2j, length(pres_labels))]) pH"
	        push!(exp_labels, label)
	
			#---normal scale---
	        #h = lines!(ax, num_df[!, xcol], num_df[!, ycol]; color = exp_colors[j])
			#---log scale---
			h = lines!(ax, num_df[!, xcol], num_df[!, ycol]; color = exp_colors[j])
	        push!(exp_plots, h)
	    end
	
	    #Theoretical (LiquidElectrolytes.jl 결과)
	    theo_colors = [RGB(0.6 - (0.2 * i/npairs), 0.5 * (1 - i/npairs), 0.3 + 0.6 * i/npairs) for i in 1:npairs]
	
	    theo_plots = Plot[]
	    theo_labels = String[]
	
	    for j in 1:npairs
	        rec = pH_vec[j]
	   		label2 = "$(pres_labels[min(2j, length(pres_labels))]) pH"
	        push!(theo_labels, label2)
	
	        j_tot = currents(rec, iohminus) .* unit_scale
			#---normal scale---
	        #h = lines!(ax, rec.voltages, j_tot; color = theo_colors[j], linestyle = :dash)
			#---log scale---
			h = lines!(ax, rec.voltages, j_tot; color = theo_colors[j], linestyle = :dash)
	        push!(theo_plots, h)
	    end
	
	    Legend(fig[1, 2], exp_plots,  exp_labels,  "Experimental"; framevisible = true)
	    Legend(fig[1, 3], theo_plots, theo_labels, "Theoretical";  framevisible = true)
	
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end

# ╔═╡ 3960e081-090e-443f-b82c-06b703686e9f
pressure_varied_sweep(elydata_Gold, sweepfun = sweep; Pvec = [0.1, 1], ispec = ico)

# ╔═╡ e210999a-7ea5-47f7-aedb-f84fe77bd653
sweep_over_L(
    elydata_Gold;
    bcond    = pnp_bcondition,
    voltages = sawtooth,
    nperiods = nperiods,
    L_values = round.(range(100, 200, length=3)),
    sweepfun = cvsweep
)

# ╔═╡ 84d1270b-8df5-4d5d-a153-da4ffdb1d283
function simulate_CO2R(grid, celldata; voltages = (-1.5:0.1:0.0) * V, kwargs...)
	kwargs 	 	= merge(solver_control, kwargs) 
    cell        = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)
	ivresult    = ivsweep(cell; voltages, store_solutions=true, kwargs...)

	cell, ivresult
end;

# ╔═╡ 11b12556-5b61-42c2-a911-4ea98a0a1e85
cell, ivresult = simulate_CO2R(grid, model)

# ╔═╡ 659091d3-60b2-4158-80e2-cd28a492e870
(~, default_index) = findmin(abs, ivresult.voltages .+ 0.9 * ufac"V");

# ╔═╡ 15fadfc2-3cf8-4fda-9aed-a79c602b1d51
plot1d(ivresult, celldata)

# ╔═╡ afb700c7-ef29-4c13-b9ea-1d40ea9534dd
let
	scale = 1 / (mol / dm^3)
	volts = ivresult.voltages[ivresult.voltages .< -0.4]
	vis = GridVisualizer(;
	                         size = (600, 400),
	                         tilte = "IV Curve",
	                         xlabel = L"\phi_{we} \, (\mathrm{V \; vs \; SHE})",
	                         ylabel = L"I / (\mathrm{mA/cm^2})",       
	                         legend = :lb,
							 yscale = :log,
		)
							 
	scalarplot!(vis,
	                volts,
	                abs.(currents(ivresult, iohminus))[ivresult.voltages .< -0.4] .* cm^2/mA;
	                color = :green,
	                clear = false,
	                linestyle = :solid,
	                label = "e⁻, we")


	
	table = readdlm("./catmap_CO2R_data/IV-Ringe-digitized.csv", ',', Float64, '\n')
	df = Dict(:voltage => table[:,1], :current => table[:,2])
	table2 = readdlm("./catmap_CO2R_data/Ringe-theorical.csv", ',', Float64, '\n')
	df2 = Dict(:voltage => table2[:,1], :current => table2[:,2])
	table3 = readdlm("./catmap_CO2R_data/Ringe-experimental.csv", ',', Float64, '\n')
	df3 = Dict(:voltage => table3[:,1], :current => table3[:,2])

	
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

	scalarplot!(vis,
						df2[:voltage],
						df2[:current],
						clear = false,
						linewidth = 0,
						markershape = :circle,
						markersize = 4,
						markevery = 1,
						color = :blue,
						label = "Ringe et. al : Theorical")
	scalarplot!(vis,
						df3[:voltage],
						df3[:current],
						clear = false,
						linewidth = 0,
						markershape = :utriangle,
						markersize = 8,
						markevery = 1,
						color = :magenta,
						label = "Ringe et. al : Experimental")
	
	    reveal(vis)
end

# ╔═╡ 2c239f3a-6335-4dde-bdcc-7bf41bc49890
conc_vs_voltage_axis_compare(ivresult; useonly_pH = false)

# ╔═╡ f8255707-2233-4e28-b542-2f3d81b31c2e
conc_vs_voltage_axis(ivresult; useonly_pH = false)

# ╔═╡ 22244e24-5b56-4933-8a09-44b601f116c3
iv_curve_axis(ivresult; cutoff=-0.4, showlegend=true)

# ╔═╡ 5caca8ea-82af-4999-93bb-a72252c456c7
function ivsweep_over_L(model;
 	L_values = round.(Int, range(80, 1500, length=6)),
    eneutral = true,
	kwargs...)
    results = Dict{Int, Any}()
    kwargs 	 	= merge(solver_control, kwargs) 

	cell, ivresult
    for L in L_values
        hmin = 1.0e-6 * μm
        hmax = 1.0    * μm
        X = ExtendableGrids.geomspace(0, L * μm, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)

        celldata = deepcopy(model)
        #celldata.eneutral = eneutral
        #celldata.tunnel   = tunnel
        #celldata.bikerman = bikerman

        reaction_kw = (model === elydata_Gold) ? (; reaction) : (;)
    #cell   = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)

        pnpcell = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata=celldata)

        results[L] = ivsweep(pnpcell; voltages, store_solutions=true, kwargs...)

    end

    return results
end

# ╔═╡ 60b410be-70f7-4053-a3db-7d777e0d3f08
if Ldependancy
	ivL = ivsweep_over_L(elydata_Gold)
end

# ╔═╡ bab42c91-2d00-463d-a921-97487e4eac67
plotcurr_over_L(ivL; species=iohminus, cutoff=-0.4, title="IV vs L (log scale)")

# ╔═╡ 9a4e01d9-f469-4427-bf4c-883adb67ae24
function pb_bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data

	if user_input_model.BC_Select == "Dirichlet"
	    ## Dirichlet ϕ=ϕ_we at Γ_we
	    boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
	elseif user_input_model.BC_Select == "Robin"
		## Robin ϕ=dϕ₀/dx 
	    boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))
	else
		## neumann ϕ=dϕ₀/dx 
	    #boundary_neumann!(f, u, bnode, species = iϕ, region = Γ_we, value = C_gap * (ϕ_we - ϕ_pzc))
	
	end

    return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ 084e2127-ea77-4894-8990-380c2e8802c7
begin
	if double_layer_curve
		#pb
		sys_pb = PBSystem(grid; bcondition = pb_bcondition,  celldata = deepcopy(model))
		result_pb = capscalc(sys_pb)
	else 
		result_pb = nothing
	end
end

# ╔═╡ 4f991d6d-3a3f-45d8-b2e0-662c5292251c
let
    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, layout = (1, 2), size = 	(650, 350))
    capsplot(vis[1, 1], result_pb, "Poisson-Boltzmann")
    capsplot(vis[1, 2], result_pnp, "Poisson-Nernst-Planck")

    reveal(vis)
end

# ╔═╡ 50ccc291-3625-4639-afd1-5209899d904e
let
	if is_Landstorfer
		f = Figure(size = (900, 900))
	    result = result_pb
	    l = 1 / length(result)
		ϕ0_pzc = 0.972
		
		ax = Axis(f[1, 1], xlabel="φ / (V vs φ_pzc)", ylabel="dlcaps / (μF / cm²)", title="CSV Plot", limits = ((-1.2, 1.2),(0, 120)))
	
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

# ╔═╡ 46d92e15-38ca-4857-8db4-60c1519523f6
steady_state_jac

# ╔═╡ 315dd351-9d68-48f1-aa7a-8f43f3dec6ac
floataside(
    md"""
    __Input Voltage Index:__ $(@bind vindex PlutoUI.Slider(1:5:length(ivresult.voltages), default=default_index))
    """,
    top = 855
)


# ╔═╡ c4876d26-e841-4e28-8303-131d4635fc23
md"""
Potential at the working electrode 
$(vshow = ivresult.voltages[vindex]; @sprintf("%+1.4f", vshow))
"""

# ╔═╡ 5dd1a1e6-7db1-479e-a684-accec53ce06a
plot1d(ivresult, celldata, vshow)

# ╔═╡ 38e147f2-7175-4975-a799-0a2e12c24368
begin
	#df_cmp = CSV.read("catmap_CO2R_data/Dist-conc.csv", DataFrame)
	plot1d(ivresult, celldata, vshow)
end

# ╔═╡ f88ecd40-1b80-4cca-a312-b23f7cfb0ad6
plot1d_makie(ivresult, celldata, vshow)

# ╔═╡ 06ae600d-3f73-49e6-858c-539079c117ab
floataside(
    md"""
    __Input Time Index:__ $(@bind it PlutoUI.Slider(1:length(pnpresult.tsol.t)-1, show_value=false))
    """,
    top = 965
)

# ╔═╡ 7454f68a-64dc-4676-b2b2-ed8fcb35d81e
floataside(
	md"""
	**Time:** $(round(pnpresult.tsol.t[it+1], digits=3)) s
	""",
	top = 865
)

# ╔═╡ 39683e98-dcbb-458b-817f-856fc6498730
floataside(
    md"""
    **Input Voltage Index:** $(@bind vindex2 PlutoUI.Slider(1:5:length(target_time), default=40))
    """,
    top = 910
)

# ╔═╡ Cell order:
# ╠═91ac9e35-71eb-4570-bef7-f63c67ce3881
# ╠═68fe205c-9dd0-441b-9f12-3ddc12ec0a0d
# ╠═a94bc4e1-506f-4e40-bfe8-1ce7e6093974
# ╠═22f2574c-abbe-4a06-89e9-14635fb30932
# ╠═bd8134d8-5a69-486e-8429-7cf810b3ccbe
# ╟─beae1479-1c0f-4a55-86e1-ad2b50174c83
# ╟─ab2184fc-0279-46d9-9ee4-88fe3e732789
# ╠═7316901c-d85d-48e9-87dc-3614ab3d81a5
# ╟─6b7cfe87-8190-40a5-8d25-e39ef8d55db5
# ╠═5a146a44-03dc-45f3-ae15-993d11c2edac
# ╠═a1b895c2-d3d1-4afd-bd85-da575effef1c
# ╠═00947475-c96e-4ecc-a1ef-5be5e3e3c864
# ╠═ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
# ╟─06f52599-7006-4a5c-ba86-0b668b6952c9
# ╟─4b64e168-5fe9-4202-9657-0d4afc237ddc
# ╟─de2c826d-6c05-47cf-b5f5-44a00ea9889c
# ╟─d8f00649-e2ed-4bdd-853f-05268f0d5353
# ╠═47b36c81-b57e-4dd0-a22f-999e4fd3ac9f
# ╠═1e877f17-0219-45f1-b640-3a25ae085dbd
# ╠═a72aa512-941c-4815-b69f-eaa9d94e6b06
# ╠═68491839-1398-4807-9a50-07c87bb349cf
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
# ╠═42efe97d-bf43-42b3-8ff0-e3d98a7597f1
# ╟─6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
# ╠═e510bce3-d33f-47bb-98d6-121eee8f2252
# ╠═848b7aeb-968f-4116-8038-b61276f02b6c
# ╠═f18dc873-1c9d-46d3-9596-92d28705e894
# ╠═12235c3c-18f2-4fc7-95ef-800f71783036
# ╟─53f12821-7d8d-4971-87fd-ad4689ec62a5
# ╟─e1e0ca0f-7f88-40f0-850e-590b25da0331
# ╠═0db74a70-af86-492c-affb-9de62ffe4455
# ╟─4f388fe0-6bc8-4a29-bccc-fa725e62c6a7
# ╠═5d179c52-43d7-4bcb-a2df-93c5806876fa
# ╠═161a810d-c05e-42ad-97ab-131059d6784a
# ╠═9d814b85-a5b6-42e5-abf4-15500bbdb717
# ╠═5f17b4f7-54d6-4ad0-9886-252854840a80
# ╠═60698d90-e67d-4d85-bc00-6935c95b5a69
# ╠═2a20d9be-6c1e-4c1f-8bb6-a7693800732d
# ╟─4f7ec19d-cd60-4c2b-a766-7557caa471c0
# ╠═952a26ce-2610-48cc-9158-eda816da3a1c
# ╠═924f8f5d-2cb0-4381-a522-509ff4c002b6
# ╠═dc203e95-7763-4b13-8408-038b933c5c9c
# ╠═9a4e01d9-f469-4427-bf4c-883adb67ae24
# ╠═36e756a9-4d9b-40ef-9e37-d86f1194cc51
# ╟─4656ee04-ae86-442f-b37c-c5563170f992
# ╠═d76d8413-c019-4728-b182-7f7cb78dede4
# ╠═084e2127-ea77-4894-8990-380c2e8802c7
# ╠═3ef57b7d-ec19-46bc-a881-0506cf5167f3
# ╠═4f991d6d-3a3f-45d8-b2e0-662c5292251c
# ╠═50ccc291-3625-4639-afd1-5209899d904e
# ╠═5df5ee46-b0d3-47a8-835b-b3f21a3cab34
# ╟─783e2058-c720-4f31-8e51-7c313813924c
# ╠═f8e6c01b-e64e-4fa2-a84a-e20f9b60287c
# ╠═f0aecbbc-3c8c-4984-8704-fd79f986beb2
# ╟─9598e2c6-521e-4f8d-82d8-a836809736f3
# ╠═f2043f2c-f3c8-4b0f-944c-7b55624dac08
# ╠═0e61a0f7-1611-4eb3-8fda-3f807a4ffca2
# ╠═8bfdf2f5-c80a-4ce0-a8e1-b315affffb5f
# ╠═b21c8394-f847-477e-ac8f-713398b81166
# ╠═44258eea-f114-4dfe-aa61-1e2cac31baa4
# ╟─5c808c71-6094-49d7-8215-e88262f34e1f
# ╟─da8390d1-47e8-451f-b12b-45b8aca7b6ec
# ╠═b4aaf070-d4ab-409a-b1e8-f5469b9f398b
# ╠═e50fe651-11d4-45ee-89dd-371a7fbc097e
# ╟─ef7212fc-a3d0-4784-b901-219204b79dc0
# ╟─b95160b5-18f7-49d9-80be-9159abd2dcd1
# ╠═9fb47b83-a853-4316-bb8d-30e65b16ef78
# ╠═39c8ef0d-aac2-4c7f-8004-4166c460ebc5
# ╟─fba53a72-5d27-4db7-8453-41200db29481
# ╟─4af42420-4630-42a8-8941-6702c16beb99
# ╟─83defc30-2532-490f-ac19-f883f06e6e2d
# ╟─57a4947d-8de3-42de-b17a-9b769ca0589d
# ╠═b2546d23-825e-4356-a640-3fd53852cdcf
# ╟─2d950a96-9404-4b78-9db2-42ae2a4c44bf
# ╟─dfd42e6e-a98e-4759-b988-52dc5a793f15
# ╠═12a4df23-3c63-41d6-bd50-d7209e423cb8
# ╟─05c8e2fa-8650-49d0-a190-b865c0fb3261
# ╠═3ba47c4f-b8f8-41e1-bada-8146062b099e
# ╠═3ec78a69-a7b3-4c32-82cb-5b863c88798a
# ╟─25eb8aa3-697e-4538-9472-ceea45fbfbd9
# ╟─11892724-1851-46f2-802d-4da45127b0af
# ╠═a42b1afb-86f2-4a11-8328-c726b614aaba
# ╟─b3649b03-25cd-4f3e-99c4-85e24ddd3d11
# ╠═fee347ff-5401-4540-a1ce-fc2e8ff0ce63
# ╟─c048e472-3983-4279-bf60-82784baa145e
# ╟─3bdaab98-c0f7-46af-86b7-d68374e8a5d0
# ╠═d38c2b43-4d8b-4be7-8d77-5a30da384541
# ╠═1f085f56-e0ee-4cb5-a37e-eb82ef3d7589
# ╟─eb920b6e-86a6-4dd6-8e66-6b7e27d81257
# ╟─79018ef0-6ab1-4420-9a52-8f8e2812fd40
# ╟─f9dade9f-8431-48a6-a2ee-2c88f178e76e
# ╟─cb9b0158-f17d-4136-8994-360f3078c7df
# ╠═c62ab378-0988-4fa5-b21d-5e1622c63c87
# ╠═52a5bbd2-0278-4d92-95f1-367797f636e6
# ╟─278dd577-2d1f-4608-aee4-f7de466cf736
# ╠═4e894347-2ce6-4c5f-a06e-7f1af1983bbc
# ╟─c23741d2-3ece-41b8-8a2f-16743425c5cd
# ╠═2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
# ╠═2420382d-227a-4063-9450-1f1726df018e
# ╟─3f30fae0-18d4-4e5f-9618-cbd9852d7857
# ╟─7dd05779-3ffc-471c-9ae9-4bb00b45b7e8
# ╠═27075c18-4da9-42f4-b5a8-d36bc7b4930d
# ╟─7b032dac-97fb-4fd7-adc7-0dbd0e34d1a5
# ╟─3d661549-a8d2-40b0-add8-b186193f90fe
# ╠═138948ff-5e55-4d9a-86f8-542e21024964
# ╟─de2baeae-eaf6-4565-9ed1-f2eb8c666839
# ╠═0607672c-9177-4717-8ddf-e07a5dd82ec4
# ╠═04790584-5822-460a-be5c-c9efb3bc26b5
# ╠═780e8faa-e346-45cb-81f1-34df2a99bc17
# ╠═8cbc4ced-7ac2-4def-97cf-ff056c1dcb4a
# ╠═0106756b-594d-4fdb-81b5-cf0739898521
# ╠═4be2d50a-7289-44a5-a4a1-2f736c466f4f
# ╠═d4fb4803-1c1b-4fd7-a782-112777f55be0
# ╠═323abbb8-6f34-4b8b-839c-e0682fed1971
# ╟─8fc7877e-c4e4-40d1-a720-7806f7dbde0a
# ╠═81c4e515-89b5-4ecf-8437-070e5a51cb4c
# ╠═58ac8edc-2432-4054-88d8-52dafe0a2a61
# ╟─fbe4aca2-6a47-4457-98bb-588a5cde0ed5
# ╠═f918dc11-80e4-4223-8f02-3d5f0a10f8e5
# ╟─69e6f136-e31c-4d71-a9a7-ea9ca6530669
# ╟─1753c20f-9b53-4120-a8c8-e2b086f46f44
# ╟─e1ef1e83-c472-4267-8450-38c65f48d3dc
# ╠═b1e64332-95a4-46a5-a45d-457c26e3fc67
# ╟─48029647-f162-459b-8824-fbf652d127f7
# ╟─4116166d-5f82-4d9b-80fb-c8035b9b6ade
# ╟─9a92d4a9-f489-4bba-9361-03cad3ea12e1
# ╟─6a11b8e7-ed7f-4972-a4d5-d713e045ee1c
# ╟─df5b1bfb-ce96-4d32-abdb-a6fcecb195a1
# ╟─e5956bb0-a33a-488d-906e-fb5a7e2473a9
# ╟─fe1e2a72-4772-4482-88da-f9e5f90e928a
# ╟─d94ec33c-3d9d-4d70-b0e1-e3d861a62821
# ╟─333492ec-9016-44c5-9059-e3cb42c05a89
# ╠═9e2b6a47-a113-4107-acdf-3901e0578898
# ╟─be26b92a-14e2-45bc-bb6f-a2664e2e3cd9
# ╟─61be3485-960b-42f8-82e6-71e213a5c9a1
# ╟─def960de-f74a-4ca8-9d95-8af4e0240b60
# ╠═3b41341e-d174-4b2f-8c19-068cb84ba571
# ╠═3f50a881-337f-4578-b0ff-a143439b0a6d
# ╠═89520d6a-7a44-41f6-92ba-3d9416ac2047
# ╠═666c55e5-f7f5-4f83-b3a3-ea6632ca5a86
# ╟─d3493ce8-85d1-4132-b3e0-4ec35ac9d36d
# ╠═f90190cc-555d-47e1-a2cb-99e5d78d4ff5
# ╠═59f05654-a8de-4e17-b2ad-60a7ac64e122
# ╠═dadf76f0-cbea-4c34-a142-41e120679674
# ╠═91242a8c-c09b-402c-a0ea-40b8e3e26ae7
# ╠═3960e081-090e-443f-b82c-06b703686e9f
# ╠═e210999a-7ea5-47f7-aedb-f84fe77bd653
# ╟─bb00b5bb-326e-47f9-a4f4-e7b4f29dd1f2
# ╟─13645622-ed9d-4843-9b83-e67313473c35
# ╟─b7cb5183-65e8-4ee8-af86-2bedd11daecc
# ╠═3a7f5ebe-aed6-4edb-84de-a99de0755453
# ╟─8c367e8f-df43-4f21-bef0-55060f36f44e
# ╟─b14d67ca-5f24-4d8a-9334-6e072e2b39eb
# ╟─deb15672-0e35-4855-b1c0-b2c0b7e78d41
# ╟─b48b4acb-ed25-4d9b-bee6-2e316ecb44c3
# ╟─0a665fc0-1230-4978-8ebd-e551c595e857
# ╟─842b074b-f808-48d8-8dc5-110ddd907f90
# ╟─31298257-d35a-4f6f-8a76-ff00d5361ced
# ╠═72269ec4-a56e-46d9-85c8-0dd8ccaf43e1
# ╠═84d1270b-8df5-4d5d-a153-da4ffdb1d283
# ╠═11b12556-5b61-42c2-a911-4ea98a0a1e85
# ╠═b976ab43-69f1-47a0-b2c6-c63e1c15cdb4
# ╠═60b410be-70f7-4053-a3db-7d777e0d3f08
# ╟─5caca8ea-82af-4999-93bb-a72252c456c7
# ╟─7a02463d-cfd9-4648-af53-f1e65d46733f
# ╠═114d2324-5289-4e44-8d77-736a9bdec365
# ╟─659091d3-60b2-4158-80e2-cd28a492e870
# ╠═c4876d26-e841-4e28-8303-131d4635fc23
# ╠═5dd1a1e6-7db1-479e-a684-accec53ce06a
# ╠═180c12b0-d410-4a5f-97bb-226e6624a39b
# ╠═15fadfc2-3cf8-4fda-9aed-a79c602b1d51
# ╟─f8b5dc8f-1f41-4600-825e-2f9653f2d925
# ╠═afb700c7-ef29-4c13-b9ea-1d40ea9534dd
# ╟─f672a256-641a-478e-b0aa-2df6e68b4d86
# ╠═bab42c91-2d00-463d-a921-97487e4eac67
# ╟─904ac4c2-50a8-4f70-8050-a0a1d4a448fa
# ╠═960ee06e-b15f-4890-912c-851691e5c1a7
# ╠═2c239f3a-6335-4dde-bdcc-7bf41bc49890
# ╠═a7537912-16d8-4312-a1a7-513695ad86de
# ╠═38e147f2-7175-4975-a799-0a2e12c24368
# ╟─c1d2305e-fb8b-4845-a414-08fff84aa9b0
# ╟─a81dd9a4-7938-4a72-b3d2-1780e8ecd536
# ╠═2ce5aa45-4aa5-4c2a-a608-f581266e55f0
# ╟─d5ab1a28-3a60-49d9-bb3e-ca589b1c79fd
# ╠═f8255707-2233-4e28-b542-2f3d81b31c2e
# ╠═c10697b7-6e67-4a5b-937e-09d97ca5b7f8
# ╠═f88ecd40-1b80-4cca-a312-b23f7cfb0ad6
# ╠═754ab149-7a65-4227-ba8b-d7d48e0092b8
# ╠═22244e24-5b56-4933-8a09-44b601f116c3
# ╠═d377af90-9b7d-4fd9-8bc8-d93e58617542
# ╟─686ac3dc-c191-4575-ba0c-d4c2551474b5
# ╠═d1ab199f-1a40-4377-bca3-7f72f3cde3a9
# ╠═32eb1122-5013-4a8e-be54-18a30c151515
# ╟─de144adb-a467-4077-8cb1-d86462f56110
# ╠═d0985ca6-fef5-4b67-9ad6-f51d84b595b4
# ╟─8ae53b8a-0fb3-4c1c-8e5f-a3782a85141c
# ╟─9d7d4d68-c9cc-4a42-a99d-ae25a1ab554c
# ╠═e5fc814f-a8e1-41ef-b81a-c3b0839a2f87
# ╠═6a9fad5b-4964-4e12-b131-8cb2628d1ab3
# ╠═d75725cd-0ef6-421f-be56-f559312e73b6
# ╠═ab0e28f4-4310-4dcc-817e-81b9e45fd501
# ╠═7454f68a-64dc-4676-b2b2-ed8fcb35d81e
# ╠═46d92e15-38ca-4857-8db4-60c1519523f6
# ╠═315dd351-9d68-48f1-aa7a-8f43f3dec6ac
# ╠═06ae600d-3f73-49e6-858c-539079c117ab
# ╠═39683e98-dcbb-458b-817f-856fc6498730
# ╠═3ac837b8-559b-41c2-8f83-1331839dcf7e
