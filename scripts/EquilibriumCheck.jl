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

# ╔═╡ aecc5e8f-1e78-4965-8f9f-4b52d850f490
begin
	using AuCO2RR
	using AuCO2RR: AuCO2RR_plots
end

# ╔═╡ 3ac837b8-559b-41c2-8f83-1331839dcf7e
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
end

# ╔═╡ bd8134d8-5a69-486e-8429-7cf810b3ccbe
Pkg.status()

# ╔═╡ 2176bc34-fc74-4532-912e-e73441b37245
isdefined(AuCO2RR, :AuCO2RR_plots)

# ╔═╡ dc90b463-7574-46e1-a7fe-d60074403747
names(AuCO2RR_plots; all=true)

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
	const vmin = -1.0
	const vmax = 1.1
	
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
		c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol/dm^3 #* 0.0001
		         # (v0/N_A)^(1/3)
		v = N_A * (a * Å)^3 #v0 * (κ + 1)
		a *= Å #8.2 * Å#8.2 * Å   
		M = M0 * v
		BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
	end
	function make_eneutral(bulk_species::Vector{BulkSpecies}
						  ;name, z, D, κ=0.0, a=8.2, v = N_A * (a * Å)^3,
						   #v0 * (κ * abs(z) + 1)
						   M=M0*v, color)
		a *= Å
		c_bulk = -mapreduce(x -> x.c_bulk * x.z, +, bulk_species)/z
		BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
	end
end;

# ╔═╡ c72ac7c7-ff6d-4aca-b6c5-61746bd146a8
use_physical_size

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

# ╔═╡ 53b4dc3e-95f0-4eee-ba1c-68c222638acd
md"""
### Boundary Condition Function
"""

# ╔═╡ 2a20d9be-6c1e-4c1f-8bb6-a7693800732d
md"""
## Double Layer Capacitance
"""

# ╔═╡ 4f7ec19d-cd60-4c2b-a766-7557caa471c0
md""" 
### System Setup
"""

# ╔═╡ 66da15be-e4e3-4592-8354-a05ef092ac86


# ╔═╡ 4656ee04-ae86-442f-b37c-c5563170f992
md"""
Run double layer curve $(@bind double_layer_curve PlutoUI.CheckBox())
"""

# ╔═╡ d76d8413-c019-4728-b182-7f7cb78dede4
md"""
### Result Plots
"""

# ╔═╡ 9598e2c6-521e-4f8d-82d8-a836809736f3
md"""
##### CSVdata
"""

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
#### General CV
"""

# ╔═╡ b95160b5-18f7-49d9-80be-9159abd2dcd1
md"""
Run general cyclic voltammetry curve $(@bind CV PlutoUI.CheckBox())
"""

# ╔═╡ 25eb8aa3-697e-4538-9472-ceea45fbfbd9
md"""
#### pH varied CV function
"""

# ╔═╡ 11892724-1851-46f2-802d-4da45127b0af
md"""
Run pH varied CV calculation $(@bind pH_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ c048e472-3983-4279-bf60-82784baa145e
md"""
#### Scan Rate CV function
"""

# ╔═╡ 3bdaab98-c0f7-46af-86b7-d68374e8a5d0
md"""
Run scan rate varied CV calculation $(@bind scan_rate_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ e0e59ef0-8b6c-4f31-8d39-c2c4bcd7f99e
md"""
#### Pressure CV function
"""

# ╔═╡ 56814250-16b2-4578-820d-2096998c84f4
md"""
Run pressure varied cyclic voltammetry $(@bind pressure_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ a2c7c4da-77cd-493f-8f98-0c86fecf271a
md"""
Run boundary layer thickness varied cyclic voltammetry $(@bind L_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ fbe4aca2-6a47-4457-98bb-588a5cde0ed5
md"""
#### Position at 0.99 Cbulk at a Given Time
"""

# ╔═╡ 1753c20f-9b53-4120-a8c8-e2b086f46f44
md"""
#### RDE_CV plots
"""

# ╔═╡ bb9ba12f-e98b-48de-b55a-d76b276ae952
md"""
Run boundary layer thickness varied cyclic voltammetry $(@bind BL_thickness_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ add42535-4311-4ac6-8f0e-b045426283e7


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

# ╔═╡ fe1e2a72-4772-4482-88da-f9e5f90e928a
md"""
#### Scan Rate CV plots
"""

# ╔═╡ d94ec33c-3d9d-4d70-b0e1-e3d861a62821
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

# ╔═╡ 333492ec-9016-44c5-9059-e3cb42c05a89
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

# ╔═╡ 9e2b6a47-a113-4107-acdf-3901e0578898
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

# ╔═╡ e4d93d39-c391-47ce-a248-6f0205761cca
md"""
Compare the simulation results: $(@bind comp PlutoUI.CheckBox(default=true))
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

# ╔═╡ d912cbca-ef9b-4319-8699-3fc7da8e73d2
md"""
### *User_input*
"""

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

# ╔═╡ ed92cece-3f89-45f5-ac17-cbc9a9abb906
sawtooth

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
			
			- **Use Physical ion size & solvation number**
			  ``κ``, ``a`` : $(Child("use_physical_size", PlutoUI.CheckBox(default=false)))
			  - **ON**  : Apply physical ion size and solvation number.
			  - **OFF** : Use the previous model setting (only potassium has a size).
			"""
        end;
        label = "Submit"
    );
    top = 350
)

# ╔═╡ ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
begin
    use_md_hydrated = user_input_ion.use_physical_size   # Bool toggle you will use

    # a: hydrated radii in water (Å), κ: MD hydration number (1st shell)
    if use_md_hydrated
        # --- hydrated radii [nm] (aqueous effective radii) ---
        a_HCO3 = 3.33   #  (Å)
        a_CO3  = 3.94   #  (Å)
        a_CO2  = 1.70   #  (Å) (often treated as vdW/effective in water)
        a_OH   = 3.00   #  (Å)
        a_H    = 2.80   #  (Å) (H3O+ effective hydrated)
        a_CO   = 1.40   #  (Å)
        a_K    = 3.31   #  (Å)

        # --- MD hydration numbers (1st shell) ---
        κ_HCO3 = 5.4     #  (MD)
        κ_CO3  = 8.5     #  (MD)
        κ_CO2  = 0.0     # neutral, typically treat as ~0 (no structured hydration number in this model)
        κ_OH   = 3.5     #  (MD)
        κ_H    = 4.0     #  (MD)
        κ_CO   = 0.0     # neutral
        κ_K    = 6.0     #   (MD)
    else
        # defalult, Stefan's Model
        a_HCO3 = 0.0;  κ_HCO3 = 0
        a_CO3  = 0.0;  κ_CO3  = 0
        a_CO2  = 0.0;  κ_CO2  = 0
        a_OH   = 0.0;  κ_OH   = 0
        a_H    = 0.0;  κ_H    = 0
        a_CO   = 0.0;  κ_CO   = 0
        a_K    = 8.2;  κ_K    = 0
    end

    const bulk = let
        bulk = [
            BulkSpecies(;name = "HCO₃⁻",
                        z = -1,
                        D = 1.185e-9,
                        c_bulk = 0.091,
                        a = a_HCO3,
                        κ = κ_HCO3,
                        color = :brown
            ),
            BulkSpecies(;name = "CO₃²⁻",
                        z = -2,
                        D = 0.923e-9,
                        c_bulk = 2.68e-6,
                        a = a_CO3,
                        κ = κ_CO3,
                        color = :violet
            ),
            BulkSpecies(;name = "CO₂",
                        z = 0,
                        D = 1.91e-9,
                        c_bulk = 0.033,
                        a = a_CO2,
                        κ = κ_CO2,
                        color = :red
            ),
            BulkSpecies(;name = "OH⁻",
                        z = -1,
                        D = 5.273e-9,
                        c_bulk = 10^(pH-14),
                        a = a_OH,
                        κ = κ_OH,
                        color = :green
            ),
            BulkSpecies(;name = "H⁺",
                        z = 1,
                        D = 9.310e-9,
                        c_bulk = 10^(-pH),
                        a = a_H,
                        κ = κ_H,
                        color = :gray
            ),
            BulkSpecies(;name="CO",
                        z = 0,
                        D = 2.23e-9,
                        c_bulk = 0.0,
                        a = a_CO,
                        κ = κ_CO,
                        color = :blue
            )
        ]

        push!(bulk, make_eneutral(bulk; name="K⁺",
                                       z = 1,
                                       D = 1.957e-9,
                                       a = a_K,
                                       κ = κ_K,
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
	top = 600
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
		σ 			= C_gap * (ϕ_we - ϕ_pzc- u[iϕ])
		local_pH 	= -log10(u[ihplus] * γ[ihplus] / (mol/dm^3))

	
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

# ╔═╡ 9a4e01d9-f469-4427-bf4c-883adb67ae24
function pb_bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data

	if user_input_model.BC_Select == "Dirichlet"
	    ## Dirichlet ϕ=ϕ_we at Γ_we
	    boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
		#boundary_dirichlet!(f, u, bnode; species = ip, region = Γ_we, value = p_we)
	elseif user_input_model.BC_Select == "Robin"
		## Robin ϕ=dϕ₀/dx 
	    boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap , C_gap * (ϕ_we - ϕ_pzc))
	else
		## neumann ϕ=dϕ₀/dx 
	    #boundary_neumann!(f, u, bnode, species = iϕ, region = Γ_we, value = C_gap * (ϕ_we - ϕ_pzc))
	
	end

    return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ 924f8f5d-2cb0-4381-a522-509ff4c002b6
begin
	model_key = user_input_model[:model_choice]
	model = model_key == "Gold_Model" ? elydata_Gold :
	        model_key == "Landstorfer_NaClO₄ model" ? elydata_NaClO₄ :
	        model_key == "Landstorfer_NaF model" ? elydata_NaF :
			model_key == "Toy model" ? elydata_toy :
	        error("Unknown model choice: $model_key")
	molarities = [0.005, 0.1]#, 0.05, 0.1, 0.5] 
end;

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

# ╔═╡ a05cf724-cd32-498e-8afb-ecbf4a1f1648
if scan_rate_varied_checkbox
	sc, saw, scresult = scanrate_varied_sweep(elydata_Gold, sawtooth, grid, pnp_bcondition, reaction; scanrates = [0.01, 0.5, 5], nperiods = user_input_cv.nperiods)
	AuCO2RR_plots.plot_scanrate_sweeps(scresult, sc; species = iohminus)
end

# ╔═╡ 6e88e1d8-1f4b-4813-8890-0cfcdc5fb967
if L_varied_checkbox
	resL = sweep_over_L_cv(
		elydata_Gold,
		[80, 100] .* μm,
		pnp_bcondition,
		sawtooth,
		reaction;
		nperiods = user_input_cv.nperiods
	)
end

# ╔═╡ 9b599280-470e-48a4-bfcd-9f3c3f2e994c
AuCO2RR_plots.plot_cv_over_L(resL, species = iohminus)

# ╔═╡ 84d1270b-8df5-4d5d-a153-da4ffdb1d283
function simulate_CO2R(grid, celldata; voltages = (-1.5:0.1:0.0) * V, kwargs...)
	kwargs 	 	= merge(solver_control, kwargs) 
    cell        = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)
	ivresult    = ivsweep(cell; voltages, store_solutions=true, kwargs...)

	cell, ivresult
end;

# ╔═╡ 60b410be-70f7-4053-a3db-7d777e0d3f08
# ╠═╡ show_logs = false
if Ldependancy
	ivL = AuCO2RR_plots.ivsweep_over_L(elydata_Gold, pnp_bcondition; voltages, solver_control)
end

# ╔═╡ bab42c91-2d00-463d-a921-97487e4eac67
if Ldependancy
	AuCO2RR_plots.plotcurr_over_L(ivL; species=iohminus, cutoff=-0.4, title="IV vs L (log scale)")
end

# ╔═╡ 7fc5e2a3-c217-4042-9a42-e66d547bef96
is_Landstorfer = model != elydata_Gold

# ╔═╡ 084e2127-ea77-4894-8990-380c2e8802c7
if double_layer_curve
	#pb
	sys_pb = PBSystem(grid; bcondition = pb_bcondition,  celldata = deepcopy(model))
	result_pb = capscalc(sys_pb, is_Landstorfer; vrange = range(vmin, vmax, length = 201))
else 
	result_pb = nothing
end

# ╔═╡ 3ef57b7d-ec19-46bc-a881-0506cf5167f3
if double_layer_curve
	#pnp
	reaction_arg = model == elydata_Gold ? (reaction) : NamedTuple()
	sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = deepcopy(model), reaction_arg)

	result_pnp = capscalc(sys_pnp, is_Landstorfer; vrange = range(vmin, vmax, length = 201))
else
	result_pnp = nothing
end

# ╔═╡ 4f991d6d-3a3f-45d8-b2e0-662c5292251c
if double_layer_curve
    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, layout = (1, 2), size = 	(650, 350))
    AuCO2RR_plots.capsplot(vis[1, 1], result_pb, "Poisson-Boltzmann"; nshow = length(result_pb[1].dlcaps), xlimits_NL = (vmin*1.1, vmax*1.1))
    AuCO2RR_plots.capsplot(vis[1, 2], result_pnp, "Poisson-Nernst-Planck"; nshow = length(result_pnp[1].dlcaps), xlimits_NL = (vmin*1.1, vmax*1.1))

    reveal(vis)
end

# ╔═╡ 2640ee7f-109d-4dcc-b8af-3f854da1a323
AuCO2RR_plots.plot_caps_comparison(; model_key, result_pb, result_pnp)

# ╔═╡ 02d12ba4-4ab3-48f6-b084-edb06cb413b1
begin
	celldata = deepcopy(model)
	pnpcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = celldata, reaction = reaction)
end

# ╔═╡ e50fe651-11d4-45ee-89dd-371a7fbc097e
function sweep(pnpdata; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    #celldata.eneutral = eneutral
	reaction_arg = model == elydata_Gold ? (; reaction) : NamedTuple()
    pnpcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = celldata, reaction = reaction)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
    )

end

# ╔═╡ 7b38e59a-d005-4cfc-ba8c-b17e7c700119
if pressure_varied_checkbox
	P_recs = pressure_varied_sweep(elydata_Gold, sweep; Pvec = [0.1, 0.2, 0.3, 0.5, 0.6, 1], ispec = 5)
	AuCO2RR_plots.plot_pressure_varied_sweep(P_recs, species = iohminus)
end

# ╔═╡ c9ae7a67-9a5f-4d9a-88c0-742f4e91fb27
	AuCO2RR_plots.plot_pressure_varied_sweep(P_recs, species = iohminus, limits = ((-0.2, 0.9), (0, 1)))


# ╔═╡ 9fb47b83-a853-4316-bb8d-30e65b16ef78
if CV
	pnpresult = sweep(model; eneutral = false, tunnel = false)
end

# ╔═╡ 7da046bf-d3b1-43a0-bdba-89b4da2f6be3
if CV
	AuCO2RR_plots.plot_time_voltage_and_dt(pnpresult, sawtooth)
end

# ╔═╡ c62ab378-0988-4fa5-b21d-5e1622c63c87
if CV
	AuCO2RR_plots.plot_cv_current(pnpresult, elydata_Gold; species = ico)
end

# ╔═╡ a64e2dc9-9be7-48b5-9d04-c448f19ed7f2
if CV
	AuCO2RR_plots.plot_conc_time_electrode(pnpresult, bulk)
end

# ╔═╡ 2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
if CV
	AuCO2RR_plots.plot_conc_profile_logx(pnpresult, bulk, X, 13)
end

# ╔═╡ 2420382d-227a-4063-9450-1f1726df018e
if CV
	path = AuCO2RR_plots.cv_conc_gif(pnpresult, bulk, X)
	LocalResource(path)
end

# ╔═╡ 3d661549-a8d2-40b0-add8-b186193f90fe
if CV
	AuCO2RR_plots.plot_cv_model_vs_koper_facets(pnpresult; species = iohminus)
end

# ╔═╡ 74c43d72-3a23-4a24-a4ae-8b18b245610a
if CV
	AuCO2RR_plots.plot_iv_with_experiment(pnpresult, ico)
end

# ╔═╡ b25f2246-0182-4d50-a606-0d81776d414f
AuCO2RR_plots.plot_conc_profile_with_delta(pnpresult, bulk, X)

# ╔═╡ 8cd25c0c-e260-4401-af12-a1def38bb7c2
if pH_varied_checkbox
	results_pH = run_pH_sweep(model, sawtooth, grid, pnp_bcondition, reaction; pH_values = [3, 10])
	AuCO2RR_plots.plot_pH_varied_sweep(results_pH; species = ico)
end

# ╔═╡ 11b12556-5b61-42c2-a911-4ea98a0a1e85
# ╠═╡ show_logs = false
cell, ivresult = simulate_CO2R(grid, model)

# ╔═╡ 58ac8edc-2432-4054-88d8-52dafe0a2a61
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

# ╔═╡ 5ccb6a73-41cc-4107-a6ee-37d906313841
(~, default_index) = findmin(abs, ivresult.voltages .+ 0.9 * ufac"V");

# ╔═╡ af333d3b-1e3a-4227-8cce-479907c11448
begin
	gifpath = AuCO2RR_plots.plot1d_movie(ivresult; bulk, grid, L, step = 3, framerate = 5)
	LocalResource(gifpath)
end

# ╔═╡ 9498845e-fa44-4d01-a7bc-33d01ec11f79
AuCO2RR_plots.plot_iv_with_ringe_refs(ivresult, species = iohminus)

# ╔═╡ 22244e24-5b56-4933-8a09-44b601f116c3
AuCO2RR_plots.iv_curve_axis(ivresult; cutoff=-0.4, showlegend=true, species = iohminus)

# ╔═╡ c6f10b66-6d06-4f2e-a7cc-780096d75785
begin 
	conc_f, conc_a = AuCO2RR_plots.conc_vs_voltage_axis_compare(ivresult; bulk, grid, compare = comp)
	conc_f
end

# ╔═╡ f8255707-2233-4e28-b542-2f3d81b31c2e
begin
    conc_out = AuCO2RR_plots.conc_vs_voltage_axis(
        ivresult;
        bulk = bulk,
        grid = grid,
        useonly_pH = false,
        showlegend = false
    )

    conc_fig = conc_out.fig
    conc_ax  = conc_out.ax
    conc_colors = conc_out.colors
    conc_species = conc_out.species

    text!(conc_ax, -1.15, 0.3, text=L"\mathrm{K^+}",
		  color=conc_colors[ikplus], fontsize=24, font="sans-bold")
	text!(conc_ax, -0.90, 0.0000002, text=L"\mathrm{H^+}", 
		  color=conc_colors[ihplus], fontsize=24, font = "sans-bold") 
	text!(conc_ax, -1.05, 5e-11, text=L"\mathrm{CO_3^{2-}}",
		  color=conc_colors[ico3], fontsize=24, font = "sans-bold") 
	text!(conc_ax, -1.0, 2.5e-7, text=L"\mathrm{HCO_3^-}", 
		  color=conc_colors[ihco3], fontsize=24, font = "sans-bold") 
	text!(conc_ax, -1.2, 5.2e-6, text=L"\mathrm{CO_2}", 
		  color=conc_colors[ico2], fontsize=24, font = "sans-bold")
	text!(conc_ax, -0.90, 8e-10, text=L"\mathrm{OH^-}",
		  color=conc_colors[iohminus], fontsize=24, font = "sans-bold") 
	text!(conc_ax, -1.15, 0.00024, text=L"\mathrm{CO}", 
	 	  color=conc_colors[ico], fontsize=24, font = "sans-bold")
    conc_fig
end

# ╔═╡ 315dd351-9d68-48f1-aa7a-8f43f3dec6ac
floataside(
    md"""
    __Input Voltage Index:__ $(@bind vindex PlutoUI.Slider(1:5:length(ivresult.voltages), default=default_index))
    """,
    top = 960
)


# ╔═╡ c4876d26-e841-4e28-8303-131d4635fc23
md"""
Potential at the working electrode 
$(vshow = ivresult.voltages[vindex]; @sprintf("%+1.4f", vshow))
"""

# ╔═╡ 15fadfc2-3cf8-4fda-9aed-a79c602b1d51
AuCO2RR_plots.plot1d(ivresult, vshow; bulk, grid, L)

# ╔═╡ ab302d08-1a6f-4553-85af-043c565b107f
begin
	cond_out = AuCO2RR_plots.plot1d_makie(ivresult, vshow; bulk=bulk, grid=grid, L=L)
	cond_fig, cond_ax, cond_colors = cond_out.fig, cond_out.ax, cond_out.colors
	
	#text!(cond_ax, 4e-11, -2.4, text=L"\mathrm{CO_2}", color=cond_colors[ico2], fontsize=24, font="sans-bold")
	cond_fig
end

# ╔═╡ 06ae600d-3f73-49e6-858c-539079c117ab
floataside(
    md"""
    __Input Time Index:__ $(@bind it PlutoUI.Slider(1:length(pnpresult.tsol.t)-1, show_value=false))
    """,
    top = 1060
)

# ╔═╡ 7454f68a-64dc-4676-b2b2-ed8fcb35d81e
floataside(
	md"""
	**Time:** $(round(pnpresult.tsol.t[it+1], digits=3)) s
	""",
	top = 1010
)

# ╔═╡ 39683e98-dcbb-458b-817f-856fc6498730
floataside(
    md"""
    **Input Voltage Index:** $(@bind vindex2 PlutoUI.Slider(1:5:length(target_time), default=40))
    """,
    top = 1110
)

# ╔═╡ Cell order:
# ╠═91ac9e35-71eb-4570-bef7-f63c67ce3881
# ╠═aecc5e8f-1e78-4965-8f9f-4b52d850f490
# ╠═bd8134d8-5a69-486e-8429-7cf810b3ccbe
# ╠═2176bc34-fc74-4532-912e-e73441b37245
# ╠═dc90b463-7574-46e1-a7fe-d60074403747
# ╟─beae1479-1c0f-4a55-86e1-ad2b50174c83
# ╟─ab2184fc-0279-46d9-9ee4-88fe3e732789
# ╠═7316901c-d85d-48e9-87dc-3614ab3d81a5
# ╟─6b7cfe87-8190-40a5-8d25-e39ef8d55db5
# ╠═5a146a44-03dc-45f3-ae15-993d11c2edac
# ╠═00947475-c96e-4ecc-a1ef-5be5e3e3c864
# ╠═c72ac7c7-ff6d-4aca-b6c5-61746bd146a8
# ╠═ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
# ╟─06f52599-7006-4a5c-ba86-0b668b6952c9
# ╟─4b64e168-5fe9-4202-9657-0d4afc237ddc
# ╟─de2c826d-6c05-47cf-b5f5-44a00ea9889c
# ╟─d8f00649-e2ed-4bdd-853f-05268f0d5353
# ╟─47b36c81-b57e-4dd0-a22f-999e4fd3ac9f
# ╟─1e877f17-0219-45f1-b640-3a25ae085dbd
# ╠═8a1047fa-e483-40d9-8904-7576f30acfb4
# ╟─8912f990-6b02-467a-bd11-92f94818b1c7
# ╟─a8157cc1-1761-4b11-a37c-9e12a9ca695e
# ╟─6b5cf93c-0df3-4a18-8786-502361736838
# ╟─d2c0642d-dfa5-4a76-bd36-ac4a735a3299
# ╟─06d45088-ab8b-4e5d-931d-b58701bf8464
# ╠═91113083-d80e-4528-be41-82d10f6860fc
# ╟─d0093605-0e35-4888-a93c-8456c698e6f0
# ╟─f0b5d356-6b97-4878-98de-bee5f380d41a
# ╟─e3eda42f-e2f3-4c10-81c4-610246ca528d
# ╟─7b87aa2a-dbaf-441c-9ad7-444abf15f664
# ╠═2b9d9bfd-d660-4b4d-8f0b-b6b5bcc0dbfa
# ╟─6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
# ╟─e510bce3-d33f-47bb-98d6-121eee8f2252
# ╟─848b7aeb-968f-4116-8038-b61276f02b6c
# ╟─f18dc873-1c9d-46d3-9596-92d28705e894
# ╟─12235c3c-18f2-4fc7-95ef-800f71783036
# ╟─53f12821-7d8d-4971-87fd-ad4689ec62a5
# ╟─9d814b85-a5b6-42e5-abf4-15500bbdb717
# ╟─e1e0ca0f-7f88-40f0-850e-590b25da0331
# ╟─0db74a70-af86-492c-affb-9de62ffe4455
# ╟─4f388fe0-6bc8-4a29-bccc-fa725e62c6a7
# ╠═5d179c52-43d7-4bcb-a2df-93c5806876fa
# ╟─161a810d-c05e-42ad-97ab-131059d6784a
# ╟─5f17b4f7-54d6-4ad0-9886-252854840a80
# ╟─53b4dc3e-95f0-4eee-ba1c-68c222638acd
# ╠═dc203e95-7763-4b13-8408-038b933c5c9c
# ╠═9a4e01d9-f469-4427-bf4c-883adb67ae24
# ╟─2a20d9be-6c1e-4c1f-8bb6-a7693800732d
# ╟─4f7ec19d-cd60-4c2b-a766-7557caa471c0
# ╠═924f8f5d-2cb0-4381-a522-509ff4c002b6
# ╠═66da15be-e4e3-4592-8354-a05ef092ac86
# ╠═084e2127-ea77-4894-8990-380c2e8802c7
# ╠═3ef57b7d-ec19-46bc-a881-0506cf5167f3
# ╟─4656ee04-ae86-442f-b37c-c5563170f992
# ╟─d76d8413-c019-4728-b182-7f7cb78dede4
# ╠═7fc5e2a3-c217-4042-9a42-e66d547bef96
# ╠═4f991d6d-3a3f-45d8-b2e0-662c5292251c
# ╠═2640ee7f-109d-4dcc-b8af-3f854da1a323
# ╟─9598e2c6-521e-4f8d-82d8-a836809736f3
# ╠═8bfdf2f5-c80a-4ce0-a8e1-b315affffb5f
# ╠═ed92cece-3f89-45f5-ac17-cbc9a9abb906
# ╟─5c808c71-6094-49d7-8215-e88262f34e1f
# ╟─da8390d1-47e8-451f-b12b-45b8aca7b6ec
# ╠═02d12ba4-4ab3-48f6-b084-edb06cb413b1
# ╠═b4aaf070-d4ab-409a-b1e8-f5469b9f398b
# ╠═e50fe651-11d4-45ee-89dd-371a7fbc097e
# ╟─ef7212fc-a3d0-4784-b901-219204b79dc0
# ╟─b95160b5-18f7-49d9-80be-9159abd2dcd1
# ╠═9fb47b83-a853-4316-bb8d-30e65b16ef78
# ╠═7da046bf-d3b1-43a0-bdba-89b4da2f6be3
# ╠═c62ab378-0988-4fa5-b21d-5e1622c63c87
# ╠═a64e2dc9-9be7-48b5-9d04-c448f19ed7f2
# ╠═2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
# ╠═2420382d-227a-4063-9450-1f1726df018e
# ╠═3d661549-a8d2-40b0-add8-b186193f90fe
# ╠═74c43d72-3a23-4a24-a4ae-8b18b245610a
# ╟─25eb8aa3-697e-4538-9472-ceea45fbfbd9
# ╟─11892724-1851-46f2-802d-4da45127b0af
# ╠═8cd25c0c-e260-4401-af12-a1def38bb7c2
# ╟─c048e472-3983-4279-bf60-82784baa145e
# ╟─3bdaab98-c0f7-46af-86b7-d68374e8a5d0
# ╠═a05cf724-cd32-498e-8afb-ecbf4a1f1648
# ╟─e0e59ef0-8b6c-4f31-8d39-c2c4bcd7f99e
# ╠═56814250-16b2-4578-820d-2096998c84f4
# ╠═7b38e59a-d005-4cfc-ba8c-b17e7c700119
# ╠═c9ae7a67-9a5f-4d9a-88c0-742f4e91fb27
# ╟─58ac8edc-2432-4054-88d8-52dafe0a2a61
# ╠═a2c7c4da-77cd-493f-8f98-0c86fecf271a
# ╠═6e88e1d8-1f4b-4813-8890-0cfcdc5fb967
# ╠═9b599280-470e-48a4-bfcd-9f3c3f2e994c
# ╟─fbe4aca2-6a47-4457-98bb-588a5cde0ed5
# ╠═b25f2246-0182-4d50-a606-0d81776d414f
# ╟─1753c20f-9b53-4120-a8c8-e2b086f46f44
# ╟─bb9ba12f-e98b-48de-b55a-d76b276ae952
# ╠═add42535-4311-4ac6-8f0e-b045426283e7
# ╠═e1ef1e83-c472-4267-8450-38c65f48d3dc
# ╠═b1e64332-95a4-46a5-a45d-457c26e3fc67
# ╠═48029647-f162-459b-8824-fbf652d127f7
# ╠═4116166d-5f82-4d9b-80fb-c8035b9b6ade
# ╟─fe1e2a72-4772-4482-88da-f9e5f90e928a
# ╟─d94ec33c-3d9d-4d70-b0e1-e3d861a62821
# ╟─333492ec-9016-44c5-9059-e3cb42c05a89
# ╟─9e2b6a47-a113-4107-acdf-3901e0578898
# ╟─be26b92a-14e2-45bc-bb6f-a2664e2e3cd9
# ╟─61be3485-960b-42f8-82e6-71e213a5c9a1
# ╟─def960de-f74a-4ca8-9d95-8af4e0240b60
# ╠═3b41341e-d174-4b2f-8c19-068cb84ba571
# ╠═3f50a881-337f-4578-b0ff-a143439b0a6d
# ╠═89520d6a-7a44-41f6-92ba-3d9416ac2047
# ╠═666c55e5-f7f5-4f83-b3a3-ea6632ca5a86
# ╟─d3493ce8-85d1-4132-b3e0-4ec35ac9d36d
# ╟─f90190cc-555d-47e1-a2cb-99e5d78d4ff5
# ╟─59f05654-a8de-4e17-b2ad-60a7ac64e122
# ╠═dadf76f0-cbea-4c34-a142-41e120679674
# ╠═91242a8c-c09b-402c-a0ea-40b8e3e26ae7
# ╟─bb00b5bb-326e-47f9-a4f4-e7b4f29dd1f2
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
# ╟─7a02463d-cfd9-4648-af53-f1e65d46733f
# ╟─114d2324-5289-4e44-8d77-736a9bdec365
# ╠═5ccb6a73-41cc-4107-a6ee-37d906313841
# ╟─c4876d26-e841-4e28-8303-131d4635fc23
# ╠═15fadfc2-3cf8-4fda-9aed-a79c602b1d51
# ╠═af333d3b-1e3a-4227-8cce-479907c11448
# ╠═ab302d08-1a6f-4553-85af-043c565b107f
# ╟─f8b5dc8f-1f41-4600-825e-2f9653f2d925
# ╠═9498845e-fa44-4d01-a7bc-33d01ec11f79
# ╠═22244e24-5b56-4933-8a09-44b601f116c3
# ╟─f672a256-641a-478e-b0aa-2df6e68b4d86
# ╠═bab42c91-2d00-463d-a921-97487e4eac67
# ╟─904ac4c2-50a8-4f70-8050-a0a1d4a448fa
# ╠═e4d93d39-c391-47ce-a248-6f0205761cca
# ╠═c6f10b66-6d06-4f2e-a7cc-780096d75785
# ╠═f8255707-2233-4e28-b542-2f3d81b31c2e
# ╟─de144adb-a467-4077-8cb1-d86462f56110
# ╠═d0985ca6-fef5-4b67-9ad6-f51d84b595b4
# ╟─8ae53b8a-0fb3-4c1c-8e5f-a3782a85141c
# ╟─9d7d4d68-c9cc-4a42-a99d-ae25a1ab554c
# ╟─d912cbca-ef9b-4319-8699-3fc7da8e73d2
# ╠═e5fc814f-a8e1-41ef-b81a-c3b0839a2f87
# ╠═6a9fad5b-4964-4e12-b131-8cb2628d1ab3
# ╠═d75725cd-0ef6-421f-be56-f559312e73b6
# ╠═ab0e28f4-4310-4dcc-817e-81b9e45fd501
# ╠═7454f68a-64dc-4676-b2b2-ed8fcb35d81e
# ╠═315dd351-9d68-48f1-aa7a-8f43f3dec6ac
# ╠═06ae600d-3f73-49e6-858c-539079c117ab
# ╠═39683e98-dcbb-458b-817f-856fc6498730
# ╟─3ac837b8-559b-41c2-8f83-1331839dcf7e
