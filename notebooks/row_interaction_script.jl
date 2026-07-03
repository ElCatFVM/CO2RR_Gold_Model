### A Pluto.jl notebook ###
# v1.0.1

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
	using Makie: to_font
	sans_font = Makie.to_font("DejaVu Sans") 
end;

# ╔═╡ aecc5e8f-1e78-4965-8f9f-4b52d850f490
begin
	using AuCO2RR
	include(joinpath(pkgdir(AuCO2RR), "plots", "AuCO2RR_plots.jl"))
	using .AuCO2RR_plots
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

# ╔═╡ f7d13047-4007-47ac-a3bb-a0b788dcd141
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
@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;

# ╔═╡ 6b7cfe87-8190-40a5-8d25-e39ef8d55db5
md"""
### Data
"""

# ╔═╡ 5a146a44-03dc-45f3-ae15-993d11c2edac
begin
	@phconstants N_A c_0 k_B e h ε_0 R
	const F = N_A * e

	const voltages = (-1.15:0.1:-0.0) * V
	const vmin = -1.0
	const vmax = 1.1
	
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
	const c̄    = 55.508mol / dm^3

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
	const c̄ 		= 55.508mol / dm^3


	
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

# ╔═╡ a4b1300f-7e8a-46ab-9efd-dba82315a966
md"""
### FVM geometry generation
"""

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

# ╔═╡ a4928999-96c7-4145-af2d-485ff13d7109
rn

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

# ╔═╡ ee9229f8-2f48-4c1f-88cf-4bcc40f6e5f6
paramsidx

# ╔═╡ fb8a9a17-ed12-4f87-9957-06e5e265fcb2
function calc_QBL_local(u, data; tolϕ = 1e-12)
    (; ip, iϕ, ε_0, pscale, ε) = data

    Δp = u[ip]
    Δϕ = u[iϕ]

    if abs(Δϕ) < tolϕ
        return zero(Δϕ)
    end
    arg = 2 * (1.0 .+ ε) * ε_0 * Δp * pscale

    if arg < 0
        error("Invalid local QBL state: Δp=$(Δp), Δϕ=$(Δϕ)")
    end

    return sign(Δϕ) * sqrt(arg)
end

# ╔═╡ 960b1e96-da59-44c1-9828-929ece1a2955
function surface_charge(u, data, boundary)
	if boundary == "Robin"
  		return C_gap * (data.ϕ_we - ϕ_pzc - u[data.iϕ])
	else
		return calc_QBL_local(u, data)
	end
end

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

# ╔═╡ dabb7ca0-ce88-47d6-9315-9192d87c244e
md"""
#### Electrolyte and Cell Data Gold model
"""

# ╔═╡ 957d729c-d310-4d67-a429-1feb7770a2fd
ElectrolyteDataElectrolyteData

# ╔═╡ 20324e33-26d0-4a95-bda3-5cb5d0fd8975
md"""
#### Electrolyte and Cell Data NaClO₄ model
"""

# ╔═╡ 848b7aeb-968f-4116-8038-b61276f02b6c
elydata_NaClO₄ = ElectrolyteData(
 		z = [-1, 1],
		κ = [15.0, 25.0],
		c_bulk = [1.0, 1.0],
		ε = 26.0,
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

# ╔═╡ 53b4dc3e-95f0-4eee-ba1c-68c222638acd
md"""
### Boundary Condition Function
"""

# ╔═╡ d8c01196-7b81-43a5-ac0c-6dac8efee3d8
md"""
#### Dirichlet Boundary Condition Function - First version
"""

# ╔═╡ 95e72a20-a621-48a6-947e-0dcf9375facf
const bcvals = [0.0, 0.0]

# ╔═╡ 5b900250-56a4-4e0b-bbc7-37bd39456312
md"""
```math
\begin{aligned}
   y_i&=  \frac{c_{i,bulk}}{\bar c}\exp\left(z_i\phi \frac{F}{RT}\right)
\end{aligned}
```
"""

# ╔═╡ 83ea329a-5d40-401d-8db1-59618c909f33
function dlcapsweep(system; nsweep = 100, δ = 1.0e-6)
    voltages = range(0, vmax, length = nsweep + 1)
    solutions = []
    dlcaps = zeros(0)
    for v in voltages
        bcvals[1] = v
        if length(solutions) == 0
            inival = VoronoiFVM.unknowns(system)
            inival .= 0
        else
            inival = solutions[end]
        end
        push!(solutions, solve(system; inival))
        # Here, we can directly access the double layer charge as
        # solution component
        Q = solutions[end][iQ, 1]
        bcvals[1] = v + δ
        solδ = solve(system; inival = solutions[end])
        Qδ = solδ[iQ, 1]
        push!(dlcaps, (Qδ[1] - Q[1]) / δ)
    end
    return voltages, dlcaps, TransientSolution(solutions, voltages)
end

# ╔═╡ dc5b1cfc-9221-45e9-a4f3-69ec4f7ebb70
md"""
---
"""

# ╔═╡ 2a20d9be-6c1e-4c1f-8bb6-a7693800732d
md"""
## Double Layer Capacitance
"""

# ╔═╡ 4f7ec19d-cd60-4c2b-a766-7557caa471c0
md""" 
### System Setup
"""

# ╔═╡ d76d8413-c019-4728-b182-7f7cb78dede4
md"""
### Result Plots
"""

# ╔═╡ 4656ee04-ae86-442f-b37c-c5563170f992
md"""
##### **Run double layer curve:** $(@bind double_layer_curve PlutoUI.CheckBox())
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
	
	Valetter_NaClO₄_100mM = CSV.read("../data/Valette_Cap_data/NaClO4_0.1M.csv", DataFrame);	
	Valette_NaClO₄_20mM = CSV.read("../data/Valette_Cap_data/NaClO4_0.02M.csv", DataFrame);
	Valetter_NaClO₄_5mM = CSV.read("../data/Valette_Cap_data/NaClO4_0.005M.csv", DataFrame);	
	
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

# ╔═╡ e50fe651-11d4-45ee-89dd-371a7fbc097e
function sweep(model, grid, bcondition, reaction, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
    celldata = deepcopy(model)
    #celldata.eneutral = eneutral
    pnpcell = PNPSystem(grid; bcondition = bcondition, celldata = celldata, reaction = reaction)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
		kwargs...
    )

end

# ╔═╡ ef7212fc-a3d0-4784-b901-219204b79dc0
md"""
#### General CV
"""

# ╔═╡ b95160b5-18f7-49d9-80be-9159abd2dcd1
md"""
##### **Run general cyclic voltammetry curve:** $(@bind CV PlutoUI.CheckBox())
"""

# ╔═╡ ad55a4c1-78b9-40d2-aec8-64ed95174e8a
md"""
Click the button below to save the Profile data as CSV files.

$(@bind export_button PlutoUI.Button("Export to CSV"))
"""

# ╔═╡ f6f26f7e-b97e-4f83-8b02-d3ff8b8ad14d
cv_sol_control = (; damp_initial = 0.1,    
				    damp_growth  = 1.1,    
			        Δu_opt = 0.03,
			        Δt_min = 3.0e-3,
			        Δt_max = 3.0e-2,
			        Δt = 3.0e-3,
			        Δt_grow = 1.2,
				 )

# ╔═╡ f519311e-d5e7-4a4d-a080-80eb92e99ea4
md"""
##### **Run `unc` curve:** $(@bind unc PlutoUI.CheckBox())
"""

# ╔═╡ 7a172418-b46a-459c-9b17-ab3f4eec8dea
md"""
##### **Run `irc` curve:** $(@bind irc PlutoUI.CheckBox())
"""

# ╔═╡ 7ff115dc-e982-4109-8db6-0bd718b13e88
md"""
##### **Run `odr` curve:** $(@bind odr PlutoUI.CheckBox())
"""

# ╔═╡ 41158868-8680-466b-a94b-9ffcd2d0ad8e
# ╠═╡ disabled = true
#=╠═╡
function cvsweep_backup_v2(
        sys, celldata;
        inival = 0,
        scanrate = 1ufac"V / s",
        vmin = -1.5ufac"V",
        vmax = 0.5ufac"V",
        nperiods = 1,
        Δt_min = 1.0e-11,
        Δu_opt = 2.0e-2,
        damp_initial = 0.5,
        kwargs...
    )
    sawtooth = SawTooth(; scanrate, vmin, vmax)
    times = [0, 2 * (vmax - vmin)] * nperiods / scanrate
    celldata.mutables.sawtooth = sawtooth
    grid = sys.vfvmsys.grid
    X = grid[XCoordinates]

    (; BP1, BP2, BP3, BP4, ε, ε0, iϕ, mutables) = celldata
    
    local_iDL = mutables.iDL
    if local_iDL == 0
        bfaceregions = grid[BFaceRegions]
        bfacenodes = grid[BFaceNodes]
        for ibf in 1:length(bfaceregions)
            if bfaceregions[ibf] == BP2  # match boundary flag marking start of Diffuse Layer
                local_iDL = bfacenodes[1, ibf]
                break
            end
        end
    end
    if local_iDL == 0
        local_iDL = 2
    end

    tff = TestFunctionFactory(sys.vfvmsys)
    tf = testfunction(tff, [BP1, BP2, BP4], [BP3])

    cvresult = CVResult()

    function pre_callback(sol, t)
        sys.vfvmsys.physics.data.ϕ_we = sawtooth(t)
    end

    function j_F(sol, t, celldata)
        (; F) = celldata
        I_boundary = VoronoiFVM.integrate(sys.vfvmsys, sys.vfvmsys.physics.breaction, sol; boundary = true)
        net_co_flux = I_boundary[ico, Γ_we]
        j_faradaic = -2 * F * net_co_flux
        return j_faradaic
    end
	
	function calc_local_surface_charge(sol_matrix, t_val)
        v_we = sawtooth(t_val)
        v_surface = sol_matrix[iϕ, 1] 
        return C_gap * (v_we - ϕ_pzc - v_surface)
    end
	
    function post(sol, oldsol, t, Δt)


		sigma_current = calc_local_surface_charge(sol, t)
        sigma_old     = calc_local_surface_charge(oldsol, t - Δt)
        j_C  		  = (sigma_current - sigma_old) / Δt

		
        push!(cvresult.i_F, j_F(sol, t, celldata))
        push!(cvresult.i_C, j_C)
        push!(cvresult.simtimes, t)
        push!(cvresult.sawtooth, sawtooth(t))
        push!(cvresult.volts, sol[iϕ, 1])
        push!(cvresult.dlvolts, sol[iϕ, 1] - sol[iϕ, local_iDL])
        push!(cvresult.Ωdrop, sol[iϕ, end] - sol[iϕ, local_iDL])
        return
    end

    function delta(sys_fvm, u, uold, t, Δt)
        return norm(u[iϕ, :] - uold[iϕ, :], Inf)
    end

    tsol = VoronoiFVM.solve(
        sys.vfvmsys;
        inival = inival,
        times = times,
        pre = pre_callback,  
        post = post,
        delta = delta,
        Δt_min = Δt_min,                        
        Δt = 1.0e-4 / scanrate,                 
        Δt_max = 5.0e-2 / scanrate,             
        damp_initial = damp_initial,            
        damp_growth = 1.1,                    
        kwargs...
    )

    cvresult.ceflow = VoronoiFVM.integrate(sys.vfvmsys, tf, tsol)

    return tsol, cvresult
end;
  ╠═╡ =#

# ╔═╡ e7e678cb-0573-4c1d-8d04-978887c0b900


# ╔═╡ 3b8bc779-69f8-4f78-b55a-df90e226bcc9
#plot_7species_contours(pnpresult_unc, X, bulk)

# ╔═╡ 22155420-3676-45f7-9bc8-8e00dfaf36be


# ╔═╡ fd5c53ac-1a27-40f9-8ae4-cfbf1aaba5d1
md"""
### CV unc`:none` Result
"""

# ╔═╡ 6f90e65d-6804-4816-9ade-677a5984e4a8
md"""
### CV irc`:pseudopotentiostat` Result
"""

# ╔═╡ 36b271e4-23d8-4752-aca2-e5294a7ad419
md"""
### CV odr`:ohmicdrop` Result
"""

# ╔═╡ 000ab285-b1bf-418b-b70c-40ff945cf2b8
md"""
### Comparesion between three mods
"""

# ╔═╡ cfb28aaa-b9cf-4906-90a8-9ca65cbc1ab1


# ╔═╡ ef1195e2-6677-4cd4-bb4b-c66332c4cfee


# ╔═╡ 9d0b6083-e49b-4b15-ad2e-a2c74f689c5e
md"""
### Comparesion with Experimental Result
"""

# ╔═╡ c048e472-3983-4279-bf60-82784baa145e
md"""
#### Scan Rate CV function
"""

# ╔═╡ 3bdaab98-c0f7-46af-86b7-d68374e8a5d0
md"""
Run scan rate varied CV calculation $(@bind scan_rate_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 6d119626-ded5-4282-bc5a-0c37985697b0
sawtooth_sr = SawTooth(vmin = -1.2, vmax = 1.2, scanup=false)

# ╔═╡ 4d398e3f-1589-4ace-b0cb-273418a36fa2
begin
	hmin_sr = 1.0e-6 	* μm
	hmax_sr = 80 * 0.02 * μm
	X_sr = ExtendableGrids.geomspace(0, 80 * μm, hmin_sr, hmax_sr)
	grid_sr = ExtendableGrids.simplexgrid(X_sr)
end

# ╔═╡ 70738a47-6476-4a87-b8a3-c9622519a1d3
	#sc_2, saw_2, scresult_2 = scanrate_varied_sweep(elydata_Gold_unc, sawtooth_sr, grid, pnp_bcondition, reaction; scanrates = [0.05, 0.5, 5], nperiods = user_input_cv.nperiods)


# ╔═╡ e0e59ef0-8b6c-4f31-8d39-c2c4bcd7f99e
md"""
#### Pressure CV function
"""

# ╔═╡ 56814250-16b2-4578-820d-2096998c84f4
md"""
Run pressure varied cyclic voltammetry $(@bind pressure_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 7b38e59a-d005-4cfc-ba8c-b17e7c700119
# ╠═╡ disabled = true
#=╠═╡
if pressure_varied_checkbox
	#P_recs_unc = pressure_varied_sweep(elydata_Gold_unc, grid, pnp_bcondition, reaction, sawtooth; Pvec = [0.1, 0.2, 0.3, 0.5, 0.6, 1], ispec = ico2)
end
  ╠═╡ =#

# ╔═╡ b48ca0c7-acfc-4d63-a3a3-ec7d46c331d1
# ╠═╡ disabled = true
#=╠═╡
if pressure_varied_checkbox
	#P_recs_irc = pressure_varied_sweep(elydata_Gold_irc, grid, pnp_bcondition, reaction, sawtooth; Pvec = [0.1, 0.2, 0.3, 0.5, 0.6, 1], ispec = ico2)
end
  ╠═╡ =#

# ╔═╡ 9c38b6f2-8f1b-49eb-bcdc-e3c792deed34
#P_recs_odr = pressure_varied_sweep(elydata_Gold_odr, grid, pnp_bcondition, reaction, sawtooth; Pvec = [0.1, 0.5, 1], ispec = ico2)

# ╔═╡ 7844c654-bf05-4b19-9556-c0f3186efded
md"""
#### Compare with Experimental Result
"""

# ╔═╡ 5a1a5d5f-821d-47b8-b045-b6554313297c
md"""
#### Check _iR_Compensation
"""

# ╔═╡ f5f9f81c-4f58-4a9b-9a6c-162d53ebb8f0
md"""
Run boundary layer thickness varied cyclic voltammetry $(@bind _iR_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ a56f239e-e5b4-4689-9d56-719241924b19
#sweep_over_L_cv(elydata_Gold_odr, [80, 150, 500, 1000] * μm , pnp_bcondition, sawtooth, reaction;  nperiods=2, eneutral=true)

# ╔═╡ 53315f32-2fef-43c7-9baf-6d33bb848c70
function plot_pressure_varied_sweep(
    P_recs;
    species=iohminus,
    fig_size=(1600, 900),
    scale=cm^2/mA,
    limits=nothing,
    # --- add experimental background (Figure_3.csv) ---
    fig3_csv::Union{Nothing,AbstractString}="../data/Langmuir_CV_data/Figure_3.csv",
    exp_title::AbstractString="Experimental",
    exp_alpha::Real=0.55,
    exp_linewidth::Real=2,
    sim_linewidth::Real=3,
)
    fig = Figure(size = fig_size)

    ax = if limits !== nothing
        Axis(fig[1, 1],
             xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
             ylabel = L"I (mA/cm^2)",
             limits = limits)
    else
        Axis(fig[1, 1],
             xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
             ylabel = L"I (mA/cm^2)")
    end

    plots  = Any[]
    labels = String[]

    # ------------------------------------------------------------
    # 1) Experimental background (Figure_3.csv)  [optional]
    # ------------------------------------------------------------
    if fig3_csv !== nothing
        try
            raw = CSV.read(fig3_csv, DataFrame; header=false)

            pres = vec(Matrix(raw[1:1, :]))
            sub  = Matrix(raw[4:end, :])

            num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
            num_df = DataFrame(num, :auto)

            npairs = size(num_df, 2) ÷ 2

            pink  = RGB(1.0, 0.7, 0.8)
            pblue = RGB(0.2, 0.5, 1.0)
            cols_exp = [RGB(pink.r + t*(pblue.r-pink.r),
                            pink.g + t*(pblue.g-pink.g),
                            pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]

            for j in 1:npairs
                xcol, ycol = 2j - 1, 2j
                lab = (j == 1) ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"

                x = num_df[!, xcol]
                y = num_df[!, ycol]

                # lighter background curves
                line = lines!(ax, x, y; color = (cols_exp[j], exp_alpha), linewidth = exp_linewidth)
                push!(plots, line)
                push!(labels, "$exp_title | $lab")
            end
        catch e
            # keep behavior: only skip for UndefVarError, otherwise rethrow
            if e isa UndefVarError
                # skip
            else
                rethrow(e)
            end
        end
    end

    # ------------------------------------------------------------
    # 2) Simulation curves (existing logic)
    # ------------------------------------------------------------
    n = length(P_recs)
    cols_sim = [RGB(1 - t, 0, t) for t in LinRange(0, 1, max(n, 1))]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        I = currents(rec, species) .* scale

        line = lines!(ax, rec.voltages, I ./ 2 ; color = cols_sim[j], linewidth = sim_linewidth)
        push!(plots, line)
        push!(labels, "Theoretical | $label")
    end

    Legend(fig[1, 2], plots, labels, "Overlay"; framevisible=true)
    return fig
end

# ╔═╡ d020a8db-2227-4a62-8c99-fa6539a261a5
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ 64fe2609-843e-49be-92b6-b462f6bf80b3
function pressure_varied_sweep(
    elydata_base, grid, bcondition, reaction, sawtooth;
    Pvec,
    ispec::Integer,
    base_value = elydata_base.c_bulk[ispec],
    scale = identity,
    sweep_kwargs...
)
    recs = Vector{Any}(undef, length(Pvec))
    for (k, p) in pairs(Pvec)
        ely = deepcopy(elydata_base)
        ely.c_bulk[ispec] = base_value * scale(p)
        pnpsys = PNPSystem(grid; celldata = ely, bcondition = bcondition, reaction = reaction)
        recs[k] = cvsweep(pnpsys; voltages = sawtooth, nperiods = 1, store_solutions = true, sweep_kwargs...)
    end
    return collect(zip(Pvec, recs))
end


# ╔═╡ b767a48f-1b20-4b85-b9a8-36d6a914c5fe
function _pressure_colname(p)
    s = @sprintf("%.4g", p)  
    s = replace(s, "." => "p", "-" => "m", "+" => "")
    return Symbol("pCO2_" * s * "_atm")
end

# ╔═╡ bb01b182-7840-4ad0-8cd6-fae57ab93173


# ╔═╡ 43e4ef6a-5c2a-4914-9972-ea93f3cb016e
md"""
### Extracted function name
"""

# ╔═╡ 9b8daa64-9e34-4da1-8136-ba5488e037c7
function cvsweep_compensated_over_L(
    model, bcondition, reaction; 
    sawtooth,
    nperiods = 1,
    L_values =[100, 500, 1000, 2000, 3000],
    f_comp = 1.0,  # 1.0 means 100% compensation (ideal overlap)
    Area = 0.000314,  
    store_solutions = true,
    solver_kwargs...,
)
    results = Dict{Int, Any}()

    for L in L_values
        @info ">>> Running Physics-Consistent Simulation: L = $(L) μm"
        
        # --- Mesh Generation ---
        hmin = 1.0e-4 * μm  
        # Keep hmax reasonable to ensure the bulk potential gradient is captured
        hmax = (L * μm) / 100.0  
        
        X = ExtendableGrids.geomspace(0, L * μm, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)
        
        celldata = deepcopy(model)
        # Ensure L is updated in celldata if the model uses it internally
        # celldata.L = L * μm 

        pnpcell = PNPSystem(
            grid;
            bcondition = bcondition,
            reaction = reaction,
            celldata = celldata
        )

        # --- Resistance Calculation ---
        # We calculate the theoretical Ohmic resistance for this specific length
        κ = calc_kappa(celldata) 
        R_bulk_theoretical = (L * μm) / (κ * Area)
        
        # R_compensated will be passed to cvsweep_COMP to 'undo' the voltage drop
        # in the recorded output.
        R_to_apply = R_bulk_theoretical * f_comp 

        @info "    Nodes: $(length(X)) | R_bulk: $(round(R_bulk_theoretical, digits=2)) Ω"

        # --- Execute Simulation ---
        results[L] = LiquidElectrolytes.cvsweep_COMP(
            pnpcell;
            voltages = sawtooth,
            nperiods = nperiods,
            Area = Area,
            R_comp = R_to_apply, # This is used in the 'post' function for V_eff calculation
            store_solutions = store_solutions,
            # Solver stability settings
            damp_initial = 0.1,   
            damp_growth = 1.2,    
            Δt_grow = 1.05,       
            max_round = 15,       
            tol_relative = 1.0e-5, 
            solver_kwargs...
        )
    end

    return results
end

# ╔═╡ c9ae7a67-9a5f-4d9a-88c0-742f4e91fb27


# ╔═╡ a2c7c4da-77cd-493f-8f98-0c86fecf271a
md"""
Run boundary layer thickness varied cyclic voltammetry $(@bind L_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 82116926-6739-4143-8cb7-29e34b8a323c
function blthickness(grd, celldata, tsol; species = 5)
    X = grd[XCoordinates]
    xbl = 0
    for it in 1:length(tsol.t)
        u = tsol[species, :, it]
        i = findlast(c -> abs(c - celldata.c_bulk[species]) > 1.0e-1, u)
        xbl = X[i]
    end
    return xbl
end

# ╔═╡ 8e2f8c2c-2c2a-46a1-90f3-8980cd6d8a51
function plot_co2_profiles(
    bulk, 
    L_values, 
    resL_COMP, 
    X_coords_dict; 
    L_varied_checkbox::Bool = true
)
    if !L_varied_checkbox
        return nothing
    end

    co2_idx = findfirst(s -> s.name == "CO₂", bulk)

    selected_Ls = L_values[1:end]
    num_panels = length(selected_Ls)
    max_L = 1e-2

    fig = Figure(size = (300, 300 * num_panels))
    
    c_min = log10(0.00001)
    c_max = log10(0.33)
    color_range = (c_min, c_max)
    
    num_levels = 12
    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    for (idx, L_val) in enumerate(selected_Ls)
        if !haskey(resL_COMP, L_val) 
            continue 
        end
        
        res = resL_COMP[L_val]
        tsol = res.tsol
        times = res.times
        X_coords = X_coords_dict[L_val]
        
        conc_matrix = [(tsol[co2_idx, ix, it] / (mol / dm^3)) 
                       for ix in 1:length(X_coords), it in 1:length(times)]
	
        ax = Axis(fig[idx, 1],
            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
            ylabel = "Time [s]",
            title = "Boundary Layer (L) = $(round(L_val/μm)) μm",
            xscale = log10,
            #limits = ((1e-10, max_L), (0, times[end])),
            xminorticksvisible = true, 
            xminorticks = IntervalsBetween(9)
        )
        
        hm = heatmap!(ax, X_coords .+ 1e-12, times, conc_matrix; 
            colorrange = color_range, 
            colormap = discrete_cmap,
            interpolate = false 
        )
	
        vlines!(ax, [L_val], color = :red, linestyle = :dash, linewidth = 2)
	
        Colorbar(fig[idx, 2], hm, label = L"\log_{10}(c_{\mathrm{CO_2}})")
    end
	
    rowgap!(fig.layout, 35)
    
    return fig
end

# ╔═╡ fc22ce7f-fd42-4f31-9792-7268ebc2928d
#unc_L_result = cvsweep_odr_over_L(elydata_Gold_unc, grid_dict, pnp_bcondition,  reaction, sawtooth; nperiods = user_input_cv.nperiods,)

# ╔═╡ 53fe324f-0024-4108-a855-66ab22952c3a
#irc_L_result = cvsweep_odr_over_L(elydata_Gold_irc, grid_dict, pnp_bcondition,  reaction, sawtooth; nperiods = user_input_cv.nperiods,)

# ╔═╡ 09d54862-2fd1-424f-aa29-d8961b2e93a4
sawtooth_src = SawTooth(vmin = -1.2, vmax = 0.8, scanup=false, scanrate = 3)

# ╔═╡ ae8d6798-a3eb-4e6e-b810-a26e81ec1438
#odr_L_result = cvsweep_odr_over_L(elydata_Gold_odr_v2, grid_dict, pnp_bcondition,  reaction, sawtooth; nperiods = 1)

# ╔═╡ d7e28d5f-16ba-4664-ba63-ec85fb29fe87
begin
	grid_dict = Dict()
    X_coords_dict = Dict()
	
	L_values = [100, 500, 1000, 1500, 2500].* μm
	
    for L_val in L_values
        hmin = 1.0e-6 * μm
        hmax = L_val * 0.02
        X = ExtendableGrids.geomspace(0, L_val, hmin, hmax)
        g = ExtendableGrids.simplexgrid(X)
        
        grid_dict[L_val] = g
        X_coords_dict[L_val] = g[Coordinates][1, :] .+ 1e-12 # include small offset
    end
	grid_dict
end

# ╔═╡ d91bcc6a-1f2d-4457-9200-1ee30953090f
function co2_log_contour(results::Dict, grid_dict, bulk;
                         scale=mol/dm^3, num_levels=24, L_val=nothing)
    
    if isnothing(L_val)
        L_val = minimum(keys(results))
    end
    
    result = results[L_val]
    grid = grid_dict[L_val]
    X = grid[Coordinates][1, :]
    
    co2_idx = findfirst(s -> s.name == "CO₂", bulk)
    times = result.tsol.t
    
    log_c_bulk = log10(result.tsol[co2_idx, end, 1] / scale)
    
    M = [log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / scale, 1e-12)) 
         for ix in 1:length(X), it in 1:length(times)]
    
    c_min = 0.0
    c_max = log10(result.tsol[co2_idx, end, 1] / scale) - log10(1e-5)
    
    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    f = Figure(size = (800, 450))
    
    ax = Axis(f[1, 1];
        xlabel = L"\text{time / s}",
        ylabel = L"x\;(\mathrm{m})",
        yscale = log10,
        yminorticksvisible = true, 
        yminorticks = IntervalsBetween(9),
        yticks = (10.0 .^ (-12:3:-6), [L"10^{-12}", L"10^{-9}", L"10^{-6}"])
    )
        
    hm = heatmap!(ax, times, X .+ 1e-12, M'; 
        colorrange = (c_min, c_max), 
        colormap = discrete_cmap,
        interpolate = false)
        
    # f[1, 2] 위치에 컬러바 배치
    cb = Colorbar(f[1, 2], hm; 
                  label = L"\log_{10}(c_{\mathrm{bulk}}) - \log_{10}(c_{\mathrm{CO_2}})", 
                  ticklabelsize = 20, 
                  labelsize = 20)
    
    #colsize!(f.layout, 2, SizeAttribute(:width, cb)) 
    
    return f
end

# ╔═╡ d91c32c8-ac55-4f77-93cb-6c5d97bcbbb2
function cvsweep_odr_over_L(
    elydata_odr, grid_dict, bcondition, reaction, sawtooth;
    nperiods = 1,
    store_solutions = true,
    solver_kwargs...
)
    results = Dict{Float64, Any}()
    L_values = sort(collect(keys(grid_dict)))

    for L in L_values
        grid = grid_dict[L]
        X    = grid[Coordinates][1, :]

        celldata = deepcopy(elydata_odr)
        celldata.Ru    = L / conductivity(celldata, celldata.c_bulk)
        celldata.x_ref = [X[end], 0, 0]

        @info ">>> :ohmicdrop sweep | L = $(L/μm) μm | Ru = $(round(celldata.Ru; digits=3)) Ω"

        pnpcell = PNPSystem(grid; bcondition = bcondition, celldata = celldata, reaction = reaction)

        results[L] = LiquidElectrolytes.cvsweep(
            pnpcell;
            voltages = sawtooth,
            nperiods = nperiods,
            store_solutions = store_solutions,
            solver_kwargs...
        )
    end

    return results
end

# ╔═╡ d9690f36-f2db-4fcb-85f0-e285e8ed5551
function electrochemistry_theme()
    Theme(
        size = (960, 540),
        Axis = (
            spinewidth = 5.5,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            xticksize = 8,
            yticksize = 8,
            xlabelsize = 25,
            ylabelsize = 25,
            xticklabelsize = 25,
            yticklabelsize = 25,
            xgridvisible = false,
            ygridvisible = false,
            xlabelpadding = 10,
            ylabelpadding = 10,
            xlabelfont = :bold,
            ylabelfont = :bold,
			yticks = LinearTicks(5),
			xticks = LinearTicks(4)
        ),
        Lines = (
            linewidth = 3, 
        )
    )
end

# ╔═╡ 2429e070-ed38-4a5a-9cdf-7818e7827021
function plot_cv_current(result, model; species=nothing, scale=cm^2/mA, color_gradient=true)
    sp = (species === nothing) ? model.cspecies[1] : species

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(540, 440)) 
        a = Axis(f[1, 1],
                 ylabel = L"I (mA/cm^2)",
                 xlabel = L"\phi (V \; \mathrm{vs}\; SHE)")
        return f, a
    end

    I = currents(result, sp) .* scale

    if color_gradient
        cols = RGBf.(range(0, 1, length=length(result.voltages)), 0.0, 0.0)
        lines!(ax, result.voltages, I; color=cols)
    else
        lines!(ax, result.voltages, I) 
    end

    return fig
end

# ╔═╡ 7bfe397e-2f49-461a-bca2-594e4bd3e607
function CV_dsp_cap_result(result, model; scale=cm^2/mA, show_diff=false)
    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(540, 440))
        a = Axis(f[1, 1],
                 ylabel = L"I \;(\mathrm{mA/cm^2})",
                 xlabel = L"\phi \;(\mathrm{V \; vs\; SHE})")
        return f, a
    end

    j_dsp = result.j_dsp .* scale      
    j_cap = result.j_cap .* scale      

    lines!(ax, result.voltages, j_dsp; color = :magenta, label = L"j_{dsp}")
    lines!(ax, result.voltages, j_cap; color = :skyblue, label = L"j_{cap}")

    if show_diff
        lines!(ax, result.voltages, j_dsp .- j_cap;
               color = :orange, linestyle = :dash, label = L"j_{dsp}-j_{cap}")
    end

    axislegend(ax; position = :rt)
    return fig
end

# ╔═╡ ad3a5236-72d0-4bbb-9fce-e74eb3ab52f0
function CV_total_current(result, model; species=nothing, scale=cm^2/mA)
    sp = (species === nothing) ? model.cspecies[1] : species

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(540, 440))
        a = Axis(f[1, 1],
                 ylabel = L"I \;(\mathrm{mA/cm^2})",
                 xlabel = L"\phi \;(\mathrm{V \; vs\; SHE})")
        return f, a
    end

    #i_F   = currents(result, sp) .* scale     
    i_cap = result.j_cap         .* scale    
    #i_tot = i_F .+ i_cap

   # lines!(ax, result.voltages, i_F;   color = :magenta,  label = L"i_F")
    lines!(ax, result.voltages, i_cap; color = :skyblue,  label = L"i_{cap}")
   # lines!(ax, result.voltages, i_tot; color = :orange,   label = L"i_{tot}")
    axislegend(ax; position = :rt)
    return fig
end

# ╔═╡ 9949de26-8fad-4dc6-9635-8aec64c3f73d
function plot_cv_total_current_tot(result, model;
                                   electrolyte,
                                   redox_species::Dict{Int,Int},
                                   co_idx,
                                   include_capacitive::Bool = true,
                                   scale = cm^2/mA,
                                   sign::Int = 1,
                                   color_F   = colorant"#F2728A",   # i_F 색
                                   color_C   = colorant"#5BA8E8",   # i_C 색
                                   mix_mode::Symbol = :mean,         # :sum 또는 :mean
                                   lw = 5)
    n_t = length(result.voltages)

    # ---- Faradaic current (total, sum over redox species) ----
    I_F = zeros(n_t)
    for (idx, n_e) in redox_species
        I_F .+= n_e .* currents(result, idx)
    end

    # ---- Capacitive current ----
    I_C = zeros(n_t)
    if include_capacitive
        if electrolyte.ircompensation == :ohmicdrop
            icc = electrolyte.icc
            node_we = 1
            I_C = [u[icc, node_we] for u in result.tsol[1:n_t]]
        else
            @warn "Capacitive current only resolved in :ohmicdrop mode " *
                  "(current mode: :$(electrolyte.ircompensation)); j_C set to 0."
        end
    end

    # ---- Scale and sign ----
    I_F_scaled = sign .* I_F          .* scale
    I_C_scaled = sign .* I_C          .* scale
    I_total    = sign .* (I_F .+ I_C) .* scale

    # ---- Sum color: RGB 합 또는 평균 ----
    cF = RGBf(color_F); cC = RGBf(color_C)
    color_tot = mix_mode === :mean ?
        RGBf((cF.r+cC.r)/2, (cF.g+cC.g)/2, (cF.b+cC.b)/2) :
        RGBf(min(cF.r+cC.r, 1f0), min(cF.g+cC.g, 1f0), min(cF.b+cC.b, 1f0))

    # ---- Plot ----
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (900, 800))
        ax_a = Axis(f[1, 2], ylabel = L"i_F \; (\mathrm{mA/cm^2})")
        ax_b = Axis(f[2, 2], ylabel = L"i_C \; (\mathrm{mA/cm^2})")
        ax_c = Axis(f[3, 2], ylabel = L"i_{tot}\; (\mathrm{mA/cm^2})",
                             xlabel = L"\phi \; (\mathrm{V \; vs \; SHE})")
        
        hidexdecorations!(ax_a; grid = false)
        hidexdecorations!(ax_b; grid = false)
        linkxaxes!(ax_a, ax_b, ax_c)
        return f, (ax_a, ax_b, ax_c)
    end
    f, (ax_a, ax_b, ax_c) = fig

    # ---- (a), (b), (c) 레이블 달기 ----
    # 플롯이 있는 f[i, 2] 내부에 Label을 얹어 좌측 상단(top-left)으로 정렬합니다.
    sub_axes = [ax_a, ax_b, ax_c]
    labels = ["(a)", "(b)", "(c)"]
    
    for i in 1:3
        Label(f[i, 1], labels[i],
              fontsize = 24,
              font = :bold,
              halign = :left,   # 왼쪽 정렬
              valign = :top,    # 상단 정렬
              padding = (15, 0, 0, 15) # 내부 여백 조절 (우, 좌, 하, 상 순서)
        )
    end

    # ---- 데이터 라인 그리기 ----
    lines!(ax_a, result.voltages, I_F_scaled; color = color_F,   linewidth = lw)
    lines!(ax_b, result.voltages, I_C_scaled; color = color_C,   linewidth = lw)
    lines!(ax_c, result.voltages, I_total;    color = color_tot, linewidth = lw)
        rowgap!(f.layout, 15) 
        
        rowsize!(f.layout, 1, Relative(0.30))
        rowsize!(f.layout, 2, Relative(0.30))
        rowsize!(f.layout, 3, Relative(0.40))
        #rowsize!(f.layout, 4, Relative(0.21)) 
    return f
end

# ╔═╡ c7f0514e-ea30-4d99-ada2-c698360ac34c
function plot_conc_time_electrode(result, bulk;
                                  nspecies=7,
                                  scale=(mol/dm^3))

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    iOH = findfirst(isequal("OH⁻"), names)
    iH  = (iOH === nothing) ? 1 : iOH # fallback
    I   = currents(result, iH) .* (cm^2/mA)

    cols = RGBf(0.3, 0.5, 1.0) 

    fig, ax_conc, ax_current = with_theme(electrochemistry_theme()) do
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

# ╔═╡ e6f43f01-15d8-4265-ae15-3673fb3cf7e3
function plot_activity_time_electrode(result, bulk, electrolyte;
                                     model_type="DGML_γ", # "DGML_γ" or "Stefan_γ"
                                     nspecies=7,
                                     scale=(mol/dm^3),
                                     ipressure=nothing)
    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    
    times = result.tsol.t
    nt    = length(times)
    
    activity_electrode = fill(NaN, nspecies, nt)
    ielectrode = 1 
    
    v0 = electrolyte.v0
    bar_c = 1.0 / v0
    RT = electrolyte.RT
    c_scale = 1.0 / scale

    for t in 1:nt
        c_all = result.tsol[:, ielectrode, t]
        
        Phi = sum(c_all[ic] * electrolyte.v[ic] for ic in 1:nspecies)
        solvent_frac = max(1.0 - Phi, eps(Float64))
        
        pnode = (ipressure !== nothing) ? result.tsol[ipressure, ielectrode, t] : 0.0

        for i in 1:nspecies
            c = c_all[i]
            v_i = electrolyte.v[i]
            size_ratio = v_i / v0
            term_conc = c / bar_c
            
            if model_type == "DGML_γ"
                term_press  = exp((1.0 - size_ratio) * pnode / (bar_c * RT))
                term_steric = solvent_frac^(-size_ratio)
                a_thermo    = term_conc * term_press * term_steric
            elseif model_type == "Stefan_γ"
                a_thermo    = term_conc * (solvent_frac^(-1.0))
            else
                a_thermo    = term_conc
            end
            
            activity_electrode[i, t] = a_thermo * (bar_c * c_scale)
        end
    end

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (960, 540))
        
        model_title = model_type == "DGML_γ" ? "Modified DGML Model (Activity)" : 
                      model_type == "Stefan_γ" ? "Stefan Model (Activity)" : "Ideal Solution"

        a = Axis(f[1, 1],
            title = model_title,
            xlabel = L"\text{time / s}",
            ylabel = L"\mathbf{a_{i,\,\mathrm{electrode}}}",
            limits = ((times[1] - (times[end] / 200), times[end] + (times[end] / 100)), (1e-12, 1e4)),
            yscale = log10
        )
        
        yt_vals = 10.0 .^ (4:-4:-12)
        yt_lbls = [L"10^{4}", L"10^{0}", L"10^{-4}", L"10^{-8}", L"10^{-12}"]
        a.yticks = (yt_vals, yt_lbls)
        
        return f, a
    end

    for i in 1:nspecies
        y = activity_electrode[i, :]
        # handle NaN and apply log-scale floor
        y_fixed = [isnan(val) || val <= eps(Float64) ? eps(Float64) : val for val in y]
        
        lines!(ax, times, y_fixed; color=colors[i])
    end

    return fig
end

# ╔═╡ 96c90e9f-9759-43d2-a13b-067d6a771cab
function CV_overlay_currents(results, model;
                             labels   = nothing,
                             species  = nothing,
                             scale    = cm^2/mA,
                             linestyles = [:solid, :dash, :dot])
    sp = (species === nothing) ? model.cspecies[1] : species
    labels = labels === nothing ? ["result $i" for i in 1:length(results)] : labels

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(620, 460))
        a = Axis(f[1, 1],
                 ylabel = L"I \;(\mathrm{mA/cm^2})",
                 xlabel = L"\phi \;(\mathrm{V \; vs\; SHE})")
        return f, a
    end

    col_F   = :magenta
    col_cap = :skyblue
    col_tot = :orange

    for (k, result) in enumerate(results)
        ls = linestyles[mod1(k, length(linestyles))]

        i_F   = currents(result, sp) .* scale
        i_cap = result.j_cap         .* scale
        i_tot = i_F .+ i_cap

        lines!(ax, result.voltages, i_F;   color = col_F,   linestyle = ls,
               label = L"i_F\;(%$(labels[k]))")
        lines!(ax, result.voltages, i_cap; color = col_cap, linestyle = ls,
               label = L"i_{cap}\;(%$(labels[k]))")
        lines!(ax, result.voltages, i_tot; color = col_tot, linestyle = ls,
               label = L"i_{tot}\;(%$(labels[k]))")
    end

    axislegend(ax; position = :rt, nbanks = 2)
    return fig
end

# ╔═╡ a2264942-aa50-4a63-823b-da26d491b040
function plot_cv_current_variedL(results, model;
                                 species = nothing,
                                 scale = cm^2/mA,
                                 color_gradient = true,   
                                 linewidth = 3.5,
                                 title = "")
    sp = (species === nothing) ? model.cspecies[1] : species

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 400))
        a = Axis(f[1, 1];
                 title  = title,
                 ylabel = L"I\;(\mathrm{mA/cm^2})",
                 xlabel = L"\phi\;(\mathrm{V \; vs\; SHE})")
        return f, a
    end

    Lkeys = sort(collect(keys(results)))
    n = length(Lkeys)

	cols = [RGBf(
	    0.6 + 0.3 * (i-1)/max(n-1,1),
	    0.8 - 0.2 * (i-1)/max(n-1,1),
	    0.7 - 0.2 * (i-1)/max(n-1,1)
	) for i in 1:n]
                                             

    for (i, L) in enumerate(Lkeys)
        I = currents(results[L], sp) .* scale ./ 2
        lines!(ax, results[L].voltages, I;
               color = cols[i], linewidth = linewidth,
               label = @sprintf("%g μm", L / μm))
    end

    axislegend(ax, "L", position = :rb, framevisible = false, legendtext = 24)
    return fig
end

# ╔═╡ 02a78c8d-d557-4769-b155-719d607418c7
function plot_cv_scanrate_grid(result_vec, model;
                               electrolyte,
                               redox_species::Dict{Int,Int},
                               co_idx,
                               scanrates = [0.05, 0.5, 5.0],
                               include_capacitive::Bool = true,
                               scale = cm^2/mA,
                               sign::Int = 1,
                               color_F   = colorant"#F2728A",
                               color_C   = colorant"#5BA8E8",
                               mix_mode::Symbol = :mean,
                               lw = 5.5)

    # 단위/라벨 헬퍼: 물리량 italic / 약어·단위 roman
    unit_i = rich("  (mA cm", superscript("−2"), ")")
    lab_iF   = rich(rich("I", font=:italic), subscript("F"),   unit_i)
    lab_iC   = rich(rich("I", font=:italic), subscript("C"),   unit_i)
    lab_itot = rich(rich("I", font=:italic), subscript("tot"), unit_i)
    lab_U    = rich(rich("U", font=:italic), "  (V vs. SHE)")

    # 1. 테마 및 축 매트릭스 생성
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 700))
        axes_matrix = Matrix{Axis}(undef, 3, 3)
        ylabels = [lab_iF, lab_iC, lab_itot]

        for col in 1:3, row in 1:3
            axes_matrix[row, col] = row == 3 ?
                Axis(f[row, col+1], xlabel = lab_U) :
                Axis(f[row, col+1])

            row < 3 && hidexdecorations!(axes_matrix[row, col]; grid = false)

            if col == 1
                axes_matrix[row, col].ylabel = ylabels[row]
                axes_matrix[row, col].ylabelpadding = 30
            else
                hideydecorations!(axes_matrix[row, col]; ticks = false, grid = false)
            end

            if row == 2
                axes_matrix[row, col].ytickformat = "{:.2f}"
                axes_matrix[row, col].yticks = LinearTicks(3)
            else
                axes_matrix[row, col].yticks = LinearTicks(4)
            end
        end

        for i in 1:3
            linkyaxes!(axes_matrix[i, 1], axes_matrix[i, 2], axes_matrix[i, 3])
            linkxaxes!(axes_matrix[1, i], axes_matrix[2, i], axes_matrix[3, i])
        end

        return f, axes_matrix
    end
    f, axes_matrix = fig

    grid_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"]
    for row in 1:3, col in 1:3
        idx = (row - 1) * 3 + col
        Label(f[row, col+1], grid_labels[idx],
              fontsize = 24, font = :bold, halign = :left, valign = :top,
              padding = (15, 0, 0, 15))
    end

    for col in 1:3
        Label(f[0, col+1],
              rich(string(scanrates[col]), "  V s", superscript("−1"));
              fontsize = 26, font = :bold, halign = :center, valign = :bottom)
    end

    cF, cC = RGBf(color_F), RGBf(color_C)
    color_tot = mix_mode === :mean ?
        RGBf((cF.r + cC.r)/2, (cF.g + cC.g)/2, (cF.b + cC.b)/2) :
        RGBf(min(cF.r + cC.r, 1f0), min(cF.g + cC.g, 1f0), min(cF.b + cC.b, 1f0))

    for col in 1:3
        res = result_vec[col]
        n_t = length(res.voltages)

        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(res, idx)
        end

        I_C = zeros(n_t)
        if include_capacitive && electrolyte.ircompensation == :ohmicdrop
            icc = electrolyte.icc
            I_C = [u[icc, 1] for u in res.tsol[1:n_t]]
        end

        I_F_scaled = sign .* I_F          .* scale
        I_C_scaled = sign .* I_C          .* scale
        I_total    = sign .* (I_F .+ I_C) .* scale

        lines!(axes_matrix[1, col], res.voltages, I_F_scaled; color = color_F,   linewidth = lw)
        lines!(axes_matrix[2, col], res.voltages, I_C_scaled; color = color_C,   linewidth = lw)
        lines!(axes_matrix[3, col], res.voltages, I_total;    color = color_tot, linewidth = lw)
    end

    rowgap!(f.layout, 15)
    colgap!(f.layout, 25)

    rowsize!(f.layout, 1, Relative(0.29))
    rowsize!(f.layout, 2, Relative(0.29))
    rowsize!(f.layout, 3, Relative(0.38))

    colsize!(f.layout, 1, Auto())
    for col in 2:4
        colsize!(f.layout, col, Relative(0.34))
    end

    return f
end

# ╔═╡ 08756476-b9d0-4bd4-be21-16cb93367dad
function plot_cv_scanrate_grid_unc(result_vec, model;
                               electrolyte,
                               redox_species::Dict{Int,Int},
                               co_idx,
                               scanrates = [0.05, 0.5, 5.0],
                               include_capacitive::Bool = true,
                               scale = cm^2/mA,
                               sign::Int = 1,
                               color_tot = "#BAC8FF",
                               lw = 5.5)

    # 라벨: 물리량 italic / 약어·단위 roman
    lab_I = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")")
    lab_U = rich(rich("U", font=:italic), "  (V vs. SHE)")

    # 1. 테마 및 축 벡터 생성 (1x3 구조)
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 350))
        axes_vec = Vector{Axis}(undef, 3)

        for col in 1:3
            axes_vec[col] = Axis(f[1, col+1], xlabel = lab_U)

            if col == 1
                axes_vec[col].ylabel = lab_I
                axes_vec[col].ylabelpadding = 30
            else
                hideydecorations!(axes_vec[col]; ticks = false, grid = false)
            end
            axes_vec[col].yticks = LinearTicks(4)
        end

        linkyaxes!(axes_vec[1], axes_vec[2], axes_vec[3])

        return f, axes_vec
    end
    f, axes_vec = fig

    grid_labels = ["(a)", "(b)", "(c)"]
    for col in 1:3
        Label(f[1, col+1], grid_labels[col],
              fontsize = 24, font = :bold, halign = :left, valign = :top,
              padding = (15, 0, 0, 15))
    end

    for col in 1:3
        Label(f[0, col+1],
              rich(string(scanrates[col]), "  V s", superscript("−1"));
              fontsize = 26, font = :bold, halign = :center, valign = :bottom)
    end

    for col in 1:3
        res = result_vec[col]
        n_t = length(res.voltages)

        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(res, idx)
        end

        I_C = zeros(n_t)
        if include_capacitive && electrolyte.ircompensation == :ohmicdrop
            icc = electrolyte.icc
            I_C = [u[icc, 1] for u in res.tsol[1:n_t]]
        end

        I_total = sign .* (I_F .+ I_C) .* scale
        lines!(axes_vec[col], res.voltages, I_total; color = color_tot, linewidth = lw)
    end

    colgap!(f.layout, 25)
    rowsize!(f.layout, 1, Relative(0.98))
    colsize!(f.layout, 1, Auto())
    for col in 2:4
        colsize!(f.layout, col, Relative(0.34))
    end
    return f
end

# ╔═╡ d3b7d864-bc1b-440a-8ddb-2096bfde0cc5
function plot_scanrate_sweeps_cv_2(
    sweep_vec, scanrates;
    electrolyte,                        # Capacitive 전류를 위해 필요
    redox_species::Dict{Int,Int},       # 모든 Faradaic 스피시즈와 전자 수 (e.g., Dict(1=>2, 2=>1))
    include_capacitive::Bool = true,
    scale = cm^2/mA,
    sign::Int = 1,
    default_lw = 2
)
    n = length(sweep_vec)
    pastel2 = cgrad([colorant"#D5F011", colorant"#16D8FF"])
    cols = [pastel2[t] for t in range(0, 1, length=max(n, 1))]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (850, 550))
        a = Axis(f[1, 1],
            ylabel = L"I_{\mathrm{total}}\;(\mathrm{mA/cm²})",
            xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})"
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        n_t = length(rec.voltages)
        
        # 1. Faradaic 전류 합산
        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(rec, idx)
        end

        # 2. Capacitive 전류 추가
        I_C = zeros(n_t)
        if include_capacitive
            if electrolyte.ircompensation == :ohmicdrop
                icc = electrolyte.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in rec.tsol[1:n_t]]
            else
                # 워닝이 너무 많이 뜨는 것을 방지하기 위해 첫 루프에서만 출력
                j == 1 && @warn "Capacitive current only resolved in :ohmicdrop mode; j_C set to 0."
            end
        end

        # 3. Total Current 계산 (부호 및 스케일 적용)
        I_total = sign .* (I_F .+ I_C) .* scale

        # 플롯 그리기
        lines!(ax, rec.voltages, I_total;
            linewidth = default_lw,
            color     = cols[j]
        )

        # 4. 텍스트 라벨링 위치 조정
        # 스캔 레이트별로 텍스트가 겹치지 않게 배치하는 기존 로직 유지
        pos_y = if j == 5
            -2.5
        elseif j == 3
            -1.5
        else
            -0.5 * j
        end

        text!(ax, "$(scanrates[j]) V/s";
            position = (-1.55, pos_y),   # 전류 크기 변화에 따라 y 위치는 조정이 필요할 수 있습니다.
            color    = cols[j],
            fontsize = 18,
            font     = :bold
        )			
    end

    return fig
end

# ╔═╡ ae50302f-02de-4e30-aa47-1baaa2c95a74
function plot_scanrate_sweeps_cv(
    sweep_vec, scanrates;
    species=iohminus,
    scale=cm^2/mA,
    default_lw = 4
)
    n = length(sweep_vec)
    pastel2 = cgrad([colorant"#D5F011", colorant"#16D8FF"])
    cols = [pastel2[t] for t in range(0, 1, length=max(n, 1))]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (850, 550))
        a = Axis(f[1, 1],
            ylabel = L"I\;(\mathrm{mA/cm²})",
            xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})"
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        I = currents(rec, species) .* scale
        lines!(ax, rec.voltages, I;
            linewidth = default_lw,
            color     = cols[j]
        )
        #text!(ax, "$(scanrates[j]) V/s";
        #    #position = (-1.2 - 0.05 * j, 0 - 1.5 * j),   # 위치는 직접 조절
        #    color    = cols[j],
        #    fontsize = 24,
        #    font     = :bold
        #)
    end

    return fig
end

# ╔═╡ 51064823-1b35-447d-95d0-2f9337e619eb
function plot_scanrate_sweeps_cap(
    sweep_vec, scanrates;
    #species=iohminus,
    scale=cm^2/mA,
    default_lw = 4                     
)
    n = length(sweep_vec)

    cols = [RGBf(0.1 + 0.6*(i/n), 
                 0.2 + 0.6*(1 - i/n), 
                 0.7 - 0.3*(i/n)) for i in 1:n]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 550)) 
        a = Axis(f[1, 1], 
            ylabel = L"I\;(\mathrm{mA/cm²})", 
            xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})"
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        color_j = cols[j]
        lw_j    = default_lw
        J_cap = rec.j_cap .* scale
        lines!(ax, rec.voltages, rec.j_cap; linewidth=lw_j, color=color_j)
    end

    if n > 1
        cgradient = cgrad(cols, categorical = true)
        
        Colorbar(fig[1, 2], 
            limits = (minimum(scanrates), maximum(scanrates)), 
            colormap = cgradient,
            label = L"\text{Scan Rate}\;(\mathrm{V/s})",
            labelsize = 20,
            ticklabelsize = 18,
            spinewidth = 2.0,
            width = 20 
        )
    end

    return fig
end

# ╔═╡ 3e1bcc1a-01f9-472c-a3eb-c332113aafbc
begin
    wanted    = ["0.1", "0.5", "1.0"]
    pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]

    raw_cv    = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
    sub    = Matrix(raw_cv[4:end, :])
    num_cv    = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num_cv, :auto)
    npairs = size(num_df, 2) ÷ 2

    n = length(wanted)
	pastel2 = cgrad([colorant"#F2728A", colorant"#5BA8E8"])
	cols1 = [pastel2[t] for t in range(0, 1, length=max(n, 1))]

    keep = findall(in(wanted), pressures[1:npairs])

    fig_cv, ax_cv = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1050, 500))
        a = Axis(f[1, 1], xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})", ylabel = L"I\;(\mathrm{mA/cm^2})")
        return f, a
    end

    # 기존 함수의 축 세팅 프레임 그대로 적용
    ax_cv.limits = ((-1.3, 0.9), (-5.5, 1.8))
    ax_cv.xticks = -1.5:0.3:1.0
    ax_cv.yticks = 1:-1:-5

    plot_objs1 = []
    labels1    = String[]

    # 3. 곡선 그리기 및 프레임 매칭
    for (k, j) in enumerate(keep)
        xcol, ycol = 2j - 1, 2j
        
        # 기존 함수의 레전드 텍스트 포맷 맞춤 ("$(p) atm")
        label_text = "$(pressures[j]) atm" 
        push!(labels1, label_text)

        # 기존 함수에서 전류 데이터를 [I ./ 2]로 스케일 조절해 그렸던 프레임 적용
        line = lines!(ax_cv, num_df[!, xcol], num_df[!, ycol] ;
                      color = cols1[k], 
                      linewidth = 5) # 시뮬레이션 선 두께 매칭
                      
        push!(plot_objs1, line)
    end

    # 4. 기존 함수의 내부 인셋 Legend 매칭 및 위치 미세조정(translate!)
    if n > 0
        leg = Legend(fig_cv[1, 1], title = "Exp", plot_objs1, labels1, L"p_{\mathrm{CO_2}}";
            framevisible = false,
            halign = :right, valign = :bottom, 
            labelsize = 20, titlesize = 23,
            padding = (0, 0, 0, 0),
            tellwidth = false, tellheight = false)
            
        translate!(leg.blockscene, -40, 40, 0)
    end

    expfig = fig_cv
end

# ╔═╡ ed4451bd-7888-4aea-8938-2a2ba85b22ad
function plot_pressure_varied_sweep_ivc(
    P_recs;
    species=iohminus,
    scale=cm^2/mA,
    limits=nothing,
    sim_linewidth::Real=5
)
    n = length(P_recs)
	pastel3 = cgrad([colorant"#FFB3BA", colorant"#A3D8FF"])  
	cols_sim = [pastel3[t] for t in range(0, 1, length=max(n, 1))]

    # 1. Figure & Axis 생성
    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1050, 500))
        a = if limits !== nothing
            Axis(f[1, 1], xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})", ylabel = L"I\;(\mathrm{mA/cm^2})", limits = limits)
        else
            Axis(f[1, 1], xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})", ylabel = L"I\;(\mathrm{mA/cm^2})")
        end
        return f, a
    end
    #ax.limits = ((-1.3, 0.9), (-1, 0.2))
    #ax.xticks = -1.5:0.3:1.0
    #ax.yticks = 0.2:-0.5:-1
    # 범례 매핑용 컨테이너 생성
    plot_objs = []
    labels    = String[]

    # 2. 곡선 그리기
    for j in 1:n
        p, rec = P_recs[j]
        I = currents(rec, species) .* scale

        label_text = "$(p) atm" 

        # lines! 인스턴스를 받아 보관
        hl = lines!(ax, rec.voltages, I ./ 2; 
            color = cols_sim[j], 
            linewidth = sim_linewidth)
            
        push!(plot_objs, hl)
        push!(labels, label_text)
    end

    # 3. 요청하신 포맷의 내부 Legend 매칭 (틀 안으로 인셋 배치)
    if n > 0
        leg = Legend(fig[1, 1], plot_objs, labels, L"p_{\mathrm{CO_2}}";
            framevisible = false,
            halign = :right, valign = :bottom, 
            labelsize = 20, titlesize = 23,
            padding = (0, 0, 0, 0),
            tellwidth = false, tellheight = false)
            
        translate!(leg.blockscene, -40, 40, 0)
    end

    return fig
end

# ╔═╡ 92417c77-5ad2-451a-9848-44c9fe1f103b
plot_pressure_varied_sweep_ivc(P_recs_unc, species = iohminus)

# ╔═╡ 6f6e779c-6ba3-49c5-9541-359ab4dbf13f
function plot_cv_current_dict(result_dict, model;
                              species = nothing,
                              scale = 1.0,
                              title = "",
                              linewidth = 3)
    sp = (species === nothing) ? model.cspecies[1] : species

    sorted_keys = sort(collect(keys(result_dict))) 
   # cols = distinct_colors(length(sorted_keys))

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (400, 300))
        a = Axis(f[1, 1];
                 title  = title,
                 ylabel = L"I\;(\mathrm{mA/cm^2})",
                 xlabel = L"\phi\;(\mathrm{V\;vs\;SHE})")
        return f, a
    end

    for (i, key) in enumerate(sorted_keys)
        result = result_dict[key]
        I = currents(result, sp) .* scale ./ 2
        lines!(ax, result.voltages, I;
              linewidth = linewidth,
               label = @sprintf("L = %g μm", key / μm))
    end

    axislegend(ax, position = :rb, framevisible = true)
    return fig
end

# ╔═╡ 842b074b-f808-48d8-8dc5-110ddd907f90
md"""
## IV Results
"""

# ╔═╡ 31298257-d35a-4f6f-8a76-ff00d5361ced
md"""
### System Setup
"""

# ╔═╡ cd310c1a-5810-4d13-9a71-96161e4c7452
begin
	    hmax_iv = 1.0 * μm 
		hmin_iv = 1.0e-6 * μm 
	    X_iv = ExtendableGrids.geomspace(0, 80 * μm, hmin_iv, hmax_iv)
	    grid_iv = ExtendableGrids.simplexgrid(X_iv)
end

# ╔═╡ 72269ec4-a56e-46d9-85c8-0dd8ccaf43e1
solver_control = (; max_round 	= 4,
					maxiters 	= 20,
              		tol_round 	= 1.0e-9,
              		verbose 	= "a",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)

# ╔═╡ 57db41d1-57c0-4eee-ba35-6f9b7e1e8263
md"""
##### **Run general polarization curve** $(@bind IV PlutoUI.CheckBox())
"""

# ╔═╡ b976ab43-69f1-47a0-b2c6-c63e1c15cdb4
md"""
##### **Boundary Layer thickness varied polarization curve:** $(@bind Ldependancy PlutoUI.CheckBox(default=false))
"""

# ╔═╡ b6f0c4bb-153c-461c-b240-dff009b77dc9


# ╔═╡ 7a02463d-cfd9-4648-af53-f1e65d46733f
md"""
### Result Plots
"""

# ╔═╡ 114d2324-5289-4e44-8d77-736a9bdec365
md"""
Show only pH: $(@bind useonly_pH PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 425a5f53-ab8b-4596-bda7-586842e13878


# ╔═╡ f2a1829d-e3d9-45c1-85a8-4ffb1564fe0f
function electrode_activity_vs_voltage(result, grid, electrolyte;
    c_ref = 1.0 * ufac"mol/dm^3",
    node_selector = :minx,
)
    tsol = LiquidElectrolytes.voltages_solutions(result)
    vgrid = LiquidElectrolytes.voltages(result)

    xcoords = grid.components[XCoordinates]
    ielectrode = node_selector === :minx ? argmin(xcoords) : argmax(xcoords)

    cspecies = electrolyte.cspecies
    ip       = LiquidElectrolytes.pressure_index(electrolyte)

    nspecies = length(cspecies)
    nv       = length(vgrid)

    γ_e = fill(NaN, nspecies, nv)
    a_e = fill(NaN, nspecies, nv)
    c_e = fill(NaN, nspecies, nv)

    for (j, U) in enumerate(vgrid)
        sol = tsol(U)
        sol === nothing && continue

        unode = view(sol, :, ielectrode)
        pnode = unode[ip]

        γ = zeros(eltype(sol), size(sol, 1))

        electrolyte.actcoeff!(γ, unode, pnode, electrolyte)

        for (k, ic) in enumerate(cspecies)
            cval = unode[ic]                   
            γval = γ[ic]                       
            aval = γval * (cval / c_ref)       

            c_e[k, j] = cval / c_ref          
            γ_e[k, j] = γval
            a_e[k, j] = aval
        end
    end

    return (
        voltages = vgrid,
        ielectrode = ielectrode,
        cspecies = cspecies,
        gamma_electrode = γ_e,
        activity_electrode = a_e,
        concentration_scaled = c_e,
    )
end

# ╔═╡ f8b5dc8f-1f41-4600-825e-2f9653f2d925
md"""
### **Polarization Curve**
"""

# ╔═╡ 94e4c398-b75e-4573-ae00-e54bdad267f6
sans = :regular

# ╔═╡ 15185eba-000c-4b6a-9692-903b6915bfd9
function iv_curve_axis(ivresult;
    cutoff = -0.4,
    showlegend = true,
    species = iohminus,
    data_dir = "../data/catmap_CO2R_data",
)
    v_all = ivresult.voltages
    mask  = v_all .< cutoff
    volts = v_all[mask]
    I_sim = abs.(currents(ivresult, species))[mask] .* (cm^2/mA)

    table2 = readdlm(joinpath(data_dir, "Ringe-theorical.csv"),    ',', Float64, '\n')
    table3 = readdlm(joinpath(data_dir, "Ringe-experimental.csv"), ',', Float64, '\n')
    df2_v = table2[:, 1]; df2_I = max.(abs.(table2[:, 2]), eps(Float64))
    df3_v = table3[:, 1]; df3_I = max.(abs.(table3[:, 2]), eps(Float64))

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (960, 540))
        a = Axis(f[1, 1];
            xlabel = L"\textbf{Voltage}\ \;\; U\ \text{(V vs. SHE)}",
            ylabel = L"\textbf{Partial Current}\ \;\; I_{CO}\ \text{(mA/cm}^2\text{)}",
            yscale = log10,
            limits = ((-1.5, -0.50), (1e-11, 1e2)),
        )
        return f, a
    end

    xt = [-1.4, -1.2, -1.0, -0.8, -0.6]
    ax.xticks = (xt, [@sprintf("%.1f", x) for x in xt])
    ax.yticks = (10.0 .^ (0:-3:-9),
                 [L"10^{0}", L"10^{-3}", L"10^{-6}", L"10^{-9}"])

    ax.xlabelfont     = :regular
    ax.ylabelfont     = :regular
    ax.xticklabelfont = :regular
    ax.yticklabelfont = :regular
    ax.xlabelsize     = 26
    ax.ylabelsize     = 26
    ax.xticklabelsize = 24
    ax.yticklabelsize = 24
    ax.spinewidth     = 5.5
    ax.xtickwidth     = 2.0
    ax.ytickwidth     = 2.0
    ax.xticksize      = 8
    ax.yticksize      = 8
    ax.xlabelpadding  = 10
    ax.ylabelpadding  = 10
    ax.xgridvisible   = false
    ax.ygridvisible   = false

    lines!(ax, volts, I_sim;
        color     = "#008b3f",
        linewidth = 8,
        label     = "LiquidElectrolyte.jl"
    )
	lines!(ax, df2_v, df2_I;
	    color     = ("#0083fe", 0.6),
	    linewidth = 10,
	    label     = "CatINT"
	)
    scatter!(ax, df3_v, df3_I;
        marker     = :circle,
        markersize = 15,
        color      = "#e52c40",
        label      = "Experiment"
    )

    if showlegend
        axislegend(ax;
            position     = :rt,
            labelsize    = 24,
            titlesize    = 24,
            labelfont    = :regular,
            titlefont    = :regular,
            framevisible = false,
            patchsize    = (40, 20)
        )
    end

    return fig
end

# ╔═╡ c5a1bc1d-a2ff-418d-b8f1-a78e06a61c6b


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

# ╔═╡ 5d0458f3-6564-43df-af01-4a4b829ce262


# ╔═╡ 11d6598c-3d6a-472a-8863-9275d5e567c6
function activity_vs_voltage_axis(
    result;
    bulk,
    grid,
    electrolyte,
    ipressure,
    model_type::String, # "Stefan(MPB)" "DGML"
    useonly_pH::Bool = false,
    showlegend::Bool = false,
)
    species  = getproperty.(bulk, :name)
    colors   = getproperty.(bulk, :color)
    nspecies = length(species)

    tsol  = LiquidElectrolytes.voltages_solutions(result)
    vgrid = result.voltages

    xcoords    = grid.components[XCoordinates]
    ielectrode = argmin(xcoords)

    scale = 1.0 / (mol / dm^3)
    nv = length(vgrid)

    activity_electrode = fill(NaN, nspecies, nv)
    gamma_electrode    = fill(NaN, nspecies, nv)
    conc_electrode     = fill(NaN, nspecies, nv)

    cspecies = electrolyte.cspecies

    for (j, v) in enumerate(vgrid)
        sol = tsol(v)
        sol === nothing && continue

        # Extract spatial unit (x=0) our unit vector and magnitude
        cnode = zeros(eltype(sol), maximum(cspecies))
        for ic in cspecies
            cnode[ic] = sol[ic, ielectrode]
        end
        pnode = sol[ipressure, ielectrode]

        v0 = electrolyte.v0
        bar_c = 1.0 / v0
        RT = electrolyte.RT
        
        Phi = sum(cnode[ic] * electrolyte.v[ic] for ic in cspecies)
        solvent_frac = max(1.0 - Phi, eps(Float64))

        for ia in 1:nspecies
            c = sol[ia, ielectrode]
            v_a = electrolyte.v[ia]
            size_ratio = v_a / v0
            
            term_conc = c / bar_c
            
            if model_type == "DMGL_γ"
				# DGML Model: Treats the electrolyte as an ideal incompressible mixture.
                # term_press: Captures the mechanical pressure penalty scaled by the specific volume difference.
                # Species with v_a > 0 are physically repelled by local pressure gradients (Barodiffusion).
                # Dimensionless species (v_a = 0) feel pure pressure-correction without steric linkage.
                term_press  = exp((1.0 - size_ratio) * pnode / (bar_c * RT))
                term_steric = solvent_frac^(-size_ratio)
                a_eff_thermo = term_conc * term_press * term_steric
                
            elseif model_type == "Stefan_γ"
				# Stefan's Model (Bikerman-Freise): Applies a global lattice-based steric penalty.
                # All species, regardless of their actual size (even v_a = 0), are subjected to the exact same penalty (1 - Φ)^-1.
                # This causes the unphysical coupled depletion of point-charge species when supporting cations overcrowd.
                a_eff_thermo = term_conc * (solvent_frac^(-1.0))
                
            else
                error("Invalid model_type. Use 'DGML' or 'Stefan'.")
            end
            
            c_scale = c * scale
            a_eff_scaled = a_eff_thermo * (bar_c * scale)

            conc_electrode[ia, j]     = c
            gamma_electrode[ia, j]    = a_eff_scaled / max(c_scale, eps(Float64))
            activity_electrode[ia, j] = a_eff_scaled
        end
    end

    fig = Figure(size=(1200, 800)) 

    ax = Axis(fig[1,1];
        xlabel = L"\mathbf{\text{U}\ \mathrm{vs.}\ \text{SHE}\ (V)}",
        ylabel = L"\mathbf{\tilde{a}_i}",
        limits = ((-1.25, -0.50), nothing),
        yscale = log10
    )

    if useonly_pH
        iH = findfirst(isequal("H⁺"), species)
        iH === nothing && error("H⁺ not found in species list.")
        
        valid_mask = .!isnan.(activity_electrode[iH, :])
        x_data = vgrid[valid_mask]
        y_data = max.(activity_electrode[iH, valid_mask], eps(Float64))
        
        lines!(ax, x_data, y_data; color=colors[iH], linewidth=5, label=species[iH])
    else
        for ia in 1:nspecies
            valid_mask = .!isnan.(activity_electrode[ia, :])
            x_data = vgrid[valid_mask]
            y_data = max.(activity_electrode[ia, valid_mask], eps(Float64))
            
            if sum(valid_mask) > 0
                lines!(ax, x_data, y_data; color=colors[ia], linewidth=5, label=species[ia])
            end
        end
    end

    showlegend && axislegend(ax, position=:rt)

    return (
        fig=fig,
        species=species,
        colors=colors,
        conc_electrode=conc_electrode,
        gamma_electrode=gamma_electrode,
        activity_electrode=activity_electrode,
        vgrid=vgrid,
    )
end

# ╔═╡ 91051ed4-9fd0-4c21-95f4-efc042060e4d
md"""
Extract the polarization Curve: $(@bind extractIV PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 43cb97f2-22ce-4940-8a17-3c04d4d2ced6
md"""
Extract the polarization Curve: $(@bind extract_conc PlutoUI.CheckBox(default=true))
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
			
			- Scan rate ``(V/s)``: $(Child("scanrate", NumberField(1e-6:1e-6:1e3; default = 0.05)))  
			- Periods: $(Child("nperiods", NumberField(1:10; default = 1)))
			- `Double64`: $(Child("double64", CheckBox()))  
			- `tunnel`: $(Child("tunnel", CheckBox()))
			- `scanup`: $(Child("scanup", CheckBox()))			
			- ``f_{comp}`` : $(Child("ircomp", NumberField(0.0:0.01:1.0; default = 0.85)))
			
			"""
		end;
	        label = "Submit"
	);
	top = 50
)

# ╔═╡ ef4d7cee-fed1-489d-8795-c2dea9361e78
sawtooth = SawTooth(
        scanrate = user_input_cv.scanrate,
        vmin     = user_input_cv.vmin, 
        vmax     = user_input_cv.vmax,
        scanup   = user_input_cv.scanup,
		vstart = 0.0; tstart = 0.0
    )

# ╔═╡ 3b07d9a6-bf83-43d0-89b6-64c112a94833
sawtooth

# ╔═╡ 9ab242ff-1c4b-4970-8db1-877674838f29
sawtooth_comp = sawtooth.vmax = 0.8

# ╔═╡ 3f030783-2e0d-48c3-8407-5d317311c021
sawtooth_v2 = SawTooth(
    scanrate = 1.0,
    vmin     = user_input_cv.vmin,
    vmax     = user_input_cv.vmax,
    scanup   = false,
    vstart   = 0.0,
    tstart   = 0.0,
)

# ╔═╡ 0a1054c5-cee7-4202-9d32-9eee6ec55265
user_input_cv.nperiods

# ╔═╡ d75725cd-0ef6-421f-be56-f559312e73b6
floataside(
    @bind user_input_ion confirm(
        PlutoUI.combine() do Child
			md"""
			###### __Parameter Set__
			
			- **Use Physical ion size & solvation number**
			  ``κ``, ``a`` : $(Child("use_physical_size", PlutoUI.CheckBox(default=true)))
			  - **ON**  : Apply physical ion size and solvation number.
			  - **OFF** : Use the previous model setting (only potassium has a size).
			"""
        end;
        label = "Submit"
    );
    top = 350
)

# ╔═╡ ab0e28f4-4310-4dcc-817e-81b9e45fd501
floataside(
    @bind user_input_model confirm(
        PlutoUI.combine() do Child
			md"""
			###### __Model Selection__  
			- Model: $(Child("model_choice", Select(["Gold_Model", "Landstorfer_NaClO₄ model", "Landstorfer_NaF model", "Toy model"])))
			---
			###### __Boundary layer thickness__   
			- ``δ``: $(Child("L", NumberField(1:80000; default = 2500))) μm  			
			---
			
			###### __Activity Coefficient__  
			- LiquidElectrolyte.Mode: $(Child("mode", Select(["DMGL_γ", "Stefan_γ"])))
			- Boundary Condition : $(Child("BC_Select", Select(["Robin", "Neumann", "Dirichlet"])))
			"""
	    end;
	    label = "Submit"
	);
	top = 600
)

# ╔═╡ 2c061d99-da09-4ccf-a4ec-73a3c9a908fe
begin
	#geometry function constant
	#Γ_we = BP1
	const Γ_we 		= 1
	const Γ_bulk 	= 2
	const L = user_input_model.L * μm
end

# ╔═╡ 05ebcce5-6904-451f-bbd6-ea4588ca05e3
begin
    #Vmax = 2 * V
    #L = user_input_model.L * μm
    hmin = 1.0e-6 	* μm
    hmax = L * 0.0070	#* μm
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ d5c6769d-35f9-461a-aae2-00122ccddd63
gridplot(grid; size = (750, 200), xlabel = "Domain Size(x)")

# ╔═╡ 303a98e5-20f5-4e59-928a-a305bf9a0b30
let
	coord = grid[Coordinates]
	x = sort(vec(coord))
	L = maximum(x)

	fig = with_theme(electrochemistry_theme()) do
	    Figure(size = (750, 320))
	end

	ax = Axis(fig[1, 1];
	    xlabel = rich("Domain size  ", rich("x", font=:italic), "  (m)"),
	    yticklabelsvisible = false, yticksvisible = false,
	    xgridvisible = false, ygridvisible = false,
	    topspinevisible = false, rightspinevisible = false,
	    leftspinevisible = false,
	)
	hideydecorations!(ax)
	ylims!(ax, -1.6, 1.6)
	xlims!(ax, -0.0002, 0.0027)
	vlines!(ax, x; ymin = 0.42, ymax = 0.58, color = (:black, 0.6), linewidth = 0.5)
	lines!(ax, [0, L], [0, 0]; color = (:firebrick, 0.4), linewidth = 1.5)

	vlines!(ax, [0.0]; ymin = 0.30, ymax = 0.70, color = :red,   linewidth = 4)
	vlines!(ax, [L];   ymin = 0.30, ymax = 0.70, color = :green, linewidth = 4)
	text!(ax, -0.0, -0.8; text = rich("electrode"), align = (:center, :top), fontsize = 18)
	text!(ax, L,   -0.8; text = rich("bulk"), align = (:center, :top), fontsize = 18)

	# --- three adjacent control volumes K, L, M ---
	i = length(x) - 20
	nodes = x[i:i+2]                     # K, L, M
	faces = [(x[i-1]+x[i])/2, (x[i]+x[i+1])/2, (x[i+1]+x[i+2])/2, (x[i+2]+x[i+3])/2]

	vlines!(ax, faces; ymin = 0.35, ymax = 0.65,
	        color = (:steelblue, 0.9), linewidth = 1.5, linestyle = :dash)
	scatter!(ax, nodes, fill(0.1, 3); color = :black, markersize = 9)

	xK, xL, xM = nodes
	dx = (faces[end] - faces[1]) * 0.3   # 라벨 가로 오프셋 (꺾인 정도)

	text!(ax, xK - dx, 1.25; text = rich("K", font=:italic),
	      align = (:center, :bottom), fontsize = 22)
	arrows!(ax, [xK - dx], [1.18], [dx], [-0.55];
	        color = :black, linewidth = 1.5, arrowsize = 11)

	text!(ax, xL, 1.25; text = rich("L", font=:italic),
	      align = (:center, :bottom), fontsize = 22)
	arrows!(ax, [xL], [1.18], [0.0], [-0.55];
	        color = :black, linewidth = 1.5, arrowsize = 11)

	text!(ax, xM + dx, 1.25; text = rich("M", font=:italic),
	      align = (:center, :bottom), fontsize = 22)
	arrows!(ax, [xM + dx], [1.18], [-dx], [-0.55];
	        color = :black, linewidth = 1.5, arrowsize = 11)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ fe0c7939-3bda-4011-8aa3-5301cd23e302
X[end]

# ╔═╡ 0fe0f681-966b-40ee-9323-da1fa8c95d80
function sweep_over_L_cv(
    elydata,
    L_values,
    bcond,
    sawtooth,
    reaction;
    nperiods::Int = 1,
    eneutral::Bool = false,
    store_solutions::Bool = true,
    sweep_kwargs...
)
    results = Dict{Float64, Any}()

    for L in L_values
        hmax = L / 100.0
        hmax = 1.0    * μm
        X    = ExtendableGrids.geomspace(0, L, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)
		R_u = L / conductivity(elydata, elydata.c_bulk)
		
        celldata = deepcopy(elydata)
        celldata.eneutral = eneutral
		celldata.Ru = R_u
        pnpcell = PNPSystem(grid; bcondition=bcond, celldata=celldata, reaction = reaction)

        results[L] = cvsweep(
            pnpcell;
            voltages = sawtooth,
            nperiods = nperiods,
            store_solutions = store_solutions,
            sweep_kwargs...
        )
    end

    return results
end

# ╔═╡ 0c083fc8-0e73-4085-93c3-606493b3c46f
L

# ╔═╡ 9d814b85-a5b6-42e5-abf4-15500bbdb717
begin
	γ_key = user_input_model.mode 
	γ_mode = γ_key == "Stefan_γ" ? Stefan_γ! : DGML_γ!
end;

# ╔═╡ ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
begin
    use_md_hydrated = user_input_ion.use_physical_size   # Bool toggle you will use

    # a: hydrated radii in water (Å), κ: MD hydration number (1st shell)
	
    if use_md_hydrated == true 

		a_HCO3 = 8.2;  κ_HCO3 = 0
        a_CO3  = 8.2;  κ_CO3  = 0
        a_CO2  = 8.2;  κ_CO2  = 0
        a_OH   = 8.2;  κ_OH   = 0
        a_H    = 8.2;  κ_H    = 0
        a_CO   = 8.2;  κ_CO   = 0
        a_K    = 8.2;  κ_K    = 0
		"""
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
						"""

	elseif use_md_hydrated == true && γ_key == "Stefan_γ"
        # Stefan's Model / Consider all effective size
        a_HCO3 = 8.5;  κ_HCO3 = 0
        a_CO3  = 9.9;  κ_CO3  = 0
        a_CO2  = 3.4;  κ_CO2  = 0
        a_OH   = 7.6;  κ_OH   = 0
        a_H    = 7.3;  κ_H    = 0
        a_CO   = 2.8;  κ_CO   = 0
        a_K    = 8.2;  κ_K    = 0
	elseif use_md_hydrated == false && γ_key == "DMGL_γ"
		# Stefan's Model / Only Potaissium has a size
        a_HCO3 = 0;  κ_HCO3 = 0
        a_CO3  = 0;  κ_CO3  = 0
        a_CO2  = 0;  κ_CO2  = 0
        a_OH   = 0;  κ_OH   = 0
        a_H    = 0;  κ_H    = 0
        a_CO   = 0;  κ_CO   = 0
        #a_K    = 3.31;  κ_K    = 6.0
        a_K    = 8.2;  κ_K    = 0.0#5.0
	else
		a_HCO3 = 0;  κ_HCO3 = 0
        a_CO3  = 0;  κ_CO3  = 0
        a_CO2  = 0;  κ_CO2  = 0
        a_OH   = 0;  κ_OH   = 0
        a_H    = 0;  κ_H    = 0
        a_CO   = 0;  κ_CO   = 0
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

# ╔═╡ 95bb0ec4-73c2-4735-b7a6-f54911881323
let
    L_values = sort(collect(keys(unc_L_result)))
    num_plots = length(L_values)
    
    result_dicts = [unc_L_result, odr_L_result]
    col_titles = ["Uncompensated", "Ohmic Drop"]
    
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (450 * 3, 250 * num_plots))
        
        co2_idx = findfirst(s -> s.name == "CO₂", bulk)
        discrete_cmap = cgrad(:jet, 24, categorical = true)
        
        ref_result = unc_L_result[L_values[1]]
        c_min = 0.0
        c_max = log10(ref_result.tsol[co2_idx, end, 1] / (mol/dm^3)) - log10(1e-5)
        
        last_hm = nothing
        
        for (j, result_dict) in enumerate(result_dicts)
            for (i, L) in enumerate(L_values)
                result = result_dict[L]
                grid = grid_dict[L]
                X = grid[Coordinates][1, :]
                times = result.tsol.t
                log_c_bulk = log10(result.tsol[co2_idx, end, 1] / (mol/dm^3))
                
                M = [log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / (mol/dm^3), 1e-12))
                     for ix in 1:length(X), it in 1:length(times)]
                
                ax = Axis(f[i, j];
                    ylabel = j == 1 ? L"x\;(\mathrm{m})" : "",
                    yscale = log10,
                    yminorticksvisible = true,
                    yminorticks = IntervalsBetween(9),
                    yticks = (10.0 .^ (-12:3:0), [L"10^{-12}", L"10^{-9}", L"10^{-6}", L"10^{-3}", L"10^{0}"]),
                    title = i == 1 ? col_titles[j] : "L = $(round(L/μm, digits=1)) μm",
                    titlesize = 24
                )
                
                # y축 decoration은 첫 번째 열만
                if j != 1
                    hideydecorations!(ax, grid = false)
                end
                
                # x축 label은 마지막 행만
                if i == num_plots
                    ax.xlabel = L"\text{time / s}"
                else
                    hidexdecorations!(ax, grid = false)
                end
                
                last_hm = heatmap!(ax, times, X .+ 1e-12, M';
                    colorrange = (c_min, c_max),
                    colormap = discrete_cmap,
                    interpolate = false)
                
                # 패널 라벨: (a), (b), ... 열 우선으로
                label_idx = (i - 1) * 3 + j
                label_char = string(Char(96 + label_idx))
                Label(f[i, j], "($label_char)", fontsize = 25, font = :bold,
                      halign = :left, valign = :top, padding = (15, 0, 0, 15),
                      tellwidth = false, tellheight = false, color = :white)
            end
        end
        
        # 공유 colorbar: 오른쪽 열 전체
        cb = Colorbar(f[1:num_plots, 4], last_hm;
                      label = L"\log_{10}(c_{\mathrm{bulk}}) - \log_{10}(c_{\mathrm{CO_2}})",
                      ticklabelsize = 20, labelsize = 20)
        
        rowgap!(f.layout, 25)
        colgap!(f.layout, 10)
        
        f
    end
    
    CairoMakie.save("co2_contours_comparison.png", fig)
    #fig
end

# ╔═╡ 8f63ec58-bc97-43ba-9bbe-e10bb50e2cfe
elydata_Gold_unc = ElectrolyteData(; 
								    ircompensation=:none,
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
									actcoeff! = γ_mode,
								    C_gap = C_gap,
								    ϕ_pzc = ϕ_pzc,
									ircompfactor = user_input_cv.ircomp,
									#redoxreaction = we_breactions
									)

# ╔═╡ cd62234e-7f4a-4d0e-aea5-989f410d8cdc
R_u = L / conductivity(elydata_Gold_unc, elydata_Gold_unc.c_bulk)

# ╔═╡ 13a0e5da-4b0e-4cf2-b43e-05f17802e48d
begin
	powlab(n) = rich("10", superscript(string(n)))
	lab_time  = rich(rich("t", font=:italic), "  (s)")

	function panel_conc_time!(fig, panel_pos, result, bulk; nspecies=7, scale=mol/dm^3, lw=4,
	                          xlabel = lab_time, legend_pos = nothing)
	    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
	                 "#C0392B","#27AE60","#2980B9"]
	    colors = sp_colors[1:nspecies]
	    names  = getproperty.(bulk, :name)
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

	function panel_time_current!(fig, panel_pos, result, model;
                             redox_species = nothing,
                             include_capacitive = true,
                             scale = cm^2/mA,
                             sign = -1,
                             color_gradient = true,
                             lw = 4,
                             xlabel = lab_time)

    ax = Axis(panel_pos; xlabel = xlabel,
              ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"))

    redox = Dict(model.cspecies[ico2] => 2)

    n_t = length(result.times)
    I_F = zeros(n_t)
    for (idx, n_e) in redox
        I_F .+= n_e .* currents(result, idx)
    end

    I_C = zeros(n_t)
    if include_capacitive
        ely = elydata_Gold_unc
        if ely.ircompensation == :ohmicdrop
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
end

# ╔═╡ ac2bef83-9ebb-4bdc-869a-44cc1e2ef9e4
function plot_7species_contours(result, X, bulk; scale=mol/dm^3, num_levels=24)
    target_species = [
        ("K⁺",    rich("K", superscript("+"))),
        ("H⁺",    rich("H", superscript("+"))),
        ("CO₂",   rich("CO", subscript("2"))),
        ("OH⁻",   rich("OH", superscript("−"))),
        ("HCO₃⁻", rich("HCO", subscript("3"), superscript("−"))),
        ("CO₃²⁻", rich("CO", subscript("3"), superscript("2−"))),
        ("CO",    rich("CO", superscript("−"))) 
    ]
    fig = Figure(size = (750, 1500)) 
    times = result.tsol.t
    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    for i in 1:7
        sp_name, sp_label = target_species[i]
        
        sp_idx = findfirst(s -> s.name == sp_name, bulk)
        if isnothing(sp_idx)
            @warn "Species '$sp_name'를 bulk에서 찾을 수 없어 이 패널은 건너뜁니다."
            continue
        end

        # 1. 데이터 가져오기 및 하한선(1e-12) 적용
        @views c_matrix = result.tsol[sp_idx, 1:length(X), 1:length(times)] ./ scale
        
        # 2. bulk 차이 대신, 순수 농도의 log10 계산
        # 완전히 0인 값들이 무한대로 터지는 것을 막기 위해 max 사용
        M = log10.(max.(c_matrix, 1e-12))

        # 혹시 모를 비정상 값 필터링
        replace!(M, Inf => -12.0, -Inf => -12.0, NaN => -12.0)

        # 3. 각 종의 데이터에 맞게 Colorbar 범위 자동 조절
        # 값이 완전히 균일한 경우(예: 평형 상태 고정)를 대비해 최소 범위를 1.0 확보
        c_min = minimum(M)
        c_max = maximum(M)
        if c_min == c_max
            c_min -= 0.5
            c_max += 0.5
        end

        # fig[i, 1] 자리에 Axis 생성
        ax = Axis(fig[i, 1];
            xlabel = (i == 7) ? "Time (s)" : "", # 맨 아래만 x축 이름 표시
            ylabel = rich(rich("x", font=:italic), "  (m)"),
            yscale = log10,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(9),
            yticks = (10.0 .^ (-12:3:-6), [powlab(-12), powlab(-9), powlab(-6)]),
            title = "$sp_name Concentration Contour" 
        )

        # 히트맵 그리기
        hm = heatmap!(ax, times, X .+ 1e-12, M';
            colorrange = (c_min, c_max),
            colormap = discrete_cmap,
            interpolate = false)

        # 컬러바 생성 (라벨을 log10(c) 형태로 변경)
        Colorbar(fig[i, 2], hm;
            label = rich("log", subscript("10"), "(", rich("c", font=:italic), 
                         subscript(sp_label), " / M)"),
            ticklabelsize = 14, labelsize = 14)
    end
    return fig
end

# ╔═╡ 358b1fba-2f7a-4a45-b2f7-e8fbbff725ee
function panel_log_contour!(panel_pos, cbar_pos, result, X, times, sp;
                            scale=mol/dm^3, num_levels=24)
    c_bulk = elydata_Gold_unc.c_bulk[sp.idx]
    log_c_bulk = log10(c_bulk / scale)

    M = [log_c_bulk - log10(max(result.tsol[sp.idx, ix, it] / scale, 1e-12))
         for ix in 1:length(X), it in 1:length(times)]
    c_min = 0.0
    c_max = log_c_bulk - log10(1e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)
    ax = Axis(panel_pos;
        xlabel = lab_time,
        ylabel = rich(rich("x", font=:italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (10.0 .^ (-12:3:-6),
                  [powlab(-12), powlab(-9), powlab(-6)]))

    hm = heatmap!(ax, times, X .+ 1e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false)

    Colorbar(cbar_pos, hm;
        label = rich("log", subscript("10"), "(", rich("c", font=:italic),
                     subscript("bulk"), ") − log", subscript("10"), "(",
                     rich("c", font=:italic), subscript(sp.label), ")"),
        ticklabelsize = 20, labelsize = 20)

    return ax, hm
end

# ╔═╡ 34857db0-530c-435b-81c7-fa4b111faddc
function QoverK(fig, panel_pos, result, bulk; scale=mol/dm^3, lw=4,
                xlabel = lab_time, legend_pos = nothing)
	
    times = result.tsol.t
    nt    = length(times)
	
    # Equilibrium Constant
    EqK_hco3 = (elydata_Gold_unc.c_bulk[ihco3] * elydata_Gold_unc.c_bulk[iohminus]) /
           		elydata_Gold_unc.c_bulk[ico3]
	EqK_co3 = (elydata_Gold_unc.c_bulk[ico2] * elydata_Gold_unc.c_bulk[iohminus]) /
           elydata_Gold_unc.c_bulk[ihco3]
	
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

# ╔═╡ 79e62b5b-7578-4c7e-80c4-90bd34168650
elydata_Gold_unc.ircompfactor

# ╔═╡ 2f0d56cb-0668-44d8-8cea-cb3b5c4036d2
function plot_combined_exp_sim_ivc(
    P_recs;
    redox_species = Dict(ico => 2),  # OH- 1개당 1 electron (CO2RR/HER 둘 다)
    include_capacitive = true,
    scale = cm^2/mA,
    sign = 1,                             # cathodic 음수 convention
    sim_linewidth::Real = 5
)
    wanted    = ["0.1", "0.5", "1.0"]
    pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]
    raw_cv = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
    sub    = Matrix(raw_cv[4:end, :])
    num_cv = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num_cv, :auto)
    npairs = size(num_df, 2) ÷ 2
    keep   = findall(in(wanted), pressures[1:npairs])
    n_exp = length(keep)
    n_sim = length(P_recs)
    pastel2  = cgrad([colorant"#F2728A", colorant"#5BA8E8"])
    cols_exp = [pastel2[t] for t in range(0, 1, length = max(n_exp, 1))]
    pastel3  = cgrad([colorant"#FFB3BA", colorant"#A3D8FF"])
    cols_sim = [pastel3[t] for t in range(0, 1, length = max(n_sim, 1))]
    fig = with_theme(electrochemistry_theme()) do
        return Figure(size = (800, 800))
    end
    ax_exp = Axis(fig[1, 1],
        ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"),
        xticklabelsvisible = false, xticksvisible = false)
    ax_sim = Axis(fig[2, 1],
        xlabel = rich(rich("ϕ", font=:italic), "  (V vs. SHE)"),
        ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"))
    linkxaxes!(ax_exp, ax_sim)
    ax_sim.xticks = -1.5:0.3:1.0
    ax_exp.limits = (nothing, (-5.5, 1.8))
    ax_exp.yticks = 1:-2:-5
  #  ax_sim.limits = ((-1.3, 0.9), (-1.0, 0.2))
   # ax_sim.yticks = 0.2:-0.5:-1.0
    # ---- exp ----
    for (k, j) in enumerate(keep)
        xcol, ycol = 2j - 1, 2j
        label_text = "$(pressures[j]) atm"
        lines!(ax_exp, num_df[!, xcol], num_df[!, ycol];
               color = cols_exp[k], linewidth = sim_linewidth)
        text!(ax_exp, label_text;
              position = (-0.79 - 0.03*k, 1.6 - 1.3*k),
              color = cols_exp[k], fontsize = 26, font = :bold)
    end
    # ---- tor ----
    for j in 1:n_sim
        p, rec = P_recs[j]
        # Faradaic: sum over redox species with electron counts
        n_t = length(rec.voltages)
        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(rec, idx)
        end
        # Capacitive (only meaningful in :ohmicdrop mode)
        I_C = zeros(n_t)
        if include_capacitive
            ely = elydata_Gold_unc   # rec 구조에 맞게 조정
            if ely.ircompensation == :ohmicdrop
                icc = ely.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in rec.tsol[1:end-1]]
            end
        end
        I = sign .* (I_F .+ I_C) .* scale .* 2
        label_text = "$(p) atm"
        lines!(ax_sim, rec.voltages, I;
               color = cols_sim[j], linewidth = sim_linewidth)
        text!(ax_sim, label_text;
              position = (-0.99 - 0.05*j, 0.7 - 0.7*j),
              color = cols_sim[j], fontsize = 26, font = :bold)
    end
    text!(ax_sim, "Theory";
          position = (-0.5, -0.8), color = :gray30,
          fontsize = 32, font = :bold)
    text!(ax_exp, "Experiment";
          position = (-0.5, -1.27), color = :gray40,
          fontsize = 32, font = :bold)
    Label(fig[1, 1, TopLeft()], "(a)";
          fontsize = 28, font = :bold, padding = (0, 5, 20, 0))
    Label(fig[2, 1, TopLeft()], "(b)";
          fontsize = 28, font = :bold, padding = (0, 5, 20, 0))
    rowgap!(fig.layout, 1, 15)
    return fig
end

# ╔═╡ 2947efab-e67e-41fa-ab91-7e3c1c09df2d
elydata_Gold_irc = ElectrolyteData(; 
								   ircompensation=:pseudopotentiostat,
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
								  	x_ref = [X[end], 0, 0],
									M     = getproperty.(bulk, :M),
								  	Γ_we  = Γ_we,
								  	Γ_bulk= Γ_bulk,
								    Ru    = R_u,
									actcoeff! = γ_mode,
								    C_gap = C_gap,
								    ϕ_pzc = ϕ_pzc,
									ircompfactor = user_input_cv.ircomp,
									#redoxreaction = we_breactions
										
									)

# ╔═╡ 196f0405-9c32-48b6-9ddb-df929fbcf8d6
elydata_Gold_irc

# ╔═╡ bd2aa378-713b-4995-aaf6-c2f3521da7a2
plot_cv_current_variedL(odr_L_result, elydata_Gold_irc; species = iohminus)

# ╔═╡ de3c2f2e-77fe-4d4f-9bbc-ae2cf8c0406f
elydata_Gold_odr = ElectrolyteData(; 
								   ircompensation=:ohmicdrop,
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
								  	x_ref = [X[end], 0, 0],
									M     = getproperty.(bulk, :M),
								  	Γ_we  = Γ_we,
								  	Γ_bulk= Γ_bulk,
								    Ru    = R_u,
									actcoeff! = γ_mode,
								    C_gap = C_gap,
								    ϕ_pzc = ϕ_pzc,
									ircompfactor = user_input_cv.ircomp,
									#redoxreaction = we_breactions
										
									)

# ╔═╡ f506e83f-fee3-4661-b30a-ff0205d0922d
function plot_cv_total_current(result, model;
                               redox_species::Dict{Int,Int},
                               include_capacitive::Bool = true,
                               scale = cm^2/mA,
                               sign::Int = 1,
                               color_gradient::Bool = true)

    n_t = length(result.voltages)

    # ---- Faradaic current: sum over all redox species ----
    # currents(result, idx) already multiplies by F; n_e applied here.
    I_F = zeros(n_t)
    for (idx, n_e) in redox_species
        I_F .+= n_e .* currents(result, idx)
    end

    # ---- Capacitive current: icc boundary species (only :ohmicdrop) ----
    I_C = zeros(n_t)
    if include_capacitive
        ely = elydata_Gold_odr
        if ely.ircompensation == :ohmicdrop
            icc = ely.icc
            node_we = 1   # working-electrode boundary node (Γ_we = 1)
            I_C = [u[icc, node_we] for u in result.tsol[1:end-1]]
        else
            @warn "Capacitive current only resolved in :ohmicdrop mode " *
                  "(current mode: :$(ely.ircompensation)); j_C set to 0."
        end
    end

    # ---- Total current with sign convention and unit scaling ----
    I_total = sign .* (I_F .+ I_C) .* scale

    # ---- Plot ----
    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (540, 440))
        a = Axis(f[1, 1],
                 ylabel = L"I \; (\mathrm{mA/cm^2})",
                 xlabel = L"\phi \; (\mathrm{V \; vs \; SHE})")
        return f, a
    end

    if color_gradient
        cols = RGBf.(range(0, 1, length = n_t), 0.0, 0.0)
        lines!(ax, result.voltages, I_total; color = cols)
    else
        lines!(ax, result.voltages, I_total)
    end

    return fig
end



# ╔═╡ 486da4cc-776c-40b8-987a-888bc7ee9614
elydata_odr_f0 = deepcopy(elydata_Gold_odr)

# ╔═╡ 3a63acb9-4d2b-4e61-a5fb-638f33ba5215
elydata_Gold_odr

# ╔═╡ 5f17b4f7-54d6-4ad0-9886-252854840a80
function activity_coefficient!(
    γ::AbstractVector,
	u,
    data,
    mode::Function;
)
    (; v, ip, pscale, p_bulk, M, M0, v0, κ, RT, nc, Mrel, tildev, v0, cspecies, rexp) = data

    if γ_mode == Stefan_γ!
        # Ringe et al. approach (volume fraction-based)
    	for ic in cspecies
        	γ[ic] = (1.0 / (1 - sum(u[i] * v[i] for i in 1:nc)))# / (mol/dm^3)))
    	end
    elseif γ_mode == DGML_γ!
        # Dreyer et al. approach
        p = u[ip] * pscale - p_bulk
        c0, barc = c0_barc(u, data)
        for ic in cspecies
		#println(c0,"\t", barc,"\t", tildev[ic],"\t", p, "\t",Mrel[ic])
			
            γ[ic] = exp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc)) #* (1 /barc)
		end

    else
		γ .= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))
	end

    return γ
end


# ╔═╡ 8a1047fa-e483-40d9-8904-7576f30acfb4
begin
	const γ_cache = DiffCache(zeros(nc), 14)
	
	function reaction(
		f, 
		u::VoronoiFVM.NodeUnknowns, 
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

		#q = zero(u[iϕ])
		#for i in data.cspecies
		#	q += data.z[i] * u[i]
		#end
		#f[iq] = u[iq] - data.F * q

		
		nothing
	end
end;

# ╔═╡ 91113083-d80e-4528-be41-82d10f6860fc
begin
	const ps_cache = DiffCache(zeros(18), 14)
	const us_cache = DiffCache(zeros(isurfaceend-isurfacestart+1), 14)
	
	function we_breactions(f, 
			u, 
			bnode, 
			data
		) where {Tval, Tv, Tc, Tp, Ti}
		(; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we, ε, cspecies) = data
				
		γ = get_tmp(γ_cache, u[ico2])
		γ_co2 	 	= activity_coefficient!(γ, u, data, γ_mode)[ico2]
		γ_co 	 	= activity_coefficient!(γ, u, data, γ_mode)[ico]
		σ 			= surface_charge(u, data, user_input_model.BC_Select)
		local_pH 	= -log10((u[ihplus] * γ[ihplus] / (mol/dm^3)))

	


		ps = get_tmp(ps_cache, u[iϕ])
		ps[paramsidx[Symbolics.rename(odesys.σ, :σ)]] = σ
		ps[paramsidx[Symbolics.rename(odesys.γCO2_aq, :γCO2_aq)]] = γ_co2 
		ps[paramsidx[Symbolics.rename(odesys.aH2O_g, :aH2O_g)]] = aH₂O 
		ps[paramsidx[Symbolics.rename(odesys.ϕ, :ϕ)]] = u[iϕ] 
		ps[paramsidx[Symbolics.rename(odesys.ϕ_we, :ϕ_we)]] = ϕ_we 
		ps[paramsidx[Symbolics.rename(odesys.local_pH, :local_pH)]] = local_pH 
		ps[paramsidx[Symbolics.rename(odesys.γCO_aq, :γCO_aq)]] = γ_co 
		ps[paramsidx[Symbolics.rename(odesys.βCOOHΔH2OΔele_t, :βCOOHΔH2OΔele_t)]] = 0.59 


		
		@views f_microkinetics!(
			f[isurfacestart:isurfaceend], 
			u[isurfacestart:isurfaceend],
			ps,
			nothing
		)

		
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
		boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = (ϕ_we - ϕ_pzc))

	elseif user_input_model.BC_Select == "Robin"
		potentialbcondition!(f,u,bnode,data,ϕ_we)
	else
		boundary_neumann!(f, u, bnode; species = ic, region, value = 0)
		potentialbcondition!(f,u,bnode,data,ϕ_we)
	 	for ic in cspecies 
			if ic == ico2 || ic == ico
	       		boundary_neumann!(f, u, bnode; species = ic, region = Γ_we, value = 0)
			else
				boundary_dirichlet!(f, u, bnode; species = ic, region = Γ_we, value = c_bulk[ic])
			end
	   end
	end
		
		
	if bnode.region == Γ_we 
		we_breactions(f, u, bnode, data)
	end


	
	return bulkbcondition(f, u, bnode, data; region = Γ_bulk)

end;

# ╔═╡ 2d714608-7993-4dec-9065-0eeb26b7d365
if unc
	pnpsys_unc = PNPSystem(grid; bcondition = pnp_bcondition, celldata = elydata_Gold_unc, reaction = reaction)
end

# ╔═╡ 9fb47b83-a853-4316-bb8d-30e65b16ef78
if unc
	pnpresult_unc = sweep(elydata_Gold_unc, grid, pnp_bcondition, reaction, sawtooth; nperiods = user_input_cv.nperiods, eneutral = true, tunnel = false)
end

# ╔═╡ 7da046bf-d3b1-43a0-bdba-89b4da2f6be3
if CV
	AuCO2RR_plots.plot_time_voltage_and_dt(pnpresult_unc, sawtooth)
end

# ╔═╡ 8510c7ea-762c-4db4-b76a-9773dc4bbedc
if CV
	path = AuCO2RR_plots.cv_conc_gif(pnpresult_unc, bulk, X; framerate=4)
	LocalResource(path)
end

# ╔═╡ 3d661549-a8d2-40b0-add8-b186193f90fe
if CV
	AuCO2RR_plots.plot_cv_model_vs_koper_facets(pnpresult_unc; species = ico)
end

# ╔═╡ 74c43d72-3a23-4a24-a4ae-8b18b245610a
if CV
	AuCO2RR_plots.plot_iv_with_experiment(pnpresult_unc, iohminus)
end

# ╔═╡ 4b67ecff-c4d4-4ac8-b38f-4e7d6b180486
if irc
	pnpsys_irc = PNPSystem(grid; bcondition = pnp_bcondition, celldata = elydata_Gold_irc, reaction = reaction)
end

# ╔═╡ 28337c8b-8687-4618-bd9d-16d94193dd5f
if irc
	pnpresult_irc = sweep(elydata_Gold_irc, grid, pnp_bcondition, reaction, sawtooth; nperiods = user_input_cv.nperiods, eneutral = true, tunnel = false, )
end

# ╔═╡ f8b13f30-4b84-41be-8968-13c79ca62e32
if irc lot_conc_time_electrode(pnpresult_irc, bulk) end

# ╔═╡ 705aab5e-b425-4da1-8c05-0b547d70b5e6
if odr 
	pnpsys_odr = PNPSystem(grid; bcondition = pnp_bcondition, celldata = elydata_Gold_odr, reaction = reaction)
end

# ╔═╡ e8e6cff2-574f-4863-8093-9a747e419cf7
if odr 
	pnpresult_odr = sweep(elydata_Gold_odr, grid, pnp_bcondition, reaction, sawtooth; nperiods = user_input_cv.nperiods, eneutral = true, tunnel = false)
end

# ╔═╡ 89087b9d-04e0-4844-be32-fa5411900fe5
begin
	ratio = pnpresult_odr.j_cap ./ pnpresult_odr.j_dsp
	extrema(ratio) 
	
end

# ╔═╡ 230eb388-5cca-45d8-ac53-e150f917f2bc
plot_cv_total_current_tot(pnpresult_odr, elydata_Gold_odr; electrolyte=elydata_Gold_odr, redox_species = Dict(iohminus => 1), co_idx = ico)

# ╔═╡ a5679e25-8eb9-4b2e-8e6b-185853bdb7c1
if odr plot_cv_total_current(pnpresult_odr, elydata_Gold_odr, redox_species=Dict(ico => 2)) end

# ╔═╡ 5cb16d2e-d1cf-43cb-821a-e7478d927dcf
if odr plot_conc_time_electrode(pnpresult_odr, bulk) end

# ╔═╡ e0c8f10c-de89-4f44-93b1-105e9f387440
blthickness(grid, elydata_Gold_odr, pnpresult_odr.tsol)

# ╔═╡ a05cf724-cd32-498e-8afb-ecbf4a1f1648
if scan_rate_varied_checkbox
	sc, saw, scresult = scanrate_varied_sweep(elydata_Gold_unc, sawtooth, grid, pnp_bcondition, reaction; scanrates = [0.05, 0.5, 5], nperiods = user_input_cv.nperiods)
end

# ╔═╡ 980a13df-b578-4296-9cf5-2c5543d6b1f8
maximum(abs.(scresult[1].voltages .- scresult_2[1].voltages))

# ╔═╡ 3c476260-2fb5-4c87-9cbe-955ad1949671
AuCO2RR_plots.plot_scanrate_sweeps_cv_2(scresult, sc; electrolyte=elydata_Gold_unc, redox_species = Dict(iohminus => 2))

# ╔═╡ 2e642590-305c-4059-b756-16227f36a71c
AuCO2RR_plots.plot_scanrate_sweeps_cv(  scresult,  sc; )

# ╔═╡ 126f65a1-c02f-4f8e-be0a-40b1baa6b697
plot_cv_scanrate_grid(scresult, elydata_Gold_unc; electrolyte = elydata_Gold_odr, redox_species=Dict(ico => 2), co_idx= ico)

# ╔═╡ 2c1153a7-49e9-4661-ba22-7e6bca236f96
AuCO2RR_plots.plot_cv_scanrate_grid_unc(scresult, elydata_Gold_unc; electrolyte = elydata_Gold_unc, redox_species=Dict(ico => 2), co_idx= ico)

# ╔═╡ 4eac37f2-b364-43ef-a57d-ce00025c07d7
AuCO2RR_plots.plot_scanrate_sweeps_cap(scresult, sc; )

# ╔═╡ 76099e9c-018b-4468-9d32-1767514ee5ec
P_recs_odr = pressure_varied_sweep(elydata_Gold_unc, grid, pnp_bcondition, reaction, sawtooth; Pvec = [0.1, 0.5, 1], ispec = ico2)

# ╔═╡ 7036f45b-b0ba-4013-935e-756523f0f646
the = plot_pressure_varied_sweep_ivc(P_recs_odr, species = iohminus)

# ╔═╡ 4c1d7bef-f890-4078-89f0-85114df11216
(expfig, the)

# ╔═╡ cdca328c-c7c6-42d9-953e-c3f5fda73e77
plot_combined_exp_sim_ivc(P_recs_odr)

# ╔═╡ 012b426e-3550-4678-a666-eeb4bd32d20e
function simulate_CO2R_dir(grid, celldata; voltages = (-1.5:0.1:0.0) * V, kwargs...)
	kwargs 	 	= merge(solver_control, kwargs) 
    cell        = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)
	ivresult    = ivsweep(cell; voltages, store_solutions=true, kwargs...)

	cell, ivresult
end;

# ╔═╡ 06ca7d68-9c76-4d1c-947a-dd64a0fe9ec3
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
	    boundary_neumann!(f, u, bnode, species = iϕ, region = Γ_we, value = C_gap * (ϕ_we - ϕ_pzc))
	
	end

    return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ 924f8f5d-2cb0-4381-a522-509ff4c002b6
begin
	model_key = user_input_model[:model_choice]
	model = model_key == "Gold_Model" ? elydata_Gold_unc :
	        model_key == "Landstorfer_NaClO₄ model" ? elydata_NaClO₄ :
	        error("Unknown model choice: $model_key")
	molarities = [0.005, 0.1]#, 0.05, 0.1, 0.5] 
end;

# ╔═╡ 165446a0-da6c-4ecc-a5ba-c96eb300af12
const c0_bulk   = c̄ - sum(model.c_bulk) # solvent bulk molar concentration

# ╔═╡ 69ed6449-5caa-441b-9d8b-e38a9fce1f1f
function pb_flux(y, u, edge, data)
    #(; Γ_we, Γ_bulk, ϕ_we, iϕ, ip, ε) = data
	return y[model.iϕ] = elydata_NaF.ε * ε_0 * (u[model.iϕ, 1] - u[model.iϕ, 2])
end

# ╔═╡ 7fc5e2a3-c217-4042-9a42-e66d547bef96
is_Landstorfer = model != elydata_Gold_unc

# ╔═╡ 084e2127-ea77-4894-8990-380c2e8802c7
if double_layer_curve
	#pb
	sys_pb = PBSystem(grid; bcondition = pb_bcondition,  celldata = deepcopy(model))
	result_pb = capscalc(sys_pb, is_Landstorfer; vrange = range(vmin, vmax, length = 201))
else 
	result_pb = nothing
end

# ╔═╡ dbd3471a-1b80-4ad5-aa2d-aa2f486af86d
if double_layer_curve 
	fig_pb, ax_pb = AuCO2RR_plots.capsplot_fixed(result_pb, "Poisson-Nernst-Planck", xlimits_L=(-1.0, 1.0), ylimits_L=(0, 200)) 
	fig_pb
end

# ╔═╡ 3ef57b7d-ec19-46bc-a881-0506cf5167f3
if double_layer_curve
	#pnp
	#reaction_arg = model == elydata_Gold ? (reaction) : NamedTuple()
	sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = deepcopy(model), reaction)

	result_pnp = capscalc(sys_pnp, is_Landstorfer; vrange = range(vmin, vmax, length = 201))
else
	result_pnp = nothing
end

# ╔═╡ e115560a-f79e-4a08-9bd2-5c11b0346827
if double_layer_curve
	fig_pnp, ax_pnp = AuCO2RR_plots.capsplot_fixed(result_pnp, "Poisson-Nernst-Planck", xlimits_L=(-1.0, 1.0), ylimits_L=(0, 200))
	fig_pnp
end

# ╔═╡ 5a0c57b6-75cc-421c-8369-6531fb4d5ef7
let
    result = pnpresult_unc
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 1260))          # 6패널이니 높이 ↑ (1050 → 1260)
        ax1, leg1 = panel_conc_time!(f, f[1, 2], result, bulk;     xlabel = "")
        ax2 = panel_time_current!(f, f[2, 2], result, model;      xlabel = "")
        ax3 = panel_time_voltage!(f, f[3, 2], result;             xlabel = "")
        ax4 = panel_time_ph!(f, f[4, 2], result;                  xlabel = "")   # (d) pH
        #ax5, hm = panel_co2_log_contour!(f, f[5, 2], f[5, 3], result, X, bulk; L_val = L)  # (e)
        ax5, leg5 = QoverK(f, f[5, 2], result, bulk; legend_pos = f[5, 3])       # (f) Q/K

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

# ╔═╡ 1b809579-668b-4efa-a4f0-7352d8dfa5d5
let
    result = pnpresult_unc
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 1260))          # 6패널이니 높이 ↑ (1050 → 1260)
        ax1, leg1 = AuCO2RR_plots.panel_conc_time!(f, f[1, 2], result, bulk;     xlabel = "")
        ax2 = AuCO2RR_plots.panel_time_current!(f, f[2, 2], result, model;      xlabel = "")
        ax3 = AuCO2RR_plots.panel_time_voltage!(f, f[3, 2], result;             xlabel = "")
        ax4 = AuCO2RR_plots.panel_time_ph!(f, f[4, 2], result;                  xlabel = "")   # (d) pH
        #ax5, hm = panel_co2_log_contour!(f, f[5, 2], f[5, 3], result, X, bulk; L_val = L)  # (e)
        ax5, leg5 = AuCO2RR_plots.QoverK(f, f[5, 2], result, bulk; legend_pos = f[5, 3])       # (f) Q/K

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

# ╔═╡ f19fc2e0-5bae-4f2f-814f-9e53cbcc5654
if unc CV_total_current(pnpresult_unc, model; ) end

# ╔═╡ e6de6e89-5dc6-4fac-a7fc-8346175283a3
if unc plot_cv_current(pnpresult_unc, model; species = ico) end

# ╔═╡ d242507d-d1bb-461f-b4fe-ab3f597d9c40
if unc plot_activity_time_electrode(pnpresult_unc, bulk, model) end

# ╔═╡ 9ada29a1-4715-4f65-bea9-fc0863997dd8
if irc CV_dsp_cap_result(pnpresult_irc, model) end

# ╔═╡ a3967d7c-2b6b-4f2c-8401-b46aeefb7025
if irc plot_cv_current(pnpresult_irc, model; species = ico) end

# ╔═╡ 701cb999-6955-4e25-bfa7-01b4c868d0d0
if irc plot_activity_time_electrode(pnpresult_irc, bulk, model) end

# ╔═╡ 4196e28e-b787-4573-bb23-2f68a4574304
if odr CV_dsp_cap_result(pnpresult_odr, model) end

# ╔═╡ d341ec4a-ac49-430a-8b32-ff0925adf11b
if odr plot_cv_current(pnpresult_odr, model; species = ico) end

# ╔═╡ 699c75b4-09ec-4349-b660-0594da78f932
if odr plot_activity_time_electrode(pnpresult_odr, bulk, model) end

# ╔═╡ 7b66e057-4ab4-4e90-b576-7b496a956d2c
if unc && irc && odr 
	CV_overlay_currents([pnpresult_unc, pnpresult_irc, pnpresult_odr], model; species = iohminus)
end

# ╔═╡ cbcd2ff5-3df2-40df-b0c5-e127f3abb964
if CV_comp
	cvL = cvsweep_compensated_over_L(model, pnp_bcondition, reaction; sawtooth)
end

# ╔═╡ 35a462d3-1f8a-4a8e-ab8e-a07e844469e0
export_L_varied_species_csv_long(
	cvL;
	species = iohminus,
	function_name = "L_drop_V"
)

# ╔═╡ 2c7ca45b-8995-4bb9-8d47-43816caa594f
fig = plot_cv_current_dict(odr_L_result, model, scale = cm^2/mA, species= iohminus)

# ╔═╡ 11b12556-5b61-42c2-a911-4ea98a0a1e85
if IV 
	cell, ivresult = simulate_CO2R(grid_iv, model)
end

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
if IV
	(~, default_index) = findmin(abs, ivresult.voltages .+ 0.9 * ufac"V");
end

# ╔═╡ af333d3b-1e3a-4227-8cce-479907c11448
if IV
	gifpath = AuCO2RR_plots.plot1d_movie(ivresult; bulk, grid, L, step = 3, framerate = 5)
	LocalResource(gifpath)
end

# ╔═╡ 9498845e-fa44-4d01-a7bc-33d01ec11f79
if IV
	AuCO2RR_plots.plot_iv_with_ringe_refs(ivresult, species = iohminus)
end

# ╔═╡ 22244e24-5b56-4933-8a09-44b601f116c3
if IV
	iv_curve_axis(ivresult; cutoff=-0.4, showlegend=true, species = iohminus)
end

# ╔═╡ c6f10b66-6d06-4f2e-a7cc-780096d75785
if IV 
	conc_f, conc_a = AuCO2RR_plots.conc_vs_voltage_axis_compare(ivresult; bulk, grid, compare = comp)
	conc_f
end

# ╔═╡ f8255707-2233-4e28-b542-2f3d81b31c2e
if IV
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



	if γ_key == "Stefan_γ"
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
	else
  		text!(conc_ax, -1.2, 0.3, text=L"\mathrm{K^+}",
			  color=conc_colors[ikplus], fontsize=24, font="sans-bold")
		text!(conc_ax, -0.90, 0.0000004, text=L"\mathrm{H^+}", 
			  color=conc_colors[ihplus], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.05, 2e-9, text=L"\mathrm{CO_3^{2-}}",
			  color=conc_colors[ico3], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.0, 7e-6, text=L"\mathrm{HCO_3^-}", 
			  color=conc_colors[ihco3], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.21, 5e-4, text=L"\mathrm{CO_2}", 
			  color=conc_colors[ico2], fontsize=24, font = "sans-bold")
		text!(conc_ax, -0.9, 7e-9, text=L"\mathrm{OH^-}",
			  color=conc_colors[iohminus], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.15, 0.003, text=L"\mathrm{CO}", 
		 	  color=conc_colors[ico], fontsize=24, font = "sans-bold")

	end
    conc_fig
end

# ╔═╡ 9076385e-d6eb-4fc7-93d8-dca97195d003
electrode_activity_vs_voltage(ivresult, grid, model)

# ╔═╡ 47d92164-456a-43a0-8ed2-fed7f9fe31a2

if IV
activity = activity_vs_voltage_axis(ivresult;
							        bulk = bulk,
							        grid = grid, 
									electrolyte = model,
									ipressure = model.ip,
									model_type = γ_key
	)
    act_fig = activity.fig
    act_colors = activity.colors
    act_species = activity.species
	act_activity = activity.activity_electrode


	if γ_key == "Stefan_γ"
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
	else
  		text!(conc_ax, -1.2, 0.3, text=L"\mathrm{K^+}",
			  color=conc_colors[ikplus], fontsize=24, font="sans-bold")
		text!(conc_ax, -0.90, 0.0000004, text=L"\mathrm{H^+}", 
			  color=conc_colors[ihplus], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.05, 2e-9, text=L"\mathrm{CO_3^{2-}}",
			  color=conc_colors[ico3], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.0, 7e-6, text=L"\mathrm{HCO_3^-}", 
			  color=conc_colors[ihco3], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.21, 5e-4, text=L"\mathrm{CO_2}", 
			  color=conc_colors[ico2], fontsize=24, font = "sans-bold")
		text!(conc_ax, -0.9, 7e-9, text=L"\mathrm{OH^-}",
			  color=conc_colors[iohminus], fontsize=24, font = "sans-bold") 
		text!(conc_ax, -1.15, 0.003, text=L"\mathrm{CO}", 
		 	  color=conc_colors[ico], fontsize=24, font = "sans-bold")

	end
	act_fig

end

# ╔═╡ 223aa524-b341-49ca-811d-44cc12d39e94
activity

# ╔═╡ 0357d55d-0476-4a36-9d7e-3fee38c30674
act_activity

# ╔═╡ 3f3d5cfd-49d8-4d29-8e19-5048b8c1a82b
activity

# ╔═╡ 26f02407-ce54-436f-9630-63c0e2d32f73
if double_layer_curve
    outdir_pb = joinpath("..", "data", "output")
    isdir(outdir_pb) || mkpath(outdir_pb)

	if user_input_ion.use_physical_size == false
		ionsize_pb = "Potassium_only"
	else
		ionsize_pb = "All_species"
	end
	
    base_pb = string("DLCap_", user_input_model.BC_Select, "_", user_input_model.mode, "_pb_", ionsize_pb)

    for i in eachindex(result_pb)
        df = DataFrame(
            Voltage = result_pb[i].voltage_range,
            Capacitance = result_pb[i].dlcaps
        )
        fname = string(base_pb, "_", i, ".csv")
        filepath = joinpath(outdir_pb, fname)

        CSV.write(filepath, df)
    end
end

# ╔═╡ bc3077b4-6816-4214-a51f-6a5c9377bb0c
if double_layer_curve
    outdir_pnp = joinpath("..", "data", "output")
    isdir(outdir_pnp) || mkpath(outdir_pnp)
	
	if user_input_ion.use_physical_size == false
		ionsize_pnp = "Potassium_only"
	else
		ionsize_pnp = "All_species"
	end
	
    base_pnp = string("DLCap_", user_input_model.BC_Select, "_", user_input_model.mode, "_pnp_", ionsize_pnp)

    for i in eachindex(result_pb)
        df = DataFrame(
            Voltage = result_pnp[i].voltage_range,
            Capacitance = result_pnp[i].dlcaps
        )
        fname = string(base_pnp, "_", i, ".csv")
        filepath = joinpath(outdir_pnp, fname)

        CSV.write(filepath, df)
    end
end

# ╔═╡ 360313f0-2dad-4f7f-8d11-1800c6d934b3
base_cv = string("CV_", user_input_model.BC_Select, "_", user_input_model.mode, ionsize_pb)

# ╔═╡ 28638585-e95c-4947-9143-9ac8d8202f80
if extractIV
    outdir_pc = joinpath("..", "data", "output")
    isdir(outdir_pc) || mkpath(outdir_pc)

	if user_input_ion.use_physical_size == false
		ionsize = "Potassium_only"
	else
		ionsize = "All_species"
	end

    base_pc = string("Polarization_Curve_", user_input_model.BC_Select, "_", user_input_model.mode, "_pnp_", ionsize)

        df = DataFrame(
            Voltage = ivresult.voltages,
            Current = currents(ivresult, iohminus)
        )
        fname = string(base_pc, ".csv")
        filepath = joinpath(outdir_pc, fname)

        CSV.write(filepath, df)
end

# ╔═╡ 7e102647-23a9-4f40-b6de-cb1938bbb23e
	base_cvname = string("ad", "_σ_", 1000 , user_input_model.BC_Select, "_", user_input_model.mode, "_pnp_", ionsize, "_Scanrate_", user_input_cv.scanrate, "_Periods_", user_input_cv.nperiods)

# ╔═╡ a34cdba3-38f6-4bc0-b26e-72c956599109
function filename(function_name) 
	σ = round(L / μm)
	base_cvname = string(function_name, "_σ_", σ , "_", user_input_model.BC_Select, "_", user_input_model.mode, "_pnp_", ionsize, "_Scanrate_", user_input_cv.scanrate, "_Periods_", user_input_cv.nperiods, "_sweep_range_", user_input_cv.vmin,"-",user_input_cv.vmax, "_cv")

	return base_cvname
end

# ╔═╡ 406fb8e5-61e6-4688-96bf-a5530f57d4fa
function export_scanrate_varied_species_csv_long(
    saws,
    scresults;
    species=iohminus,
    outdir::AbstractString = "../data/output",
    function_name::AbstractString = "scanrate_varied",
    current_scale = cm^2 / mA,
    scanrate_scale = 1.0,
    scanrate_get = saw -> saw.scanrate,   # <- how to extract scan rate from saw
)
    @assert length(saws) == length(scresults)

    scanrates = Float64[]
    voltages  = Float64[]
    values    = Float64[]

    for (saw, rec) in zip(saws, scresults)
        sr = Float64(scanrate_get(saw) * scanrate_scale)

        V = rec.voltages
        Y = currents(rec, species) .* current_scale

        @assert length(V) == length(Y)

        append!(scanrates, fill(sr, length(V)))
        append!(voltages, V)
        append!(values, Y)
    end

    df = DataFrame(
        ScanRate = scanrates,
        Voltage  = voltages,
        Value    = values,
    )

    fname = filename(function_name)
    outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    return df
end

# ╔═╡ 5ccb0682-09de-4ae3-95e6-a8403d540d9b
if scan_rate_varied_checkbox
	export_scanrate_varied_species_csv_long(
	    saw,
	    scresult;
	    species=iohminus,
	    function_name="scanrate_varied_iohminus"
	)
end

# ╔═╡ 9e44f14b-a799-4ca5-8641-a3726780a4fe
function export_cv_profile_csv(
    rec;
    species = iohminus,
    outdir::AbstractString = "../data/output",
    function_name::AbstractString = "cv_profile",
    current_scale = cm^2 / mA,
)
    V = rec.voltages
    I = currents(rec, species) .* current_scale

    @assert length(V) == length(I)

    df = DataFrame(
        Voltage = V,
        Current = I,
    )

    fname = filename(function_name)
    outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    return df
end

# ╔═╡ fd2fe768-020c-4751-bde0-8f75843580f7
if CV
	export_cv_profile_csv(pnpresult)
end

# ╔═╡ e6dca43d-69b5-4d35-903c-6742b16a4715
function export_pressure_varied_species_csv_long(
    P_recs;
    species=iohminus,
    outdir::AbstractString = "../data/output",
    function_name::AbstractString = "pressure_varied",	
    scale=cm^2/mA,
)
    pressures = Float64[]
    voltages  = Float64[]
    values    = Float64[]

    for (p, rec) in P_recs
        V = rec.voltages
        Y = currents(rec, species) .* scale

        @assert length(V) == length(Y)

        append!(pressures, fill(Float64(p), length(V)))
        append!(voltages, V)
        append!(values, Y)
    end

    df = DataFrame(
        Pressure = pressures,
        Voltage  = voltages,
        Value    = values,
    )
	fname = filename(function_name)
	outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    return df
end

# ╔═╡ 4231590b-18c4-4af8-b798-364c529e47d8
if pressure_varied_checkbox
	df_out = export_pressure_varied_species_csv_long(
	    P_recs;
	    species=iohminus,
	    scale=cm^2/mA,
	)
end

# ╔═╡ f000e6ad-0224-43f4-bdf5-1c0c4308160f
if extractIV
	outdir_ac = joinpath("..", "data", "output")
	isdir(outdir_ac) || mkpath(outdir_ac)

	if user_input_ion.use_physical_size == false
	    ionsize_ac = "Potassium_only"
	else
	    ionsize_ac = "All_species"
	end

	base_ac = string(
	    "Activity_Curve_",
	    user_input_model.BC_Select, "_",
	    user_input_model.mode, "_pnp_",
	    ionsize_ac
	)
	
	Vol_ac = activity.vgrid
	Act_ac   = activity.activity_electrode
	
	spc_ac = ["K⁺", "H⁺", "HCO₃⁻", "CO₃²⁻", "CO₂", "OH⁻", "CO"]
	
	@assert size(Act_ac,1) == length(spc_ac)
	@assert size(Act_ac,2) == length(Vol_ac)
	
	df_ac = DataFrame(Voltage = Vol_ac)
	
	for i in 1:length(spc_ac)
	    df_ac[!, spc_ac[i]] = vec(Act_ac[i, :])
	end
	
	fname_ac = string(base_ac, ".csv")
	filepath_ac = joinpath(outdir_ac, fname_ac)
	
	CSV.write(filepath_ac, df_ac; bom=true)

end

# ╔═╡ d52f9cb8-9a48-43de-ab75-92380f7aca2e
if extract_conc	
    outdir_conc = joinpath("..", "data", "output")
    isdir(outdir_conc) || mkpath(outdir_conc)
	
	Volts  = conc_out.vgrid               
	Concent = conc_out.conc_electrode        
	sp = String.(conc_out.species)    
	
	@assert size(Concent, 1) == length(sp) == 7
	@assert size(Concent, 2) == length(Volts)
	
	df_conc = DataFrame(Voltage = Volts)
	
	for i in 1:length(sp)
	    df_conc[!, sp[i]] = vec(Concent[i, :])
	end
	base_conc = string("Concentration_", user_input_model.BC_Select, "_", user_input_model.mode,"_",ionsize ,".csv")

    filepath_conc = joinpath(outdir_conc, base_conc)

	CSV.write(filepath_conc, df_conc; bom=true)  


	
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
if IV
	AuCO2RR_plots.plot1d(ivresult, vshow; bulk, grid, L)
end

# ╔═╡ ab302d08-1a6f-4553-85af-043c565b107f
if IV
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
# ╠═f7d13047-4007-47ac-a3bb-a0b788dcd141
# ╟─beae1479-1c0f-4a55-86e1-ad2b50174c83
# ╠═ab2184fc-0279-46d9-9ee4-88fe3e732789
# ╠═7316901c-d85d-48e9-87dc-3614ab3d81a5
# ╠═6b7cfe87-8190-40a5-8d25-e39ef8d55db5
# ╠═5a146a44-03dc-45f3-ae15-993d11c2edac
# ╠═00947475-c96e-4ecc-a1ef-5be5e3e3c864
# ╠═165446a0-da6c-4ecc-a5ba-c96eb300af12
# ╟─a4b1300f-7e8a-46ab-9efd-dba82315a966
# ╠═05ebcce5-6904-451f-bbd6-ea4588ca05e3
# ╠═2c061d99-da09-4ccf-a4ec-73a3c9a908fe
# ╠═d5c6769d-35f9-461a-aae2-00122ccddd63
# ╠═303a98e5-20f5-4e59-928a-a305bf9a0b30
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
# ╟─a4928999-96c7-4145-af2d-485ff13d7109
# ╠═6b5cf93c-0df3-4a18-8786-502361736838
# ╟─d2c0642d-dfa5-4a76-bd36-ac4a735a3299
# ╟─06d45088-ab8b-4e5d-931d-b58701bf8464
# ╠═91113083-d80e-4528-be41-82d10f6860fc
# ╠═ee9229f8-2f48-4c1f-88cf-4bcc40f6e5f6
# ╠═960b1e96-da59-44c1-9828-929ece1a2955
# ╠═fb8a9a17-ed12-4f87-9957-06e5e265fcb2
# ╟─d0093605-0e35-4888-a93c-8456c698e6f0
# ╟─f0b5d356-6b97-4878-98de-bee5f380d41a
# ╟─e3eda42f-e2f3-4c10-81c4-610246ca528d
# ╟─7b87aa2a-dbaf-441c-9ad7-444abf15f664
# ╟─6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
# ╟─dabb7ca0-ce88-47d6-9315-9192d87c244e
# ╠═957d729c-d310-4d67-a429-1feb7770a2fd
# ╠═8f63ec58-bc97-43ba-9bbe-e10bb50e2cfe
# ╠═2947efab-e67e-41fa-ab91-7e3c1c09df2d
# ╠═de3c2f2e-77fe-4d4f-9bbc-ae2cf8c0406f
# ╠═196f0405-9c32-48b6-9ddb-df929fbcf8d6
# ╠═cd62234e-7f4a-4d0e-aea5-989f410d8cdc
# ╠═fe0c7939-3bda-4011-8aa3-5301cd23e302
# ╟─20324e33-26d0-4a95-bda3-5cb5d0fd8975
# ╠═848b7aeb-968f-4116-8038-b61276f02b6c
# ╟─53f12821-7d8d-4971-87fd-ad4689ec62a5
# ╠═9d814b85-a5b6-42e5-abf4-15500bbdb717
# ╟─e1e0ca0f-7f88-40f0-850e-590b25da0331
# ╠═0db74a70-af86-492c-affb-9de62ffe4455
# ╟─4f388fe0-6bc8-4a29-bccc-fa725e62c6a7
# ╠═5d179c52-43d7-4bcb-a2df-93c5806876fa
# ╠═5f17b4f7-54d6-4ad0-9886-252854840a80
# ╠═53b4dc3e-95f0-4eee-ba1c-68c222638acd
# ╠═dc203e95-7763-4b13-8408-038b933c5c9c
# ╠═9a4e01d9-f469-4427-bf4c-883adb67ae24
# ╟─d8c01196-7b81-43a5-ac0c-6dac8efee3d8
# ╠═69ed6449-5caa-441b-9d8b-e38a9fce1f1f
# ╠═95e72a20-a621-48a6-947e-0dcf9375facf
# ╟─5b900250-56a4-4e0b-bbc7-37bd39456312
# ╠═83ea329a-5d40-401d-8db1-59618c909f33
# ╟─dc5b1cfc-9221-45e9-a4f3-69ec4f7ebb70
# ╟─2a20d9be-6c1e-4c1f-8bb6-a7693800732d
# ╟─4f7ec19d-cd60-4c2b-a766-7557caa471c0
# ╠═924f8f5d-2cb0-4381-a522-509ff4c002b6
# ╠═084e2127-ea77-4894-8990-380c2e8802c7
# ╠═3ef57b7d-ec19-46bc-a881-0506cf5167f3
# ╠═26f02407-ce54-436f-9630-63c0e2d32f73
# ╠═bc3077b4-6816-4214-a51f-6a5c9377bb0c
# ╟─d76d8413-c019-4728-b182-7f7cb78dede4
# ╠═4656ee04-ae86-442f-b37c-c5563170f992
# ╠═7fc5e2a3-c217-4042-9a42-e66d547bef96
# ╠═e115560a-f79e-4a08-9bd2-5c11b0346827
# ╠═dbd3471a-1b80-4ad5-aa2d-aa2f486af86d
# ╟─9598e2c6-521e-4f8d-82d8-a836809736f3
# ╠═8bfdf2f5-c80a-4ce0-a8e1-b315affffb5f
# ╟─5c808c71-6094-49d7-8215-e88262f34e1f
# ╟─da8390d1-47e8-451f-b12b-45b8aca7b6ec
# ╠═ef4d7cee-fed1-489d-8795-c2dea9361e78
# ╠═e50fe651-11d4-45ee-89dd-371a7fbc097e
# ╟─ef7212fc-a3d0-4784-b901-219204b79dc0
# ╠═b95160b5-18f7-49d9-80be-9159abd2dcd1
# ╟─ad55a4c1-78b9-40d2-aec8-64ed95174e8a
# ╠═f6f26f7e-b97e-4f83-8b02-d3ff8b8ad14d
# ╠═f519311e-d5e7-4a4d-a080-80eb92e99ea4
# ╠═3f030783-2e0d-48c3-8407-5d317311c021
# ╠═2d714608-7993-4dec-9065-0eeb26b7d365
# ╠═9fb47b83-a853-4316-bb8d-30e65b16ef78
# ╠═7a172418-b46a-459c-9b17-ab3f4eec8dea
# ╠═4b67ecff-c4d4-4ac8-b38f-4e7d6b180486
# ╠═28337c8b-8687-4618-bd9d-16d94193dd5f
# ╠═7ff115dc-e982-4109-8db6-0bd718b13e88
# ╠═705aab5e-b425-4da1-8c05-0b547d70b5e6
# ╠═e8e6cff2-574f-4863-8093-9a747e419cf7
# ╠═7da046bf-d3b1-43a0-bdba-89b4da2f6be3
# ╠═89087b9d-04e0-4844-be32-fa5411900fe5
# ╟─41158868-8680-466b-a94b-9ffcd2d0ad8e
# ╠═2429e070-ed38-4a5a-9cdf-7818e7827021
# ╠═7bfe397e-2f49-461a-bca2-594e4bd3e607
# ╠═ad3a5236-72d0-4bbb-9fce-e74eb3ab52f0
# ╠═e7e678cb-0573-4c1d-8d04-978887c0b900
# ╠═13a0e5da-4b0e-4cf2-b43e-05f17802e48d
# ╠═3b8bc779-69f8-4f78-b55a-df90e226bcc9
# ╠═ac2bef83-9ebb-4bdc-869a-44cc1e2ef9e4
# ╠═22155420-3676-45f7-9bc8-8e00dfaf36be
# ╠═358b1fba-2f7a-4a45-b2f7-e8fbbff725ee
# ╠═34857db0-530c-435b-81c7-fa4b111faddc
# ╠═5a0c57b6-75cc-421c-8369-6531fb4d5ef7
# ╠═79e62b5b-7578-4c7e-80c4-90bd34168650
# ╠═1b809579-668b-4efa-a4f0-7352d8dfa5d5
# ╠═f506e83f-fee3-4661-b30a-ff0205d0922d
# ╠═230eb388-5cca-45d8-ac53-e150f917f2bc
# ╠═9949de26-8fad-4dc6-9635-8aec64c3f73d
# ╠═fd5c53ac-1a27-40f9-8ae4-cfbf1aaba5d1
# ╠═f19fc2e0-5bae-4f2f-814f-9e53cbcc5654
# ╠═e6de6e89-5dc6-4fac-a7fc-8346175283a3
# ╠═d242507d-d1bb-461f-b4fe-ab3f597d9c40
# ╠═6f90e65d-6804-4816-9ade-677a5984e4a8
# ╠═9ada29a1-4715-4f65-bea9-fc0863997dd8
# ╠═a3967d7c-2b6b-4f2c-8401-b46aeefb7025
# ╠═f8b13f30-4b84-41be-8968-13c79ca62e32
# ╠═701cb999-6955-4e25-bfa7-01b4c868d0d0
# ╠═36b271e4-23d8-4752-aca2-e5294a7ad419
# ╠═4196e28e-b787-4573-bb23-2f68a4574304
# ╠═d341ec4a-ac49-430a-8b32-ff0925adf11b
# ╠═a5679e25-8eb9-4b2e-8e6b-185853bdb7c1
# ╠═5cb16d2e-d1cf-43cb-821a-e7478d927dcf
# ╠═699c75b4-09ec-4349-b660-0594da78f932
# ╟─c7f0514e-ea30-4d99-ada2-c698360ac34c
# ╠═e6f43f01-15d8-4265-ae15-3673fb3cf7e3
# ╠═000ab285-b1bf-418b-b70c-40ff945cf2b8
# ╠═cfb28aaa-b9cf-4906-90a8-9ca65cbc1ab1
# ╠═8510c7ea-762c-4db4-b76a-9773dc4bbedc
# ╠═7b66e057-4ab4-4e90-b576-7b496a956d2c
# ╠═96c90e9f-9759-43d2-a13b-067d6a771cab
# ╠═bd2aa378-713b-4995-aaf6-c2f3521da7a2
# ╠═a2264942-aa50-4a63-823b-da26d491b040
# ╠═ef1195e2-6677-4cd4-bb4b-c66332c4cfee
# ╠═9d0b6083-e49b-4b15-ad2e-a2c74f689c5e
# ╠═3d661549-a8d2-40b0-add8-b186193f90fe
# ╠═74c43d72-3a23-4a24-a4ae-8b18b245610a
# ╠═3b07d9a6-bf83-43d0-89b6-64c112a94833
# ╠═fd2fe768-020c-4751-bde0-8f75843580f7
# ╟─c048e472-3983-4279-bf60-82784baa145e
# ╟─3bdaab98-c0f7-46af-86b7-d68374e8a5d0
# ╠═6d119626-ded5-4282-bc5a-0c37985697b0
# ╠═486da4cc-776c-40b8-987a-888bc7ee9614
# ╠═980a13df-b578-4296-9cf5-2c5543d6b1f8
# ╠═4d398e3f-1589-4ace-b0cb-273418a36fa2
# ╠═a05cf724-cd32-498e-8afb-ecbf4a1f1648
# ╠═70738a47-6476-4a87-b8a3-c9622519a1d3
# ╠═3c476260-2fb5-4c87-9cbe-955ad1949671
# ╠═2e642590-305c-4059-b756-16227f36a71c
# ╠═126f65a1-c02f-4f8e-be0a-40b1baa6b697
# ╠═02a78c8d-d557-4769-b155-719d607418c7
# ╠═2c1153a7-49e9-4661-ba22-7e6bca236f96
# ╠═08756476-b9d0-4bd4-be21-16cb93367dad
# ╠═4eac37f2-b364-43ef-a57d-ce00025c07d7
# ╠═d3b7d864-bc1b-440a-8ddb-2096bfde0cc5
# ╠═ae50302f-02de-4e30-aa47-1baaa2c95a74
# ╠═51064823-1b35-447d-95d0-2f9337e619eb
# ╠═e0e59ef0-8b6c-4f31-8d39-c2c4bcd7f99e
# ╟─56814250-16b2-4578-820d-2096998c84f4
# ╠═7b38e59a-d005-4cfc-ba8c-b17e7c700119
# ╠═b48ca0c7-acfc-4d63-a3a3-ec7d46c331d1
# ╠═3a63acb9-4d2b-4e61-a5fb-638f33ba5215
# ╠═76099e9c-018b-4468-9d32-1767514ee5ec
# ╠═92417c77-5ad2-451a-9848-44c9fe1f103b
# ╠═3e1bcc1a-01f9-472c-a3eb-c332113aafbc
# ╠═4c1d7bef-f890-4078-89f0-85114df11216
# ╠═7036f45b-b0ba-4013-935e-756523f0f646
# ╠═9ab242ff-1c4b-4970-8db1-877674838f29
# ╠═9c38b6f2-8f1b-49eb-bcdc-e3c792deed34
# ╠═ed4451bd-7888-4aea-8938-2a2ba85b22ad
# ╠═cdca328c-c7c6-42d9-953e-c3f5fda73e77
# ╠═2f0d56cb-0668-44d8-8cea-cb3b5c4036d2
# ╠═7844c654-bf05-4b19-9556-c0f3186efded
# ╟─5a1a5d5f-821d-47b8-b045-b6554313297c
# ╟─f5f9f81c-4f58-4a9b-9a6c-162d53ebb8f0
# ╟─a56f239e-e5b4-4689-9d56-719241924b19
# ╟─0fe0f681-966b-40ee-9323-da1fa8c95d80
# ╟─406fb8e5-61e6-4688-96bf-a5530f57d4fa
# ╟─53315f32-2fef-43c7-9baf-6d33bb848c70
# ╠═d020a8db-2227-4a62-8c99-fa6539a261a5
# ╠═64fe2609-843e-49be-92b6-b462f6bf80b3
# ╠═b767a48f-1b20-4b85-b9a8-36d6a914c5fe
# ╠═bb01b182-7840-4ad0-8cd6-fae57ab93173
# ╠═0a1054c5-cee7-4202-9d32-9eee6ec55265
# ╠═43e4ef6a-5c2a-4914-9972-ea93f3cb016e
# ╠═9e44f14b-a799-4ca5-8641-a3726780a4fe
# ╠═7e102647-23a9-4f40-b6de-cb1938bbb23e
# ╠═5ccb0682-09de-4ae3-95e6-a8403d540d9b
# ╟─a34cdba3-38f6-4bc0-b26e-72c956599109
# ╟─e6dca43d-69b5-4d35-903c-6742b16a4715
# ╠═35a462d3-1f8a-4a8e-ab8e-a07e844469e0
# ╠═cbcd2ff5-3df2-40df-b0c5-e127f3abb964
# ╟─9b8daa64-9e34-4da1-8136-ba5488e037c7
# ╠═360313f0-2dad-4f7f-8d11-1800c6d934b3
# ╠═4231590b-18c4-4af8-b798-364c529e47d8
# ╠═c9ae7a67-9a5f-4d9a-88c0-742f4e91fb27
# ╟─58ac8edc-2432-4054-88d8-52dafe0a2a61
# ╠═a2c7c4da-77cd-493f-8f98-0c86fecf271a
# ╠═2c7ca45b-8995-4bb9-8d47-43816caa594f
# ╠═e0c8f10c-de89-4f44-93b1-105e9f387440
# ╠═82116926-6739-4143-8cb7-29e34b8a323c
# ╠═6f6e779c-6ba3-49c5-9541-359ab4dbf13f
# ╠═8e2f8c2c-2c2a-46a1-90f3-8980cd6d8a51
# ╠═95bb0ec4-73c2-4735-b7a6-f54911881323
# ╠═fc22ce7f-fd42-4f31-9792-7268ebc2928d
# ╠═53fe324f-0024-4108-a855-66ab22952c3a
# ╠═09d54862-2fd1-424f-aa29-d8961b2e93a4
# ╠═ae8d6798-a3eb-4e6e-b810-a26e81ec1438
# ╠═d7e28d5f-16ba-4664-ba63-ec85fb29fe87
# ╠═d91bcc6a-1f2d-4457-9200-1ee30953090f
# ╠═d91c32c8-ac55-4f77-93cb-6c5d97bcbbb2
# ╠═d9690f36-f2db-4fcb-85f0-e285e8ed5551
# ╟─842b074b-f808-48d8-8dc5-110ddd907f90
# ╟─31298257-d35a-4f6f-8a76-ff00d5361ced
# ╠═cd310c1a-5810-4d13-9a71-96161e4c7452
# ╠═0c083fc8-0e73-4085-93c3-606493b3c46f
# ╠═72269ec4-a56e-46d9-85c8-0dd8ccaf43e1
# ╠═012b426e-3550-4678-a666-eeb4bd32d20e
# ╠═06ca7d68-9c76-4d1c-947a-dd64a0fe9ec3
# ╟─57db41d1-57c0-4eee-ba35-6f9b7e1e8263
# ╠═11b12556-5b61-42c2-a911-4ea98a0a1e85
# ╟─b976ab43-69f1-47a0-b2c6-c63e1c15cdb4
# ╠═60b410be-70f7-4053-a3db-7d777e0d3f08
# ╠═b6f0c4bb-153c-461c-b240-dff009b77dc9
# ╟─7a02463d-cfd9-4648-af53-f1e65d46733f
# ╟─114d2324-5289-4e44-8d77-736a9bdec365
# ╠═5ccb6a73-41cc-4107-a6ee-37d906313841
# ╟─c4876d26-e841-4e28-8303-131d4635fc23
# ╠═15fadfc2-3cf8-4fda-9aed-a79c602b1d51
# ╠═af333d3b-1e3a-4227-8cce-479907c11448
# ╠═425a5f53-ab8b-4596-bda7-586842e13878
# ╠═ab302d08-1a6f-4553-85af-043c565b107f
# ╠═9076385e-d6eb-4fc7-93d8-dca97195d003
# ╠═f2a1829d-e3d9-45c1-85a8-4ffb1564fe0f
# ╟─f8b5dc8f-1f41-4600-825e-2f9653f2d925
# ╠═9498845e-fa44-4d01-a7bc-33d01ec11f79
# ╠═94e4c398-b75e-4573-ae00-e54bdad267f6
# ╠═22244e24-5b56-4933-8a09-44b601f116c3
# ╠═15185eba-000c-4b6a-9692-903b6915bfd9
# ╠═c5a1bc1d-a2ff-418d-b8f1-a78e06a61c6b
# ╟─f672a256-641a-478e-b0aa-2df6e68b4d86
# ╠═bab42c91-2d00-463d-a921-97487e4eac67
# ╟─904ac4c2-50a8-4f70-8050-a0a1d4a448fa
# ╟─e4d93d39-c391-47ce-a248-6f0205761cca
# ╠═c6f10b66-6d06-4f2e-a7cc-780096d75785
# ╠═5d0458f3-6564-43df-af01-4a4b829ce262
# ╠═f8255707-2233-4e28-b542-2f3d81b31c2e
# ╠═47d92164-456a-43a0-8ed2-fed7f9fe31a2
# ╠═223aa524-b341-49ca-811d-44cc12d39e94
# ╠═11d6598c-3d6a-472a-8863-9275d5e567c6
# ╠═91051ed4-9fd0-4c21-95f4-efc042060e4d
# ╠═28638585-e95c-4947-9143-9ac8d8202f80
# ╠═f000e6ad-0224-43f4-bdf5-1c0c4308160f
# ╠═0357d55d-0476-4a36-9d7e-3fee38c30674
# ╠═3f3d5cfd-49d8-4d29-8e19-5048b8c1a82b
# ╠═43cb97f2-22ce-4940-8a17-3c04d4d2ced6
# ╠═d52f9cb8-9a48-43de-ab75-92380f7aca2e
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
