### A Pluto.jl notebook ###
# v0.20.25

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

# ╔═╡ 5a27d95a-21d2-4e8a-a2fd-895eed32b113


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

# ╔═╡ 9b7350ca-e176-4ebc-8a67-00d6737ab1a2


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

# ╔═╡ 960b1e96-da59-44c1-9828-929ece1a2955
function surface_charge(u, data, boundary)
    # vanilla LE 에서 mode 별 σ 가정 분리:
    #  :none → Robin BC 가 Stern relation enforce → σ = C_gap·(ϕ_we−ϕ_pzc−u[iϕ])
    #  :pseudopotentiostat / :ohmicdrop → library 가 Stern 없이 u[iϕ,1]≈ϕ_we 로 고정
    #    Stern 식 쓰면 σ 가 상수가 되어 microkinetics 가 applied 변화에 응답 안 함.
    #    Levey 식 Stern-less framework 에서 σ 는 applied overpotential 에만 의존.
    if data.ircompensation == :none
        return C_gap * (data.ϕ_we - ϕ_pzc - u[data.iϕ])
    else
        return C_gap * (data.ϕ_we - ϕ_pzc)
    end
end

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
    # Excluded-volume 자체-제동: c0 → 0⁺ 또는 음수면 ratio → ∞ → γ → ∞ (강한 제동)
    ratio = barc / max(c0, eps())
    barc_safe = max(barc, eps())
    for ic in cspecies
        γ[ic] = rexp(tildev[ic] * p / RT) * ratio^Mrel[ic] * (1 / (v0 * barc_safe))
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
    denom = max(1 - sum(c[i] * v[i] for i in 1:nc), eps())
    for ic in cspecies
        γ[ic] = 1.0 / denom
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

# ╔═╡ 02d12ba4-4ab3-48f6-b084-edb06cb413b1
# ╠═╡ disabled = true
#=╠═╡
begin
	celldata = deepcopy(model)
	pnpcell = PNPSystem(grid, celldata; bcondition = pnp_bcondition, celldata = celldata, reaction = reaction)
end
  ╠═╡ =#

# ╔═╡ ef7212fc-a3d0-4784-b901-219204b79dc0
md"""
#### General CV
"""

# ╔═╡ a72db110-37c3-4f4d-b37e-fe613c1ae827
md"""
##### **Run general cyclic voltammetry curve_with iR comp:** $(@bind CV_comp PlutoUI.CheckBox())
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
cv_sol_control = (; damp_initial = 0.01,    
				    damp_growth  = 1.1,    
			        Δu_opt = 0.03,
			        Δt_min = 3.0e-3,
			        Δt_max = 3.0e-2,
			        Δt = 3.0e-3,
			        Δt_grow = 1.2,
				 )

# ╔═╡ 28337c8b-8687-4618-bd9d-16d94193dd5f
# ╠═╡ disabled = true
#=╠═╡
pnpresult_irc = sweep(elydata_Gold_irc, grid, pnp_bcondition, reaction, sawtooth; nperiods = user_input_cv.nperiods, eneutral = false, tunnel   = false, Δt_min   = 1.0e-4)
  ╠═╡ =#

# ╔═╡ e8e6cff2-574f-4863-8093-9a747e419cf7
# ╠═╡ disabled = true
#=╠═╡
pnpresult_odr = sweep(elydata_Gold_odr, grid, pnp_bcondition, reaction, sawtooth; nperiods = user_input_cv.nperiods, eneutral = false, tunnel = false)
  ╠═╡ =#

# ╔═╡ 5b3ef1cb-c063-42b2-a0af-1e1a36d1e180


# ╔═╡ cdf52b70-94ad-45db-b82d-f1268cead86e
function plot_conc_time_electrode(result, bulk;
                                  nspecies=7,
                                  fig_size=(850, 500), 
                                  scale=(mol/dm^3))

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = fig_size)
    
    ax_conc = Axis(fig[1, 1],
                   xlabel = L"time / s",
                   ylabel = L"\text{Electrode Concentration}\;\;c_{i}^\ddagger / (\mathrm{mol/dm^3})",
                   limits = ((times[1]-(times[end] / 200), times[end] + (times[end] / 100)), (1e-12, 1e4)),
                   yscale = log10,
                   rightspinevisible = false) 
	cols 	= 	RGBf.(0.5, 0.5, 1)
    I 		= currents(result, iohminus) .* cm^2/mA
	ax_conc.xlabelsize = 25
    ax_current = Axis(fig[1, 1],
				ylabel = L"\text{Current density} / (\mathrm{mA/cm^{2}})",                      yaxisposition = :right,
                      ygridvisible = false, 
                      rightspinecolor = cols, 
                      ylabelcolor = cols,
                      yticklabelcolor = cols)

    linkxaxes!(ax_conc, ax_current)
	#ax_conc.xlabelsize = 25
	#ax_conc.ylabelsize = 25
	#ax_current.ylabelsize = 25

    for i in 1:nspecies
        y = conc[i, :]
        y_fixed = map(c -> (c > 1e-128 ? c : 1e-128), y)
        lines!(ax_conc, times, y_fixed; color=colors[i], label=string(names[i]))
    end

    lines!(ax_current, result.times, I ./ 2; color=cols, linewidth = 4, label="Current")

    Legend(fig[1, 2], ax_conc; labelsize=10, backgroundcolor=RGBA(1, 1, 1, 0.5))
    
    return fig
end

# ╔═╡ 2d421be5-9287-4e8d-9a7e-d63c95d9f370
md"""
t/s: $@bind time PlutoUI.Slider(range(tsol_irc.t[begin], tsol_irc.t[end], length = 101), show_value = true)
"""

# ╔═╡ c62b2553-e42f-4002-b206-9af0d8655c3d
#plottsol(
#    grid, celldata_Gold_unc, tsol_unc; species = celldata_Gold_unc.iO,
#	limits=(0,0.035)
#)

# ╔═╡ 089843dc-e9b0-4d43-afb7-d3477af39587
#plottsol(
#    grid, celldata_Gold_irc, tsol_irc; species = celldata_Gold_irc.iO,#
#	limits=(0,0.035)
#)

# ╔═╡ cdffb77b-89de-4aba-87b6-21e29bac0ae0
#plottsol(
#    grid, celldata_Gold_odr, tsol_odr; species = celldata_Gold_irc.iO,
#	limits=(0,0.035)
#)

# ╔═╡ 1efbe0ce-699f-4b6a-ab82-7b12d88e150b


# ╔═╡ 6d3dc6b1-82ea-448a-af86-ed6001ed5446


# ╔═╡ adfc0c3b-10d0-4f99-8b42-6da9a26de219


# ╔═╡ 239abbd9-928d-4f33-a9a1-517c986d4202
figscale=0.6

# ╔═╡ 6027ad72-7d72-47c2-9025-f00af8380c7f
cmax=0.005

# ╔═╡ fd5c53ac-1a27-40f9-8ae4-cfbf1aaba5d1
md"""
### CV unc`:=none` Result
"""

# ╔═╡ e6f43f01-15d8-4265-ae15-3673fb3cf7e3
function plot_activity_time_electrode(result, bulk, electrolyte;
    model_type="DGML_γ", # "DGML_γ" 또는 "Stefan_γ"
    nspecies=7,
    fig_size=(800, 500),
    scale=(mol/dm^3),
    ipressure=nothing 
)
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

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
               xlabel = L"time / s",
               ylabel = L"a_{i,\,\mathrm{electrode}} \text{ (Activity)}",
               limits = ((times[1]-(times[end] / 200), times[end] + (times[end] / 100)), (1e-12, 1e4)),
               yscale = log10)

    for i in 1:nspecies
        y = activity_electrode[i, :]
        y_fixed = map(val -> (isnan(val) || val <= 1e-128 ? 1e-128 : val), y)
        
        lines!(ax, times, y_fixed; color=colors[i], linewidth=2, label=string(names[i]))
    end

    Legend(fig[1, 2], ax; labelsize=10, backgroundcolor=RGBA(1, 1, 1, 0.5))
    return fig
end

# ╔═╡ 9d0b6083-e49b-4b15-ad2e-a2c74f689c5e


# ╔═╡ ec296961-151f-4219-a7f2-fb34e7facaa5


# ╔═╡ f1035efd-3e58-4c62-8799-a8750bdf137e


# ╔═╡ 25eb8aa3-697e-4538-9472-ceea45fbfbd9
md"""
#### pH varied CV function
"""

# ╔═╡ 11892724-1851-46f2-802d-4da45127b0af
md"""
Run pH varied CV calculation $(@bind pH_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ b2201bee-c72e-4498-bde0-7772dd2928c2
function plot_pH_varied_sweep(
    pH_recs;
    species = iohminus,
    fig_size = (1600, 900),
    scale = cm^2 / mA,
    legend_title = "Theoretical",
)

    fig = Figure(size = fig_size)
    ax = Axis(
        fig[1, 1],
        xlabel = L"\phi \; (\mathrm{V\ vs\ SHE})",
        ylabel = L"I \; (\mathrm{mA/cm^2})",
    )

    n = length(pH_recs)
    cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]

    plots  = Any[]
    labels = String[]

    for (j, item) in pairs(pH_recs)

        # ---- support both:
        # 1) old format: (pH, rec)
        # 2) new format: (pH=..., cH=..., cOH=..., c_bulk=..., record=...)
        pH, rec = if item isa NamedTuple
            if haskey(item, :record) && haskey(item, :pH)
                (item.pH, item.record)
            else
                error("NamedTuple input must contain at least :pH and :record")
            end
        elseif item isa Tuple
            if length(item) >= 2
                (item[1], item[2])
            else
                error("Tuple input must have at least 2 elements: (pH, rec)")
            end
        else
            error("Unsupported element type in pH_recs: $(typeof(item))")
        end

        ivres = hasproperty(rec, :ivresult) ? getproperty(rec, :ivresult) : rec
        I = currents(ivres, species) .* scale

        plt = lines!(ax, ivres.voltages, I; color = cols[j], linewidth = 3)
        push!(plots, plt)
        push!(labels, "pH = $(pH)")
    end

    Legend(fig[1, 2], plots, labels, legend_title; framevisible = true)

    return fig
end

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

# ╔═╡ a2c7c4da-77cd-493f-8f98-0c86fecf271a
md"""
Run boundary layer thickness varied cyclic voltammetry $(@bind L_varied_checkbox PlutoUI.CheckBox())
"""

# ╔═╡ 0db921e1-80ab-4d79-9035-df2dc33c0c3c
function sweep_COMP_single(L_val, pnpdata, grid, model, pnp_bcondition, sawtooth, reaction; f_comp = 1.0, eneutral = true, tunnel = false, bikerman = true, nperiods = 1)
    
    Area = 3.14e-6
    celldata = deepcopy(model)
    
    pnpcell = PNPSystem(grid, celldata; bcondition = pnp_bcondition, celldata = celldata, reaction = reaction)
    
    R_total = calc_R_comp(L_val, celldata; f_comp = f_comp, A = Area)
    R_comp = f_comp * R_total

    result = LiquidElectrolytes.cvsweep_COMP(
        pnpcell;
        Area = Area,
        voltages = sawtooth,
        nperiods = nperiods,
        R_comp = R_comp,
        store_solutions = true
    )
    return result
end

# ╔═╡ d7e28d5f-16ba-4664-ba63-ec85fb29fe87
begin
	grid_dict = Dict()
    X_coords_dict = Dict()
	
	L_values = [100, 500, 1000].* μm
	
    for L_val in L_values
        hmin = 1.0e-6 * μm
        hmax = 1.0    * μm
        X = ExtendableGrids.geomspace(0, L_val, hmin, hmax)
        g = ExtendableGrids.simplexgrid(X)
        
        grid_dict[L_val] = g
        X_coords_dict[L_val] = g[Coordinates][1, :] .+ 1e-12 # offset 포함
    end
	grid_dict
end

# ╔═╡ 6e88e1d8-1f4b-4813-8890-0cfcdc5fb967
if L_varied_checkbox
#	resL = sweep_over_L_cv(
#		elydata_Gold,
#		L_values,
#		pnp_bcondition,
#		sawtooth,
#		reaction; 
#		nperiods = user_input_cv.nperiods
#	)
end

# ╔═╡ 68e2cd9e-e8c0-491a-91c5-b30b16ed3f20


# ╔═╡ fbe4aca2-6a47-4457-98bb-588a5cde0ed5
md"""
#### Position at 0.99 Cbulk at a Given Time
"""

# ╔═╡ b1e64332-95a4-46a5-a45d-457c26e3fc67
const target_time = 20

# ╔═╡ 48029647-f162-459b-8824-fbf652d127f7
let
    try
        L_values = sort(collect(keys(RDE_result)))

        conc_vals = Float64[]

        for L in L_values
            rec = RDE_result[L]             
            result_L = rec.result           

            times = result_L.tsol.t
            idx   = argmin(abs.(times .- target_time))

            c_ico2 = result_L.tsol.u[idx][ico2, 1] / (mol / dm^3) 
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
              		verbose 	= "n",
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

# ╔═╡ b4aaf070-d4ab-409a-b1e8-f5469b9f398b
begin 
	sawtooth = SawTooth(
        scanrate = user_input_cv.scanrate,
        vmin = user_input_cv.vmin , vmax = user_input_cv.vmax,
		scanup = user_input_cv.scanup
    )
	    const nperiods = user_input_cv.nperiods
end

# ╔═╡ ed92cece-3f89-45f5-ac17-cbc9a9abb906
sawtooth

# ╔═╡ 87855181-ffbb-4a15-9852-0ded03746f18
sawtooth

# ╔═╡ 7da046bf-d3b1-43a0-bdba-89b4da2f6be3
if CV
	AuCO2RR_plots.plot_time_voltage_and_dt(pnpresult, sawtooth)
end

# ╔═╡ bf3a3756-59a9-4767-a713-034ef37ab55f
let
	    fig = Figure(size = (200,200))
		axV = Axis(fig[1, 1],
        	       xlabel = "Time (s)",
            	   ylabel = "Voltage")
		lines!(axV, pnpresult.times, sawtooth.(pnpresult.times))
	fig
end

# ╔═╡ fdbaa3a0-6b2b-45ad-8618-a63a7f720049
sawtooth.(pnpresult.times)

# ╔═╡ 0a1054c5-cee7-4202-9d32-9eee6ec55265
user_input_cv.nperiods

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

# ╔═╡ ab0e28f4-4310-4dcc-817e-81b9e45fd501
floataside(
    @bind user_input_model confirm(
        PlutoUI.combine() do Child
			md"""
			###### __Model Selection__  
			- Model: $(Child("model_choice", Select(["Gold_Model", "Landstorfer_NaClO₄ model", "Landstorfer_NaF model", "Toy model"])))
			---
			###### __Boundary layer thickness__   
			- ``δ``: $(Child("L", NumberField(1:80000; default = 100))) μm  			
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

# ╔═╡ 2d3e4cc2-7823-4fd3-b512-4d042f15cc68
function makegrid(celldata, L; hmin = 1.0e-6 * μm, hmax = 1.0 * μm)
	
    "Position of boundary BP1 (electrode surface)."
    x_BP1 = 0.0 * μm

    "Position of boundary BP2 (inner Helmholtz plane)."
    x_BP2 = 0.00029 * μm

    "Position of boundary BP3 (plane of electron transfer / outer Helmholtz plane)."
    x_BP3 = 0.00059 * μm

    "Approximate outer edge of the electric double layer."
    x_DL = user_input_model.L

    "Integer label for Domain 1."
    Domain1::Int = 1

    "Integer label for Domain 2."
    Domain2::Int = 2

    "Integer label for Domain 3."
    Domain3::Int = 3

    "Integer label for boundary point BP1."
    BP1::Int = 1

    "Integer label for boundary point BP2."
    BP2::Int = 2

    "Integer label for boundary point BP3."
    BP3::Int = 3

    "Integer label for boundary point BP4."
    BP4::Int = 4

    X_raw = ExtendableGrids.geomspace(0.0, L, hmin, hmax)
    
    X = sort(unique([X_raw; x_BP1; x_BP2; x_BP3; L]))

    grd = simplexgrid(X)
    
    safe_tol = 1.0e-12 

    cellmask!(grd, [x_BP1], [x_BP2], Domain1, tol = safe_tol)
    cellmask!(grd, [x_BP2], [x_BP3], Domain2, tol = safe_tol)
    cellmask!(grd, [x_BP3], [L], Domain3, tol = safe_tol)

    bfacemask!(grd, [x_BP1], [x_BP1], BP1, tol = safe_tol)
    bfacemask!(grd, [x_BP2], [x_BP2], BP2, tol = safe_tol)
    bfacemask!(grd, [x_BP3], [x_BP3], BP3, tol = safe_tol)
    bfacemask!(grd, [L], [L], BP4, tol = safe_tol)
    
    return X, grd

end


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
    hmax = 1.0	* μm 
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)
	grid_CV    = simplexgrid(geomspace(0, L, 1e-5*μm, hmax))
	#X, grid = makegrid(elydata_Gold_unc, L)
end

# ╔═╡ d5c6769d-35f9-461a-aae2-00122ccddd63
gridplot(grid; size = (750, 200))

# ╔═╡ 50dc3a93-6eeb-4f6f-ac42-65d8d597614c
blthickness(grid, celldata_Gold_unc, tsol_unc; species = 1)/ufac"nm"

# ╔═╡ ea02ccff-20ea-426d-b102-e93a30e8fdb4
blthickness(grid, celldata_Gold_irc, tsol_irc; species = 1)/ufac"nm"

# ╔═╡ ce8dd32b-2ea3-4a42-9436-69b06e751984
blthickness(grid, celldata_Gold_odr, tsol_odr; species = ico)/ufac"nm"

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
        a_HCO3 = 1.56   #  (Å)
        a_CO3  = 1.78   #  (Å)
        a_CO2  = 1.65   #  (Å) (often treated as vdW/effective in water)
        a_OH   = 1.40   #  (Å)
        a_H    = 0.01   #  (Å) (H3O+ effective hydrated)
        a_CO   = 1.88   #  (Å)
        a_K    = 1.38   #  (Å)

        # --- MD hydration numbers (1st shell) ---
        κ_HCO3 = 2.7     #  (MD)
        κ_CO3  = 4.25     #  (MD)
        κ_CO2  = 0.0     # neutral, typically treat as ~0 (no structured hydration number in this model)
        κ_OH   = 1.5     #  (MD)
        κ_H    = 2.0     #  (MD)
        κ_CO   = 0.0     # neutral
        κ_K    = 3.0     #   (MD)
		"""
	#elseif use_md_hydrated == true && γ_key == "Stefan_γ"
    #    # Stefan's Model / Consider all effective size
    #    a_HCO3 = 8.5;  κ_HCO3 = 0
    #    a_CO3  = 9.9;  κ_CO3  = 0
    #    a_CO2  = 3.4;  κ_CO2  = 0
    #    a_OH   = 7.6;  κ_OH   = 0
    #    a_H    = 7.3;  κ_H    = 0
    #    a_CO   = 2.8;  κ_CO   = 0
    #    a_K    = 8.2;  κ_K    = 0
		
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

# ╔═╡ e23eaa0a-a4ec-4c05-8526-15df6ad0eec1
    κ_bulk = F^2 / (R*T) * sum((getproperty.(bulk, :z) .^ 2) .* getproperty.(bulk, :D) .* getproperty.(bulk, :c_bulk))


# ╔═╡ a64e2dc9-9be7-48b5-9d04-c448f19ed7f2
if CV
	plot_conc_time_electrode(pnpresult, bulk)
end

# ╔═╡ b25f2246-0182-4d50-a606-0d81776d414f
if CV
	c_b = AuCO2RR_plots.plot_conc_profile_with_delta(pnpresult, bulk, X)
end

# ╔═╡ bf991ac4-8ff2-4a49-9f1e-16ff60b12493
c_b[1]

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


# ╔═╡ 8f63ec58-bc97-43ba-9bbe-e10bb50e2cfe
elydata_unc = ElectrolyteData(; ircompensation=:none,
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
									x_ref = [grid[Coordinates][end], 0.0, 0.0],   # f_ref ∈ (0,1)
								  	Γ_bulk= Γ_bulk,
								   #ϕ_pzc = ϕ_pzc,
									actcoeff! = γ_mode,
									ircompfactor = user_input_cv.ircomp,
									#redoxreaction = we_breactions
									)

# ╔═╡ e50fe651-11d4-45ee-89dd-371a7fbc097e
function sweep(model, grid, bcondition, reaction, sawtooth; nperiods = 1, eneutral = true, tunnel = false, kwargs...)
    celldata = deepcopy(model)
	#kwargs 	 	= merge(solver_control, kwargs)
    #celldata.eneutral = eneutral
	reaction_arg = model == elydata_unc ? (; reaction) : NamedTuple()
    pnpcell = PNPSystem(grid, celldata; bcondition = bcondition, reaction = reaction)
    return result = LiquidElectrolytes.cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
		kwargs...
    )

end

# ╔═╡ 7b38e59a-d005-4cfc-ba8c-b17e7c700119
if pressure_varied_checkbox
	P_recs = pressure_varied_sweep(elydata_Gold, sweep; Pvec = [0.1, 0.2, 0.3, 0.5, 0.6, 1], ispec = 5)
	AuCO2RR_plots.plot_pressure_varied_sweep(P_recs, species = iohminus)
end

# ╔═╡ c9ae7a67-9a5f-4d9a-88c0-742f4e91fb27
if pressure_varied_checkbox
	AuCO2RR_plots.plot_pressure_varied_sweep(P_recs, species = iohminus, limits = ((-0.6, 0.9), (0, 1)))
end

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
        denom = max(1 - sum(u[i] * v[i] for i in 1:nc), eps())
    	for ic in cspecies
        	γ[ic] = 1.0 / denom
    	end
    elseif γ_mode == DGML_γ!
        # Dreyer et al. approach
        p = u[ip] * pscale - p_bulk
        c0, barc = c0_barc(u, data)
        ratio = barc / max(c0, eps())     # c0→0⁺ 또는 음수 ⇒ ratio→∞ ⇒ γ→∞ (제동)
        barc_safe = max(barc, eps())
        for ic in cspecies
            γ[ic] = exp(tildev[ic] * p / RT) * ratio^Mrel[ic] * (1 / (v0 * barc_safe))
		end

    else
        denom = max(1 - v[ikplus] * u[ikplus] / (mol/dm^3), eps())
		γ .= 1.0 / denom
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
		# Excluded-volume overflow 시 H⁺ activity 가 폭주해 pH 가 -∞ 로 튀는 걸 방지
		# microkinetics 의 exp(rate(pH)) 가 NaN 으로 cascading 되는 것을 막기 위한 clamp
		local_pH 	= -log10(u[ihplus] * γ[ihplus] / (mol/dm^3))

	


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

# ╔═╡ 2947efab-e67e-41fa-ab91-7e3c1c09df2d
elydata_Gold_irc = ElectrolyteData(; ircompensation=:pseudopotentiostat,
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
									x_ref=[L * (1 - user_input_cv.ircomp), 0, 0],   
								   Ru=L/κ_bulk,
								  	Γ_bulk= Γ_bulk,
									actcoeff! = γ_mode,
									ircompfactor = user_input_cv.ircomp,
								    redoxreaction = we_breactions

									)

# ╔═╡ de3c2f2e-77fe-4d4f-9bbc-ae2cf8c0406f
elydata_Gold_odr = ElectrolyteData(; ircompensation=:ohmicdrop,
	                               	nc = size(bulk)[1],
									na    = na,
									z     = getproperty.(bulk, :z),
								  	D     = getproperty.(bulk, :D),
								  	T     = T,
								  	eneutral=false,
								  	κ     = getproperty.(bulk, :κ),
		                            c_bulk= getproperty.(bulk, :c_bulk),
								    v0 	  = v0,
								    Ru = grid[Coordinates][end] / (F^2 / (R*T) * sum((getproperty.(bulk, :z) .^ 2) .* getproperty.(bulk, :D) .* getproperty.(bulk, :c_bulk))),
									v     = getproperty.(bulk, :v),
									M0 	  = M0,
								    x_ref = [grid[Coordinates][end], 0.0, 0.0],
									M     = getproperty.(bulk, :M),
								  	Γ_we  = Γ_we,
								  	Γ_bulk= Γ_bulk,
									actcoeff! = γ_mode,
									ircompfactor = user_input_cv.ircomp,
									redoxreaction = we_breactions
									);

# ╔═╡ dc203e95-7763-4b13-8408-038b933c5c9c
function pnp_bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, iϕ, ϕ_we, ircompensation) = data
    if ircompensation == :none
        boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * ϕ_we)  # ϕ_we 는 이미 PZC-shifted
    end
    # :pseudopotentiostat / :ohmicdrop 는 library generic op 가 iϕ@1 처리
    if bnode.region == Γ_we
        we_breactions(f, u, bnode, data)
    end
    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
end


# ╔═╡ 9fb47b83-a853-4316-bb8d-30e65b16ef78
pnpresult_unc = cvsweep(PNPSystem(grid; bcondition=pnp_bcondition, celldata=elydata_unc, reaction=reaction); voltages=sawtooth)

# ╔═╡ f668c334-ba69-4ab9-b719-b4078305e68b
pnpresult_unc.times[end], pnpresult_unc.voltages[end], sawtooth(pnpresult_unc.times[end])


# ╔═╡ 2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
AuCO2RR_plots.plot_conc_profile_logx(pnpresult_unc, bulk, X, 80)

# ╔═╡ 2420382d-227a-4063-9450-1f1726df018e
begin
	path = AuCO2RR_plots.cv_conc_gif(pnpresult_unc, bulk, X; framerate=4)
		LocalResource(path)
end

# ╔═╡ 3d661549-a8d2-40b0-add8-b186193f90fe
	AuCO2RR_plots.plot_cv_model_vs_koper_facets(pnpresult_unc; species = iohminus)


# ╔═╡ 74c43d72-3a23-4a24-a4ae-8b18b245610a
	AuCO2RR_plots.plot_iv_with_experiment(pnpresult_unc, iohminus)


# ╔═╡ 017d8b3a-5923-411d-bfbd-31e3167e5ac0
pnpresult_irc = cvsweep(PNPSystem(grid; bcondition=pnp_bcondition, celldata=elydata_Gold_irc, reaction=reaction); voltages=sawtooth)

# ╔═╡ d93cd922-d021-40b1-8a4d-c9f2bfa4869a
pnpresult_irc.times[end], pnpresult_irc.voltages[end], sawtooth(pnpresult_irc.times[end])

# ╔═╡ 9e01784a-c575-41e4-b82a-6c65d78004c1
pnpresult_odr = cvsweep(PNPSystem(grid; bcondition=pnp_bcondition, celldata=elydata_Gold_odr, reaction=reaction); voltages=sawtooth)

# ╔═╡ ac192132-1778-471c-9e8d-958f3bfe400a
(unc_v0 = pnpresult_unc.voltages[1], unc_vend = pnpresult_unc.voltages[end], 
 irc_v0 = pnpresult_irc.voltages[1], irc_vend = pnpresult_irc.voltages[end],
 odr_v0 = pnpresult_odr.voltages[1], odr_vend = pnpresult_odr.voltages[end],
 sawtooth_at_end = sawtooth(pnpresult_unc.times[end]))


# ╔═╡ 606cc901-0c24-4588-be13-d05f55b6907f
let
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="ϕ_we (applied, V)", ylabel="I (mA/cm²)")
    for (rec, lbl, col) in [(pnpresult_unc, "unc", :black),
                            (pnpresult_irc, "irc", :red),
                            (pnpresult_odr, "odr", :blue)]
        x = sawtooth.(rec.times)
        I = currents(rec, iohminus) .* cm^2/mA
        lines!(ax, x, I; color=col, label=lbl)
    end
    axislegend(ax)
    fig
end


# ╔═╡ 3a07864e-9579-4fd4-a62a-d68472b9f9c6
pnpresult_odr.times[end], pnpresult_odr.voltages[end], sawtooth(pnpresult_odr.times[end])

# ╔═╡ 311660ce-0168-47eb-8d23-00a8e750fd84
function simulate(
        grid, 
        celldata,
		reaction; 
        Δu_opt = 2.0e-2,
        damp_initial = 0.5,
        Δt_min = 1.0e-11,
        max_round = 3,
        tol_round = 1.0e-8,
        verbose = "",
        handle_exceptions = true,
        vmin = -0.5ufac"V",
        vmax = 0.1ufac"V",
        nperiods = 1,
        scanrate = 1.0ufac"V / s",
        kwargs...
    )
	celldata = deepcopy(celldata)
    sys = PNPSystem(grid; bcondition = pnp_bcondition, celldata = celldata, reaction = reaction)
    
    initial_guess_vector = LiquidElectrolytes.unknowns(sys)
    
    inisol = VoronoiFVM.solve(
        sys.vfvmsys; 
        inival = initial_guess_vector, 
        damp_initial = 0.1, 
        time = 0.0
    )
    
    return result = LiquidElectrolytes.cvsweep(
        sys;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
    )
    
    return tsol, cvresult
end

# ╔═╡ b65575ce-7c5e-4a70-bb50-3bd8054516aa
function run(L_values; ircompensation = :pseudopotentiostat)
    tsols = []
    cvresults = []
    grids = []
    celldatas = []

    for L in L_values
        if ircompensation == :pseudopotentiostat
            celldata = celldata_Gold_irc
        elseif ircompensation == :none
            celldata = celldata_Gold_unc
        else
            celldata = celldata_Gold_odr
        end
        

        grd = makegrid(celldata)
        
        tsol, cvresult = simulate(
            grd, 
            elydata, 
            celldata; 
            vmin = user_input_cv.vmin,
            vmax = user_input_cv.vmax,
            nperiods = user_input_cv.nperiods,
            scanrate = user_input_cv.scanrate,
            Δu_opt = 2.0e-2,
            damp_initial = 0.5,
            Δt_min = 1.0e-11,
            max_round = 3,
            tol_round = 1.0e-9,
            handle_exceptions = true
        )
        
        push!(cvresults, cvresult)
        push!(tsols, tsol)
        push!(grids, grd)
        push!(celldatas, celldata)
    end
    
    return cvresults, tsols, grids, celldatas
end;

# ╔═╡ a05cf724-cd32-498e-8afb-ecbf4a1f1648
if scan_rate_varied_checkbox
	sc, saw, scresult = scanrate_varied_sweep(elydata_Gold, sawtooth, grid, pnp_bcondition, reaction; scanrates = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5], nperiods = user_input_cv.nperiods)
	AuCO2RR_plots.plot_scanrate_sweeps(scresult, sc; species = iohminus)
end

# ╔═╡ d91c32c8-ac55-4f77-93cb-6c5d97bcbbb2
if L_varied_checkbox
    @info "Starting COMP CV sweeps for $(length(L_values)) cases..."
    
    resL_COMP = Dict(L_val => begin
        @info "Running for L = $(round(L_val/μm)) μm..."
        
        X = ExtendableGrids.geomspace(0, L_val, 1e-6*μm, 1.0*μm)
        grid = ExtendableGrids.simplexgrid(X)
        
        sweep_COMP_single(
            L_val, elydata_Gold, grid, elydata_Gold, 
            pnp_bcondition, sawtooth, reaction; 
            nperiods = user_input_cv.nperiods
        )
    end for L_val in L_values)

    @info "All sweeps completed!"
end

# ╔═╡ ee83fb10-b4a1-4084-a644-0541500fc9cb
let
	if L_varied_checkbox
	    co2_idx = findfirst(s -> s.name == "CO₂" || s.name == "CO2", bulk)
	    selected_Ls = L_values[1:end]
	    num_panels = length(selected_Ls)
	
	    # 1. Figure Setup
	    fig = Figure(size = (1100, 380 * num_panels))
	    
	    # View limits (Interface focus)
	    max_L_view = 1e-2 
	    x_limits = (1e-12, max_L_view)
	    y_limits = (-6, 0.5) 
	
		colors = :Pastel1_7
		
	    for (idx, L_val) in enumerate(selected_Ls)
	        if !haskey(resL_COMP, L_val) continue end
	        
	        res = resL_COMP[L_val]
	        tsol = res.tsol
	        times = res.times
	        nt = length(times)
	        # --- GRID HANDLING ---
			X_coords = X_coords_dict[L_val]
	
	        # --- Axis Setup ---
	        ax = Axis(fig[idx, 1],
	            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
	            ylabel = L"\log_{10}(c_{\mathrm{CO_2}})",
	            xscale = log10,
	            limits = (x_limits, y_limits),
	            title = "Boundary Layer (L) = $(round(L_val/μm)) μm",
	            xgridvisible = false, ygridvisible = false
	        )
	
	        # 2. Time-step selection (e.g., 15 lines per panel)
			num_points = 100 
	        t_indices = round.(Int, range(1, nt, length=num_points))
	        
	        # --- Color setup using your idea ---
	        time_cols = resample_cmap(colors, length(t_indices))
	        
	        local line_obj
	        
	        for (i, ti) in enumerate(t_indices)
	            c_profile = [log10(max(tsol[co2_idx, ix, ti] / (mol / dm^3), 1e-25)) for ix in 1:size(tsol, 2)]
	            
	            line_obj = lines!(ax, X_coords, c_profile; 
	                color = time_cols[i], 
	                linewidth = 2.5
	            )
	        end
	
	        # 3. Add Colorbar to the right [idx, 2]
			Colorbar(fig[idx, 2], 
			    limits = ((times[1]), (times[end])), # Use first and last elements
	            colormap = colors,
	            label = "Time [s]",
	            # Ensure ticks is also a range of numbers, not vectors
	            ticks = range(times[1], times[end], length=5) 
	        )
	
	        # Bulk Boundary & Reference
	        vlines!(ax, [L_val], color = :red, linestyle = :dash, linewidth = 2)
	        hlines!(ax, [log10(max(tsol[co2_idx, end, 1] / (mol / dm^3), 1e-25))], 
	                color = :black, linestyle = :dot, alpha = 0.5)
	    end
	
	    rowgap!(fig.layout, 35)
	    colgap!(fig.layout, 25) 
	    fig
	end
end

# ╔═╡ 65e8592d-a3f8-4f38-9c7e-8085a5c28b73
let
	if L_varied_checkbox
	    # 1. CO2 Index 확인
	    co2_idx = findfirst(s ->  s.name == "CO₂", bulk)
	    selected_Ls = L_values[1:end]
	    num_panels = length(selected_Ls)
	
	    max_L = 1e-2
	
	    fig = Figure(size = (1100, 380 * num_panels))
	    
	    # --- COLOR RANGE 조절 ---
	    c_min = log10(0.00001)
	    c_max = log10(0.33)
	    color_range = (c_min, c_max)
	
	    for (idx, L_val) in enumerate(selected_Ls)
	        if !haskey(resL_COMP, L_val) continue end
			num_levels = 4 
			discrete_cmap = resample_cmap(:Pastel1_7, num_levels)        
	        res = resL_COMP[L_val]
	        tsol = res.tsol
	        times = res.times
	        # Grid Handling: Consistent with each L's specific mesh
	        # --- GRID HANDLING ---
			X_coords = X_coords_dict[L_val]
	        
	        conc_matrix = [log10(max(tsol[co2_idx, ix, it] / (mol / dm^3), 1e-25)) 
	                       for ix in 1:length(X_coords), it in 1:length(times)]
	
	        # --- Axis Setup ---
	        ax = Axis(fig[idx, 1],
	            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
	            ylabel = "Time [s]",
	            title = "Boundary Layer (L) = $(round(L_val/μm)) μm",
	            xscale = log10,
	            # 모든 패널의 X축 범위를 max_L로 고정하여 스케일 통일
	            limits = ((1e-10, max_L), (0, 80)),
	            xminorticksvisible = true, xminorticks = IntervalsBetween(9)
	        )
			time_cols = 
	        # Heatmap
	        hm = heatmap!(ax, X_coords .+ 1e-12, times, conc_matrix; 
	            #limits = ((1e-10, max_L), (nothing, nothing)),
	            colorrange = color_range, 
	            colormap = Reverse(discrete_cmap),
	            interpolate = false
	        )
	
	        # 물리적 경계(Bulk Boundary)를 시각적으로 표시 (선택 사항)
	        vlines!(ax, [L_val], color = :red, linestyle = :dash, linewidth = 2)
	
	        Colorbar(fig[idx, 2], hm, label = L"\log_{10}(c_{\mathrm{CO_2}})")
	    end
	
	    rowgap!(fig.layout, 35)
	    display(fig)
	end
end

# ╔═╡ 48bfb18c-d8af-453e-b93a-e0fc5748d293
if L_varied_checkbox
	for L_val in L_values
	    if haskey(resL_COMP, L_val)
	        res_data = resL_COMP[L_val]
	        
	        num_steps = size(res_data.tsol)
	        
	        time_points = res_data.times[1:size(res_data.voltages)]
	        oh_surface  = [res_data.tsol[iohminus, 1, w] / (mol / dm^3) for w in 1:num_steps]
	        applied_V   = res_data.voltages[1:num_steps]
	        df = DataFrame(
	            Time_s = time_points,
	            Voltage_V = applied_V,
	            OH_Surface_M = oh_surface
	        )
	        
	        L_int = Int(L_val / μm)
	        filename = "concentration_L_$(L_int)um.csv"
	        
	        CSV.write(filename, df)
	        @info "Saved: $filename"
	    else
	        @warn "Key L = $(L_val) not found in resL_COMP"
	    end
	end
end

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
	model = model_key == "Gold_Model" ? elydata_unc :
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
is_Landstorfer = model != elydata_unc

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
	reaction_arg = model == elydata_Gold ? (reaction) : NamedTuple()
	sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = deepcopy(model), reaction_arg)

	result_pnp = capscalc(sys_pnp, is_Landstorfer; vrange = range(vmin, vmax, length = 201))
else
	result_pnp = nothing
end

# ╔═╡ e115560a-f79e-4a08-9bd2-5c11b0346827
if double_layer_curve
	fig_pnp, ax_pnp = AuCO2RR_plots.capsplot_fixed(result_pnp, "Poisson-Nernst-Planck", xlimits_L=(-1.0, 1.0), ylimits_L=(0, 200))
	fig_pnp
end

# ╔═╡ e6de6e89-5dc6-4fac-a7fc-8346175283a3
plot_cv_current(pnpresult_irc, model; species = ico)

# ╔═╡ d242507d-d1bb-461f-b4fe-ab3f597d9c40
plot_activity_time_electrode(pnpresult_odr, bulk, model)

# ╔═╡ 8cd25c0c-e260-4401-af12-a1def38bb7c2
if pH_varied_checkbox
	results_pH = run_pH_sweep(model, sawtooth, grid, pnp_bcondition, reaction; counter_index = 1, eneutral = true, pH_values = [2, 4, 7, 10])
end

# ╔═╡ 65df72b8-7954-4d7f-8ecd-ad168a303347
plot_pH_varied_sweep(results_pH; species = iohminus)

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

# ╔═╡ 7d1bd699-e398-484d-abab-37f86884f0ed
if L_varied_checkbox
    # 1. Selection of L values
    selected_Ls = L_values[1:end]

    # Meta-data
    species_names = getproperty.(bulk, :name)
    species_colors = getproperty.(bulk, :color)  
    nspecies = length(model.cspecies)

    num_panels = length(selected_Ls)
    fig = Figure(size = (800, 450 * num_panels)) 
    
    # Global Plotting Limits
    x_limits = (-1e-4, 6e-3) 
    y_limits_conc = (-12, 4) 
    y_limits_curr = (-1, 1) 
    
    axes_conc = []
    
    for (idx, L_val) in enumerate(selected_Ls)
        if !haskey(resL_COMP, L_val)
            @warn "Data for L = $(L_val/μm) μm not found. Skipping."
            continue
        end

        res = resL_COMP[L_val]
        tsol = res.tsol

        # --- FIX: GET ONLY THE INTEGER INDEX ---
		
        #target_V = -0.7
		#ti = findmin(abs.(res.voltages .- target_V))[2]
		target_time = 60
		ti = findmin(abs.(res.times .- target_time))[2]
		
		
        # Adding at the end is mandatory to get the Integer index
        
        L_label = round(L_val / μm, digits=0)
        actual_V = round(res.voltages[ti], digits=2)
		actual_T = round(res.times[ti], digits=2) 

		
        # --- GRID HANDLING ---
		
        # Ensure X_coords is extracted from the solution's system to match nx_current
		X_coords = X_coords_dict[L_val]
		nx_current = size(tsol, 2)
        
        # --- Axis 1: Concentration (Left Y-Axis) ---
        ax_conc = Axis(fig[idx, 1],
            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
            ylabel = L"\log_{10}(c_i)",
            #xscale = log10,
            limits = (x_limits, y_limits_conc),
			title  = "L = $L_label μm | t ≈ $actual_T s (ϕ_eff ≈ $actual_V V)",      
			xgridvisible = false,
			ygridvisible = false,
            spinewidth = 4,
            ylabelcolor = :blue, yticklabelcolor = :blue
        )
        push!(axes_conc, ax_conc)
        
        # --- Axis 2: Current Density (Right Y-Axis) ---
        ax_curr = Axis(fig[idx, 1],
            ylabel = L"j~(\mathrm{mA/cm^2})",
            #xscale = log10,
            yaxisposition = :right,
            limits = (x_limits, y_limits_curr),
            ylabelcolor = :skyblue, yticklabelcolor = :skyblue,
            ygridvisible = false,
            spinewidth = 4
        )
        hidespines!(ax_curr, :l, :t, :b)
        hidexdecorations!(ax_curr)

        # 2. Plotting Concentration Profiles
        for i in 1:nspecies
            # ti is now a single Integer, indexing will work
            c_profile = tsol[i, 1:nx_current, ti]
            y_log = map(c -> log10(max(c, 1e-25)), c_profile)
            
            lines!(ax_conc, X_coords[1:nx_current] .+ 1e-12, y_log; 
                   color = species_colors[i], linewidth = 3.5, label = species_names[i])
        end
        
        # 3. Current Reference
		phi_profile = tsol[model.iϕ, 1:nx_current, ti]
        lines!(ax_curr, X_coords[1:nx_current] .+ 1e-12, (phi_profile); 
               color = :skyblue, linewidth = 5, label = "Potential")
        
        hlines!(ax_curr, 0, color = (:black, 0.4), linestyle = :dash, linewidth = 2)
        
        #text!(ax_curr, 1e-6, j_val; text = "j = $(round(j_val, digits=2)) mA/cm²", 
		#	  align = (:left, :bottom), color = :red, fontsize = 20, font = :bold)

        if idx == 1
            axislegend(ax_conc; position = :lb, labelsize = 16, framevisible = false)
        end
    end
    
    rowgap!(fig.layout, 20)
    
    fig
end

# ╔═╡ 001209ff-074c-4bf4-bbd5-d8c2e88bfe8f
let
	if L_varied_checkbox
	    # 1. Selection of L values
	    target_indices = [1, 5, 6, 7, 8, length(L_values)]
	    selected_Ls = L_values[1:end]
	
	    # Meta-data
	    species_names = getproperty.(bulk, :name)
	    species_colors = getproperty.(bulk, :color)  
	    nspecies = length(model.cspecies)
	
	    num_panels = length(selected_Ls)
	    fig = Figure(size = (800, 450 * num_panels)) 
	    
	    # Global Plotting Limits
	    x_limits = (1e-12, 1e-2) 
	    y_limits_conc = (-12, 4) 
	    y_limits_curr = (-1, 1) 
	    
	    axes_conc = []
	    
	    for (idx, L_val) in enumerate(selected_Ls)
	        if !haskey(resL_COMP, L_val)
	            @warn "Data for L = $(L_val/μm) μm not found. Skipping."
	            continue
	        end
	
	        res = resL_COMP[L_val]
	        tsol = res.tsol
	
	        # --- FIX: GET ONLY THE INTEGER INDEX ---
			
	        #target_V = -0.9
			#ti = findmin(abs.(res.voltages .- target_V))[2]
			target_time = 25
			ti = findmin(abs.(res.times .- target_time))[2]
			
			
	        # Adding at the end is mandatory to get the Integer index
	        
	        L_label = round(L_val / μm, digits=0)
	        actual_V = round(res.voltages[ti], digits=2)
			actual_T = round(res.times[ti], digits=2) 
			
	        # Ensure X_coords is extracted from the solution's system to match nx_current
			X_coords = X_coords_dict[L_val]
	        nx_current = size(tsol, 2)
	        
	        # --- Axis 1: Concentration (Left Y-Axis) ---
	        ax_conc = Axis(fig[idx, 1],
	            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
	            ylabel = L"\log_{10}(c_i)",
	            xscale = log10,
	            limits = (x_limits, y_limits_conc),
				title  = "L = $L_label μm | t ≈ $actual_T s (ϕ_eff ≈ $actual_V V)",      
				xgridvisible = false,
				ygridvisible = false,
	            spinewidth = 4,
	            ylabelcolor = :blue, yticklabelcolor = :blue
	        )
	        push!(axes_conc, ax_conc)
	        
	        # --- Axis 2: Current Density (Right Y-Axis) ---
	        ax_curr = Axis(fig[idx, 1],
	            ylabel = L"j~(\mathrm{mA/cm^2})",
	            xscale = log10,
	            yaxisposition = :right,
	            limits = (x_limits, y_limits_curr),
	            ylabelcolor = :skyblue, yticklabelcolor = :skyblue,
	            ygridvisible = false,
	            spinewidth = 4
	        )
	        hidespines!(ax_curr, :l, :t, :b)
	        hidexdecorations!(ax_curr)
	
	        # 2. Plotting Concentration Profiles
	        for i in 1:nspecies
	            # ti is now a single Integer, indexing will work
	            c_profile = tsol[i, 1:nx_current, ti]
	            y_log = map(c -> log10(max(c, 1e-25)), c_profile)
	            
	            lines!(ax_conc, X_coords[1:nx_current] .+ 1e-12, y_log; 
	                   color = species_colors[i], linewidth = 3.5, label = species_names[i])
	        end
	        
	        # 3. Current Reference
			phi_profile = tsol[model.iϕ, 1:nx_current, ti]
	        lines!(ax_curr, X_coords[1:nx_current] .+ 1e-12, (phi_profile); 
	               color = :skyblue, linewidth = 6, label = "Potential")
	        
	        hlines!(ax_curr, 0, color = (:black, 0.4), linestyle = :dash, linewidth = 2)
	        
	        #text!(ax_curr, 1e-6, j_val; text = "j = $(round(j_val, digits=2)) mA/cm²", 
			#	  align = (:left, :bottom), color = :red, fontsize = 20, font = :bold)
	
	        if idx == 1
	            axislegend(ax_conc; position = :lb, labelsize = 16, framevisible = false)
	        end
	    end
	    
	    rowgap!(fig.layout, 20)
	    
	    fig
	end
end

# ╔═╡ d0bcf4f4-00d0-4408-8007-2d04933f7a0b
let
	if L_varied_checkbox
	    # 1. Selection of L values (Keep the same order)
	    target_indices = [1, 5, 6, 7, 8, length(L_values)]
	    selected_Ls = L_values[1:end]
	    num_panels = length(selected_Ls)
	
	    # 2. Meta-data
	    species_names = getproperty.(bulk, :name)
	    species_colors = getproperty.(bulk, :color)  
	    nspecies = length(model.cspecies)
	    
	    fig = Figure(size = (1000, 380 * num_panels))
	    
	    # --- CRITICAL: FIXED X-AXIS BASED ON THE LARGEST L ---
	    # Setting limits to the largest L in the set (the last one)
	    max_L = selected_Ls[end]
	    x_limits = (1e-12, max_L)
	    
	    y_limits_conc = (-12, 4)
	    y_limits_phi = (-1.5, 0.5)
	
	    panels_data = []
	
	    for (idx, L_val) in enumerate(selected_Ls)
	        res = resL_COMP[L_val]
	        nx_local = size(res.tsol, 2)
	        
	        # Grid Handling: Consistent with each L's specific mesh
			X_coords = X_coords_dict[L_val]
	
	        # --- Axis Setup ---
	        ax_conc = Axis(fig[idx, 1],
	            ylabel = L"\log_{10}(c_i)", xscale = log10,
	            limits = (x_limits, y_limits_conc), # Fixed X-axis for all
	            xgridvisible = false, ygridvisible = false, spinewidth = 3,
	            ylabelcolor = :blue, yticklabelcolor = :blue
	        )
	        ax_phi = Axis(fig[idx, 1],
	            ylabel = L"\Phi_{\mathrm{sol}}~(\mathrm{V})", xscale = log10, yaxisposition = :right,
	            limits = (x_limits, y_limits_phi),
	            ylabelcolor = :skyblue, yticklabelcolor = :skyblue,
	            ygridvisible = false, spinewidth = 3
	        )
	        hidespines!(ax_phi, :l, :t, :b)
	        hidexdecorations!(ax_phi)
	
	        # Observables
	        obs_concs = [Observable(fill(NaN, nx_local)) for _ in 1:nspecies]
	        obs_phi = Observable(fill(NaN, nx_local))
	        
	        # Plotting
	        for i in 1:nspecies
	            lines!(ax_conc, X_coords, obs_concs[i] ; 
	                   color = species_colors[i], linewidth = 3, label = species_names[i])
	        end
	        lines!(ax_phi, X_coords, obs_phi; color = :skyblue, linewidth = 5)
	        hlines!(ax_phi, 0, color = (:black, 0.4), linestyle = :dash)
	
	        if idx == 1
	            axislegend(ax_conc; position = :lb, framevisible = false, labelsize = 14)
	        end
	        vlines!(ax_conc, [L_val], color = :red, linestyle = :dash, linewidth = 2)
	
	        push!(panels_data, (concs=obs_concs, phi=obs_phi, axis=ax_conc, L=round(L_val/μm)))
	    end
	
	    # 3. GIF Recording (Slower framerate and synced steps)
	    # Finding the shortest common time steps
	    min_nt = minimum([size(resL_COMP[L].tsol, 3) for L in selected_Ls])
	    t_indices = 1:5:min_nt 
	
	    println("Recording GIF... Slower speed set via framerate=8.")
	
	    record(fig, "__L_varied_dynamic_profiles_EN.gif", t_indices; framerate = 8) do ti
	        for (idx, L_val) in enumerate(selected_Ls)
	            res = resL_COMP[L_val]
	            tsol = res.tsol
	            nx_local = size(tsol, 2)
	            
	            # Update Concentrations
	            for i in 1:nspecies
	                c_raw = tsol[i, 1:nx_local, ti] / (mol / dm^3) 
	                panels_data[idx].concs[i][] = map(c -> log10(max(c, 1e-25)), c_raw)
	            end
	            
	            # Update Potential
	            p_raw = tsol[model.iϕ, 1:nx_local, ti]
	            panels_data[idx].phi[] = p_raw
	
	
				
	            # Update Titles (English)
	            V_now = round(res.voltages[ti-1], digits=2)
	            T_now = round(res.times[ti], digits=2)
	            panels_data[idx].axis.title = "Thickness (L) = $(panels_data[idx].L) μm | Time = $T_now s (V_eff = $V_now V)"
	        end
	    end
	end
end

# ╔═╡ 8932c3d0-e103-4173-8d85-aa281ef58f53
let
	if L_varied_checkbox
	    # 1. Identify CO2 Index
	    # Replace "CO2" with the exact string name used in your 'bulk' species setup
	    co2_idx = findfirst(s -> s.name == "CO₂", bulk)
	    if isnothing(co2_idx)
	        error("CO2 species not found in bulk model. Please check the name.")
	    end
	    co2_color = bulk[co2_idx].color
	    co2_name = bulk[co2_idx].name
	
	    # 2. Selection of L values
	    target_indices = [1, 5, 6, 7, 8, length(L_values)]
	    selected_Ls = L_values[1:end]
	    num_panels = length(selected_Ls)
	
	    fig = Figure(size = (1000, 380 * num_panels))
	    
	    # Fixed X-axis based on the largest L for comparison
	    max_L = selected_Ls[end]
	    x_limits = (1e-12, max_L)
	    
	    # Specific Y-limits for CO2 (Adjust based on your bulk concentration)
	    y_limits_conc = (-8, 2) 
	    y_limits_phi = (-1.5, 0.5)
	
	    panels_data = []
	
	    for (idx, L_val) in enumerate(selected_Ls)
	        res = resL_COMP[L_val]
	        nx_local = size(res.tsol, 2)
	
			# Grid Handling: Consistent with each L's specific mesh
	        # --- GRID HANDLING ---
			X_coords = X_coords_dict[L_val]
	
	        # --- Axis Setup ---
	        ax_conc = Axis(fig[idx, 1],
	            ylabel = L"\log_{10}(c_{\mathrm{CO_2}})", xscale = log10,
	            limits = (x_limits, y_limits_conc),
	            xgridvisible = false, ygridvisible = false, spinewidth = 3,
	            ylabelcolor = co2_color, yticklabelcolor = co2_color
	        )
	        ax_phi = Axis(fig[idx, 1],
	            ylabel = L"\Phi_{\mathrm{sol}}~(\mathrm{V})", xscale = log10,
				yaxisposition = :right,
	            limits = (x_limits, y_limits_phi),
	            ylabelcolor = :skyblue, yticklabelcolor = :skyblue,
	            ygridvisible = false, spinewidth = 3
	        )
	        hidespines!(ax_phi, :l, :t, :b)
	        hidexdecorations!(ax_phi)
	
	        # Observables for CO2 and Potential
	        obs_co2 = Observable(fill(NaN, nx_local))
	        obs_phi = Observable(fill(NaN, nx_local))
	        
	        # Plotting
	        # Dash line for initial bulk concentration reference
	        hlines!(ax_conc, [log10(max(res.tsol[co2_idx, end, 1] / (mol / dm^3), 1e-25))], 
	                color = co2_color, linestyle = :dot, alpha = 0.5)
	        
	        lines!(ax_conc, X_coords, obs_co2; color = co2_color, linewidth = 4, label = co2_name)
	        lines!(ax_phi, X_coords, obs_phi; color = :skyblue, linewidth = 5, label = "Potential")
	        
	        hlines!(ax_phi, 0, color = (:black, 0.4), linestyle = :dash)
	
	        if idx == 1
	            axislegend(ax_conc; position = :lb, framevisible = false)
	        end
	
	        push!(panels_data, (co2=obs_co2, phi=obs_phi, axis=ax_conc, L=round(L_val/μm)))
	        vlines!(ax_conc, [L_val], color = :red, linestyle = :dash, linewidth = 2)
	
	    end
	
	    # 3. Slower GIF Recording
	    min_nt = minimum([size(resL_COMP[L].tsol, 3) for L in selected_Ls])
	    t_indices = 1:5:min_nt 
	
	    record(fig, "CO2_Evolution_L_Varied.gif", t_indices; framerate = 8) do ti
	        for (idx, L_val) in enumerate(selected_Ls)
	            res = resL_COMP[L_val]
	            tsol = res.tsol
	            nx_local = size(tsol, 2)
	            
	            # Update CO2 only
	            c_raw = tsol[co2_idx, 1:nx_local, ti] / (mol / dm^3)
	            panels_data[idx].co2[] = map(c -> log10(max(c, 1e-25)), c_raw)
	            
	            # Update Potential
	            p_raw = tsol[model.iϕ, 1:nx_local, ti]
	            panels_data[idx].phi[] = p_raw
	            
	            # English Title
	            V_now = round(res.voltages[ti], digits=2)
	            T_now = round(res.times[ti], digits=2)
	            panels_data[idx].axis.title = "Thickness L = $(panels_data[idx].L) μm | Time = $T_now s (V_eff = $V_now V)"
	        end
	    end
	end
end

# ╔═╡ 11b12556-5b61-42c2-a911-4ea98a0a1e85
if IV 
	cell, ivresult = simulate_CO2R(grid, model)
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
	AuCO2RR_plots.iv_curve_axis(ivresult; cutoff=-0.4, showlegend=true, species = iohminus)
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

# ╔═╡ 9739f8a5-f97e-4013-997f-d65cbf357ae5
if CV_comp
	compresult = sweep_COMP(model; eneutral = false, tunnel = false)
	export_cv_profile_csv(compresult; function_name = "cv_profile_drop_")
end

# ╔═╡ fd2fe768-020c-4751-bde0-8f75843580f7
if CV
	export_cv_profile_csv(pnpresult)
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
    scanrate_get = saw -> saw.scanrate,   # <- saw에서 scan rate 꺼내는 방식
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
        Y = currents(rec, species) .* scale   # <- current 추출

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

# ╔═╡ 5cb3fb17-6ce9-4230-9a3f-4df9f49293da
function export_L_varied_COMPENSATED_conc_csv(
    res_dict;                 
    species_list = [iohminus],    
    outdir::AbstractString = "../data/output",
    function_name::AbstractString = "L_varied_compensated_total",
    current_scale = cm^2 / mA,    
    L_scale = 1e6                 # Unit conversion for L (m -> μm)
)
    Ls       = Float64[]
    voltages = Float64[] 
    values   = Float64[]

    for L in sort(collect(keys(res_dict)))
        res = res_dict[L]
        
        I_sim = sum(currents(res, sp) for sp in species_list)
        
        
        Y_scaled = I_sim .* current_scale 
        
        @assert length(V_eff) == length(I_sim)

        append!(Ls, fill(Float64(L) * L_scale, length(V_eff)))  
        append!(voltages, V_eff) 
        append!(values, Y_scaled)
    end

    df = DataFrame(
        BoundaryLayerThickness = Ls,
        Voltage = voltages, 
        Value   = values,
    )
    # Generate the file path and save the DataFrame as a CSV file
    fname = filename(function_name)
    outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    @info "Data successfully exported! (Total current of all specified species, internal iR compensation applied)"
    return df
end

# ╔═╡ 2266c928-276f-4cec-8fee-eca23cc0e1fd
function export_L_varied_COMPENSATED_csv(
    res_dict;                     
    species_list,                 # Array of species indices, e.g., [iohminus, ihplus]
    outdir::AbstractString = "../data/output",
    function_name::AbstractString = "L_varied_compensated_total",
    current_scale = cm^2 / mA,    
    L_scale = 1e6                 
)
    # Initialize a DataFrame to hold everything
    final_df = DataFrame()

    # Sort L values to ensure the CSV is ordered
    sorted_keys = sort(collect(keys(res_dict)))

    for L in sorted_keys
        res = res_dict[L]
        num_steps = length(res.voltages)
        # 1. Basic Columns
        current_L_df = DataFrame(
            L_um = fill(Float64(L) * L_scale, num_steps),
            Voltage_V = res.voltages,
            Current_mAcm2 = sum(currents(res, sp) for sp in species_list) .* current_scale,
			Times_s = res.times
        )
		species = getproperty.(bulk, :name)        
        # 2. Dynamic Concentration Columns
        # This creates a column for each species provided in species_list
        for sp in species_list
            # We use the species index/name as part of the column header
            col_name = "$(species[sp])_Surface_M" 
            current_L_df[!, col_name] = [res.tsol[sp, 1, w] / (mol / dm^3) for w in 1:num_steps]
        end

        # Append this L's data to the master DataFrame
        append!(final_df, current_L_df)
    end

    # File Export Logic
	#function_name = 
    outfile = joinpath(outdir, filename(function_name) * ".csv")
    mkpath(dirname(outfile))
    
    CSV.write(outfile, final_df)

    @info "Combined Data Exported to: $outfile"
    return final_df
end

# ╔═╡ 6035d86f-c7f0-4fd8-b79b-83e405665f59
if L_varied_checkbox
	df_L_comp = export_L_varied_COMPENSATED_csv(
	    resL_COMP; function_name = "L_varied_compensated_timestep_", species_list = model.cspecies
	)
end

# ╔═╡ e3cb8c7c-ec7e-4848-8e41-21bcc106d4e1
if L_varied_checkbox
	df_L_uncomp = export_L_varied_COMPENSATED_csv(
	    resL; function_name = "L_varied_uncompensated_", species_list = model.cspecies
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
# ╠═5a27d95a-21d2-4e8a-a2fd-895eed32b113
# ╠═bd8134d8-5a69-486e-8429-7cf810b3ccbe
# ╠═2176bc34-fc74-4532-912e-e73441b37245
# ╠═f7d13047-4007-47ac-a3bb-a0b788dcd141
# ╠═beae1479-1c0f-4a55-86e1-ad2b50174c83
# ╠═9b7350ca-e176-4ebc-8a67-00d6737ab1a2
# ╟─ab2184fc-0279-46d9-9ee4-88fe3e732789
# ╠═7316901c-d85d-48e9-87dc-3614ab3d81a5
# ╟─6b7cfe87-8190-40a5-8d25-e39ef8d55db5
# ╠═5a146a44-03dc-45f3-ae15-993d11c2edac
# ╠═00947475-c96e-4ecc-a1ef-5be5e3e3c864
# ╠═165446a0-da6c-4ecc-a5ba-c96eb300af12
# ╟─a4b1300f-7e8a-46ab-9efd-dba82315a966
# ╠═2d3e4cc2-7823-4fd3-b512-4d042f15cc68
# ╠═05ebcce5-6904-451f-bbd6-ea4588ca05e3
# ╠═2c061d99-da09-4ccf-a4ec-73a3c9a908fe
# ╠═d5c6769d-35f9-461a-aae2-00122ccddd63
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
# ╠═a4928999-96c7-4145-af2d-485ff13d7109
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
# ╠═8f63ec58-bc97-43ba-9bbe-e10bb50e2cfe
# ╠═e23eaa0a-a4ec-4c05-8526-15df6ad0eec1
# ╠═2947efab-e67e-41fa-ab91-7e3c1c09df2d
# ╠═de3c2f2e-77fe-4d4f-9bbc-ae2cf8c0406f
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
# ╠═ed92cece-3f89-45f5-ac17-cbc9a9abb906
# ╟─5c808c71-6094-49d7-8215-e88262f34e1f
# ╟─da8390d1-47e8-451f-b12b-45b8aca7b6ec
# ╠═02d12ba4-4ab3-48f6-b084-edb06cb413b1
# ╠═b4aaf070-d4ab-409a-b1e8-f5469b9f398b
# ╠═87855181-ffbb-4a15-9852-0ded03746f18
# ╠═e50fe651-11d4-45ee-89dd-371a7fbc097e
# ╟─ef7212fc-a3d0-4784-b901-219204b79dc0
# ╠═a72db110-37c3-4f4d-b37e-fe613c1ae827
# ╠═9739f8a5-f97e-4013-997f-d65cbf357ae5
# ╠═b95160b5-18f7-49d9-80be-9159abd2dcd1
# ╟─ad55a4c1-78b9-40d2-aec8-64ed95174e8a
# ╠═f6f26f7e-b97e-4f83-8b02-d3ff8b8ad14d
# ╠═9fb47b83-a853-4316-bb8d-30e65b16ef78
# ╠═017d8b3a-5923-411d-bfbd-31e3167e5ac0
# ╠═9e01784a-c575-41e4-b82a-6c65d78004c1
# ╠═28337c8b-8687-4618-bd9d-16d94193dd5f
# ╠═e8e6cff2-574f-4863-8093-9a747e419cf7
# ╠═7da046bf-d3b1-43a0-bdba-89b4da2f6be3
# ╠═5b3ef1cb-c063-42b2-a0af-1e1a36d1e180
# ╠═bf3a3756-59a9-4767-a713-034ef37ab55f
# ╠═fdbaa3a0-6b2b-45ad-8618-a63a7f720049
# ╠═a64e2dc9-9be7-48b5-9d04-c448f19ed7f2
# ╠═cdf52b70-94ad-45db-b82d-f1268cead86e
# ╠═311660ce-0168-47eb-8d23-00a8e750fd84
# ╟─2d421be5-9287-4e8d-9a7e-d63c95d9f370
# ╠═c62b2553-e42f-4002-b206-9af0d8655c3d
# ╠═089843dc-e9b0-4d43-afb7-d3477af39587
# ╠═cdffb77b-89de-4aba-87b6-21e29bac0ae0
# ╟─b65575ce-7c5e-4a70-bb50-3bd8054516aa
# ╠═1efbe0ce-699f-4b6a-ab82-7b12d88e150b
# ╠═6d3dc6b1-82ea-448a-af86-ed6001ed5446
# ╠═adfc0c3b-10d0-4f99-8b42-6da9a26de219
# ╠═239abbd9-928d-4f33-a9a1-517c986d4202
# ╠═6027ad72-7d72-47c2-9025-f00af8380c7f
# ╠═50dc3a93-6eeb-4f6f-ac42-65d8d597614c
# ╠═ea02ccff-20ea-426d-b102-e93a30e8fdb4
# ╠═ce8dd32b-2ea3-4a42-9436-69b06e751984
# ╠═e6de6e89-5dc6-4fac-a7fc-8346175283a3
# ╠═ac192132-1778-471c-9e8d-958f3bfe400a
# ╠═606cc901-0c24-4588-be13-d05f55b6907f
# ╠═fd5c53ac-1a27-40f9-8ae4-cfbf1aaba5d1
# ╠═d242507d-d1bb-461f-b4fe-ab3f597d9c40
# ╠═f668c334-ba69-4ab9-b719-b4078305e68b
# ╠═d93cd922-d021-40b1-8a4d-c9f2bfa4869a
# ╠═3a07864e-9579-4fd4-a62a-d68472b9f9c6
# ╟─e6f43f01-15d8-4265-ae15-3673fb3cf7e3
# ╠═2754c3f8-c22b-4389-8aab-a6ab93a9ca9c
# ╠═9d0b6083-e49b-4b15-ad2e-a2c74f689c5e
# ╠═2420382d-227a-4063-9450-1f1726df018e
# ╠═3d661549-a8d2-40b0-add8-b186193f90fe
# ╠═ec296961-151f-4219-a7f2-fb34e7facaa5
# ╠═74c43d72-3a23-4a24-a4ae-8b18b245610a
# ╠═fd2fe768-020c-4751-bde0-8f75843580f7
# ╠═9e44f14b-a799-4ca5-8641-a3726780a4fe
# ╠═f1035efd-3e58-4c62-8799-a8750bdf137e
# ╟─25eb8aa3-697e-4538-9472-ceea45fbfbd9
# ╟─11892724-1851-46f2-802d-4da45127b0af
# ╠═8cd25c0c-e260-4401-af12-a1def38bb7c2
# ╠═65df72b8-7954-4d7f-8ecd-ad168a303347
# ╠═b2201bee-c72e-4498-bde0-7772dd2928c2
# ╟─c048e472-3983-4279-bf60-82784baa145e
# ╟─3bdaab98-c0f7-46af-86b7-d68374e8a5d0
# ╠═a05cf724-cd32-498e-8afb-ecbf4a1f1648
# ╠═5ccb0682-09de-4ae3-95e6-a8403d540d9b
# ╠═406fb8e5-61e6-4688-96bf-a5530f57d4fa
# ╟─e0e59ef0-8b6c-4f31-8d39-c2c4bcd7f99e
# ╠═56814250-16b2-4578-820d-2096998c84f4
# ╠═7b38e59a-d005-4cfc-ba8c-b17e7c700119
# ╠═b767a48f-1b20-4b85-b9a8-36d6a914c5fe
# ╠═bb01b182-7840-4ad0-8cd6-fae57ab93173
# ╠═0a1054c5-cee7-4202-9d32-9eee6ec55265
# ╠═43e4ef6a-5c2a-4914-9972-ea93f3cb016e
# ╠═7e102647-23a9-4f40-b6de-cb1938bbb23e
# ╠═a34cdba3-38f6-4bc0-b26e-72c956599109
# ╠═e6dca43d-69b5-4d35-903c-6742b16a4715
# ╠═35a462d3-1f8a-4a8e-ab8e-a07e844469e0
# ╠═cbcd2ff5-3df2-40df-b0c5-e127f3abb964
# ╟─9b8daa64-9e34-4da1-8136-ba5488e037c7
# ╠═360313f0-2dad-4f7f-8d11-1800c6d934b3
# ╠═4231590b-18c4-4af8-b798-364c529e47d8
# ╠═c9ae7a67-9a5f-4d9a-88c0-742f4e91fb27
# ╟─58ac8edc-2432-4054-88d8-52dafe0a2a61
# ╟─a2c7c4da-77cd-493f-8f98-0c86fecf271a
# ╟─0db921e1-80ab-4d79-9035-df2dc33c0c3c
# ╟─d7e28d5f-16ba-4664-ba63-ec85fb29fe87
# ╟─7d1bd699-e398-484d-abab-37f86884f0ed
# ╟─001209ff-074c-4bf4-bbd5-d8c2e88bfe8f
# ╟─d0bcf4f4-00d0-4408-8007-2d04933f7a0b
# ╟─ee83fb10-b4a1-4084-a644-0541500fc9cb
# ╟─65e8592d-a3f8-4f38-9c7e-8085a5c28b73
# ╟─8932c3d0-e103-4173-8d85-aa281ef58f53
# ╠═d91c32c8-ac55-4f77-93cb-6c5d97bcbbb2
# ╠═48bfb18c-d8af-453e-b93a-e0fc5748d293
# ╠═5cb3fb17-6ce9-4230-9a3f-4df9f49293da
# ╠═6e88e1d8-1f4b-4813-8890-0cfcdc5fb967
# ╠═6035d86f-c7f0-4fd8-b79b-83e405665f59
# ╠═e3cb8c7c-ec7e-4848-8e41-21bcc106d4e1
# ╠═68e2cd9e-e8c0-491a-91c5-b30b16ed3f20
# ╠═2266c928-276f-4cec-8fee-eca23cc0e1fd
# ╟─fbe4aca2-6a47-4457-98bb-588a5cde0ed5
# ╠═b25f2246-0182-4d50-a606-0d81776d414f
# ╠═bf991ac4-8ff2-4a49-9f1e-16ff60b12493
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
# ╠═012b426e-3550-4678-a666-eeb4bd32d20e
# ╠═06ca7d68-9c76-4d1c-947a-dd64a0fe9ec3
# ╟─57db41d1-57c0-4eee-ba35-6f9b7e1e8263
# ╠═11b12556-5b61-42c2-a911-4ea98a0a1e85
# ╟─b976ab43-69f1-47a0-b2c6-c63e1c15cdb4
# ╠═60b410be-70f7-4053-a3db-7d777e0d3f08
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
# ╠═22244e24-5b56-4933-8a09-44b601f116c3
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
