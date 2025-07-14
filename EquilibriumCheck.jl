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

# ╔═╡ 2b901eca-db3b-4ad2-b0ee-e031854c57fa
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using Pkg
	Pkg.activate(@__DIR__) 
	using Revise
	using CairoMakie
	using GridVisualize
	using Colors
	if isdefined(Main,:PlutoRunner)
   		default_plotter!(CairoMakie)
 		CairoMakie.activate!(type="svg")
    end
end
  ╠═╡ =#

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using PlutoUI, Latexify, HypertextLiteral
    using VoronoiFVM
    using ExtendableGrids
    using LinearAlgebra
    using NLsolve
    using Unitful
    using LessUnitful
	using DelimitedFiles
	using PreallocationTools
    using LinearSolve, ExtendableSparse
	using CatmapInterface 
    using LessUnitful.MoreUnitful
    using LiquidElectrolytes
    using Test
	using Printf
	using Catalyst
end

# ╔═╡ b95ff168-68dc-4172-be97-df5362be6c48
begin
	using CSV, DataFrames
end

# ╔═╡ ecde8a95-f660-4b53-9296-cd4d4022c95f
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
end

# ╔═╡ ef660f6f-9de3-4896-a65e-13c60df5de1e
md"""
# Double layer capcacitance comparison
"""

# ╔═╡ 8bf52bab-4830-4936-b3c9-78c1ccd0406d
#=╠═╡
Pkg.status("LiquidElectrolytes")
  ╠═╡ =#

# ╔═╡ 4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="",depth=5)
  ╠═╡ =#

# ╔═╡ 852d9c74-b2aa-49c0-9bd0-0ccb1afcc520
pkgdir(LiquidElectrolytes)

# ╔═╡ ac27c318-9a00-4287-bd36-3a97d65b5459
md"""
## 8. Tests
"""

# ╔═╡ fe48d05b-99bd-48b4-a044-4dd8e8d18b5d
begin
    SI(x) = Float64(Unitful.ustrip(Unitful.upreferred(1 * x)))
    const V = SI(Unitful.V)
    const eV = SI(Unitful.eV)
    const nm = SI(Unitful.nm)
    const cm = SI(Unitful.cm)
    const μF = SI(Unitful.μF)
end


# ╔═╡ b3b55993-cb87-43c6-a084-22fa33f02f8c
begin
	@unitfactors mol dm m s K μm bar Pa μA Å;
	@phconstants N_A c_0 k_B e h
end

# ╔═╡ 5fe96d0b-7bd0-4183-901d-727e966d434b
begin
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
	
	const scheme 	= :μex

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

# ╔═╡ a23eece5-8e94-4c0b-b487-e742a37e714e
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
		c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol/dm^3 * 0.01
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

# ╔═╡ 59855587-c6c2-4af6-a713-0b710cf2b0fe
begin
	const at = 0.0
	const κt = 0.0
	const ak = 8.2
	const κk = 0.0
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
							κ = κt, 
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
							κ = κt,  
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

# ╔═╡ 3be02c97-5c28-4370-97c1-e3f9faaba62a
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

# ╔═╡ 595715e5-f108-4167-b104-ac7c6f652e48
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
							   	#actcoeff! = DGL_gamma!
							   )

# ╔═╡ 53ae411f-42e8-41b3-ad38-e0189a001acf
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

# ╔═╡ 6b69df37-8754-457c-93cc-e9bacfa47d9a
elydata_NaClO₄ = ElectrolyteData(
 		z = [-1, 1],
		κ = [8.0, 8.0],
		#v = [25.0, 25.0],
		c_bulk = [0.5, 0.5],
		ε = 26.0,
		#vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
)

# ╔═╡ 045573a4-aff5-4262-aca2-ead598a37bf2
elydata_Default = ElectrolyteData(
	 	z = [-1, 1],
)

# ╔═╡ 25bafe0e-f2fc-4a5c-828e-89fdccbc250c
md"""
### Reaction
"""

# ╔═╡ 6d1f0f93-876a-4ec9-9763-f1abcb36cc07
md"""
#### Boundary Reaction
"""

# ╔═╡ 0bf42bde-0c30-4df2-8fa0-82e493add078
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

# ╔═╡ c11fdb45-b46e-40b6-b5a7-9c115aa5fee5
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

# ╔═╡ e9f29b23-0f46-4515-9ba3-06cd25c1d741
rn

# ╔═╡ 949d126e-d863-40f4-b902-8b3b4c97b1b5
md"""
#### Buffer Reaction
"""

# ╔═╡ 5ddce46c-22a5-427a-8cb7-f45f216cefe0
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

# ╔═╡ e2ab9bb4-a1b5-4d49-a760-013ad0fc4c68
md"""
#### Reaction Rates
"""

# ╔═╡ 3a9940b9-c7e1-484c-bbf2-14c0c69d685b
md"""
##### Electrode Reaction
"""

# ╔═╡ 12cbfb8b-edb6-4335-8d80-0d6fe0eb9d3a
begin
	const ps_cache = DiffCache(zeros(8), 13)
	const us_cache = DiffCache(zeros(isurfaceend-isurfacestart+1), 13)
	
	function we_breactions(f, 
			u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
			bnode, 
			data
		) where {Tval, Tv, Tc, Tp, Ti}
		(; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we) = data
		
		γ_co2 	= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))
		γ_co 	= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))
		σ 			= C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)
		local_pH 	= -log10((u[ihplus] / (mol/dm^3)))

		ps = get_tmp(ps_cache, u[iϕ])
		for (p, default_value) in odesys.defaults
			ps[paramsidx[p]] = default_value
			println(default_value)
		end
		
		ps = get_tmp(ps_cache, u[iϕ])
		ps[paramsidx[Symbolics.rename(odesys.σ, :σ)]] = σ
		ps[paramsidx[Symbolics.rename(odesys.γCO2_aq, :γCO2_aq)]] = γ_co2 
		ps[paramsidx[Symbolics.rename(odesys.aH2O_g, :aH2O_g)]] = aH₂O 
		ps[paramsidx[Symbolics.rename(odesys.ϕ, :ϕ)]] = u[iϕ] 
		ps[paramsidx[Symbolics.rename(odesys.ϕ_we, :ϕ_we)]] = ϕ_we 
		ps[paramsidx[Symbolics.rename(odesys.local_pH, :local_pH)]] = local_pH 
		ps[paramsidx[Symbolics.rename(odesys.γCO_aq, :γCO_aq)]] = γ_co 
		ps[paramsidx[Symbolics.rename(odesys.βCOOHΔH2OΔele_t, :βCOOHΔH2OΔele_t)]] = 0.59 

		if bnode.region == Γ_we && size(f,1) ≥ isurfaceend
			@views f_microkinetics!(
				f[isurfacestart:isurfaceend], 
				u[isurfacestart:isurfaceend],
				ps,
				nothing
			)
		elseif bnode.region == Γ_we
			nothing
		end
		
		# conversion from turnover frequency (appropriate for change in coverage) to production rate (per unit area) (approprite for change in concentration) by S = number of free catalyst sites in mole per unit area
		f[ico2] *= S
		f[iohminus] *= S
		f[ico] *= S
		f[ikplus] *= S
	end
end

# ╔═╡ 1c94b9ba-429d-44fc-887e-e028ba070cc9
md"""
##### Bulk Reaction
"""

# ╔═╡ 1e52766d-12a9-46cd-be96-5c13a046944f
begin
	const γ_cache = DiffCache(zeros(nc), 12)
	
	function reaction(
		f, 
		u::VoronoiFVM.NodeUnknowns{Tv, Tc, Tp, Ti}, 
		node, 
		data
	) where {Tv, Tc, Tp, Ti}  
		
		(; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, pscale, p_bulk) = data

		# compute activity coefficients according to the approach in Ringe et al.
		γ = get_tmp(γ_cache, u[ico2])

		γ .= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))

		@views f_buffer!(
			f[ibufferstart:ibufferend], 
			u[ibufferstart:ibufferend],
			γ[ibufferstart:ibufferend],
			nothing
		)
		nothing
	end
end;

# ╔═╡ 442fe098-497b-404f-80a0-880bc95d5e02
inival = VoronoiFVM.unknowns(sys_sy, inival = 0);

# ╔═╡ c381803b-daad-4778-8d79-5abcecbce9ee
molarities = [0.005, 0.01, 0.02, 0.04, 0.1, 0.5, 1] 

# ╔═╡ 6df5efba-c6bc-4bc9-ae86-514e8171c8f1
md"""
#### γ(activity coefficients) function
"""

# ╔═╡ 237a3b9e-60bc-48c8-b5d1-ea954e65c401
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

# ╔═╡ 72f17b01-34cc-4196-99ae-a64af3f864c4
function DGML_gamma!(γ, c, p, electrolyte)
	
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
    for ic in cspecies
        γ[ic] = rexp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc))
    end
    return nothing
end

# ╔═╡ 38061646-9c66-4f9c-a0b5-5090dc62f8fe
md"""
#### Algebraic pressure equation
"""

# ╔═╡ e114ec0d-13d3-4455-b1c9-d1c5d76671d9
md"""
#### Pressure poisson problem
"""

# ╔═╡ 9b1dc273-9938-43a0-ac10-1928a80f89d8
md"""
#### Poisson Nernst-Planck from LiquidElectrolytes
"""

# ╔═╡ 98464285-2bd4-4631-8c4f-8790fe15cb93
md"""
#### Poisson-Boltzmann from LiquidElectrolytes
"""

# ╔═╡ a8e26e1a-a9ac-4d51-b09c-7951acd4b4b7
function pb_bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data
	
    ## Dirichlet ϕ=ϕ_we at Γ_we
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_bulk, value = data.ϕ_bulk)
    boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)

	## Robin ϕ=dϕ₀/dx
	#boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))


    return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ e5dbe1db-3b69-49ed-9f32-e65579069c47
md"""
### Sweep Functions
"""

# ╔═╡ 2fb75c2a-c877-4c78-abf8-6f89706a58fe
md"""
## 9. Tests Results
"""

# ╔═╡ 289d2c59-e920-47fe-b9ad-cb0a33ef0c9c
md"""
### Model Select
"""

# ╔═╡ c75a852d-e3b8-46e5-bcdd-5c41aef36c64
@bind model_choice Select(["Gold_Model", "Landstorfer_NaClO₄ model", "Landstorfer_NaF model"])

# ╔═╡ 791ccb34-e761-4e65-a9ef-95eac5395376
model = model_choice == "Landstorfer_NaClO₄ model" ? elydata_NaClO₄ : 
        model_choice == "Gold_Model" ? elydata_Gold : 
        elydata_NaF

# ╔═╡ 267a7233-e34a-4c4a-8627-a49bb45ccd25
is_Landstorfer = model != elydata_Gold

# ╔═╡ 70e1a34b-9041-4151-91aa-4dd7907a5b13
function capscalc(sys, molarities)
    result = []
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
		        volts = voltages(r)
		        caps = r.dlcaps
		    end
		    cdl0 = dlcap0(data)
		    @info "elapsed=$(t)"
		    push!(result, (voltages = volts, dlcaps = caps, cdl0 = cdl0, molarity = 	molarities[imol]))
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
	        #volts = voltages(r)
	        caps = r.dlcaps
	    end
	    cdl0 = dlcap0(data)
	    @info "elapsed=$(t)"
	    push!(result, (voltages = volts, dlcaps = caps, cdl0 = cdl0))
	end
    return result
end


# ╔═╡ 398b3511-4f7c-4436-9fe8-8edd76e3e0e7
result_sy = capscalc(sys_sy, molarities)

# ╔═╡ ca3bd6ba-1b3d-42c7-b008-8012b06368e4
result_pp = capscalc(sys_sy, molarities)

# ╔═╡ 53cdf6d7-a025-49e0-af7b-cc0838cfb422
function pnp_bcondition(f, u, bnode, data::ElectrolyteData)
	(; Γ_we, Γ_bulk, ϕ_we, iϕ) = data
	
    ## Dirichlet ϕ=ϕ_we at Γ_we
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
	
	## Robin ϕ=dϕ₀/dx
	#boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))
	
	if model == elydata_Gold && bnode.region == Γ_we
		we_breactions(f, u, bnode, data)
	end
		
	return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ 1fcfbad3-2fad-4eee-a1a3-031dc29c9083
function resultcompare(r1, r2; tol = 1.0e-3)
    for i in 1:length(r1)
        for f in fieldnames(typeof(r1[i]))
            if !isapprox(r1[i][f], r2[i][f]; rtol = tol)
                return false
            end
        end
    end
    return true
end

# ╔═╡ 0b6f33b9-41d4-48fd-8026-8a3bddcc1989
md"""
### Result plot

Compare with Fig 4.2 of [Fuhrmann (2015)](https://dx.doi.org/10.1016/j.cpc.2015.06.004)
"""

# ╔═╡ a7677fc3-a83d-4fab-8d00-b5c2454056c9
md"""
#### datafiles
"""

# ╔═╡ d18fe756-b0b9-44d7-8872-6b7812108c16
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

# ╔═╡ 4fe01b26-ad9f-44b8-8900-916a5aeb5ad6
md"""
#### Plotting Functions
"""

# ╔═╡ a22a5421-05bf-484f-a2d3-91a06a0c6476
# ╠═╡ skip_as_script = true
#=╠═╡
function capsplot(vis, result, title)
	if is_Landstorfer
    	hmol = 1 / length(result)
    	for imol in 1:length(result)
        	c = RGB(imol * hmol, 0, 1-imol * hmol)
        	scalarplot!(
            vis, result[imol].voltages, result[imol].dlcaps / (μF / cm^2),
            color = c, clear = false, label = "$(result[imol].molarity)M", 						markershape = :none, title = title,  xlabel = "φ / (V vs φ_pzc)", ylabel 			= "dlcaps / (μF / cm²)"
        )
        scalarplot!(
            vis, [0], [result[imol].cdl0] / (μF / cm^2),
            clear = false, markershape = :circle, markersize = 8, label = ""
        )
    	end
	else
		scalarplot!(
            vis, result[1].voltages, result[1].dlcaps / (μF / cm^2), limits=					(-1, 100), xlimits=(-1.1, 1.1), color = :green, clear = false, label = 				"$title", title = title, markershape = :none, yscale=10, xlabel = "φ / (V vs φ_pzc)", 				ylabel = "dlcaps / (μF / cm²)")
		scalarplot!(
            vis, [0], [result[1].cdl0] / (μF / cm^2),
            clear = false, markershape = :circle, markersize = 8, label = ""
        )
	end
    return vis
end;
  ╠═╡ =#

# ╔═╡ a4a01dcb-8c02-441b-ba25-8e8c062d7d58
md"""
### Compare EDL Plots
"""

# ╔═╡ 7cb51a9b-6357-4ec6-af6a-9113dea61661
md"""
#### Plotting Functions
"""

# ╔═╡ 87f2b4c4-b163-4ae2-86b6-0266dff1da19
function capsplot_v(vis, result_named)
    color = [:magenta, :blue]
    for (i, (name, res)) in enumerate(result_named)
        scalarplot!(
            vis, res[1].voltages, res[1].dlcaps / (μF / cm^2);
            limits = (-1, 100), xlimits = (-1.1, 1.1),
            color = color[i], clear = false, label = name,
            markershape = :none, yscale = 10,
            xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)"
        )
    end
    return vis
end;

# ╔═╡ 4c1f6b31-ce09-4fba-b827-460e8a0d7e1a
md"""
### κ(Solvation Number) Plots
"""

# ╔═╡ 0e734e72-fc3a-48c7-b1d5-c0a380768eec
function caps(sys)
	dls = dlcapsweep(sys, voltages = range(-1, 1, length = 401))	
	return dls 
end

# ╔═╡ fae68c38-be85-4718-8ee6-f900150e2b9a
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


# ╔═╡ 43a3d5f4-16ce-4402-bb95-d759b07bd573
md"""
### CV Result
"""

# ╔═╡ 1f7971ad-80cc-4bc1-a2f7-a912186537f7
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


# ╔═╡ 05703d80-3299-4692-9d23-f44c2f371f7b
md"""
#### Plotting Functions
"""

# ╔═╡ f88ca8d4-ecbb-4eee-b4af-5854fbd16e33
md"""
### Concentration Plots
"""

# ╔═╡ 7ab1de97-81f8-4e1d-b585-6cda82139959
solver_control = (; max_round 	= 4,
					maxiters 	= 20,
              		tol_round 	= 1.0e-9,
              		verbose 	= "a",
              		reltol 		= 1.0e-8,
              		tol_mono 	= 1.0e-10)

# ╔═╡ 07c4aa03-7be0-483e-a108-0eb3faeb7d88
function simulate_CO2R(grid, celldata; voltages = (-1.5:0.1:0.0) * V, kwargs...)
    kwargs 	 	= merge(solver_control, kwargs) 
    cell        = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)
	ivresult    = ivsweep(cell; voltages = (-1.5:0.1:0.0), store_solutions=true, kwargs...)
	cell, ivresult
end;

# ╔═╡ 4056a626-86d3-4884-ac54-886dbb23e6d8
md"""
#### Plotting Functions
"""

# ╔═╡ a29478aa-139d-4367-bf05-a2a82bd44163
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

# ╔═╡ 7ee0629e-a1de-456a-9627-956f92806ed6
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
            	width: 400px;
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

# ╔═╡ e32e89bd-005e-4af3-a6f3-8ebd7730471c
floataside(
    @bind guidata confirm(
        #! format: off
        PlutoUI.combine() do Child
md"""
  __User Data__	``\quad O + ne^- \leftrightharpoons R``
 - ``z_R``: $(Child("zR", NumberField(-2:2;default=-1))) 
   ``\quad n:`` $(Child("n", NumberField(0:2;default=1)))
   ``\quad κ:`` $(Child("κ", NumberField(0:10;default=10)))
 - ``M/(mol/L)`` $(Child("fgmol",NumberField(0.1:0.1:4;default=0.1)))
   ``M_{bg}/(mol/L)`` $(Child("bgmol",NumberField(0.1:0.1:4;default=2)))
 - scanrate/``(V/s)``: $(Child("scanrate", TextField(6;default="0.1")))
   nperiods: $(Child("nperiods", NumberField(1:10;default=1)))
 - ``L/μm``: $(Child("L", TextField(10;default="80")))
 - Double64: $(Child("double64",CheckBox()))
 - tunnel: $(Child("tunnel", CheckBox()))
   ``\quad β/cm^{-1}``: $(Child("β", TextField(10;default="1.0e8")))
"""
        end,
        #! format: on
        label = "Submit"
    );
    top = 50
)

# ╔═╡ a629e8a1-b1d7-42d8-8c17-43475785218e
begin
    Vmax = 2 * V

    L = parse(Float64, guidata.L) * μm

    hmin = 1.0e-5 	* μm

    hmax = 1.0 		* μm 

    X = ExtendableGrids.geomspace(0, L, hmin, hmax)

    grid = ExtendableGrids.simplexgrid(X)
end;

# ╔═╡ cf646a34-bd94-49af-8f8e-ec06446e18ca
begin
	reaction_arg = model == elydata_Gold ? (; reaction) : NamedTuple()
	sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = model, reaction_arg...)
end

# ╔═╡ 966ed6ab-d6fa-43f1-9ddb-45eb024d949c
result_pnp = capscalc(sys_pnp, molarities)

# ╔═╡ b8608787-6f71-44e1-90c9-bad6544bc4c0
sys_pb = PBSystem(grid; celldata = deepcopy(model), bcondition = pb_bcondition)

# ╔═╡ f2ba0e8a-4a9f-4b98-85b8-d54c71fd3616
result_pb = capscalc(sys_pb, molarities)

# ╔═╡ 85856abf-ee16-424a-ac06-97f76e32e444
# ╠═╡ skip_as_script = true
#=╠═╡
let
    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, layout = (2, 2), size = 	(650, 650))

    capsplot(vis[1, 1], result_sy, "Algebraic pressure")
    capsplot(vis[1, 2], result_pp, "Pressure Poisson")
    capsplot(vis[2, 1], result_pb, "Poisson-Boltzmann")
    capsplot(vis[2, 2], result_pnp, "Poisson-Nernst-Planck")

    reveal(vis)
end
  ╠═╡ =#

# ╔═╡ c4c62b30-6e5b-40ba-b922-4ed40d04f1ea
#=╠═╡
let
	if is_Landstorfer
		f = Figure()
	    result = result_pb
	    l = 1 / length(result)
		ϕ0_pzc = 0.972
		
		ax = Axis(f[1, 1], xlabel="φ / (V vs φ_pzc)", ylabel="dlcaps / (μF / cm²)", title="CSV Plot")
	
		if model_choice == "Landstorfer_NaClO₄ model"
			Low_c0 = lines!(ax, Landstorfer_NaClO₄_5mM.voltages .+ ϕ0_pzc, Landstorfer_NaClO₄_5mM.dlcaps, color = :darkblue, linestyle = :dash)
			High_c0 = lines!(ax, Landstorfer_NaClO₄_100mM.voltages .+ ϕ0_pzc, Landstorfer_NaClO₄_100mM.dlcaps, color = :red, linestyle = :dash)
		else	
			Low_c0 = lines!(ax, Landstorfer_NaF_5mM.voltages .+ ϕ0_pzc, Landstorfer_NaF_5mM.dlcaps, color = :darkblue, linestyle = :dash)
			High_c0 = lines!(ax, Landstorfer_NaF_100mM.voltages .+ ϕ0_pzc, Landstorfer_NaF_100mM.dlcaps, color = :red, linestyle = :dash)
		end
		
	    k = []  
	    legend_labels = [model_choice*"\t 5mM", model_choice*"\t 100mM"]  
		
	    for i in 1:length(result)
	        c = RGB(i * l, 0.0, 1 - i * l)
			push!(k, lines!(ax, result_pb[i].voltages, result_pb[i].dlcaps / (μF / cm^2), color=c, label="Result $i"))
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
  ╠═╡ =#

# ╔═╡ e181c648-7f4e-473a-92ed-6fde8c177202
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 88d38a68-1f8a-425a-bbae-90355a2213d0
sys_Default = PBSystem(grid; celldata = deepcopy(elydata_Default), bcondition = pb_bcondition)

# ╔═╡ 88f8d0a0-e5bb-4e1f-89ed-425bcc9a93b1
cell, result = simulate_CO2R(grid, elydata_Gold; voltages)

# ╔═╡ 171e028a-dffc-4a87-adf8-311942ebb1d6
begin 
	sawtooth = SawTooth(
        scanrate = parse(Float64, guidata.scanrate),
        vmin = -1.0, vmax = 1.0
    )
	    const nperiods = guidata.nperiods

end

# ╔═╡ 2259f2f2-f83a-4195-86ee-1fe7ddefd2f3
function sweep(pnpdata; eneutral = true, tunnel = false, bikerman = true)
    celldata = deepcopy(pnpdata)
    celldata.eneutral = eneutral
    pnpcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = model, reaction_arg...)
    return result = cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
    )

end

# ╔═╡ bcfc1095-b478-4f52-84ed-65b9513129e7
# ╠═╡ show_logs = false
pnpresult = sweep(elydata_Gold; eneutral = false, tunnel = false)

# ╔═╡ e1d3786d-d80b-46e9-8450-c9b25cdffc5e
let
    fig = Figure(size = (600, 200))
    ax = Axis(fig[1, 1], yscale = log10)
    T = pnpresult.times
    #lines!(ax,T, voltages.(T))
    lines!(ax, T[2:end], T[2:end] - T[1:(end - 1)])
    fig
end

# ╔═╡ e8af7132-3b5d-4cc6-860d-1951822bede4
# ╠═╡ show_logs = false
nnpresult = sweep(elydata_Gold; eneutral = true, tunnel = false)

# ╔═╡ 13bad9f1-33fa-4d14-a85a-0a0b76b81db3
let
    fig = Figure(size = (650, 400))
    ax = Axis(fig[1, 1])
    lines!(
        ax, voltages(pnpresult), currents(pnpresult, ihplus, electrode = :we),
        color = RGBf.(range(0.1, 1, length(voltages(pnpresult))), 0.0, 0.0)
    )
    lines!(ax, voltages(nnpresult), currents(nnpresult, ico2), color = :gray)
    #ylims!(-0.0001, 0.0001)
    fig
end

# ╔═╡ 84c0c523-e5df-4399-93b3-25c7fb7558d8
let
    fig = Figure()
    ax = Axis(fig[1, 1])
    T = 0:1.0e-3:10
    lines!(ax, T, sawtooth.(T))
    fig
end

# ╔═╡ ba2d77f3-2f8d-43af-a4ac-aebc551189e8
floataside(
    md"""
    Show only pH: $(@bind useonly_pH PlutoUI.CheckBox(default=false))

    """, top = 330
)

# ╔═╡ 7ea1c62c-a606-426c-94c0-f3d5778207f6
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

# ╔═╡ 558b7dcf-4e50-4f1f-a1d5-dc6d1c611ea5
plot1d(result, elydata_Gold)

# ╔═╡ c5107d7a-383f-40bf-bdb7-5bc8a6d79681
floataside(
    md"""
    __Time:__ $(@bind it PlutoUI.Slider(1:length(pnpresult.tsol.t)-1, show_value=false))
    """, top = 395
)


# ╔═╡ 7896e772-4390-4653-afcf-d7abb8063598
let
    ZR = zstr(guidata.zR)
    ZO = zstr(guidata.zR + guidata.n)

	species = getproperty.(bulk, :name)
	colors = getproperty.(bulk, :color)	
	
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
    ylims!(ax1, 1.0e-10, 1.0e2)

    lines!(ax1, XX, ru[ico2, 2:end], color = bulk.colors[2], linestyle = :solid, label = L"O^{%$(ZO)}, pnp")
	#lines!(ax1, XX, ru[ikplus, 2:end], color = :red, linestyle = :solid, label = 		L"O^{%$(ZO)}, pnp")
    lines!(ax1, XX, ru[ico, 2:end], color = :blue, linestyle = :solid, label = 			L"R^{%$(ZR)}, pnp")
    lines!(ax1, XX, ru[ihco3, 2:end], color = :lightgreen, linestyle = :solid, label 	= L"S^-, pnp")
    lines!(ax1, XX, ru[ihplus, 2:end], color = :orange, linestyle = :solid, label = 	L"X^+, pnp")

    Legend(
        fig[1, 2], ax1; labelsize = 10,
        backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)

    )

    fig
end

# ╔═╡ d6f9f77a-948b-4583-86a5-ef91be3477ee
html"""<hr>"""

# ╔═╡ Cell order:
# ╟─ef660f6f-9de3-4896-a65e-13c60df5de1e
# ╠═2b901eca-db3b-4ad2-b0ee-e031854c57fa
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═8bf52bab-4830-4936-b3c9-78c1ccd0406d
# ╠═b95ff168-68dc-4172-be97-df5362be6c48
# ╟─4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╠═852d9c74-b2aa-49c0-9bd0-0ccb1afcc520
# ╟─ac27c318-9a00-4287-bd36-3a97d65b5459
# ╠═fe48d05b-99bd-48b4-a044-4dd8e8d18b5d
# ╠═b3b55993-cb87-43c6-a084-22fa33f02f8c
# ╠═5fe96d0b-7bd0-4183-901d-727e966d434b
# ╠═a23eece5-8e94-4c0b-b487-e742a37e714e
# ╠═59855587-c6c2-4af6-a713-0b710cf2b0fe
# ╟─3be02c97-5c28-4370-97c1-e3f9faaba62a
# ╠═595715e5-f108-4167-b104-ac7c6f652e48
# ╠═53ae411f-42e8-41b3-ad38-e0189a001acf
# ╠═6b69df37-8754-457c-93cc-e9bacfa47d9a
# ╠═045573a4-aff5-4262-aca2-ead598a37bf2
# ╟─25bafe0e-f2fc-4a5c-828e-89fdccbc250c
# ╟─6d1f0f93-876a-4ec9-9763-f1abcb36cc07
# ╟─0bf42bde-0c30-4df2-8fa0-82e493add078
# ╟─e9f29b23-0f46-4515-9ba3-06cd25c1d741
# ╠═c11fdb45-b46e-40b6-b5a7-9c115aa5fee5
# ╟─949d126e-d863-40f4-b902-8b3b4c97b1b5
# ╟─5ddce46c-22a5-427a-8cb7-f45f216cefe0
# ╟─e2ab9bb4-a1b5-4d49-a760-013ad0fc4c68
# ╟─3a9940b9-c7e1-484c-bbf2-14c0c69d685b
# ╠═12cbfb8b-edb6-4335-8d80-0d6fe0eb9d3a
# ╟─1c94b9ba-429d-44fc-887e-e028ba070cc9
# ╠═1e52766d-12a9-46cd-be96-5c13a046944f
# ╠═a629e8a1-b1d7-42d8-8c17-43475785218e
# ╠═442fe098-497b-404f-80a0-880bc95d5e02
# ╠═c381803b-daad-4778-8d79-5abcecbce9ee
# ╠═70e1a34b-9041-4151-91aa-4dd7907a5b13
# ╟─6df5efba-c6bc-4bc9-ae86-514e8171c8f1
# ╟─237a3b9e-60bc-48c8-b5d1-ea954e65c401
# ╠═72f17b01-34cc-4196-99ae-a64af3f864c4
# ╟─38061646-9c66-4f9c-a0b5-5090dc62f8fe
# ╠═267a7233-e34a-4c4a-8627-a49bb45ccd25
# ╠═398b3511-4f7c-4436-9fe8-8edd76e3e0e7
# ╟─e114ec0d-13d3-4455-b1c9-d1c5d76671d9
# ╠═ca3bd6ba-1b3d-42c7-b008-8012b06368e4
# ╟─9b1dc273-9938-43a0-ac10-1928a80f89d8
# ╠═53cdf6d7-a025-49e0-af7b-cc0838cfb422
# ╠═cf646a34-bd94-49af-8f8e-ec06446e18ca
# ╠═966ed6ab-d6fa-43f1-9ddb-45eb024d949c
# ╟─98464285-2bd4-4631-8c4f-8790fe15cb93
# ╠═a8e26e1a-a9ac-4d51-b09c-7951acd4b4b7
# ╠═b8608787-6f71-44e1-90c9-bad6544bc4c0
# ╠═88d38a68-1f8a-425a-bbae-90355a2213d0
# ╠═f2ba0e8a-4a9f-4b98-85b8-d54c71fd3616
# ╟─e5dbe1db-3b69-49ed-9f32-e65579069c47
# ╠═171e028a-dffc-4a87-adf8-311942ebb1d6
# ╠═2259f2f2-f83a-4195-86ee-1fe7ddefd2f3
# ╟─2fb75c2a-c877-4c78-abf8-6f89706a58fe
# ╟─289d2c59-e920-47fe-b9ad-cb0a33ef0c9c
# ╟─c75a852d-e3b8-46e5-bcdd-5c41aef36c64
# ╟─791ccb34-e761-4e65-a9ef-95eac5395376
# ╟─1fcfbad3-2fad-4eee-a1a3-031dc29c9083
# ╟─0b6f33b9-41d4-48fd-8026-8a3bddcc1989
# ╠═85856abf-ee16-424a-ac06-97f76e32e444
# ╟─a7677fc3-a83d-4fab-8d00-b5c2454056c9
# ╟─d18fe756-b0b9-44d7-8872-6b7812108c16
# ╟─4fe01b26-ad9f-44b8-8900-916a5aeb5ad6
# ╠═a22a5421-05bf-484f-a2d3-91a06a0c6476
# ╟─a4a01dcb-8c02-441b-ba25-8e8c062d7d58
# ╟─c4c62b30-6e5b-40ba-b922-4ed40d04f1ea
# ╠═7cb51a9b-6357-4ec6-af6a-9113dea61661
# ╠═87f2b4c4-b163-4ae2-86b6-0266dff1da19
# ╟─4c1f6b31-ce09-4fba-b827-460e8a0d7e1a
# ╟─0e734e72-fc3a-48c7-b1d5-c0a380768eec
# ╟─fae68c38-be85-4718-8ee6-f900150e2b9a
# ╟─e181c648-7f4e-473a-92ed-6fde8c177202
# ╟─43a3d5f4-16ce-4402-bb95-d759b07bd573
# ╟─bcfc1095-b478-4f52-84ed-65b9513129e7
# ╟─e8af7132-3b5d-4cc6-860d-1951822bede4
# ╠═84c0c523-e5df-4399-93b3-25c7fb7558d8
# ╠═e1d3786d-d80b-46e9-8450-c9b25cdffc5e
# ╠═13bad9f1-33fa-4d14-a85a-0a0b76b81db3
# ╠═7896e772-4390-4653-afcf-d7abb8063598
# ╠═1f7971ad-80cc-4bc1-a2f7-a912186537f7
# ╟─05703d80-3299-4692-9d23-f44c2f371f7b
# ╟─f88ca8d4-ecbb-4eee-b4af-5854fbd16e33
# ╠═7ab1de97-81f8-4e1d-b585-6cda82139959
# ╠═07c4aa03-7be0-483e-a108-0eb3faeb7d88
# ╠═88f8d0a0-e5bb-4e1f-89ed-425bcc9a93b1
# ╠═558b7dcf-4e50-4f1f-a1d5-dc6d1c611ea5
# ╠═4056a626-86d3-4884-ac54-886dbb23e6d8
# ╠═7ea1c62c-a606-426c-94c0-f3d5778207f6
# ╟─a29478aa-139d-4367-bf05-a2a82bd44163
# ╟─7ee0629e-a1de-456a-9627-956f92806ed6
# ╟─e32e89bd-005e-4af3-a6f3-8ebd7730471c
# ╟─ba2d77f3-2f8d-43af-a4ac-aebc551189e8
# ╟─c5107d7a-383f-40bf-bdb7-5bc8a6d79681
# ╟─d6f9f77a-948b-4583-86a5-ef91be3477ee
# ╟─ecde8a95-f660-4b53-9296-cd4d4022c95f
