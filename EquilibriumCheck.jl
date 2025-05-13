### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

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
end
  ╠═╡ =#

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using PlutoUI, Latexify
    using VoronoiFVM
    using ExtendableGrids
    using LinearAlgebra
    using NLsolve
    using Unitful
    using LessUnitful
	using DelimitedFiles
	using PreallocationTools
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

# ╔═╡ ef660f6f-9de3-4896-a65e-13c60df5de1e
md"""
# Double layer capcacitance comparison
"""

# ╔═╡ 4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="",depth=5)
  ╠═╡ =#

# ╔═╡ 852d9c74-b2aa-49c0-9bd0-0ccb1afcc520
pkgdir(LiquidElectrolytes)

# ╔═╡ 462f512b-9b92-442b-bae6-b4aa168b32ee
# ╠═╡ disabled = true
#=╠═╡
begin
	mylog=RLog()
	LiquidElectrolytes.rlog(x::Number)=mylog(x)
end
  ╠═╡ =#

# ╔═╡ 920b7d84-56c6-4958-aed9-fc67ba0c43f6
md"""
## 1. Intro

This code implements the model described in
[Müller, R., Fuhrmann, J., & Landstorfer, M. (2020). Modeling polycrystalline electrode-electrolyte interfaces: The differential capacitance. Journal of The Electrochemical Society, 167(10), 106512](https://iopscience.iop.org/article/10.1149/1945-7111/ab9cca/meta)


The code is part of the LiquidElectrolytes.jl package.

Equation numbers refer to the paper.


Concentrations are given in number densities, and calculations are done in mole fractions.
"""

# ╔═╡ 87ac16f4-a4fc-4205-8fb9-e5459517e1b8
md"""
If not stated otherwise, all calculations and calculation results are in coherent SI units.
"""

# ╔═╡ 7d77ad32-3df6-4243-8bad-b8df4126e6ea
md"""
## 2. Model data
"""

# ╔═╡ 4cabef42-d9f9-43fe-988e-7b54462dc775
md"""
#### EquilibriumData
"""

# ╔═╡ ba2428ac-fdf7-4eae-acf5-dc0f20e876ee
begin
	@unitfactors mol dm m s K μm bar Pa μA Å;
	@phconstants N_A c_0 k_B e h
end

# ╔═╡ 30c6a176-935b-423f-9447-86f78746322f
md"""
#### debyelength(data)

```math
L_{Debye}=\sqrt{ \frac{(1+χ)ε_0k_BT}{e^2n_E}}
```
"""

# ╔═╡ f3049938-2637-401d-9411-4d7be07c19ca
md"""
#### set_molarity!(data,M)
"""

# ╔═╡ a21545da-3b53-47af-b0c4-f253b37dc84f
md"""

#### dlcap0(data)
Double layer capacitance at $φ=0$
```math
C_{dl,0}=\sqrt{\frac{2(1+χ) ε_0e^2 n_E}{k_BT}}
```
"""

# ╔═╡ 5a210961-19fc-40be-a5f6-033a80f1414d
md"""
Check with Bard/Faulkner: the value must be $(22.8u"μF/cm^2")
"""

# ╔═╡ 9b57f6ed-02f8-48ba-afa2-0766fe8c0c4c
md"""
#### Species indices
"""

# ╔═╡ 5fed71ec-35fb-4804-99ff-e1eaf18fac1b
begin
    const iφ = 1
    const ip = 2
    const iA = 1
    const iC = 2
end;

# ╔═╡ 5eca37ba-f858-45fb-a66a-3795327dfd18
md"""
## 3. Model equations
"""

# ╔═╡ a26cf11b-0ce1-4c1d-a64d-1917178ff676
md"""
### Mole fractions
Equilibrium expression for mole fractions (``α≥0``) (16)
```math
y_α(φ,p)=y_α^E\exp\left(\frac{-z_αe}{k_BT}(φ- φ^E)-\frac{v_α}{k_BT}(p-p^E)\right)
```
"""

# ╔═╡ cdd1d359-08fa-45a1-a857-e19f2adefcab
md"""
#### y_α(φ,p,α,data)

Ion molar fractions
"""

# ╔═╡ 188f67d8-2ae8-474c-8e58-68b8b4fde02e
function y_α(φ, p, α, data)
    η_φ = data.z[α] * data.e * (φ - data.E_ref)
    η_p = data.v[α] * (p * data.pscale - data.p_ref)
    return data.y_E[α] * exp(-(η_φ + η_p) / (data.kT))
end;

# ╔═╡ f70eed13-a6c2-4d54-9f30-113367afaf7d
md"""
#### y0(p,data)

Solvent molar fraction
"""

# ╔═╡ d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
y0(p, data) = data.y0_E * exp(-data.v0 * (p * data.pscale - data.p_ref) / (data.kT));

# ╔═╡ f6f004a6-d71b-4813-a363-9f51dc37e42a
md"""
### Poisson equation
Poisson equation (32a)

```math
-∇⋅(1+χ)ε_0∇φ = q(φ,p)
```
"""

# ╔═╡ 3810cc88-07f1-4741-853f-331e71c87923
md"""
#### poisson_flux!(f,u,edge,data)

VoronoiFVM flux function for left hand side of Poisson equation
"""

# ╔═╡ 0e2d20a1-5f26-4263-9a91-3b40b2c2996a
function poisson_flux!(f, u, edge, data)
    return f[iφ] = (1.0 + data.χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
end;

# ╔═╡ 824c610b-6e5e-48a3-be37-19104f52d1d9
md"""
#### Space charge expression
"""

# ╔═╡ 2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
md"""
Solvated ion volumes:
```math
	v_α=(1+κ_α)v_0
```
"""

# ╔═╡ b1e062c6-f245-4edc-aa02-871e2c776998
md"""
Incompressibility condition (14)
```math
\begin{aligned}
1&=∑\limits_α v_αn_α= n ∑\limits_α v_α y_α\\
n&=\frac1{∑\limits_α v_α y_α}
\end{aligned}
```
"""

# ╔═╡ c4cc940c-74aa-45f8-a2fa-6016d7c3c145
md"""
Space charge
```math
\begin{aligned}
q(φ,p)&=e∑\limits_α z_αn_α = ne∑\limits_α z_αy_α\\
      &=e\frac{∑\limits_α z_αy_α(\phi,p)}{∑\limits_α v_α y_α(\phi,p)}
\end{aligned}
```
"""

# ╔═╡ b07246b8-aec5-4161-8879-8cefb350aced
function spacecharge(φ, p, data)
    y = y0(p, data)
    sumyz = zero(eltype(p))
    sumyv = data.v0 * y
    for α in 1:(data.N)
        y = y_α(φ, p, α, data)
        sumyz += data.z[α] * y
        sumyv += data.v[α] * y
    end
    return data.e * sumyz / sumyv
end

# ╔═╡ b41838bb-3d5b-499c-9eb5-137c252ae366
md"""
#### Sum of mole fractions
"""

# ╔═╡ a468f43a-aa20-45dc-9c21-77f5adf2d700
function ysum(φ, p, data)
    sumy = y0(p, data)
    for α in 1:(data.N)
        sumy += y_α(φ, p, α, data)
    end
    return sumy
end

# ╔═╡ 978bf1d3-4758-4d01-b1e5-8aed1db9024f
md"""
#### spacecharge\_and\_ysum!(f,u,node,data)

VoronoiFVM reaction function. This assumes that terms are on the left hand side.
In addition to the space charge it calculates the residuum of  the 
definition of ``y_\alpha`` (32b):
```math
∑_α y_α(φ,p)=1
```
However, direct usage of this equation leads to slow convergence of Newton's method.
So we use

```math
\log\left(∑_α y_α(φ,p)\right)=0
```

instead.
"""

# ╔═╡ 13fc2859-496e-4f6e-8b22-36d9d55768b8
md"""
#### update_derived!(data)

Update derived data in data record.

Calculate bulk mole fractions from incompressibiltiy:
```math
\begin{aligned}
∑\limits_αv_αn_α^E&=1\\
n_0^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)\\
n^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)+ ∑\limits_{α>0}n_α^E\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}(v_α-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}((1+ κ_α)v_0-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}κ_αv_0n_α^E\right)\\
   &=\frac1{v_0}-∑\limits_{α>0}κ_αn_α^E\\
y_α^E&=\frac{n_α^E}{n^E}
\end{aligned}
```
"""

# ╔═╡ 32db42f3-5084-4908-9b53-59291b6133c5
function derived(κ, v0, n_E, T)
    c0 = 1 / v0
    barc = 0.0
    v = (1.0 .+ κ) .* v0
    N = length(κ)
    for α in 1:N
        barc += n_E[α]
        c0 -= n_E[α] * (1 + κ[α])
    end
    barc += c0
    y_E = n_E / barc
    y0_E = c0 / barc
    U_T = ph"k_B" * T / ph"e"
    return (; v, y_E, y0_E, U_T)
end;

# ╔═╡ 0d825f88-cd67-4368-90b3-29f316b72e6e
begin
    """
    	EquilibriumData
    Data structure containg data for equilibrum calculations
    """
    Base.@kwdef mutable struct EquilibriumData
        N::Int64 = 2                     # number of ionic species
        T::Float64 = 298.15 * ufac"K"        # temperature
        kT::Float64 = ph"k_B" * T             # temperature
        p_ref::Float64 = 1.0e5 * ufac"Pa"        # referece pressure
        pscale::Float64 = 1.0 * ufac"GPa"         # pressure scaling nparameter
        E_ref::Float64 = 0.0 * ufac"V"           # reference voltage
        n0_ref::Float64 = 55.508 * ph"N_A" / ufac"dm^3"  # solvent molarity
        v0::Float64 = 1 / n0_ref              # solvent molecule volume
        χ::Float64 = 15                    # dielectric susceptibility
        z::Vector{Int} = [-1, 1]                # ion charge numbers
        κ::Vector{Int} = [10, 10]               # ion solvation numbers
        molarity::Float64 = 0.1 * ph"N_A" / ufac"dm^3"
        n_E::Vector{Float64} = [molarity, molarity]  # bulk ion number densities
        μ_e::Vector{Float64} = [0.0]             # grain facet electron chemical potential

        e::Float64 = ph"e"
        ε_0::Float64 = ph"ε_0"

        v::Vector{Float64} = derived(κ, v0, n_E, T).v   # ion volumes
        y_E::Vector{Float64} = derived(κ, v0, n_E, T).y_E # bulk ion mole fractions
        y0_E::Float64 = derived(κ, v0, n_E, T).y0_E       # bulk solvent mole fraction
        U_T::Float64 = derived(κ, v0, n_E, T).U_T     # Temperature voltage k_BT/e0
    end

    function EquilibriumData(electrolyte::AbstractElectrolyteData)
        return EquilibriumData(;
            N = electrolyte.nc,
            T = electrolyte.T,
            p_ref = electrolyte.p_bulk,
            pscale = electrolyte.pscale,
            E_ref = electrolyte.ϕ_bulk,
            n0_ref = ph"N_A" / electrolyte.v0,
            χ = electrolyte.ε - 1.0,
            z = -electrolyte.z,
            κ = electrolyte.κ,
            molarity = ph"N_A" * electrolyte.c_bulk[1],
            n_E = ph"N_A" * electrolyte.c_bulk
        )
    end
end

# ╔═╡ 00e536dc-34aa-4a1a-93de-4eb3f5e0a348
LiquidElectrolytes.debyelength(data::EquilibriumData) = sqrt((1 + data.χ) * data.ε_0 * data.kT / (ph"e"^2 * data.n_E[1]))

# ╔═╡ 1065b3e0-60bf-497c-b7fb-c5a065737f77
# ╠═╡ skip_as_script = true
#=╠═╡
debyelength(EquilibriumData(molarity=0.01ph"N_A"/ufac"dm^3"))|>u"nm"
  ╠═╡ =#

# ╔═╡ 5d6340c4-2ddd-429b-a60b-3de5570a7398
function set_molarity!(data::EquilibriumData, M_E)
    n_E = M_E * ph"N_A" / ufac"dm^3"
    data.molarity = n_E
    return data.n_E = fill(n_E, data.N)
end

# ╔═╡ 1d22b09e-99c1-4026-9505-07bdffc98582
LiquidElectrolytes.dlcap0(data::EquilibriumData) = sqrt(2 * (1 + data.χ) * ph"ε_0" * ph"e"^2 * data.n_E[1] / (ph"k_B" * data.T));

# ╔═╡ fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╠═╡ skip_as_script = true
#=╠═╡
let
    data=EquilibriumData()
    set_molarity!(data,0.01)
    data.χ=78.49-1
    cdl0=dlcap0(data)|>u"μF/cm^2"
    @assert cdl0 ≈ 22.84669184882525u"μF/cm^2"
end
  ╠═╡ =#

# ╔═╡ 3d9a47b8-2754-4a21-84a4-39cbeab12286
begin
    function update_derived!(data::EquilibriumData)
        (; κ, v0, n_E, T) = data
        return data.v, data.y_E, data.y0_E, data.U_T = derived(κ, v0, n_E, T)
    end
    update_derived!(::ElectrolyteData) = nothing
end

# ╔═╡ b1e333c0-cdaa-4242-b71d-b54ff71aef83
let
    data = EquilibriumData()
    set_molarity!(data, 0.01)
    update_derived!(data)
    sumyz = 0.0
    sumyv = data.y0_E * data.v0
    sumy = data.y0_E
    for α in 1:data.N
        v = (1.0 + data.κ[α]) * data.v0
        sumyz += data.y_E[α] * data.z[α]
        sumyv += data.y_E[α] * v
        sumy += data.y_E[α]
    end
    @assert sumy ≈ 1.0
end


# ╔═╡ 243d27b5-a1b8-4127-beec-d5643ad07855
md"""
### Bulk boundary condition
From (32d):
```math
∇ φ\to 0
```
we choose homogeneous Neumann boundary conditions 
```math
\partial_n φ=0
```

which do not need any implementation.
"""

# ╔═╡ 005289e8-6979-49fe-b20f-66afd207baea
md"""
### Electrode boundary condition
See (32c) !!! bug in paper

```math
φ|_{Σ_i}=\frac{1}{e} \mu_{e,i}- (E-E^{ref})
```
"""

# ╔═╡ cbd3fbab-e95a-41d1-98c2-3cd8aec9ce18
md"""
#### φ_Σ(ifacet,data,E)

Calculate potential boundary value for each facet from applied voltage `E`.
"""

# ╔═╡ 0c5ed337-9310-417d-a1f6-7d69dd8c377b
φ_Σ(ifacet, data, E) = data.μ_e[ifacet] / data.e - (E - data.E_ref);

# ╔═╡ 0bbd9482-d17d-4027-8eec-450807cff792
md"""
## 4. System setup and solution
"""

# ╔═╡ 04f5584c-14af-4b68-9bcc-7f36b545bef7
md"""
#### create\_equilibrium\_system(grid,data)

Create equlibrium system, enable species and apply zero voltage
"""

# ╔═╡ c8822d32-affe-473e-8dbf-84aa83b3580c
md"""
#### apply_voltage!(sys,E)

Apply voltage `E` to system.
"""

# ╔═╡ d885ac23-ddfa-495c-b93b-54032c8a5c1f
function apply_voltage!(sys, E)
    data = sys.physics.data
    nbc = num_bfaceregions(sys.grid)
    nfacets = length(data.μ_e)
    @assert nbc > nfacets
    for ifacet in 1:nfacets
        boundary_dirichlet!(sys, iφ, ifacet, φ_Σ(ifacet, data, E))
    end
    return sys
end;

# ╔═╡ 93428d11-a3dc-4e29-ae6d-48ba37082c74
md"""
## 5. Postprocessing
"""

# ╔═╡ 7020a6f3-f49d-4fa3-bae2-2a6dad8a1fcd
md"""
#### calc_φ(sol,sys)

Obtain electrostatic potential from solution
"""

# ╔═╡ f3279037-01ed-4596-8e5a-86afe4c02c5f
calc_φ(sol, sys) = sol[iφ, :];

# ╔═╡ c5c8e124-be7e-4d06-ba23-dd72a88e4a18
md"""
#### calc_p(sol,sys)

Obtain pressure from solution
"""

# ╔═╡ 2afd54ca-4240-4f07-b38a-242ba0485b45
calc_p(sol, sys) = sol[ip, :] * sys.physics.data.pscale;

# ╔═╡ 55bd7b9a-a191-4a0b-9c6b-13733be5023e
md"""
#### c_num!(c,φ,p, data)
Calculate number concentration at discretization node
```math
	n_α=ny_α
```
"""

# ╔═╡ 3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
function c_num!(c, φ, p, data)
    y = y0(p, data)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        c[α] = y_α(φ, p, α, data)
        sumyv += c[α] * data.v[α]
    end
    return c ./= sumyv
end;

# ╔═╡ 97c5942c-8eb4-4b5c-8951-87ac0c9f396d
function c0_num!(c, φ, p, data)
    y = y0(p, data)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        c[α] = y_α(φ, p, α, data)
        sumyv += c[α] * data.v[α]
    end
    return y / sumyv
end;

# ╔═╡ 0c54efd0-f279-4dc6-8b00-ba092dd13f44
md"""
#### calc_cnum(sol,sys)

Obtain ion number densities from system
"""

# ╔═╡ 800dfed8-9f29-4138-96f8-e8bf1f2f00e6
function calc_cnum(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    nnodes = num_nodes(grid)
    conc = zeros(data.N, nnodes)
    for i in 1:nnodes
        @views c_num!(conc[:, i], sol[iφ, i], sol[ip, i], data)
    end
    return conc
end;

# ╔═╡ 24910762-7d56-446b-a758-d8e830fe9a09
function calc_c0num(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    nnodes = num_nodes(grid)
    c0 = zeros(nnodes)
    conc = zeros(data.N)
    for i in 1:nnodes
        @views c0[i] = c0_num!(conc, sol[iφ, i], sol[ip, i], data)
    end
    return c0
end;

# ╔═╡ 9fe3ca93-c051-426e-8b9a-cc59f59319ad
md"""
#### calc_cmol(sol,sys)

Obtain ion  molarities (molar densities in mol/L)  from system
"""

# ╔═╡ 2ee34d76-7238-46c2-94d1-a40d8b017af6
calc_cmol(sol, sys) = calc_cnum(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

# ╔═╡ 79cc671b-ef6e-42da-8641-61e43f221cb1
calc_c0mol(sol, sys) = calc_c0num(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

# ╔═╡ f4b2f509-0769-4df7-956e-e8bfc9ccd89a
md"""
#### calc_QBL(sol,sys)
Obtain boundary layer charge (35)
```math
Q^{BL}(E)=-\frac{1}{Σ} ∫_{Ω_E} q dx
```
"""

# ╔═╡ 65955950-2879-4b8c-bf73-d63e07d2ad96
md"""
Poly surface charge (34)

```math
Q_s(E)= ∑_i s_i \hat Q_s(E-E^{ref} + \frac1{e}\mu_{e,i})
```
"""

# ╔═╡ d8f80c62-b2d6-456f-9650-e8102e968673
md"""
Single surface charge (26)

```math
\hat Q_s = \frac{∑_\alpha z_αey_{s,α} +\sum_α \sum_β ν_{αβ}z_αey_{s,β}}{a_V^{ref}+ ∑_{α} a_α^{ref}z_αey_{s,α}} 
```

We assume that there are no surface reactions, so we assume ``\hat Q_s=0`` and ``Q_s=0``.
"""

# ╔═╡ 77f913ea-f89f-48f6-9dd2-e7cd0b6150b6
md"""
## 6. Solution scenarios
"""

# ╔═╡ bb6ef288-373f-4944-bc85-37ab327dc4d5
md"""
#### dlcapsweep_equi(sys)

Calculate double layer capacitance. Return vector of voltages `V` and vector of double layer capacitances `C`.
"""

# ╔═╡ 7a607454-7b75-4313-920a-2dbdad258015
md"""
## 7. The pressure Poisson equation
"""

# ╔═╡ 9cb8324c-896f-40f8-baa8-b7d47a93e9f5
md"""
An alternative possibility to handle the pressure has been introduced in 

[J. Fuhrmann, “Comparison and numerical treatment of generalised Nernst–Planck models,” Computer Physics Communications, vol. 196, pp. 166–178, 2015.](https://dx.doi.org/10.1016/j.cpc.2015.06.004).

Starting with the momentum balance in mechanical equilibrium
```math
	\nabla p = -q\nabla \varphi
```
by taking the divergence on both sides of the equation, one derives the pressure Poisson problem
```math
\begin{aligned}
	-\Delta p &= \nabla\cdot q\nabla \varphi & \text{in}\; \Omega\\
      p&=p_{bulk} & \text{on}\; \Gamma_{bulk}\\
	(\nabla p + q\nabla \varphi)\cdot \vec n &=0 & \text{on}\; \partial\Omega\setminus\Gamma_{bulk}\\
\end{aligned}
```
"""

# ╔═╡ 003a5c0b-17c7-4407-ad23-21c0ac000fd4
md"""
The bulk Dirichlet boundary condition for the pressure is necessary to make the solution unique. It is reasonable to set the ``\varphi`` to a bulk value at ``\Gamma_{bulk}`` as well, and to calculate ``p_{bulk}`` from the molar fraction sum constraint.
"""

# ╔═╡ e1c13f1e-5b67-464b-967b-25e3a93e33d9
function spacecharge!(f, u, node, data)
    φ = u[iφ]
    p = u[ip]
    return f[iφ] = -spacecharge(u[iφ], u[ip], data)
end;

# ╔═╡ 64e47917-9c61-4d64-a6a1-c6e8c7b28c59
function poisson_and_p_flux!(f, u, edge, data)
    f[iφ] = (1.0 + data.χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
    q1 = spacecharge(u[iφ, 1], u[ip, 1], data)
    q2 = spacecharge(u[iφ, 2], u[ip, 2], data)
    return f[ip] = (u[ip, 1] - u[ip, 2]) + (u[iφ, 1] - u[iφ, 2]) * (q1 + q2) / (2 * data.pscale)
end;

# ╔═╡ 48670f54-d303-4c3a-a191-06e6592a2e0a
function ysum(sys, sol)
    data = sys.physics.data
    n = size(sol, 2)
    sumy = zeros(n)
    for i in 1:n
        sumy[i] = ysum(sol[iφ, i], sol[ip, i], data)
    end
    return sumy
end

# ╔═╡ 042a452a-1130-4a56-a1b9-b2674803e445
function spacecharge_and_ysum!(f, u, node, data)
    φ = u[iφ]
    p = u[ip]
    f[iφ] = -spacecharge(φ, p, data)
    return f[ip] = log(ysum(φ, p, data)) # this behaves much better with Newton's method
end;

# ╔═╡ 6e3dbf34-c1c9-460b-9546-b2ee8ee99d68
function create_equilibrium_system(
        grid,
        data::EquilibriumData = EquilibriumData();
        Γ_bulk = 0
    )
    update_derived!(data)
    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_flux!,
        reaction = spacecharge_and_ysum!,
        species = [iφ, ip]
    )
    if Γ_bulk > 0
        boundary_dirichlet!(sys, iφ, Γ_bulk, 0.0)
    end
    return apply_voltage!(sys, 0)
end;

# ╔═╡ 49466829-9459-4dc8-85cc-c67460e290d2
calc_QBL(sol, sys) = VoronoiFVM.integrate(sys, spacecharge_and_ysum!, sol)[iφ, 1]

# ╔═╡ 77f49da5-ffd2-4148-93a6-f45382ba6d91
function dlcapsweep_equi(
        sys; vmax = 2 * ufac"V", nsteps = 21, δV = 1.0e-3 * ufac"V",
        molarity = nothing,
        verbose = false
    )

    if !isnothing(molarity)
        error("The molarity kwarg of dlcapsweep_equie has been removed. Pass the molarity information with set_molarity!.")
    end


    data = sys.physics.data
    update_derived!(data)
    apply_voltage!(sys, 0)

    c = VoronoiFVM.NewtonControl()
    #	c.damp_growth=1.1
    c.verbose = verbose
    c.tol_round = 1.0e-10
    c.max_round = 3
    c.damp_initial = 0.01
    c.damp_growth = 2

    inival = solve(sys; inival = 0, control = c)
    vstep = vmax / (nsteps - 1)

    c.damp_initial = 1

    function rundlcap(dir)
        volts = zeros(0)
        caps = zeros(0)
        volt = 0.0
        sol = inival
        for iv in 1:nsteps
            apply_voltage!(sys, volt)
            c.damp_initial = 1
            sol = solve(sys; inival = sol, control = c)

            Q = calc_QBL(sol, sys)
            apply_voltage!(sys, volt + dir * δV)
            c.damp_initial = 1
            sol = solve(sys; inival = sol, control = c)
            Qδ = calc_QBL(sol, sys)
            push!(caps, (Q - Qδ) / (dir * δV))
            push!(volts, volt)
            volt += dir * vstep
        end
        return volts, caps
    end
    Vf, Cf = rundlcap(1)
    Vr, Cr = rundlcap(-1)
    return vcat(reverse(Vr), Vf), vcat(reverse(Cr), Cf)
end

# ╔═╡ 7bf3a130-3b47-428e-916f-4a0ec1237844
function create_equilibrium_pp_system(
        grid,
        data::EquilibriumData = EquilibriumData();
        Γ_bulk = 0
    )
    update_derived!(data)

    sys = VoronoiFVM.System(
        grid; data = data, flux = poisson_and_p_flux!,
        reaction = spacecharge!, species = [iφ, ip]
    )
    if Γ_bulk > 0
        logysum!(y, p) = y[1] = log(ysum(0.0, p[1], data))
        res = nlsolve(logysum!, [0.0]; autodiff = :forward, method = :newton, xtol = 1.0e-10, ftol = 1.0e-20)
        boundary_dirichlet!(sys, iφ, Γ_bulk, 0.0)
        boundary_dirichlet!(sys, ip, Γ_bulk, res.zero[1])
    end
    return apply_voltage!(sys, 0)
end;

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


# ╔═╡ 5fe96d0b-7bd0-4183-901d-727e966d434b
begin
	#const voltages = (-1.0:0.01:1.0) * V
	
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

# ╔═╡ b1a470b7-17af-456f-8e58-c64cd697f0bd
# ╠═╡ disabled = true
#=╠═╡
begin
	const bulk = let 
		bulk = [
		BulkSpecies(;name="HCO₃⁻", z=-1, D=1.185e-9, c_bulk=0.091, color=:brown),
		BulkSpecies(;name="CO₃²⁻", z=-2, D=0.923e-9, c_bulk=2.68e-5, color=:violet),
		BulkSpecies(;name="CO₂", z=0, D=1.91e-9, c_bulk=0.033, a=0.0, color=:red),
		BulkSpecies(;name="OH⁻", z=-1, D=5.273e-9, c_bulk=10^(pH-14), color=:green),
		BulkSpecies(;name="H⁺", z=1, D=9.310e-9, c_bulk=10^(-pH), a=0.0, color=:gray),
		BulkSpecies(;name="CO", z=0, D=2.23e-9, c_bulk=0.0, a=0.0, color=:blue)
		]
		push!(bulk, make_eneutral(bulk;name="K⁺",z=1, D=1.957e-9, a=8.2, color=:orange))
		sort(bulk, by=x->species_dict[x.name])
	end
end;
  ╠═╡ =#

# ╔═╡ 1afdcbff-29d9-4e09-b791-54c4ce55a30d
md"""
#### Size of BulkSpecies(Tempo)
"""

# ╔═╡ 59855587-c6c2-4af6-a713-0b710cf2b0fe
begin
	const at = 8.2
	const κt = 8.0
	const ak = 8.2
	const κk = 8.0
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
elydata = ElectrolyteData(;
                               nc = size(bulk)[1],
							   z     = getproperty.(bulk, :z),
							   D     = getproperty.(bulk, :D),
							   eneutral=false,
							   c_bulk = getproperty.(bulk, :c_bulk),
				    		   v0 	  = v0,
							   κ     = getproperty.(bulk, :κ),
							   v     = getproperty.(bulk, :v),
							   M0 	  = M0,
							   M     = getproperty.(bulk, :M),
							   #scheme = :act,
                               Γ_we = 1,
                               Γ_bulk = 2)

# ╔═╡ 25bafe0e-f2fc-4a5c-828e-89fdccbc250c
md"""
### Reaction
"""

# ╔═╡ 6d1f0f93-876a-4ec9-9763-f1abcb36cc07
md"""
#### Boundary Reaction
"""

# ╔═╡ c11fdb45-b46e-40b6-b5a7-9c115aa5fee5
begin
	catmap_params 		= CatmapInterface.parse_catmap_input("catmap_CO2R_template.mkm")
	rn 					= create_reaction_network(catmap_params)
	odesys 				= convert(ODESystem, rn; combinatoric_ratelaws=false)
	odesys 				= CatmapInterface.liquidize(odesys, catmap_params)
	vars 				= states(odesys)
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
		ps[paramsidx[odesys.σ]] = σ
		ps[paramsidx[odesys.γCO2_aq]] = γ_co2
		ps[paramsidx[odesys.aH2O_g]] = aH₂O
		ps[paramsidx[odesys.ϕ]] = u[iϕ]
		ps[paramsidx[odesys.ϕ_we]] = ϕ_we
		ps[paramsidx[odesys.local_pH]] = local_pH
		ps[paramsidx[odesys.γCO_aq]] = γ_co
		ps[paramsidx[odesys.βCOOHΔH2OΔele_t]] = 0.59

	    #println[1.0 / (1 - v[ikplus] * u[ikplus] / (mol/dm^3))]
		#@show size(f)

		if bnode.region == Γ_we && size(f,1) ≥ isurfaceend
			@views f_microkinetics!(
				f[isurfacestart:isurfaceend], 
				u[isurfacestart:isurfaceend],
				ps,
				nothing
			)
		#elseif bnode.region == Γ_we
		#	nothing
		end
		
		#f = 0
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

# ╔═╡ 00464966-2b1e-455c-a3a1-2af61c6649b7
dlcap_exact = 0.22846691848825248

# ╔═╡ c53791b6-6c12-483a-910f-6183149fac80
# ╠═╡ disabled = true
#=╠═╡
@test dlcap0(elydata) ≈ dlcap_exact
  ╠═╡ =#

# ╔═╡ 05334798-a072-41ae-b23e-f884baadb071
begin
    equidata = EquilibriumData()
    set_molarity!(equidata, 0.01)
    equidata.χ = 78.49 - 1
end

# ╔═╡ ddb3e60b-8571-465f-acf3-2403fb884363
@test dlcap0(equidata) |> unitfactor ≈ dlcap_exact

# ╔═╡ 512b631e-93ec-42fc-8416-8d10ca97f23d
# ╠═╡ disabled = true
#=╠═╡
@test dlcap0(EquilibriumData(elydata)) ≈ dlcap_exact
  ╠═╡ =#

# ╔═╡ a629e8a1-b1d7-42d8-8c17-43475785218e
begin
    Vmax = 2 * V

    L = 20nm

    hmin = 0.05 * nm

    hmax = 0.5 * nm

    X = ExtendableGrids.geomspace(0, L, hmin, hmax)

    grid = ExtendableGrids.simplexgrid(X)

    data = EquilibriumData(elydata)
end;

# ╔═╡ cdb7e8a1-dcdf-4e7a-9ecf-121f51b485c3
sys_sy = create_equilibrium_system(grid, data)

# ╔═╡ 31a1f686-f0b6-430a-83af-187df411b293
sys_pp = create_equilibrium_pp_system(grid, data, Γ_bulk = 2)

# ╔═╡ 442fe098-497b-404f-80a0-880bc95d5e02
inival = unknowns(sys_sy, inival = 0);

# ╔═╡ 1c0145d5-76b1-48c1-8852-de1a2668285a
molarities = [elydata.c_bulk,elydata.c_bulk*0.1, elydata.c_bulk*0.01, elydata.c_bulk*0.001]

# ╔═╡ 699eaadb-80f4-4230-9054-27df7c224c99
# ╠═╡ disabled = true
#=╠═╡
function capscalc(sys)
    result = []
    for imol in 1:length(molarities)
        if isa(sys.physics.data, EquilibriumData)
            set_molarity!(sys.physics.data, molarities[imol])
            t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
        else
            #sys.physics.data.c_bulk[1] .= molarities[imol] * ufac"mol/dm^3"
            t = @elapsed r = dlcapsweep(
                sys,
                voltages = range(-1, 1, length = 201),
				#reaction = reaction ###buffer reaction
            )
            volts = voltages(r)
            caps = r.dlcaps
        end
        cdl0 = dlcap0(sys.physics.data)
        @info "elapsed=$(t)"
        push!(
            result,
            (
                voltages = volts,
                dlcaps = caps,
                cdl0 = cdl0,
                molarity = molarities[imol],
            )
        )
    end
    return result
end
  ╠═╡ =#

# ╔═╡ 70e1a34b-9041-4151-91aa-4dd7907a5b13
function capscalc(sys)
    result = []
    #for imol in 1:length(molarities)
        if isa(sys.physics.data, EquilibriumData)
            #set_molarity!(sys.physics.data, molarities[imol])
            t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
        else
            #sys.physics.data.c_bulk[1] .= molarities[imol] * ufac"mol/dm^3"
            t = @elapsed r = dlcapsweep(
                sys,
                voltages = range(-1, 1, length = 201),
            )
            volts = voltages(r)
            caps = r.dlcaps
        end
        cdl0 = dlcap0(sys.physics.data)
        @info "elapsed=$(t)"
        push!(
            result,
            (
                voltages = volts,
                dlcaps = caps,
                cdl0 = cdl0,
                #molarity = molarities[imol],
            )
        )
    #end
    return result
end

# ╔═╡ 38061646-9c66-4f9c-a0b5-5090dc62f8fe
md"""
#### Algebraic pressure equation
"""

# ╔═╡ 398b3511-4f7c-4436-9fe8-8edd76e3e0e7
result_sy = capscalc(sys_sy)

# ╔═╡ e114ec0d-13d3-4455-b1c9-d1c5d76671d9
md"""
#### Pressure poisson problem
"""

# ╔═╡ ca3bd6ba-1b3d-42c7-b008-8012b06368e4
result_pp = capscalc(sys_pp)

# ╔═╡ 9b1dc273-9938-43a0-ac10-1928a80f89d8
md"""
#### Poisson Nernst-Planck from LiquidElectrolytes
"""

# ╔═╡ 53cdf6d7-a025-49e0-af7b-cc0838cfb422
function pnp_bcondition(f, u, bnode, data::ElectrolyteData)
    (; iϕ, Γ_we, ϕ_we) = data
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
	#if bnode.region == Γ_we
	#	we_breactions(f, u, bnode, data)
	#end
	#boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))
    return bulkbcondition(f, u, bnode, data)
end

# ╔═╡ cf646a34-bd94-49af-8f8e-ec06446e18ca
sys_pnp = PNPSystem(grid; bcondition = pnp_bcondition, celldata = elydata, reaction=reaction)

# ╔═╡ 966ed6ab-d6fa-43f1-9ddb-45eb024d949c
result_pnp = capscalc(sys_pnp)

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

	#boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))
    boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)

    return bulkbcondition(f, u, bnode, data)
	#return boundary_robin!(f, u, bnode, iϕ, Γ_bulk, C_gap, C_gap * (ϕ_we - ϕ_pzc))
end

# ╔═╡ 88d38a68-1f8a-425a-bbae-90355a2213d0
sys_pb = PBSystem(grid; celldata = deepcopy(elydata), bcondition = pb_bcondition)

# ╔═╡ f2ba0e8a-4a9f-4b98-85b8-d54c71fd3616
result_pb = capscalc(sys_pb)

# ╔═╡ 289d2c59-e920-47fe-b9ad-cb0a33ef0c9c
md"""
#### Result comparison
"""

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

# ╔═╡ 25a183d9-c6a2-4ac3-a283-a02e4e9231dd
@test resultcompare(result_pp, result_sy; tol = 5.0e-2)

# ╔═╡ 6d1d8ae2-6a9e-48c1-a545-1f7354125bf0
@test resultcompare(result_pb, result_pp; tol = 5.0e-3)

# ╔═╡ b43c5e74-5010-4870-a058-d3ad2c1ed548
@test resultcompare(result_pp, result_pnp; tol = 5.0e-1)

# ╔═╡ ee76e884-86e6-45f6-bbb2-c8e73daa5883
@test resultcompare(result_pb, result_pnp; tol = 5.0e-1)

# ╔═╡ 0b6f33b9-41d4-48fd-8026-8a3bddcc1989
md"""
#### Result plot

Compare with Fig 4.2 of [Fuhrmann (2015)](https://dx.doi.org/10.1016/j.cpc.2015.06.004)
"""

# ╔═╡ a22a5421-05bf-484f-a2d3-91a06a0c6476
# ╠═╡ skip_as_script = true
#=╠═╡
function capsplot(vis, result, title)
    hmol = 1 / length(result)
    for imol in 1:length(result)
        c = RGB(1 - imol * hmol, 0, imol * hmol)
		try
        	scalarplot!(
            vis, result[imol].voltages, (result[imol].dlcaps / (μF / cm^2)), limits=(-1, 100), xlimits=(-1.1, 1.1), bg = :transparent, color = c, clear = false, label = " ", markershape = :none, title = title, yscale=10, xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)"
        )
		catch 
			continue
		end
        scalarplot!(
            vis, [0], [result[imol].cdl0] / (μF / cm^2),
            clear = false, markershape = :circle, markersize = 8, label = "",
        )
    end
    return vis
end
  ╠═╡ =#

# ╔═╡ 85856abf-ee16-424a-ac06-97f76e32e444
# ╠═╡ skip_as_script = true
#=╠═╡
let
    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, layout = (2, 2), size = (650, 650), backgroundcolor = :transparent)

    capsplot(vis[1, 1], result_sy, "Algebraic pressure")
    capsplot(vis[1, 2], result_pp, "Pressure Poisson")
    capsplot(vis[2, 1], result_pb, "Poisson-Boltzmann")
    capsplot(vis[2, 2], result_pnp, "Poisson-Nernst-Planck")

    reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 87f2b4c4-b163-4ae2-86b6-0266dff1da19
#=╠═╡
function capsplot_v(vis, result)
    hmol = 1 / length(result)
	color = [:red, :green, :blue, :lightgray]
    for i in 1:length(result)
		scalarplot!(
            vis, result[i][1].voltages, result[i][1].dlcaps / (μF / cm^2), limits=(-1, 100), xlimits=(-1.1, 1.1), bg = :transparent, color = color[i], clear = false, label = "result", markershape = :none, yscale=10, xlabel = "φ / (V vs φ_pzc)", ylabel = "dlcaps / (μF / cm²)")
    end
    return vis
end
  ╠═╡ =#

# ╔═╡ 970871ac-a5f2-4e30-9420-489acdbe79f9
 sym = Symbol("result_" * "pp")

# ╔═╡ a90686f1-5f57-43e9-b22a-c8725775864f
#=╠═╡
let
    vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, size = (650, 650), backgroundcolor = :transparent)
    plots = []
    results = [
		#result_pp, 
		result_pb, 
		result_pnp, 
		#result_sy
		]

    for i in 1:length(results)
        result = results[i]
        try
            # 비어있는 결과가 아니면 push
            if !(isempty(result))
                push!(plots, result)
            end
        catch
            # 비어있거나 예외 발생 시 무시
            continue
        end
    end

    capsplot_v(vis, plots)
    reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 7afb1a46-2675-4b70-be39-5100fe2c2274
[result_sy, result_pp, result_pb, result_pnp][4][1]

# ╔═╡ b757aee6-49dc-4ead-abbf-2c60ddfad5a6
function save_dlcaps_to_csv(result, filepath)
    # DLCapSweepResult 객체 또는 객체가 담긴 벡터 처리
    volt = result[1].voltages
    cdl_vals = result[1].dlcaps / (μF / cm^2)
    df = DataFrame(voltage = volt, cdl = cdl_vals)
    CSV.write(filepath, df)
    println("Results saved to CSV: ", filepath)
end

# ╔═╡ 576a3999-c408-44b3-851a-d02c9290a374
save_dlcaps_to_csv(result_pb, "size_modified_pb")

# ╔═╡ 4c1f6b31-ce09-4fba-b827-460e8a0d7e1a
md"""
#### (Tempo) Solvation Number Plots
"""

# ╔═╡ 0e734e72-fc3a-48c7-b1d5-c0a380768eec
function caps(sys)
	dls = dlcapsweep(sys, voltages = range(-1, 1, length = 401))	
	return dls 
end

# ╔═╡ fae68c38-be85-4718-8ee6-f900150e2b9a
#=╠═╡
function capsplot_κ(vis, sys; n::Int=23)
    color = [RGB((i/n), 0.5, i/n) for i in 1:n]
    dls = LiquidElectrolytes.DLCapSweepResult[]
    κ_values = Float64[]

    sys = deepcopy(sys)
    κ = 1.0
    sys.physics.data.κ .= κ

    for j in 1:n
        try
            result = caps(sys)
            push!(dls, result)
            push!(κ_values, κ)

            scalarplot!(
                vis,
                result.voltages,
                result.dlcaps / (μF / cm^2),
                limits = (-1, 100),
                xlimits = (-1.1, 1.1),
                bg = :transparent,
                color = color[j],
                clear = false,
                label = "κ = $κ",
                markershape = :none,
                yscale = 10,
                xlabel = "φ / (V vs φ_pzc)",
                ylabel = "dlcaps / (μF / cm²)",
                backgroundcolor = :gray,
            )
        catch e
            @warn "caps 실패 at κ=$κ" exception=e
        end
        κ += 2.0
        sys.physics.data.κ .= κ
    end
	sys.physics.data.κ .= κt # Default Value
end
  ╠═╡ =#

# ╔═╡ 228c8672-03cc-4dc0-ba0c-74ee8928fb11
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ e181c648-7f4e-473a-92ed-6fde8c177202
#=╠═╡
let
	vis = GridVisualizer(Plotter = CairoMakie, legend = :lt, size = (650, 650), background = :gray)
	capsplot_κ(vis, sys_pb)
    reveal(vis)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─ef660f6f-9de3-4896-a65e-13c60df5de1e
# ╠═2b901eca-db3b-4ad2-b0ee-e031854c57fa
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═b95ff168-68dc-4172-be97-df5362be6c48
# ╟─4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╠═852d9c74-b2aa-49c0-9bd0-0ccb1afcc520
# ╠═462f512b-9b92-442b-bae6-b4aa168b32ee
# ╟─920b7d84-56c6-4958-aed9-fc67ba0c43f6
# ╟─87ac16f4-a4fc-4205-8fb9-e5459517e1b8
# ╟─7d77ad32-3df6-4243-8bad-b8df4126e6ea
# ╟─4cabef42-d9f9-43fe-988e-7b54462dc775
# ╠═0d825f88-cd67-4368-90b3-29f316b72e6e
# ╠═ba2428ac-fdf7-4eae-acf5-dc0f20e876ee
# ╟─30c6a176-935b-423f-9447-86f78746322f
# ╠═00e536dc-34aa-4a1a-93de-4eb3f5e0a348
# ╠═1065b3e0-60bf-497c-b7fb-c5a065737f77
# ╟─f3049938-2637-401d-9411-4d7be07c19ca
# ╠═5d6340c4-2ddd-429b-a60b-3de5570a7398
# ╟─a21545da-3b53-47af-b0c4-f253b37dc84f
# ╠═1d22b09e-99c1-4026-9505-07bdffc98582
# ╟─5a210961-19fc-40be-a5f6-033a80f1414d
# ╠═fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╟─9b57f6ed-02f8-48ba-afa2-0766fe8c0c4c
# ╠═5fed71ec-35fb-4804-99ff-e1eaf18fac1b
# ╟─5eca37ba-f858-45fb-a66a-3795327dfd18
# ╟─a26cf11b-0ce1-4c1d-a64d-1917178ff676
# ╟─cdd1d359-08fa-45a1-a857-e19f2adefcab
# ╠═188f67d8-2ae8-474c-8e58-68b8b4fde02e
# ╟─f70eed13-a6c2-4d54-9f30-113367afaf7d
# ╠═d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
# ╟─f6f004a6-d71b-4813-a363-9f51dc37e42a
# ╟─3810cc88-07f1-4741-853f-331e71c87923
# ╠═0e2d20a1-5f26-4263-9a91-3b40b2c2996a
# ╟─824c610b-6e5e-48a3-be37-19104f52d1d9
# ╟─2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
# ╟─b1e062c6-f245-4edc-aa02-871e2c776998
# ╟─c4cc940c-74aa-45f8-a2fa-6016d7c3c145
# ╠═b07246b8-aec5-4161-8879-8cefb350aced
# ╟─b41838bb-3d5b-499c-9eb5-137c252ae366
# ╠═a468f43a-aa20-45dc-9c21-77f5adf2d700
# ╟─978bf1d3-4758-4d01-b1e5-8aed1db9024f
# ╠═042a452a-1130-4a56-a1b9-b2674803e445
# ╟─13fc2859-496e-4f6e-8b22-36d9d55768b8
# ╠═32db42f3-5084-4908-9b53-59291b6133c5
# ╠═3d9a47b8-2754-4a21-84a4-39cbeab12286
# ╠═b1e333c0-cdaa-4242-b71d-b54ff71aef83
# ╟─243d27b5-a1b8-4127-beec-d5643ad07855
# ╟─005289e8-6979-49fe-b20f-66afd207baea
# ╟─cbd3fbab-e95a-41d1-98c2-3cd8aec9ce18
# ╠═0c5ed337-9310-417d-a1f6-7d69dd8c377b
# ╟─0bbd9482-d17d-4027-8eec-450807cff792
# ╟─04f5584c-14af-4b68-9bcc-7f36b545bef7
# ╠═6e3dbf34-c1c9-460b-9546-b2ee8ee99d68
# ╟─c8822d32-affe-473e-8dbf-84aa83b3580c
# ╠═d885ac23-ddfa-495c-b93b-54032c8a5c1f
# ╟─93428d11-a3dc-4e29-ae6d-48ba37082c74
# ╟─7020a6f3-f49d-4fa3-bae2-2a6dad8a1fcd
# ╠═f3279037-01ed-4596-8e5a-86afe4c02c5f
# ╟─c5c8e124-be7e-4d06-ba23-dd72a88e4a18
# ╠═2afd54ca-4240-4f07-b38a-242ba0485b45
# ╟─55bd7b9a-a191-4a0b-9c6b-13733be5023e
# ╠═3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
# ╠═97c5942c-8eb4-4b5c-8951-87ac0c9f396d
# ╟─0c54efd0-f279-4dc6-8b00-ba092dd13f44
# ╠═800dfed8-9f29-4138-96f8-e8bf1f2f00e6
# ╠═24910762-7d56-446b-a758-d8e830fe9a09
# ╟─9fe3ca93-c051-426e-8b9a-cc59f59319ad
# ╠═2ee34d76-7238-46c2-94d1-a40d8b017af6
# ╠═79cc671b-ef6e-42da-8641-61e43f221cb1
# ╟─f4b2f509-0769-4df7-956e-e8bfc9ccd89a
# ╠═49466829-9459-4dc8-85cc-c67460e290d2
# ╟─65955950-2879-4b8c-bf73-d63e07d2ad96
# ╟─d8f80c62-b2d6-456f-9650-e8102e968673
# ╟─77f913ea-f89f-48f6-9dd2-e7cd0b6150b6
# ╟─bb6ef288-373f-4944-bc85-37ab327dc4d5
# ╠═77f49da5-ffd2-4148-93a6-f45382ba6d91
# ╟─7a607454-7b75-4313-920a-2dbdad258015
# ╟─9cb8324c-896f-40f8-baa8-b7d47a93e9f5
# ╟─003a5c0b-17c7-4407-ad23-21c0ac000fd4
# ╠═e1c13f1e-5b67-464b-967b-25e3a93e33d9
# ╠═64e47917-9c61-4d64-a6a1-c6e8c7b28c59
# ╠═7bf3a130-3b47-428e-916f-4a0ec1237844
# ╠═48670f54-d303-4c3a-a191-06e6592a2e0a
# ╟─ac27c318-9a00-4287-bd36-3a97d65b5459
# ╠═fe48d05b-99bd-48b4-a044-4dd8e8d18b5d
# ╠═5fe96d0b-7bd0-4183-901d-727e966d434b
# ╠═a23eece5-8e94-4c0b-b487-e742a37e714e
# ╠═b1a470b7-17af-456f-8e58-c64cd697f0bd
# ╠═1afdcbff-29d9-4e09-b791-54c4ce55a30d
# ╠═59855587-c6c2-4af6-a713-0b710cf2b0fe
# ╟─3be02c97-5c28-4370-97c1-e3f9faaba62a
# ╠═595715e5-f108-4167-b104-ac7c6f652e48
# ╟─25bafe0e-f2fc-4a5c-828e-89fdccbc250c
# ╟─6d1f0f93-876a-4ec9-9763-f1abcb36cc07
# ╟─e9f29b23-0f46-4515-9ba3-06cd25c1d741
# ╠═c11fdb45-b46e-40b6-b5a7-9c115aa5fee5
# ╟─949d126e-d863-40f4-b902-8b3b4c97b1b5
# ╟─5ddce46c-22a5-427a-8cb7-f45f216cefe0
# ╟─e2ab9bb4-a1b5-4d49-a760-013ad0fc4c68
# ╟─3a9940b9-c7e1-484c-bbf2-14c0c69d685b
# ╠═12cbfb8b-edb6-4335-8d80-0d6fe0eb9d3a
# ╟─1c94b9ba-429d-44fc-887e-e028ba070cc9
# ╠═1e52766d-12a9-46cd-be96-5c13a046944f
# ╟─00464966-2b1e-455c-a3a1-2af61c6649b7
# ╠═c53791b6-6c12-483a-910f-6183149fac80
# ╟─05334798-a072-41ae-b23e-f884baadb071
# ╠═ddb3e60b-8571-465f-acf3-2403fb884363
# ╠═512b631e-93ec-42fc-8416-8d10ca97f23d
# ╠═a629e8a1-b1d7-42d8-8c17-43475785218e
# ╠═cdb7e8a1-dcdf-4e7a-9ecf-121f51b485c3
# ╠═31a1f686-f0b6-430a-83af-187df411b293
# ╠═442fe098-497b-404f-80a0-880bc95d5e02
# ╠═1c0145d5-76b1-48c1-8852-de1a2668285a
# ╠═699eaadb-80f4-4230-9054-27df7c224c99
# ╠═70e1a34b-9041-4151-91aa-4dd7907a5b13
# ╟─38061646-9c66-4f9c-a0b5-5090dc62f8fe
# ╠═398b3511-4f7c-4436-9fe8-8edd76e3e0e7
# ╟─e114ec0d-13d3-4455-b1c9-d1c5d76671d9
# ╠═ca3bd6ba-1b3d-42c7-b008-8012b06368e4
# ╟─9b1dc273-9938-43a0-ac10-1928a80f89d8
# ╠═53cdf6d7-a025-49e0-af7b-cc0838cfb422
# ╠═cf646a34-bd94-49af-8f8e-ec06446e18ca
# ╠═966ed6ab-d6fa-43f1-9ddb-45eb024d949c
# ╟─98464285-2bd4-4631-8c4f-8790fe15cb93
# ╠═a8e26e1a-a9ac-4d51-b09c-7951acd4b4b7
# ╠═88d38a68-1f8a-425a-bbae-90355a2213d0
# ╠═f2ba0e8a-4a9f-4b98-85b8-d54c71fd3616
# ╟─289d2c59-e920-47fe-b9ad-cb0a33ef0c9c
# ╠═1fcfbad3-2fad-4eee-a1a3-031dc29c9083
# ╠═25a183d9-c6a2-4ac3-a283-a02e4e9231dd
# ╠═6d1d8ae2-6a9e-48c1-a545-1f7354125bf0
# ╠═b43c5e74-5010-4870-a058-d3ad2c1ed548
# ╠═ee76e884-86e6-45f6-bbb2-c8e73daa5883
# ╟─0b6f33b9-41d4-48fd-8026-8a3bddcc1989
# ╠═a22a5421-05bf-484f-a2d3-91a06a0c6476
# ╠═85856abf-ee16-424a-ac06-97f76e32e444
# ╠═87f2b4c4-b163-4ae2-86b6-0266dff1da19
# ╠═970871ac-a5f2-4e30-9420-489acdbe79f9
# ╠═a90686f1-5f57-43e9-b22a-c8725775864f
# ╠═7afb1a46-2675-4b70-be39-5100fe2c2274
# ╠═b757aee6-49dc-4ead-abbf-2c60ddfad5a6
# ╠═576a3999-c408-44b3-851a-d02c9290a374
# ╠═4c1f6b31-ce09-4fba-b827-460e8a0d7e1a
# ╠═0e734e72-fc3a-48c7-b1d5-c0a380768eec
# ╠═fae68c38-be85-4718-8ee6-f900150e2b9a
# ╠═228c8672-03cc-4dc0-ba0c-74ee8928fb11
# ╠═e181c648-7f4e-473a-92ed-6fde8c177202
