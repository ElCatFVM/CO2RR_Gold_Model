module GoldModel
using LiquidElectrolytes
using VoronoiFVM
using ExtendableGrids
using LessUnitful
using CatmapInterface
using PreallocationTools
using Catalyst
using DelimitedFiles
using DrWatson

Base.@kwdef struct SpeciesLayout
    # species not involved in reactions
    ikplus::Int = 1
    # species involved in buffer reactions but not in surface reactions
    ibufferstart::Int = 2
    ihplus::Int = 2
    ihco3::Int = 3
    ico3::Int = 4
    # species involved in buffer reactions and surface reactions
    isurfacestart::Int = 5
    ico2::Int = 5
    iohminus::Int = 6
    ibufferend::Int = 6
    # species involved in surface reactions but not in buffer reactions
    ico::Int = 7
    nc::Int = 7
    ico_t::Int = 8
    icooh_t::Int = 9
    ico2_t::Int = 10
    isurfaceend::Int = 10
    na::Int = 3 # CO_t, CO2_t, COOH_t
    Γ_we::Int = 1
    Γ_bulk::Int = 2
end


Base.@kwdef struct ReactionData
    # bulk constants
    pH::Float64 = 6.8
    T::Float64 = 298.0 * ufac"K"
    Hcp_CO::Float64 = 9.7e-6 * ufac"mol / (m^3 * Pa)"
    Hcp_CO2::Float64 = 3.3e-4 * ufac"mol / (m^3 * Pa)"

    ## reaction rate constants for bulk reactions
    ### CO2 + OH- <=> HCO3-
    kbe1::Float64 = 4.44e7 / ufac"(mol / dm^3)"
    kbf1::Float64 = 5.93e3 / ufac"(mol / dm^3) / s"
    kbr1::Float64 = kbf1 / kbe1
    ### HCO3- + OH- <=> CO3-- + H2O
    kbe2::Float64 = 4.66e3 / ufac"(mol / dm^3)"
    kbf2::Float64 = 1.0e8 / ufac"(mol / dm^3) / s"
    kbr2::Float64 = kbf2 / kbe2
    ### CO2 + H20 <=> HCO3- + H+
    kae1::Float64 = 4.44e-7 * ufac"(mol / dm^3)"
    kaf1::Float64 = 3.7e-2 / ufac"s"
    kar1::Float64 = kaf1 / kae1
    ### HCO3- <=> CO3-- + H+
    kae2::Float64 = 4.66e-5 / ufac"(mol / dm^3)"
    kaf2::Float64 = 59.44e3 / ufac"(mol / dm^3) / s"
    kar2::Float64 = kaf2 / kae2
    ### autoprotolyse
    kwe::Float64 = 1.0e-14 * ufac"(mol / dm^3)^2"
    kwf::Float64 = 2.4e-5 * ufac"(mol / dm^3) / s"
    kwr::Float64 = kwf / kwe
    aH₂O::Float64 = 1.0 #* mol/dm^3
    c̄::Float64 = 55.508 * ufac"mol / dm^3"

    M0::Float64 = 18.0153 * ufac"g/mol"
    # surface nts
    S::Float64 = 9.61e-5 / ph"N_A" * (1.0e10)^2 * ufac"mol / m^2"
    C_gap::Float64 = 20 * ufac"μF / cm^2"
    ϕ_pzc::Float64 = 0.16 * ufac"V"
    v0::Float64 = ph"N_A" * (8.2 * 1.0e-1 * ufac"nm")^3 #18.048 * ufac"cm^3/mol"# # 1 / (55.4 * ufac"M") #
    xref::Float64=10* ufac"nm" # reference electrode position
end

function make_species_dict(specieslayout)
    return Dict(
        "K⁺" => specieslayout.ikplus,
        "HCO₃⁻" => specieslayout.ihco3,
        "CO₃²⁻" => specieslayout.ico3,
        "CO₂" => specieslayout.ico2,
        "OH⁻" => specieslayout.iohminus,
        "H⁺" => specieslayout.ihplus,
        "CO" => specieslayout.ico,
        "CO_t" => specieslayout.ico_t,
        "COOH_t" => specieslayout.icooh_t,
        "CO2_t" => specieslayout.ico2_t,
    )
end

function make_species_dict_catmap(specieslayout)
    return Dict(
        "OH_g" => specieslayout.iohminus,
        "CO2_aq" => specieslayout.ico2,
        "CO_aq" => specieslayout.ico,
        "CO_t" => specieslayout.ico_t,
        "COOH_t" => specieslayout.icooh_t,
        "CO2_t" => specieslayout.ico2_t,
    )
end

elydata_NaClO₄() = ElectrolyteData(
    z = [-1, 1],
    κ = [15.0, 25.0],
    c_bulk = [1.0, 1.0],
    ε = 26.0,
)

elydata_NaF() = ElectrolyteData(
    z = [-1, 1],
    κ = [25.0, 25.0],
    #v = [25.0, 25.0],
    c_bulk = [0.5, 0.5],
    #vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
    #vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
    v0 = 18.048 * ufac"cm^3" / ufac"mol",
)

elydata_Au(bulk, γ, specieslayout, reactiondata, ircompensation) = ElectrolyteData(;
    nc = size(bulk)[1],
    na = specieslayout.na,
    z = getproperty.(bulk, :z),
    D = getproperty.(bulk, :D),
    T = reactiondata.T,
    eneutral = false,
    κ = getproperty.(bulk, :κ),
    c_bulk = getproperty.(bulk, :c_bulk),
    v0 = reactiondata.v0,
    v = getproperty.(bulk, :v),
    M0 = reactiondata.M0,
    M = getproperty.(bulk, :M),
    Γ_we = specieslayout.Γ_we,
    Γ_bulk = specieslayout.Γ_bulk,
    C_gap = reactiondata.C_gap,
    ϕ_pzc = reactiondata.ϕ_pzc,
    ircompensation = ircompensation,
    xref=[reactiondata.xref],
    actcoeff! = γ,
)

function DGML_γ!(γ, c, p, electrolyte)
    return LiquidElectrolytes.DGML_gamma!(γ, c, p, electrolyte)
end

function Stefan_γ!(γ, c, p, electrolyte)
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp, c_bulk, v, nc) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
    for ic in cspecies
        γ[ic] = (1.0 / (1 - sum(c[i] * v[i] for i in 1:nc))) # (mol/dm^3)))
    end
    return γ
end

function Potassium_γ!(γ, c, p, electrolyte)
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp, c_bulk, v, nc) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
    γ .= 1 # ????
    ikplus = 1 # need to find a way to provide this otherwise.
    # Probably it will suffice to pass v=0 for all species except of K+
    γ[ikplus] = 1.0 / (1 - v[ikplus] * c[ikplus]) # / (mol/dm^3))
    return γ
end


Base.@kwdef struct BulkSpecies
    name::String = ""
    z::Int = 0
    D::Float64 = 1.0e-9 * ufac"m^2/s"
    c_bulk::Union{Nothing, Float64} = nothing
    κ::Float64 = 0
    a::Float64 = 1.0 # Å
    v::Float64 = v = ph"N_A" * (a * 1.0e-1 * ufac"nm")^3
    M::Float64 = ReactionData().M0 * v
    color::String = "#000000"

end

function make_eneutral(
        bulk_species::Vector{BulkSpecies}
        ; name, z, D, κ = 0.0, a = 8.2, v = ph"N_A" * (a * 1.0e-1 * ufac"nm")^3,
        #v0 * (κ * abs(z) + 1)
        M = ReactionData().M0 * v, color
    )
    a *= 1.0e-1 * ufac"nm"
    c_bulk = -mapreduce(x -> x.c_bulk * x.z, +, bulk_species) / z
    return BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
end

function buffer_system(reactiondata)
    @variables t
    @species HCO₃⁻(t), CO₃²⁻(t), CO₂(t), OH⁻(t), H⁺(t)
    @parameters γHCO₃⁻ γCO₃²⁻ γCO₂ γOH⁻ γH⁺

    (; kbf1, kbr1, kbf2, kbr2, kaf1, kar1, kaf2, kar2, kwf, kwr, aH₂O) = reactiondata
    buffer_rn = @reaction_network buffer begin
        ($kbf1 * γCO₂ * γOH⁻, $kbr1 * γHCO₃⁻), CO₂ + OH⁻ <--> HCO₃⁻
        ($kbf2 * γHCO₃⁻ * γOH⁻, $kbr2 * $aH₂O * γCO₃²⁻), HCO₃⁻ + OH⁻ <--> CO₃²⁻
        ($kaf1 * $aH₂O * γCO₂, $kar1 * γHCO₃⁻ * γH⁺), CO₂ <--> HCO₃⁻ + H⁺
        ($kaf2 * γHCO₃⁻, $kar2 * γCO₃²⁻ * γH⁺), HCO₃⁻ <--> CO₃²⁻ + H⁺
        ($kwf * $aH₂O, $kwr * γH⁺ * γOH⁻), ∅ <--> H⁺ + OH⁻
    end

    f_buffer! = CatmapInterface.generate_function(
        buffer_rn;
        dvs = [H⁺, HCO₃⁻, CO₃²⁻, CO₂, OH⁻],
        ps = [γH⁺, γHCO₃⁻, γCO₃²⁻, γCO₂, γOH⁻]
    )
    return f_buffer!
end


"""
A microkinetic modeling approach is taken:

The reaction mechanism for the \$CO_2\$ reduction is divided into four elementary reactions at the electrode surface:

1. Adsorption of \$CO_2\$ molecules at the oxygen atoms
\${CO_2}_{(aq)} + * \rightleftharpoons {CO_2 *}_{(ad)}\$

2. First proton-coupled electron transfer
\${CO_2*}_{(ad)} + H_2O_{(l)} + e^- \rightleftharpoons COOH*_{(ad)} + OH^-_{(aq)}\$

3. Second proton-coupled electron transfer
\$COOH*_{(ad)} + e^- \rightleftharpoons CO*_{(ad)} + OH^{-}_{(aq)}\$
with the transition state: \$*CO-OH^{TS}\$

4. Desorption of \$CO\$
\$*CO_{(ad)} \rightleftharpoons CO_{(aq)} + *\$
"""
function surface_reaction(specieslayout)
    catmap_params = CatmapInterface.parse_catmap_input(datadir("catmap_CO2R_data", "catmap_CO2R_template.mkm"))
    rn = create_reaction_network(catmap_params; symbolic_formation_energies = false)
    odesys = convert(ODESystem, rn; combinatoric_ratelaws = false)
    odesys = CatmapInterface.liquidize(odesys, catmap_params)
    vars = Catalyst.unknowns(odesys)
    species_dict_catmap = make_species_dict_catmap(specieslayout)
    f_microkinetics! = CatmapInterface.generate_function(
        odesys;
        dvs = sort(vars, by = x -> species_dict_catmap[string(operation(x))])
    )
    paramsidx = paramsmap(odesys)
    return f_microkinetics!, paramsidx, odesys

end


function calc_QBL_local(u, data; tolϕ = 1.0e-12)
    (; ip, iϕ, ε_0, pscale, ε) = data

    Δp = u[ip] * pscale + 0.1
    Δϕ = u[iϕ]

    arg = 2 * (1.0 .+ ε) * ε_0 * Δp
    arg = abs(arg)
    if arg < 0
        error("Invalid local QBL state: Δp=$(Δp), Δϕ=$(Δϕ)")
        arg = 0
    end

    Q = sign(Δϕ) * sqrt(arg)
    return Q
end


function create_model(;
        γ_select = "DGML",
        use_md_hydrated = true,
        BC_model = :Robin,
        model = :Gold,
        ircompensation = NoIRCompensation(),
        specieslayout = SpeciesLayout(),
        reactiondata = ReactionData()
    )

    if γ_select == "DGML"
        γ = DGML_γ!
    elseif γ_select == "Stefan"
        γ = Stefan_γ!
    elseif γ_select == "Potassium"
        γ = Potassium_γ!
    else
        error("undefined γ model: $(γ_select)")
    end


    if use_md_hydrated == true && γ == DGML_γ!
        println("hydrated_DGML")
        # --- hydrated radii [nm] (aqueous effective radii) ---
        a_HCO3 = 3.33   #  (Å)
        a_CO3 = 3.94   #  (Å)
        a_CO2 = 1.7   #  (Å) (often treated as vdW/effective in water)
        a_OH = 3.0   #  (Å)
        a_H = 2.8   #  (Å) (H3O+ effective hydrated)
        a_CO = 1.4   #  (Å)
        a_K = 3.31   #  (Å)

        # --- MD hydration numbers (1st shell) ---
        κ_HCO3 = 5.4     #  (MD)
        κ_CO3 = 8.5     #  (MD)
        κ_CO2 = 0.0     # neutral, typically treat as ~0 (no structured hydration number in this model)
        κ_OH = 3.5     #  (MD)
        κ_H = 4.0     #  (MD)
        κ_CO = 0.0     # neutral
        κ_K = 6.0     #   (MD)
    elseif use_md_hydrated == true && γ == Stefan_γ!
        println("hydrated_Stefan")
        # Stefan's Model / Consider all effective size
        a_HCO3 = 8.5
        a_CO3 = 9.9
        a_CO2 = 3.4
        a_OH = 7.6
        a_H = 7.3
        a_CO = 2.8
        a_K = 8.2

        κ_HCO3 = 0
        κ_CO3 = 0
        κ_CO2 = 0
        κ_OH = 0
        κ_H = 0
        κ_CO = 0
        κ_K = 0
    elseif use_md_hydrated == false && γ == DGML_γ!
        println("nonhydrated_DGML")
        # --- hydrated radii [nm] (aqueous effective radii) ---
        a_HCO3 = 3.33   #  (Å)
        a_CO3 = 3.94   #  (Å)
        a_CO2 = 1.7   #  (Å) (often treated as vdW/effective in water)
        a_OH = 3.0   #  (Å)
        a_H = 2.8   #  (Å) (H3O+ effective hydrated)
        a_CO = 1.4   #  (Å)
        a_K = 3.31   #  (Å)

        # a_HCO3 = 8.5
        # a_CO3 = 9.9
        # a_CO2 = 3.4
        # a_OH = 7.6
        # a_H = 7.3
        # a_CO = 2.8
        # a_K = 8.2

        κ_HCO3 = 0
        κ_CO3 = 0
        κ_CO2 = 0
        κ_OH = 0
        κ_H = 0
        κ_CO = 0
        κ_K = 0

    elseif use_md_hydrated == false && γ == Potassium_γ!
        println("Potassium")
        fac = 0.8
        a_HCO3 = 0;  κ_HCO3 = 0
        a_CO3 = 0;  κ_CO3 = 0
        a_CO2 = 0;  κ_CO2 = 0
        a_OH = 0;  κ_OH = 0
        a_H = 0;  κ_H = 0
        a_CO = 0;  κ_CO = 0
        a_K = 8.2;  κ_K = 0
    elseif use_md_hydrated == false && γ == Stefan_γ!
        println("Same_size_MPB_Stefan")
        a_HCO3 = 8.2      ;  κ_HCO3 = 0
        a_CO3 = 8.2      ;  κ_CO3 = 0
        a_CO2 = 8.2      ;  κ_CO2 = 0
        a_OH = 8.2      ;  κ_OH = 0
        a_H = 8.2      ;  κ_H = 0
        a_CO = 8.2      ;  κ_CO = 0
        a_K = 8.2      ;  κ_K = 0
    else
        error("undefined case:  use_md_hydrated = $(use_md_hydrated), γ=$(γ)")
    end

    bulk = [
        BulkSpecies(;
            name = "HCO₃⁻",
            z = -1,
            D = 1.185e-9,
            c_bulk = 0.091 * ufac"mol / dm^3",
            a = a_HCO3,
            κ = κ_HCO3,
            color = "#7B5C3E"
        ),
        BulkSpecies(;
            name = "CO₃²⁻",
            z = -2,
            D = 0.923e-9,
            c_bulk = 2.68e-6 * ufac"mol / dm^3",
            a = a_CO3,
            κ = κ_CO3,
            color = "#222222"
        ),
        BulkSpecies(;
            name = "CO₂",
            z = 0,
            D = 1.91e-9,
            c_bulk = 0.033 * ufac"mol / dm^3",
            a = a_CO2,
            κ = κ_CO2,
            color = "#C0392B"
        ),
        BulkSpecies(;
            name = "OH⁻",
            z = -1,
            D = 5.273e-9,
            c_bulk = 10^(reactiondata.pH - 14) * ufac"mol / dm^3",
            a = a_OH,
            κ = κ_OH,
            color = "#27AE60"
        ),
        BulkSpecies(;
            name = "H⁺",
            z = 1,
            D = 9.31e-9,
            c_bulk = 10^(-reactiondata.pH) * ufac"mol / dm^3",
            a = a_H,
            κ = κ_H,
            color = "#888888"
        ),
        BulkSpecies(;
            name = "CO",
            z = 0,
            D = 2.23e-9,
            c_bulk = 0.0 * ufac"mol / dm^3",
            a = a_CO,
            κ = κ_CO,
            color = "#2980B9"
        ),
    ]

    push!(
        bulk, make_eneutral(
            bulk; name = "K⁺",
            z = 1,
            D = 1.957e-9,
            a = a_K,
            κ = κ_K,
            color = "#E07B39"
        )
    )
    species_dict = make_species_dict(specieslayout)

    bulk = sort(bulk, by = x -> species_dict[x.name])
    bulknames = [ b.name for b in bulk]
    bulkcolors = [ b.color for b in bulk]
    @show bulknames
    @show bulkcolors

    γ_cache = DiffCache(zeros(specieslayout.nc), 14)

    f_buffer! = buffer_system(reactiondata)

    function reaction(
            f,
            u,
            node,
            data
        )

        (; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, pscale, p_bulk, actcoeff!) = data


        γ = get_tmp(γ_cache, u[1])
        # compute activity coefficients according to the approach in Ringe et al.
        actcoeff!(γ, u, u[ip], data)
        (; ibufferstart, ibufferend) = specieslayout
        @views f_buffer!(
            f[ibufferstart:ibufferend],
            u[ibufferstart:ibufferend],
            γ[ibufferstart:ibufferend],
            nothing
        )
        return nothing
    end

    ps_cache = DiffCache(zeros(18), 14)
    us_cache = DiffCache(zeros(specieslayout.isurfaceend - specieslayout.isurfacestart + 1), 14)
    f_microkinetics!, paramsidx, odesys = surface_reaction(specieslayout)

    function we_breactions(
            f,
            u,
            bnode,
            data
        )
        (; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we, iϕ_we, ε, cspecies, actcoeff!) = data

        γ = get_tmp(γ_cache, u[1])
        actcoeff!(γ, u, u[ip], data)

        γ_co2 = γ[specieslayout.ico2]
        γ_co = γ[specieslayout.ico]

        if isactive(data.ircompensation)
            ϕ_we_set = u[iϕ_we]
        else
            ϕ_we_set = ϕ_we
        end
        #ρ = F .* sum(z[k] .* u[k] for k in 1:cspecies)

        if BC_model == :Robin
            σ = reactiondata.C_gap * (ϕ_we_set - reactiondata.ϕ_pzc - u[data.iϕ])
        else
            σ = calc_QBL_local(u, data)
        end


        #        σ = surface_charge(u, data, user_input_model.BC_Select)
        local_pH = -log10(u[specieslayout.ihplus] * γ[specieslayout.ihplus] / ufac"(mol / dm^3)")

        #a_ohminus = u[specieslayout.iohminus] * γ[specieslayout.iohminus] / ufac"(mol / dm^3)"
        #local_pH  = 14 + log10(a_ohminus)


        #for (p, default_value) in odesys.defaults
        #	ps[paramsidx[p]] = default_value
        #	println(default_value)
        #end

        ps = get_tmp(ps_cache, u[iϕ])
        ps[paramsidx[Symbolics.rename(odesys.σ, :σ)]] = σ
        ps[paramsidx[Symbolics.rename(odesys.γCO2_aq, :γCO2_aq)]] = γ_co2
        ps[paramsidx[Symbolics.rename(odesys.aH2O_g, :aH2O_g)]] = reactiondata.aH₂O
        ps[paramsidx[Symbolics.rename(odesys.ϕ, :ϕ)]] = u[iϕ]
        ps[paramsidx[Symbolics.rename(odesys.ϕ_we, :ϕ_we)]] = ϕ_we_set
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
            f[specieslayout.isurfacestart:specieslayout.isurfaceend],
            u[specieslayout.isurfacestart:specieslayout.isurfaceend],
            ps,
            nothing
        )
        #if ϕ_we == -0.9
        #	println(u[isurfacestart:isurfaceend])
        #end

        #            global coeff = u[isurfacestart:isurfaceend]
        # conversion from turnover frequency (appropriate for change in coverage) to
        # production rate (per unit area) (approprite for change in concentration)
        # by S = number of free catalyst sites in mole per unit area
        S = reactiondata.S
        f[specieslayout.ico2] *= S
        f[specieslayout.ico] *= S
        #f[specieslayout.iohminus] *= S
        # f is a sink. r_oh > 0 → OH⁻ consumed (anodic); r_oh < 0 → OH⁻ produced (cathodic).
        # Neither H⁺ (1.6e-7 M) nor OH⁻ (6.3e-8 M) can sustain the ~1e-4 mol/m²/s proton
        # turnover of this step — H₂O (55.5 M) is the actual reservoir. So always book the
        # stoichiometry on the species being *produced*; the dilute ions re-equilibrate
        # afterwards through the buffer network. Branch-free (max/min instead of `if`) so the
        # sparsity tracer can walk both paths.
        r_oh  = f[specieslayout.iohminus] * S
        r_pos = max(r_oh, zero(r_oh))     # anodic  part → H⁺ produced
        r_neg = min(r_oh, zero(r_oh))     # cathodic part → OH⁻ produced

        f[specieslayout.iohminus] = r_neg
        f[specieslayout.ihplus]  -= r_pos
        
        return
    end


    function pb_bcondition(f, u, bnode, data)
        (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data

        if BC_model == :Dirichlet
            ## Dirichlet ϕ=ϕ_we at Γ_we
            boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
            #boundary_dirichlet!(f, u, bnode; species = ip, region = Γ_we, value = p_we)
        elseif BC_model == :Robin
            ## Robin ϕ=dϕ₀/dx
            boundary_robin!(f, u, bnode, iϕ, Γ_we, reactiondata.C_gap, reactiondata.C_gap * (ϕ_we - reactiondata.ϕ_pzc))
        else
            ## neumann ϕ=dϕ₀/dx
            #boundary_neumann!(f, u, bnode, species = iϕ, region = Γ_we, value = C_gap * (ϕ_we - ϕ_pzc))

        end

        return bulkbcondition(f, u, bnode, data)
    end

    function pnp_bcondition(
            f,
            u,
            bnode,
            data
        )
        (; Γ_we, Γ_bulk, ϕ_we, iϕ, ϕ_bulk, ip, p_bulk, c_bulk, cspecies, ircompensation) = data

        if bnode.region == Γ_we
            # --- Working-electrode potential BC ---
            # For :pseudopotentiostat / :ohmicdrop, LiquidElectrolytes' generic
            if !isactive(ircompensation)
                if BC_model == :Dirichlet
                    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = (ϕ_we - reactiondata.ϕ_pzc))
                elseif BC_model == :Robin
                    potentialbcondition!(f, u, bnode, data, ϕ_we)
                end
            end

            # --- Surface (faradaic) reactions ---
            if !isa(ircompensation,OhmicDropEstimation) && model == :Gold
                we_breactions(f, u, bnode, data)
            end
        end

        return bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    end
if model == :Gold
    if isa(ircompensation, OhmicDropEstimation)
        ircompensation=copy(ircompensation; redoxreaction=we_breactions)
    end
    elydata = elydata_Au(bulk, γ, specieslayout, reactiondata, ircompensation)
    elseif model == :Landstorfer_NaClO₄
        elydata = elydata_NaClO₄
    elseif model == :Landstorfer_NaF
        elydata = elydata_NaF
    end
return (bcondition = pnp_bcondition, reaction = reaction, elydata = elydata, bulk = bulk, bulknames = bulknames, bulkcolors = bulkcolors, species_dict = species_dict)
end

end
