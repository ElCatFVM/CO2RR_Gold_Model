using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

module sweeps
    using AuCO2RR
    using AuCO2RR: AuCO2RR_plots
    using LessUnitful
    using LiquidElectrolytes
    using ExtendableGrids
    using VoronoiFVM
    using CatmapInterface
    using Catalyst
    using DrWatson
    using PreallocationTools
    using DataFrames, CSV

    sweepcomparedir(args...) = datadir("sweepcompare", args...)

    @unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å nm mm;
    @phconstants N_A c_0 k_B e h ε_0 R
    const F = N_A * e

    const Γ_we = 1
    const Γ_bulk = 2

    # bulk constants
    const pH = 6.8
    const T = 298.0 * K
    const Hcp_CO = 9.7e-6 * mol / (m^3 * Pa)
    const Hcp_CO2 = 3.3e-4 * mol / (m^3 * Pa)
    # species not involved in reactions
    const ikplus = 1
    # species involved in buffer reactions but not in surface reactions
    const ibufferstart = 2
    const ihplus = 2
    const ihco3 = 3
    const ico3 = 4
    # species involved in buffer reactions and surface reactions
    const isurfacestart = 5
    const ico2 = 5
    const iohminus = 6
    const ibufferend = 6
    # species involved in surface reactions but not in buffer reactions
    const ico = 7
    const nc = 7
    ## reaction rate constants for bulk reactions
    ### CO2 + OH- <=> HCO3-
    const kbe1 = 4.44e7 / (mol / dm^3)
    const kbf1 = 5.93e3 / (mol / dm^3) / s
    const kbr1 = kbf1 / kbe1
    ### HCO3- + OH- <=> CO3-- + H2O
    const kbe2 = 4.66e3 / (mol / dm^3)
    const kbf2 = 1.0e8 / (mol / dm^3) / s
    const kbr2 = kbf2 / kbe2
    ### CO2 + H20 <=> HCO3- + H+
    const kae1 = 4.44e-7 * (mol / dm^3)
    const kaf1 = 3.7e-2 / s
    const kar1 = kaf1 / kae1
    ### HCO3- <=> CO3-- + H+
    const kae2 = 4.66e-5 / (mol / dm^3)
    const kaf2 = 59.44e3 / (mol / dm^3) / s
    const kar2 = kaf2 / kae2
    ### autoprotolyse
    const kwe = 1.0e-14 * (mol / dm^3)^2
    const kwf = 2.4e-5 * (mol / dm^3) / s
    const kwr = kwf / kwe
    const aH₂O = 1.0 #* mol/dm^3
    const c̄ = 55.508mol / dm^3

    # surface constants
    const S = 9.61e-5 / N_A * (1.0e10)^2 * mol / m^2
    const C_gap = 20 * μF / cm^2
    const ϕ_pzc = 0.16 * V
    const ico_t = 8
    const icooh_t = 9
    const ico2_t = 10
    const isurfaceend = 10
    const na = 3 # CO_t, CO2_t, COOH_t
    const M0 = 18.0153 * ufac"g/mol"
    const v0 = N_A * (8.2 * Å)^3 #18.048 * ufac"cm^3/mol"# # 1 / (55.4 * ufac"M") #
    const c̄ = 55.508mol / dm^3


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

    const elydata_NaClO₄ = ElectrolyteData(
        z = [-1, 1],
        κ = [15.0, 25.0],
        c_bulk = [1.0, 1.0],
        ε = 26.0,
    )

    elydata_NaF = ElectrolyteData(
        z = [-1, 1],
        κ = [25.0, 25.0],
        #v = [25.0, 25.0],
        c_bulk = [0.5, 0.5],
        #vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
        #vrel = [1.7973*10e-5*46, 1.7973*10e-5*46]
        v0 = 18.048 * ufac"cm^3" / ufac"mol",
    )

    elydata_Au(bulk, γ) = ElectrolyteData(;
        nc = size(bulk)[1],
        na = na,
        z = getproperty.(bulk, :z),
        D = getproperty.(bulk, :D),
        T = T,
        eneutral = false,
        κ = getproperty.(bulk, :κ),
        c_bulk = getproperty.(bulk, :c_bulk),
        v0 = v0,
        v = getproperty.(bulk, :v),
        M0 = M0,
        M = getproperty.(bulk, :M),
        Γ_we = Γ_we,
        Γ_bulk = Γ_bulk,
        actcoeff! = γ,
    )


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

    function BulkSpecies(; name, z, c_bulk = nothing, a, D, κ = nothing, color)
        D *= m^2 / s
        c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol / dm^3 #* 0.0001
        # (v0/N_A)^(1/3)
        v = N_A * (a * Å)^3 #v0 * (κ + 1)
        a *= Å #8.2 * Å#8.2 * Å
        M = M0 * v
        return BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
    end

    function make_eneutral(
            bulk_species::Vector{BulkSpecies}
            ; name, z, D, κ = 0.0, a = 8.2, v = N_A * (a * Å)^3,
            #v0 * (κ * abs(z) + 1)
            M = M0 * v, color
        )
        a *= Å
        c_bulk = -mapreduce(x -> x.c_bulk * x.z, +, bulk_species) / z
        return BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
    end

    function buffer_system()
        @variables t
        @species HCO₃⁻(t), CO₃²⁻(t), CO₂(t), OH⁻(t), H⁺(t)
        @parameters γHCO₃⁻ γCO₃²⁻ γCO₂ γOH⁻ γH⁺
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
    function surface_reaction()
        catmap_params = CatmapInterface.parse_catmap_input(datadir("catmap_CO2R_data", "catmap_CO2R_template.mkm"))
        rn = create_reaction_network(catmap_params; symbolic_formation_energies = false)
        odesys = convert(ODESystem, rn; combinatoric_ratelaws = false)
        odesys = CatmapInterface.liquidize(odesys, catmap_params)
        vars = Catalyst.unknowns(odesys)
        f_microkinetics! = CatmapInterface.generate_function(
            odesys;
            dvs = sort(vars, by = x -> species_dict_catmap[string(operation(x))])
        )
        paramsidx = paramsmap(odesys)
        return f_microkinetics!, paramsidx, odesys

    end

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
        γ[ikplus] = 1.0 / (1 - v[ikplus] * c[ikplus]) # / (mol/dm^3))
        return γ
    end


    function calc_QBL_local(u, data; tolϕ = 1.0e-12)
        (; ip, iϕ, ε_0, pscale, ε) = data

        Δp = u[ip]
        Δϕ = u[iϕ]

        arg = 2 * (1.0 .+ ε) * ε_0 * Δp * pscale

        if arg < 0
            error("Invalid local QBL state: Δp=$(Δp), Δϕ=$(Δϕ)")
        end

        Q = sign(Δϕ) * sqrt(arg)
        return Q
    end

    function create_model(;
            γ = DGML_γ!,
            use_md_hydrated = true,
            BC_model = :Robin,
            model = :Gold
        )

        if use_md_hydrated == true && γ == "DMGL_γ"
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
        elseif use_md_hydrated == true && γ == "Stefan_γ"
            # Stefan's Model / Consider all effective size
            a_HCO3 = 8.5;  κ_HCO3 = 0
            a_CO3 = 9.9;  κ_CO3 = 0
            a_CO2 = 3.4;  κ_CO2 = 0
            a_OH = 7.6;  κ_OH = 0
            a_H = 7.3;  κ_H = 0
            a_CO = 2.8;  κ_CO = 0
            a_K = 8.2;  κ_K = 0
        elseif use_md_hydrated == false && γ == "DMGL_γ"
            # Stefan's Model / Only Potaissium has a size
            a_HCO3 = 0;  κ_HCO3 = 0
            a_CO3 = 0;  κ_CO3 = 0
            a_CO2 = 0;  κ_CO2 = 0
            a_OH = 0;  κ_OH = 0
            a_H = 0;  κ_H = 0
            a_CO = 0;  κ_CO = 0
            #a_K    = 3.31;  κ_K    = 6.0
            a_K = 8.2;  κ_K = 0.0 #5.0
        else
            a_HCO3 = 0;  κ_HCO3 = 0
            a_CO3 = 0;  κ_CO3 = 0
            a_CO2 = 0;  κ_CO2 = 0
            a_OH = 0;  κ_OH = 0
            a_H = 0;  κ_H = 0
            a_CO = 0;  κ_CO = 0
            a_K = 8.2;  κ_K = 0
        end

        bulk = [
            BulkSpecies(;
                name = "HCO₃⁻",
                z = -1,
                D = 1.185e-9,
                c_bulk = 0.091,
                a = a_HCO3,
                κ = κ_HCO3,
                color = :brown
            ),
            BulkSpecies(;
                name = "CO₃²⁻",
                z = -2,
                D = 0.923e-9,
                c_bulk = 2.68e-6,
                a = a_CO3,
                κ = κ_CO3,
                color = :violet
            ),
            BulkSpecies(;
                name = "CO₂",
                z = 0,
                D = 1.91e-9,
                c_bulk = 0.033,
                a = a_CO2,
                κ = κ_CO2,
                color = :red
            ),
            BulkSpecies(;
                name = "OH⁻",
                z = -1,
                D = 5.273e-9,
                c_bulk = 10^(pH - 14),
                a = a_OH,
                κ = κ_OH,
                color = :green
            ),
            BulkSpecies(;
                name = "H⁺",
                z = 1,
                D = 9.31e-9,
                c_bulk = 10^(-pH),
                a = a_H,
                κ = κ_H,
                color = :gray
            ),
            BulkSpecies(;
                name = "CO",
                z = 0,
                D = 2.23e-9,
                c_bulk = 0.0,
                a = a_CO,
                κ = κ_CO,
                color = :blue
            ),
        ]

        push!(
            bulk, make_eneutral(
                bulk; name = "K⁺",
                z = 1,
                D = 1.957e-9,
                a = a_K,
                κ = κ_K,
                color = :orange
            )
        )
        bulk = sort(bulk, by = x -> species_dict[x.name])
        bulknames = [ b.name for b in bulk]
        bulkcolors = [ b.color for b in bulk]
        @show bulknames
        @show bulkcolors

        if model == :Gold
            elydata = elydata_Au(bulk, γ)
        elseif model == :Landstorfer_NaClO₄
            elydata = elydata_NaClO₄
        elseif model == :Landstorfer_NaF
            elydata = elydata_NaF
        end

        γ_cache = DiffCache(zeros(nc), 14)

        f_buffer! = buffer_system()

        function reaction(
                f,
                u,
                node,
                data
            )

            (; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, pscale, p_bulk, actcoeff!) = data


            γ = get_tmp(γ_cache, u[ico2])
            # compute activity coefficients according to the approach in Ringe et al.
            actcoeff!(γ, u, u[ip], data)
            @views f_buffer!(
                f[ibufferstart:ibufferend],
                u[ibufferstart:ibufferend],
                γ[ibufferstart:ibufferend],
                nothing
            )
            return nothing
        end

        ps_cache = DiffCache(zeros(18), 14)
        us_cache = DiffCache(zeros(isurfaceend - isurfacestart + 1), 14)
        f_microkinetics!, paramsidx, odesys = surface_reaction()

        function we_breactions(
                f,
                u,
                bnode,
                data
            )
            (; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we, ε, cspecies, actcoeff!) = data

            γ = get_tmp(γ_cache, u[ico2])
            actcoeff!(γ, u, u[ip], data)

            γ_co2 = γ[ico2]
            γ_co = γ[ico]

            #ρ = F .* sum(z[k] .* u[k] for k in 1:cspecies)

            if BC_model == :Robin
                σ = C_gap * (data.ϕ_we - ϕ_pzc - u[data.iϕ])
            else
                σ = calc_QBL_local(u, data)
            end


            #        σ = surface_charge(u, data, user_input_model.BC_Select)
            local_pH = -log10(u[ihplus] * γ[ihplus] / (mol / dm^3))


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
            return f[ikplus] *= S
        end


        function pb_bcondition(f, u, bnode, data)
            (; Γ_we, Γ_bulk, ϕ_we, iϕ, ip) = data

            if BC_model == :Dirichlet
                ## Dirichlet ϕ=ϕ_we at Γ_we
                boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
                #boundary_dirichlet!(f, u, bnode; species = ip, region = Γ_we, value = p_we)
            elseif BC_model == :Robin
                ## Robin ϕ=dϕ₀/dx
                boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))
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
            (; Γ_we, Γ_bulk, ϕ_we, iϕ, ϕ_bulk, ip, p_bulk, c_bulk, cspecies) = data

            if BC_model == :Dirichlet
                boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = (ϕ_we - ϕ_pzc))
            elseif BC_model == :Robin
                boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))
            end


            if bnode.region == Γ_we && model == :Gold
                we_breactions(f, u, bnode, data)
            end


            return bulkbcondition(f, u, bnode, data; region = Γ_bulk)
        end
        return (bcondition = pnp_bcondition, reaction = reaction, elydata = elydata, bulknames = bulknames)
    end

    function sweep_IV(;
        voltages = (-1.5:0.1:0.0) * V,
        actcoeff = :DGML_γ!,
        hydrated = true,
        bcmodel = :Robin,
        model = :Gold,
        kwargs...
    )
        solver_control = (;
            max_round = 4,
            maxiters = 20,
            tol_round = 1.0e-9,
            verbose = "",
            reltol = 1.0e-8,
            tol_mono = 1.0e-10,
            damp_initial = 0.1,
            damp_growth = 2,
        )

        kwargs = merge(solver_control, kwargs)
        L = 80 * μm
        hmin = 1.0e-6 * μm
        hmax = 1.0 * μm
        X = ExtendableGrids.geomspace(0, L, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)

        γ = getproperty(sweeps, actcoeff)

        modeldata = create_model(
            γ = γ,
            use_md_hydrated = hydrated,
            BC_model = bcmodel,
            model = model
        )

        cell = PNPSystem(
            grid;
            modeldata.bcondition,
            reaction = modeldata.reaction,
            celldata = deepcopy(modeldata.elydata)
        )

        ivresult = ivsweep(cell; voltages, store_solutions = true, kwargs...)

        electrolyte = modeldata.elydata
        cspecies = electrolyte.cspecies
        ip = LiquidElectrolytes.pressure_index(electrolyte)
        nspecies = length(cspecies)
        nv = length(ivresult.voltages)
        tsol = LiquidElectrolytes.voltages_solutions(ivresult)
        ielectrode = 1

        scale = 1.0 / (mol / dm^3)

        γ_e = fill(NaN, nspecies, nv)
        a_e = fill(NaN, nspecies, nv)
        c_e = fill(NaN, nspecies, nv)

        model_type = actcoeff == :Stefan_γ! ? "Stefan_γ" : "DMGL_γ"

        for j in 1:nv
            sol = tsol.u[j]  

            unode = view(sol, :, ielectrode)
            pnode = unode[ip]

            v0 = electrolyte.v0
            bar_c = 1.0 / v0
            RT = electrolyte.RT

            Phi = sum(unode[ic] * electrolyte.v[ic] for ic in cspecies)
            solvent_frac = max(1.0 - Phi, eps(Float64))

            for (i, ic) in enumerate(cspecies)
                c = unode[ic]
                v_a = electrolyte.v[ic]
                size_ratio = v_a / v0

                term_conc = c / bar_c

                if model_type == "DMGL_γ"
                    term_press = exp((1.0 - size_ratio) * pnode / (bar_c * RT))
                    term_steric = solvent_frac^(-size_ratio)
                    a_eff_thermo = term_conc * term_press * term_steric

                elseif model_type == "Stefan_γ"
                    a_eff_thermo = term_conc * solvent_frac^(-1.0)

                else
                    error("Invalid model_type. Use 'DMGL_γ' or 'Stefan_γ'.")
                end

                c_scale = c * scale
                a_eff_scaled = a_eff_thermo * (bar_c * scale)

                c_e[i, j] = c_scale
                γ_e[i, j] = a_eff_scaled / max(c_scale, eps(Float64))
                a_e[i, j] = a_eff_scaled
            end
        end

        df = DataFrame(
            Voltage = ivresult.voltages,
            Current = currents(ivresult, iohminus)
        )

        specnames = modeldata.bulknames
        for i in 1:length(specnames)
            df[!, "a_" * specnames[i]] = vec(a_e[i, :])
            df[!, "c_" * specnames[i]] = vec(c_e[i, :])
            df[!, "γ_" * specnames[i]] = vec(γ_e[i, :])
        end
        @show df[!, "c_CO"]
        params = (actcoeff = actcoeff, hydrated = hydrated, bcmodel = bcmodel, model = model)
        fname = savename("sweep_iv", params, "csv")
        CSV.write(sweepcomparedir(fname), df)

        return cell, ivresult
    end


    function main()
        sweep_IV(actcoeff = :DGML_γ!, hydrated = true, bcmodel = :Robin, model = :Gold)
        sweep_IV(actcoeff = :Stefan_γ!, hydrated = true, bcmodel = :Robin, model = :Gold)
        sweep_IV(actcoeff = :Potassium_γ!, hydrated = true, bcmodel = :Robin, model = :Gold)
        return nothing
    end

end


if !isinteractive()
    sweeps.main()
end
