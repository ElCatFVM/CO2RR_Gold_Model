sweepcomparedir(args...) = datadir("sweepcompare", args...)


# CSV generating sweeps
function ivsweep_csv(;
        specieslayout = GoldModel.SpeciesLayout(),
        reactiondata = GoldModel.ReactionData(),
        voltages = range(-1.25, reactiondata.ϕ_pzc, length = 31) * ufac"V",
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
        damp_initial = 0.5,
        damp_growth = 2,
    )

    kwargs = merge(solver_control, kwargs)
    L = 80 * ufac"μm"
    hmin = 1.0e-6 * ufac"μm"
    hmax = 1.0 * ufac"μm"
    X = ExtendableGrids.geomspace(0, L, hmin, hmax)
    grid = ExtendableGrids.simplexgrid(X)

    γ = getproperty(GoldModel, actcoeff)

    modeldata = GoldModel.create_model(
        γ = γ,
        use_md_hydrated = hydrated,
        BC_model = bcmodel,
        model = model
    )

    cell = PNPSystem(grid; modeldata.bcondition, reaction = modeldata.reaction, celldata = deepcopy(modeldata.elydata))
    ivresult = ivsweep(cell; voltages, pzc = reactiondata.ϕ_pzc, store_solutions = true, kwargs...)

    electrolyte = modeldata.elydata
    cspecies = electrolyte.cspecies
    ip = LiquidElectrolytes.pressure_index(electrolyte)
    iϕ = LiquidElectrolytes.voltage_index(electrolyte)
    nspecies = length(cspecies)
    nv = length(ivresult.voltages)
    tsol = LiquidElectrolytes.voltages_solutions(ivresult)
    ielectrode = 1

    c_ref = 1.0 * ufac"mol/dm^3"
    γ_e = zeros(nspecies, nv)
    a_e = zeros(nspecies, nv)
    c_e = zeros(nspecies, nv)
    ϕ = zeros(nv)
    p = zeros(nv)
    for iv in 1:nv
        sol = tsol.u[iv]
        unode = view(sol, :, ielectrode)
        pnode = unode[ip]
        ϕnode = unode[iϕ]
        γ = zeros(size(sol, 1))
        electrolyte.actcoeff!(γ, unode, pnode, electrolyte)
        for ic in cspecies
            cval = unode[ic]
            γval = γ[ic]
            aval = γval * cval

            c_e[ic, iv] = cval / c_ref
            γ_e[ic, iv] = γval
            a_e[ic, iv] = aval
        end
        ϕ[iv] = ϕnode
        p[iv] = pnode
    end

    df = DataFrame(
        Voltage = ivresult.voltages,
        Current = currents(ivresult, specieslayout.iohminus)
    )

    specnames = modeldata.bulknames
    df[!, "ϕ"] = ϕ
    df[!, "p"] = p
    for i in 1:length(specnames)
        df[!, "a_" * specnames[i]] = vec(a_e[i, :])
        df[!, "c_" * specnames[i]] = vec(c_e[i, :])
        df[!, "γ_" * specnames[i]] = vec(γ_e[i, :])
    end
    params = (actcoeff = actcoeff, hydrated = hydrated, bcmodel = bcmodel, model = model)
    fname = savename("sweep_iv", params, "csv")
    CSV.write(sweepcomparedir(fname), df)
    return cell, ivresult
end;
