sweepcomparedir(args...) = datadir("sweepcompare", args...)


# CSV generating sweeps
function ivsweep_csv(;
        specieslayout = GoldModel.SpeciesLayout(),
        reactiondata = GoldModel.ReactionData(),
        voltages = range(-1.25, reactiondata.ϕ_pzc, length = 31) * ufac"V",
        actcoeff = :DGML,
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


    modeldata = GoldModel.create_model(
        γ_select = string(actcoeff),
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


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 6: CSV export utilities.
# `filename` used the notebook globals `L`, `user_input_model`,
# `user_input_cv`, `ionsize`; here those are passed via keyword arguments.
# The export functions forward them through a `meta` NamedTuple, e.g.
#   meta = (L=L, BC_Select="Robin", mode="DMGL_γ", ionsize="All_species",
#           scanrate=0.05, nperiods=1, vmin=-1.2, vmax=0.8)
# `species` defaults to 6 (= iohminus index for the Gold CO2RR model).
# =====================================================================

# ── moved from cell a34cdba3 (filename)  [notebook globals → keyword args] ──
function filename(
        function_name;
        L, BC_Select, mode, ionsize, scanrate, nperiods, vmin, vmax
    )
    σ = round(L / ufac"μm")
    base_cvname = string(
        function_name, "_σ_", σ, "_", BC_Select, "_", mode,
        "_pnp_", ionsize, "_Scanrate_", scanrate,
        "_Periods_", nperiods, "_sweep_range_", vmin, "-", vmax, "_cv"
    )
    return base_cvname
end

# ── moved from cell 406fb8e5 (export_scanrate_varied_species_csv_long) ──
function export_scanrate_varied_species_csv_long(
        saws,
        scresults;
        meta,
        species = 6,  # iohminus index
        outdir::AbstractString = "../data/output",
        function_name::AbstractString = "scanrate_varied",
        current_scale = ufac"cm^2" / ufac"mA",
        scanrate_scale = 1.0,
        scanrate_get = saw -> saw.scanrate,
    )
    @assert length(saws) == length(scresults)

    scanrates = Float64[]
    voltages = Float64[]
    values = Float64[]

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
        Voltage = voltages,
        Value = values,
    )

    fname = filename(function_name; meta...)
    outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    return df
end

# ── moved from cell 9e44f14b (export_cv_profile_csv) ──
function export_cv_profile_csv(
        rec;
        meta,
        species = 6,  # iohminus index
        outdir::AbstractString = "../data/output",
        function_name::AbstractString = "cv_profile",
        current_scale = ufac"cm^2" / ufac"mA",
    )
    V = rec.voltages
    I = currents(rec, species) .* current_scale

    @assert length(V) == length(I)

    df = DataFrame(
        Voltage = V,
        Current = I,
    )

    fname = filename(function_name; meta...)
    outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    return df
end

# ── moved from cell e6dca43d (export_pressure_varied_species_csv_long) ──
function export_pressure_varied_species_csv_long(
        P_recs;
        meta,
        species = 6,  # iohminus index
        outdir::AbstractString = "../data/output",
        function_name::AbstractString = "pressure_varied",
        scale = ufac"cm^2" / ufac"mA",
    )
    pressures = Float64[]
    voltages = Float64[]
    values = Float64[]

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
        Voltage = voltages,
        Value = values,
    )

    fname = filename(function_name; meta...)
    outfile = joinpath(outdir, fname * ".csv")
    mkpath(dirname(outfile))
    CSV.write(outfile, df)

    return df
end
