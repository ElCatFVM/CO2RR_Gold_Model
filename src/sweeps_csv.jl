"Path into `data/sweepcompare`, where the exported sweep CSVs live."
sweepcomparedir(args...) = datadir("sweepcompare", args...)


"""
    ivsweep_csv(; specieslayout, reactiondata, voltages, actcoeff, hydrated, bcmodel, model, kwargs...)

Run an IV sweep and write voltage and current to a CSV under [`sweepcomparedir`](@ref).

The filename encodes the model options, so repeated runs with different `actcoeff` /
`hydrated` / `bcmodel` settings land in separate files — see [`filename`](@ref).

!!! warning "The current column is the OH⁻ flux"
    OH⁻ takes part in the buffer network, so its boundary flux is not purely faradaic, and
    it is also the species whose stoichiometry the electrode reaction re-routes. The
    plotting layer uses CO instead; this export has not followed.
"""
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
"""
    filename(function_name; L, BC_Select, mode, ionsize, scanrate, nperiods, vmin, vmax)

Result filename encoding the sweep parameters, so runs that differ only in model options
do not overwrite one another. `L` is converted to μm.
"""
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
"""
    export_scanrate_varied_species_csv_long(saws, scresults; meta, species, outdir, ...)

Long-format CSV of one species' concentration and current against scan rate: one row per
(scan rate, sample) rather than one column per scan rate.
"""
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

# ==========================================================================
# Publication CSVs.
#
# The four files behind the polarization/activity/concentration/capacitance
# figure of `plot_publication_notebook.jl`. They used to be written by cells of
# `notebooks/row_interaction_script.jl`, which is not part of the package, so
# the export half of that pipeline could not be reached from `using AuCO2RR`
# and the naming rule lived only in a notebook.
#
# The names are reproduced exactly, including their inconsistencies: the
# capacitance and polarization files carry a solver tag (`pb` / `pnp`) and the
# concentration file does not. Changing that would orphan every file already in
# `data/output`.
# ==========================================================================

"""
    publication_csv_name(kind; bc, mode, ionsize, solver = nothing, index = nothing)

Filename stem for a publication CSV, without the extension.

`kind` is `"DLCap"`, `"Polarization_Curve"`, `"Activity_Curve"` or `"Concentration"`;
`bc` is the boundary-condition tag (`"Robin"`), `mode` the activity model (`"Stefan_γ"`),
`ionsize` the radius set (`"Same_Size"`, `"All_species"`, `"Potassium_only"`). `solver`
inserts `pb` or `pnp`; `index` appends a molarity index.

    publication_csv_name("Polarization_Curve"; bc="Robin", mode="Stefan_γ",
                         ionsize="Same_Size", solver="pnp")
    # "Polarization_Curve_Robin_Stefan_γ_pnp_Same_Size"
"""
function publication_csv_name(
        kind::AbstractString;
        bc::AbstractString,
        mode::AbstractString,
        ionsize::AbstractString,
        solver = nothing,
        index = nothing,
    )
    parts = [kind, bc, mode]
    solver === nothing || push!(parts, string(solver))
    push!(parts, ionsize)
    index === nothing || push!(parts, string(index))
    return join(parts, "_")
end

"""
    export_dlcap_csv(results; bc, mode, ionsize, solver, outdir)

One CSV per molarity of a `capscalc` result: `Voltage`, `Capacitance`.

`results` is the vector `capscalc` returns; each entry supplies `voltage_range` and
`dlcaps`. Files are numbered from 1 in the order of that vector, which is the order of the
molarities passed to `capscalc` — the index in the filename carries no other meaning, so
record which molarity each one was.

Returns the paths written.
"""
function export_dlcap_csv(
        results;
        bc::AbstractString,
        mode::AbstractString,
        ionsize::AbstractString,
        solver::AbstractString = "pnp",
        outdir::AbstractString = joinpath("..", "data", "output"),
    )
    mkpath(outdir)
    base = publication_csv_name("DLCap"; bc, mode, ionsize, solver)
    return map(enumerate(results)) do (i, r)
        df = DataFrame(Voltage = r.voltage_range, Capacitance = r.dlcaps)
        path = joinpath(outdir, string(base, "_", i, ".csv"))
        CSV.write(path, df)
        path
    end
end

"""
    export_polarization_csv(ivresult; species, scale, bc, mode, ionsize, outdir)

Polarization curve as `Voltage`, `Current`.

`scale` defaults to 1, i.e. SI, because that is what the existing files hold; pass
`cm^2 / mA` for mA cm⁻².

!!! note "Species"
    Defaults to OH⁻, matching the files already in `data/output`. That is valid for a
    cathodic sweep — the OH⁻ and CO fluxes agree exactly there — but reports **zero** on an
    anodic branch, because the electrode reaction books the anodic half of the proton
    stoichiometry on H⁺. Pass `species = ico` with the current doubled if a sweep ever
    crosses into oxidation.
"""
function export_polarization_csv(
        ivresult;
        species = iohminus,
        scale = 1.0,
        bc::AbstractString,
        mode::AbstractString,
        ionsize::AbstractString,
        solver::AbstractString = "pnp",
        outdir::AbstractString = joinpath("..", "data", "output"),
    )
    mkpath(outdir)
    df = DataFrame(
        Voltage = ivresult.voltages,
        Current = currents(ivresult, species) .* scale,
    )
    path = joinpath(
        outdir,
        publication_csv_name("Polarization_Curve"; bc, mode, ionsize, solver) * ".csv",
    )
    CSV.write(path, df)
    return path
end

"""
    export_activity_csv(activity; bc, mode, ionsize, outdir)

Surface activities against potential: `Voltage` plus one column per species.

`activity` is the record `AuCO2RR_plots.activity_vs_voltage_axis` returns — its `vgrid`,
`activity_electrode` and `species` fields. Taking the species names from the record rather
than from a list written out at the call site is what keeps the columns aligned with the
rows of the matrix.
"""
function export_activity_csv(
        activity;
        bc::AbstractString,
        mode::AbstractString,
        ionsize::AbstractString,
        solver::AbstractString = "pnp",
        outdir::AbstractString = joinpath("..", "data", "output"),
    )
    path = joinpath(
        outdir,
        publication_csv_name("Activity_Curve"; bc, mode, ionsize, solver) * ".csv",
    )
    return _write_species_vs_voltage(
        path, activity.vgrid, activity.activity_electrode, activity.species, outdir
    )
end

"""
    export_concentration_csv(conc_out; bc, mode, ionsize, outdir)

Surface concentrations against potential: `Voltage` plus one column per species.

`conc_out` is the record `AuCO2RR_plots.conc_vs_voltage_axis` returns. Note the filename
carries **no** solver tag, unlike the other three — that is how the existing files are named.
"""
function export_concentration_csv(
        conc_out;
        bc::AbstractString,
        mode::AbstractString,
        ionsize::AbstractString,
        outdir::AbstractString = joinpath("..", "data", "output"),
    )
    path = joinpath(
        outdir, publication_csv_name("Concentration"; bc, mode, ionsize) * ".csv"
    )
    return _write_species_vs_voltage(
        path, conc_out.vgrid, conc_out.conc_electrode, conc_out.species, outdir
    )
end

"""
    _write_species_vs_voltage(path, vgrid, M, species, outdir)

Write a `Voltage` column plus one column per species from a `species × voltage` matrix.

The shape assertions are the point: a silently transposed matrix would write a file of the
wrong width with plausible-looking numbers, and nothing downstream would notice.
`bom = true` so the species names survive being opened in Excel.
"""
function _write_species_vs_voltage(path, vgrid, M, species, outdir)
    size(M, 1) == length(species) || error(
        "matrix has $(size(M, 1)) rows but $(length(species)) species names"
    )
    size(M, 2) == length(vgrid) || error(
        "matrix has $(size(M, 2)) columns but $(length(vgrid)) voltages"
    )
    mkpath(outdir)
    df = DataFrame(Voltage = collect(vgrid))
    for (i, name) in enumerate(species)
        df[!, String(name)] = vec(M[i, :])
    end
    CSV.write(path, df; bom = true)
    return path
end

# ── moved from cell 9e44f14b (export_cv_profile_csv) ──
"""
    export_cv_profile_csv(...)

CSV of one voltammogram's spatial profile — concentration against position at the stored
times.
"""
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
"""
    export_pressure_varied_species_csv_long(...)

Long-format CSV of one species' concentration and current against CO₂ partial pressure.
Column names come from [`_pressure_colname`](@ref).
"""
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
