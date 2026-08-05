"""
sweep_over_L_cv(elydata, L_values, bcond, sawtooth, reaction; ...)

Sweep CV over boundary-layer thickness L (μm).
For each L: build grid -> PNPSystem -> run cvsweep with a fresh deepcopy of `sawtooth`.

Returns Dict{Float64, Any} mapping L -> CVSweepResult (or whatever cvsweep returns).
"""
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
        hmin = 1.0e-6 * μm
        hmax = 1.0 * μm
        X = ExtendableGrids.geomspace(0, L, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)

        celldata = deepcopy(elydata)
        celldata.eneutral = eneutral

        pnpcell = PNPSystem(grid; bcondition = bcond, celldata = celldata, reaction = reaction)

        saw = SawTooth(
            scanrate = sawtooth.scanrate,
            vmin = sawtooth.vmin,
            vmax = sawtooth.vmax,
            scanup = sawtooth.scanup
        )

        results[L] = cvsweep(
            pnpcell;
            voltages = saw,
            nperiods = nperiods,
            store_solutions = store_solutions,
            sweep_kwargs...
        )
    end

    return results
end


"""
    scanrate_varied_sweep(elydata, sawtooth, grid, bcond, reaction; scanrates, ...)

CV at each scan rate in `scanrates`. Returns `(scanrates, saws, sweeps)`.

!!! warning "Shadowed"
    `iv.jl` defines this name too and is included later, so **that** definition is the one
    in scope. This one also rebuilds the `SawTooth` without `vstart`/`tstart` and does not
    forward `sweep_kwargs` to `cvsweep`.
"""
function scanrate_varied_sweep(
        elydata,
        sawtooth,
        grid,
        bcond,
        reaction;
        nperiods = 1,
        scanrates = [0.002, 0.003, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 10.0],
        eneutral = false,
        sweep_kwargs...
    )
    sweep_vec = Vector{Any}(undef, length(scanrates))
    saws = Vector{Any}(undef, length(scanrates))

    for (i, sr) in pairs(scanrates)
        saw = SawTooth(
            scanrate = sr,
            vmin = sawtooth.vmin,
            vmax = sawtooth.vmax,
            scanup = sawtooth.scanup
        )
        saws[i] = saw
        celldata = deepcopy(elydata)
        celldata.eneutral = eneutral

        pnpcell = PNPSystem(grid; bcondition = bcond, celldata = celldata, reaction = reaction)
        sweep_vec[i] = cvsweep(pnpcell; voltages = saw, nperiods, store_solutions = true)
    end

    return (scanrates = scanrates, saws = saws, sweeps = sweep_vec)
end


"""
    run_pH_sweep(elydata, sawtooth, grid, bcond, reaction; pH_values, eneutral, counter_index, ...)

CV at each bulk pH in `pH_values`.

Sets H⁺ and OH⁻ from the pH and the water constant `Kw`. With `eneutral = true` the species
at `counter_index` is adjusted to restore bulk electroneutrality, and the run errors rather
than continuing if that would need a negative concentration. Returns one named tuple per pH
carrying the concentrations used and the sweep record.
"""
function run_pH_sweep(
        elydata,
        sawtooth,
        grid,
        bcond,
        reaction;
        nperiods = 1,
        pH_values = [3.0, 4.0, 5.0, 6.0, 6.8, 7.0, 8.0, 9.0],
        hplus_index::Int = 2,
        ohminus_index::Int = 6,
        eneutral::Bool = false,
        counter_index::Union{Nothing, Int} = nothing,
        Kw::Float64 = 1.0e-14,
        store_solutions::Bool = true,
        sweep_kwargs...
    )
    npH = length(pH_values)
    results = Vector{NamedTuple}(undef, npH)

    for (k, pH) in pairs(pH_values)
        ely = deepcopy(elydata)

        # --- basic checks
        nsp = length(ely.c_bulk)
        @assert 1 <= hplus_index <= nsp "hplus_index out of bounds"
        @assert 1 <= ohminus_index <= nsp "ohminus_index out of bounds"
        @assert hplus_index != ohminus_index "hplus_index and ohminus_index must be different"

        # NOTE:
        # Replace `ely.z` with your actual charge-number field if needed
        charges = ely.z
        @assert length(charges) == nsp "length(charges) != length(c_bulk)"

        # --- set H+ / OH- from pH
        cH = 10.0^(-pH)
        cOH = Kw / cH

        ely.c_bulk[hplus_index] = cH
        ely.c_bulk[ohminus_index] = cOH

        # --- enforce electroneutrality by adjusting one counter ion
        if eneutral
            counter_index === nothing &&
                error("eneutral=true requires counter_index")

            @assert 1 <= counter_index <= nsp "counter_index out of bounds"
            @assert counter_index != hplus_index "counter_index must differ from hplus_index"
            @assert counter_index != ohminus_index "counter_index must differ from ohminus_index"

            zc = charges[counter_index]
            abs(zc) > eps(Float64) ||
                error("counter species must have nonzero charge")

            q_except = zero(eltype(ely.c_bulk))
            for i in eachindex(ely.c_bulk)
                if i != counter_index
                    q_except += charges[i] * ely.c_bulk[i]
                end
            end

            c_counter_new = -q_except / zc

            if c_counter_new < 0
                error(
                    "electroneutral correction produced negative concentration " *
                        "for counter_index=$counter_index at pH=$pH: $c_counter_new"
                )
            end

            ely.c_bulk[counter_index] = c_counter_new
        end

        # --- rebuild system
        pnpcell = PNPSystem(
            grid;
            bcondition = bcond,
            celldata = ely,
            reaction = reaction,
        )

        rec = cvsweep(
            pnpcell;
            voltages = sawtooth,
            nperiods = nperiods,
            store_solutions = store_solutions,
            sweep_kwargs...
        )

        results[k] = (
            pH = pH,
            cH = cH,
            cOH = cOH,
            c_bulk = copy(ely.c_bulk),
            record = rec,
        )
    end

    return results
end


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 4: extra CV sweep drivers.
# =====================================================================

# ── moved from cell e50fe651 (sweep)  [dropped dead `reaction_arg` line that
#    referenced the notebook global `elydata_Gold_unc`; it was unused] ──
"""
    sweep(model, grid, bcondition, reaction, sawtooth; nperiods = 1, kwargs...)

One cyclic voltammogram. Copies the electrolyte so the caller's model is untouched, builds
a `PNPSystem` and runs `cvsweep`. `kwargs` go straight to the solver.
"""
function sweep(model, grid, bcondition, reaction, sawtooth; nperiods = 1, eneutral = true, tunnel = false, bikerman = true, kwargs...)
    celldata = deepcopy(model)
    #celldata.eneutral = eneutral
    pnpcell = PNPSystem(grid; bcondition = bcondition, celldata = celldata, reaction = reaction)
    return result = LiquidElectrolytes.cvsweep(
        pnpcell;
        voltages = sawtooth,
        nperiods,
        store_solutions = true,
        kwargs...
    )

end

# ── moved from cell 9b8daa64 (cvsweep_compensated_over_L) ──
"""
    cvsweep_compensated_over_L(model, bcondition, reaction; sawtooth, L_values, f_comp, Area, ...)

CV over boundary-layer thickness with the ohmic drop compensated.

Builds a grid per `L`, computes the bulk resistance from the conductivity and `Area`, and
passes `f_comp × R_bulk` to `cvsweep_COMP`. `f_comp = 1` is full compensation. Returns
`Dict(L => result)` keyed by `L` in μm.
"""
function cvsweep_compensated_over_L(
        model, bcondition, reaction;
        sawtooth,
        nperiods = 1,
        L_values = [100, 500, 1000, 2000, 3000],
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

        pnpcell = PNPSystem(
            grid;
            bcondition = bcondition,
            reaction = reaction,
            celldata = celldata
        )

        # --- Resistance Calculation ---
        κ = calc_kappa(celldata)
        R_bulk_theoretical = (L * μm) / (κ * Area)

        # R_compensated will be passed to cvsweep_COMP to 'undo' the voltage drop
        R_to_apply = R_bulk_theoretical * f_comp

        @info "    Nodes: $(length(X)) | R_bulk: $(round(R_bulk_theoretical, digits = 2)) Ω"

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

# ── moved from cell d91c32c8 (cvsweep_odr_over_L) ──
"""
    cvsweep_odr_over_L(elydata_odr, grid_dict, bcondition, reaction, sawtooth; ...)

CV over boundary-layer thickness on pre-built grids.

Takes `grid_dict = Dict(L => grid)` rather than building meshes itself. When the
electrolyte carries an `OhmicDropEstimation`, the uncompensated resistance is rescaled per
`L` as `Ru = L / κ`. Returns `Dict(L => result)`.
"""
function cvsweep_odr_over_L(
        elydata_odr, grid_dict, bcondition, reaction, sawtooth;
        nperiods = 1,
        store_solutions = true,
        unknown_storage = :dense,
        solver_kwargs...
    )
    results = Dict{Float64, Any}()
    L_values = sort(collect(keys(grid_dict)))

    for L in L_values
        grid = grid_dict[L]
        if isa(elydata_odr.ircompensation, OhmicDropEstimation)
            Ru=L/conductivity(elydata_odr, elydata_odr.c_bulk)
            celldata=copy(elydata_odr;
                          ircompensation=copy(elydata_odr.ircompensation;Ru=Ru)
                          )
        else
            celldata=copy(elydata_odr)
        end
        @info ">>> :$(celldata.ircompensation) sweep | L = $(L / μm) μm"

        pnpcell = PNPSystem(grid; bcondition, celldata, reaction, unknown_storage)
        @time results[L] = LiquidElectrolytes.cvsweep(
            pnpcell;
            voltages = sawtooth,
            nperiods,
            store_solutions,
            solver_kwargs...
        )
    end

    return results
end

# ── moved from cell 82116926 (blthickness), fixed ──

"""
    blthickness(grid, celldata, tsol; species = 5, atol = 1.0e-1)

Estimate the diffusion boundary-layer thickness of `species` (default 5 = CO₂) from a
transient solution.

Scans every time snapshot in `tsol` and finds the outermost grid node where the
concentration deviates from the bulk value `celldata.c_bulk[species]` by more than
`atol` (mol/m³). Returns the **maximum** such position over all time steps (SI units,
i.e. meters); 0 if the profile never deviates.

Set `atol` relative to the species: the default 0.1 mol/m³ suits CO₂
(`c_bulk ≈ 33 mol/m³`), but trace species (OH⁻, H⁺) need a much smaller `atol`.
"""
function blthickness(grid, celldata, tsol; species = 5, atol = 1.0e-1)
    X = grid[XCoordinates]
    cb = celldata.c_bulk[species]
    xbl = zero(eltype(X))
    for it in 1:length(tsol.t)
        u = tsol[species, :, it]
        i = findlast(c -> abs(c - cb) > atol, u)
        i === nothing && continue          # profile already at bulk everywhere
        xbl = max(xbl, X[i])
    end
    return xbl
end

"""
    blthickness_t(grid, celldata, tsol; species = 5, atol = 1.0e-1)

Time-resolved boundary-layer thickness: same criterion as [`blthickness`](@ref) but
evaluated at every stored time. Returns `(times, δ)`.
"""
function blthickness_t(grid, celldata, tsol; species = 5, atol = 1.0e-1)
    X = grid[XCoordinates]
    cb = celldata.c_bulk[species]
    t = tsol.t
    δ = zeros(eltype(X), length(t))
    for it in eachindex(t)
        u = tsol[species, :, it]
        i = findlast(c -> abs(c - cb) > atol, u)
        δ[it] = i === nothing ? zero(eltype(X)) : X[i]
    end
    return t, δ
end


# module?
