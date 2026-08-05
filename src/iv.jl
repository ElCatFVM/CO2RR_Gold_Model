

"""
    pressure_varied_sweep(elydata_base, sweepfun; Pvec, ispec, base_value, scale, kwargs...)

Sweep one species' bulk concentration over a list of pressures.

For each `p` in `Pvec` it sets `c_bulk[ispec] = base_value * scale(p)` — Henry's law when
`scale` is `identity` and `base_value` is the saturated concentration — and calls
`sweepfun(ely; kwargs...)`. Returns `(p, result)` pairs, the shape the pressure plots
expect.

`sweepfun` decides whether this is a CV or an IV study: pass a closure that builds a
`PNPSystem` and calls `cvsweep` to get pressure-resolved voltammograms.
"""
function pressure_varied_sweep(
    elydata_base, sweepfun;
    Pvec,
    ispec::Integer,
    base_value = elydata_base.c_bulk[ispec],
    scale = identity,
    sweep_kwargs...
)
    recs = Vector{Any}(undef, length(Pvec))
    for (k, p) in pairs(Pvec)
        ely = deepcopy(elydata_base)
        ely.c_bulk[ispec] = base_value .* scale(p)
        recs[k] = sweepfun(ely; sweep_kwargs...)
    end
    return collect(zip(Pvec, recs))
end


"""
    sweep_over_L_cv(model; L_values, bcond, sweepfun, eneutral, reaction, store_solutions)

Sweep over boundary-layer thickness with a caller-supplied `sweepfun`, building a grid per
`L` (μm). Returns `Dict(L => result)`.

Overrides the `cv.jl` definition of the same name — `iv.jl` is included later.
"""
function sweep_over_L_cv(
    model;
    L_values,
    bcond,
    sweepfun,
    eneutral = true,
    reaction = nothing,
    store_solutions = true,
)
    results = Dict{Float64, Any}()

    for L in L_values
        hmin = 1.0e-6 * μm
        hmax = 1.0    * μm
        X    = ExtendableGrids.geomspace(0, L * μm, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)

        celldata = deepcopy(model)
        celldata.eneutral = eneutral

        reaction_kw = (reaction === nothing) ? NamedTuple() : (; reaction)

        pnpcell = PNPSystem(grid; bcondition=bcond, celldata=celldata, reaction_kw...)

        results[L] = sweepfun(
            pnpcell;
           # voltages = voltages,
           # nperiods = nperiods,
            store_solutions = store_solutions,
        )
    end

    return results
end



"""
    scanrate_varied_sweep(elydata, user_input_cv; scanrates, sweepfun, sweep_kwargs...)

Run `sweepfun` once per scan rate in `scanrates`.

Overrides the `cv.jl` definition of the same name — `iv.jl` is included later.
"""
function scanrate_varied_sweep(
    elydata,
    user_input_cv;
    scanrates = [0.002, 0.003, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 10.0],
    eneutral = false,
    tunnel   = false,
    bikerman = true,
    sweepfun = sweep2,
    sweep_kwargs...
)
    sweep_vec = Vector{Any}(undef, length(scanrates))
    saws      = Vector{Any}(undef, length(scanrates))

    for (i, sr) in pairs(scanrates)
        saw = SawTooth(
            scanrate = sr,
            vmin     = user_input_cv.vmin,
            vmax     = user_input_cv.vmax,
            scanup   = user_input_cv.scanup
        )
        saws[i] = saw
        sweep_vec[i] = sweepfun(elydata, saw; eneutral=eneutral, tunnel=tunnel, bikerman=bikerman, sweep_kwargs...)
    end

    return (scanrates=scanrates, saws=saws, sweeps=sweep_vec)
end



"""
    ivsweep_over_L(model, bcondition; voltages, L_values, solver_control, kwargs...)

Steady-state IV curve at each boundary-layer thickness in `L_values` (μm), building a grid
per `L`. Returns `Dict(L => ivresult)`.
"""
function ivsweep_over_L(model, bcondition;
	voltages,
 	L_values = round.(Int, range(80, 1500, length=6)),
    eneutral = true,
	solver_control,
	kwargs...)
    results = Dict{Int, Any}()
    kwargs 	 	= merge(solver_control, kwargs) 

	#cell, ivresult
    for L in L_values
        hmin = 1.0e-6 * μm
        hmax = 1.0    * μm
        X = ExtendableGrids.geomspace(0, L * μm, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)

        celldata = deepcopy(model)
        #celldata.eneutral = eneutral
        #celldata.tunnel   = tunnel
        #celldata.bikerman = bikerman

#        reaction_kw = (model === elydata_Gold) ? (; reaction) : (;)
    #cell   = PNPSystem(grid; bcondition=pnp_bcondition, reaction=reaction, celldata)

        pnpcell = PNPSystem(grid; bcondition=bcondition, reaction=reaction, celldata=celldata)

        results[L] = ivsweep(pnpcell; voltages, store_solutions=true, kwargs...)

    end

    return results
end


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 5: pressure column-name helper + IV simulate drivers.
# simulate_CO2R/_dir referenced the notebook globals `pnp_bcondition`,
# `reaction`, `solver_control` in their bodies; here they are passed as
# arguments so the function works inside the package.
# =====================================================================

# ── moved from cell b767a48f (_pressure_colname) ──
"""
    _pressure_colname(p)

CSV column name for a CO₂ partial pressure, e.g. `0.25` → `:pCO2_0p25_atm`. Decimal points
become `p` and minus signs `m` so the result is a valid identifier.
"""
function _pressure_colname(p)
    s = @sprintf("%.4g", p)
    s = replace(s, "." => "p", "-" => "m", "+" => "")
    return Symbol("pCO2_" * s * "_atm")
end

# ── moved from cell 06ca7d68 (simulate_CO2R)  [bcondition/reaction/solver_control now args] ──
"""
    simulate_CO2R(grid, celldata; bcondition, reaction, solver_control, voltages, kwargs...)

One IV sweep on a given grid and electrolyte. Returns `(cell, ivresult)`.
"""
function simulate_CO2R(grid, celldata;
        bcondition,
        reaction,
        solver_control = (;),
        voltages = (-1.5:0.1:0.0) * V,
        kwargs...)
    kwargs   = merge(solver_control, kwargs)
    cell     = PNPSystem(grid; bcondition = bcondition, reaction = reaction, celldata)
    ivresult = ivsweep(cell; voltages, store_solutions = true, kwargs...)
    return cell, ivresult
end

# ── moved from cell 012b426e (simulate_CO2R_dir)  [identical body to simulate_CO2R] ──
"""
    simulate_CO2R_dir(grid, celldata; bcondition, reaction, solver_control, voltages, kwargs...)

As [`simulate_CO2R`](@ref); kept as a separate entry point for the directory-scanning
scripts. The bodies are currently identical.
"""
function simulate_CO2R_dir(grid, celldata;
        bcondition,
        reaction,
        solver_control = (;),
        voltages = (-1.5:0.1:0.0) * V,
        kwargs...)
    kwargs   = merge(solver_control, kwargs)
    cell     = PNPSystem(grid; bcondition = bcondition, reaction = reaction, celldata)
    ivresult = ivsweep(cell; voltages, store_solutions = true, kwargs...)
    return cell, ivresult
end


 # module?
