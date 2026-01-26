

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


function sweep_over_L(
    model;
    L_values,
    bcond,
    voltages,
    nperiods,
    sweepfun,                
    eneutral = true,
    tunnel   = false,
    bikerman = true,
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
            voltages = voltages,
            nperiods = nperiods,
            store_solutions = store_solutions,
        )
    end

    return results
end



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


function run_pH_sweep(
    elydata_base;
    pH_values = [3.0, 4.0, 5.0, 6.0, 6.8, 7.0, 8.0, 9.0],
    hplus_index::Int = 2,
    ohminus_index::Int = 6,
    eneutral::Bool = true,
    tunnel::Bool = false,
)
    results = Any[]

    for pH in pH_values
        ely = deepcopy(elydata_base)
        ely.c_bulk[hplus_index]  = 10.0^(-pH)
        ely.c_bulk[ohminus_index] = 10.0^(-14 + pH)

        rec = sweep(ely; eneutral=eneutral, tunnel=tunnel)
        push!(results, rec)
    end

    return (pH_values = pH_values, results = results)
end


 # module?
