


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
        hmax = 1.0    * μm
        X    = ExtendableGrids.geomspace(0, L, hmin, hmax)
        grid = ExtendableGrids.simplexgrid(X)

        celldata = deepcopy(elydata)
        celldata.eneutral = eneutral

        pnpcell = PNPSystem(grid; bcondition=bcond, celldata=celldata)

        saw = SawTooth(
            scanrate = sawtooth.scanrate,
            vmin     = sawtooth.vmin,
            vmax     = sawtooth.vmax,
            scanup   = sawtooth.scanup
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
    saws      = Vector{Any}(undef, length(scanrates))

    for (i, sr) in pairs(scanrates)
        saw = SawTooth(
            scanrate = sr,
            vmin     = sawtooth.vmin,
            vmax     = sawtooth.vmax,
            scanup   = sawtooth.scanup
        )
        saws[i] = saw
        celldata = deepcopy(elydata)
        celldata.eneutral = eneutral    

        pnpcell = PNPSystem(grid; bcondition=bcond, celldata=celldata, reaction = reaction)
        sweep_vec[i] = cvsweep(pnpcell; voltages = saw, nperiods, store_solutions = true)
    end

    return (scanrates=scanrates, saws=saws, sweeps=sweep_vec)
end


function run_pH_sweep(
    elydata, sawtooth, grid, bcond, reaction; 
    nperiods = 1,
    pH_values = [3.0, 4.0, 5.0, 6.0, 6.8, 7.0, 8.0, 9.0],
    hplus_index::Int = 2,
    ohminus_index::Int = 6,
    eneutral = false,
    sweep_kwargs...
)
    recs = Vector{Any}(undef, length(pH_values))

    for (k, pH) in pairs(pH_values)
        ely = deepcopy(elydata)
        ely.c_bulk[hplus_index]  = 10.0^(-pH)
        ely.c_bulk[ohminus_index] = 10.0^(-14 + pH)

        pnpcell = PNPSystem(grid; bcondition=bcond, celldata=ely, reaction = reaction)
        recs[k] = cvsweep(pnpcell; voltages = sawtooth, nperiods, store_solutions = true)
    end

    return collect(zip(pH_values, recs))
end



 # module?
