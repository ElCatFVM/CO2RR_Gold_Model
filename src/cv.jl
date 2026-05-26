


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

        pnpcell = PNPSystem(grid; bcondition=bcond, celldata=celldata, reaction = reaction)

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
    counter_index::Union{Nothing,Int} = nothing,
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
        cH  = 10.0^(-pH)
        cOH = Kw / cH

        ely.c_bulk[hplus_index]  = cH
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



 # module?
