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
