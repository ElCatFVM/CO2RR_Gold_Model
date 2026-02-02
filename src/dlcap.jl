function capscalc(sys, Au_Check::Bool; molarities= molarities = [0.005, 0.05, 0.5], comb = [0.01, 0.1, 1], vrange = range(-1, 1, length = 201))
    result = []
    

    if Au_Check
        for imol in 1:length(molarities)
            if !isa(sys, AbstractElectrochemicalSystem)
                data = sys.physics.data
                set_molarity!(data, molarities[imol])
                t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
            else
                data = sys.vfvmsys.physics.data
                data.c_bulk .= molarities[imol] * ufac"mol/dm^3"
                t = @elapsed r = dlcapsweep(sys, voltages = vrange)
                volts = r.voltages
                caps  = r.dlcaps
            end
            cdl0 = dlcap0(data)
            @info "elapsed=$(t)"
            push!(result, (voltage_range=volts, dlcaps=caps, cdl0=cdl0, molarity=molarities[imol]))
        end

    else
        # Scale c_bulk by factors in `comb` to study concentration dependence of C_gap
        if !isa(sys, AbstractElectrochemicalSystem)
            data = sys.physics.data
            base_c_bulk = deepcopy(data.c_bulk)

            for f in comb
                data.c_bulk .= base_c_bulk .* f
                t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)

                cdl0 = dlcap0(data)
                @info "elapsed=$(t)"
                push!(result, (voltage_range=volts, dlcaps=caps, cdl0=cdl0, comb=f))
            end
        else
            data = sys.vfvmsys.physics.data
            base_c_bulk = deepcopy(data.c_bulk)

            for f in comb
                data.c_bulk .= base_c_bulk .* f
                t = @elapsed r = dlcapsweep(sys, voltages = vrange)

                volts = r.voltages
                caps  = r.dlcaps
                cdl0  = dlcap0(data)

                @info "elapsed=$(t)"
                push!(result, (voltage_range=volts, dlcaps=caps, cdl0=cdl0, comb=f))
            end
        end
    end

    return result
end