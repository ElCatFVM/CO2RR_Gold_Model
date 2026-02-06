function capscalc(sys, Au_Check::Bool; molarities= molarities = [0.005, 0.05, 0.5], vrange = range(-1, 1, length = 201))
    result = []
    
    # Index Sequence: [1:K+, 2:H+, 3:HCO3, 4:CO3, 5:CO2, 6:OH, 7:CO]
    # Constants used: K1=4.44e-7, K2=4.66e-11, Kw=1.0e-14, CO2=0.033
    conc_map = Dict(
        0.01 => [0.01, 1.465e-6, 0.01000, 3.181e-7, 0.033, 6.826e-9, 0.0],
        #0.03 => [0.03, 4.885e-7, 0.02999, 2.861e-6, 0.033, 2.047e-8, 0.0],
        #0.05 => [0.05, 2.931e-7, 0.04998, 7.946e-6, 0.033, 3.411e-8, 0.0],
        0.1  => [0.1,  1.466e-7, 0.09994, 3.176e-5, 0.033, 6.821e-8, 0.0],
        0.3  => [0.3,  4.893e-8, 0.29943, 2.852e-4, 0.033, 2.044e-7, 0.0],
        #0.5  => [0.5,  2.940e-8, 0.49842, 7.901e-4, 0.033, 3.402e-7, 0.0]
    )

    molarities = [
        0.01, 
        #0.03, 
        #0.05, 
        0.1, 
        0.3, 
        #0.5
    ]
    

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
      data = isa(sys, AbstractElectrochemicalSystem) ? sys.vfvmsys.physics.data : sys.physics.data
        
        for c_total in molarities
            c_vals = conc_map[c_total]
            

            for i in 1:length(c_vals)
                data.c_bulk[i] = c_vals[i] * ufac"mol/dm^3"
            end
            
            if isa(sys, AbstractElectrochemicalSystem)
                t = @elapsed r = dlcapsweep(sys, voltages = vrange)
                volts, caps = r.voltages, r.dlcaps
            else
                t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
            end
            
            cdl0 = dlcap0(data)
            @info "Molarity $(c_total)M calculated in $(round(t, digits=2))s"
            
            push!(result, (voltage_range=volts, dlcaps=caps, cdl0=cdl0, molarity=c_total))
        end
    end

    return result
end