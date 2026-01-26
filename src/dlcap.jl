


function capscalc(sys; molarities =molarities)
    result = []
	vrange = range(-1, 1, length = 201)
	if is_Landstorfer
		for imol in 1:length(molarities)
		    if !isa(sys, AbstractElectrochemicalSystem)
		     	data = sys.physics.data
		        set_molarity!(data, molarities[imol])
		      	t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 		101)
		    else
		        data = sys.vfvmsys.physics.data
		        data.c_bulk .= molarities[imol] * ufac"mol/dm^3"
		    	t = @elapsed r = dlcapsweep(
		                sys,
		     	        voltages = range(-1, 1, length = 201)
		        )
		        volts = vrange
		        caps = r.dlcaps
		    end
		    cdl0 = dlcap0(data)
		    @info "elapsed=$(t)"
		    push!(result, (voltage_range = volts, dlcaps = caps, cdl0 = cdl0, molarity = 	molarities[imol]))
		end
	else
	   	if !isa(sys, AbstractElectrochemicalSystem)
	   		data = sys.physics.data
	        t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
		else
	        data = sys.vfvmsys.physics.data
	        t = @elapsed r = dlcapsweep(
	    	        sys,
	                voltages = vrange
	        )
			volts = r.voltages
			caps = r.dlcaps
	    end
	    cdl0 = dlcap0(data)
	    @info "elapsed=$(t)"
	    push!(result, (voltage_range = volts, dlcaps = caps, cdl0 = cdl0))
	end
    return result
end
