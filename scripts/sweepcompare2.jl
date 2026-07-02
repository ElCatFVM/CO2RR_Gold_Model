using AuCO2RR




function main()
    ivsweep_csv(actcoeff = :DGML_γ!, hydrated = false, bcmodel = :Robin, model = :Gold)
    ivsweep_csv(actcoeff = :Stefan_γ!, hydrated = true, bcmodel = :Robin, model = :Gold)
    ivsweep_csv(actcoeff = :Potassium_γ!, hydrated = false, bcmodel = :Robin, model = :Gold)
    return nothing
end


if !isinteractive()
    main()
end
