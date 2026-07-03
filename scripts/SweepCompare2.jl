"""
    SweepCompare2

Script which uses the shared AuCO2RR package project.

Used as third step after SweepCompare1 during  refactoring.
"""
module SweepCompare2
    using AuCO2RR

    function main()
        ivsweep_csv(actcoeff = :DGML, hydrated = false, bcmodel = :Robin, model = :Gold)
        ivsweep_csv(actcoeff = :Stefan, hydrated = true, bcmodel = :Robin, model = :Gold)
        ivsweep_csv(actcoeff = :Potassium, hydrated = false, bcmodel = :Robin, model = :Gold)
        return nothing
    end

end
abspath(PROGRAM_FILE) == @__FILE__() &&  SweepCompare2.main(ARGS)
