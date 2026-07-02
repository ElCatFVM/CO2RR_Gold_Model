# plots/cvplot.jl

using CairoMakie
using CSV
using DataFrames




"""
    plot_time_voltage_and_dt(pnpresult, sawtooth;
                             fig_size=(700, 450),
                             dt_yscale=log10,
                             use_times_field=true,
                             t_step=1.0)

Create a 2-row figure:
- Top: time vs voltage `sawtooth.(T)`
- Bottom: time-step sizes `Δt` (log-scale by default)

Arguments
- `pnpresult`: result object containing time information (expects `pnpresult.times` or `pnpresult.tsol.t`)
- `sawtooth`: callable voltage protocol (e.g., `SawTooth` object) supporting `sawtooth(t)`
- `use_times_field`: if true, use `pnpresult.times`; otherwise use `pnpresult.tsol.t`
- `t_step`: sampling step used only when `use_times_field=false` and we build `T = 0:t_step:tend`

Returns
- `fig`
"""
function plot_time_voltage_and_dt(pnpresult, sawtooth;
                                  fig_size=(700, 450),
                                  dt_yscale=log10,
                                  use_times_field=true,
                                  t_step=1.0)

    fig = Figure(size = fig_size)

    axV = Axis(fig[1, 1],
               xlabel = "Time (s)",
               ylabel = "Voltage")

    axdt = Axis(fig[2, 1],
                xlabel = "Time (s)",
                ylabel = "Δt (s)",
                yscale = dt_yscale)

    # Choose time vector source.
    if use_times_field && hasproperty(pnpresult, :times)
        T = pnpresult.times
    else
        tend = pnpresult.tsol.t[end]
        T = 0:t_step:tend
    end

    # Plot time-voltage.
    lines!(axV, T, sawtooth.(T))

    # Plot time-step sizes (Δt).
    if length(T) >= 2
        dt = T[2:end] .- T[1:end-1]
        lines!(axdt, T[2:end], dt)
    end

    return fig
end



"""
    plot_cv_current(result, model; species=nothing, fig_size=(650, 400), scale=cm^2/mA)

Plot `currents(result, species) .* scale` versus `result.voltages`.
If `species` is `nothing`, use `model.cspecies[1]`.
Returns `fig`.
"""
function plot_cv_current(result, m;
                         species=nothing,
                         fig_size=(650, 400),
                         scale=cm^2/mA,
                         color_gradient=true)
    model = m.elydata

    sp = (species === nothing) ? model.cspecies[1] : species

    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1],
              ylabel = lab_current,
              xlabel = lab_voltage)

    I = currents(result, sp) .* scale

    if color_gradient
        cols = RGBf.(range(0, 1, length(result.voltages)), 0.0, 0.0)
        lines!(ax, result.voltages, I; color=cols)
    else
        lines!(ax, result.voltages, I)
    end

    return fig
end


"""
    plot_conc_time_electrode(result, bulk; nspecies=7, fig_size=(650, 400), scale=(mol/dm^3))

Plot electrode-adjacent concentrations `result.tsol[i, 1, t] / scale` versus time for `i=1:nspecies`
(using `log10` y-scale; nonpositive values are shown as `NaN`). Returns `fig`.
"""
function plot_conc_time_electrode(result, m;
                                  nspecies=7,
                                  fig_size=(800, 500),
                                  scale=(mol/dm^3)) # if the unit is uncertain, test with 1.0 first
    bulk = m.bulk

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
               xlabel = lab_time,
               ylabel = rich(rich("c", font = :italic), subscript("i, electrode"),
                             "  (mol/dm", superscript("3"), ")"),
               limits = ((times[1]-(times[end] / 200), times[end] + (times[end] / 100)), (1e-12, 1e4)),
               yscale = log10,
               yticks = (10.0 .^ (4:-4:-12),
                         [powlab(4), powlab(0), powlab(-4), powlab(-8), powlab(-12)]))

    for i in 1:nspecies
        y = conc[i, :]
        
        y_fixed = map(c -> (c > 1e-128 ? c : 1e-128), y)
        
        lines!(ax, times, y_fixed; color=colors[i], label=string(names[i]))
    end


    Legend(fig[1, 2], ax; labelsize=10, backgroundcolor=RGBA(1, 1, 1, 0.5))
    return fig
end

"""
    plot_conc_profile_logx(result, bulk, X;
                           t_index=1,
                           fig_size=(650, 400),
                           x_limits=(1e-12, 1e-3),
                           y_limits=(-14, 1),
                           x_offset=1e-14,
                           scale=(mol/dm^3))

Plot log10(concentration) profiles versus distance from the electrode (log10 x-scale) at a given time index.
Returns `fig`.
"""
function plot_conc_profile_logx(result, m, X, t_index;
                                fig_size=(650, 400),
                                x_limits=(1e-12, 1e-3),
                                y_limits=(-14, 4), # range in terms of log10(c)
                                x_offset=1e-14,
                                scale=1.0)
    bulk = m.bulk

    tsol = result.tsol
    nvar, nx, nt = size(tsol)

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    nspecies = min(nvar, length(names), length(colors))

    nplot = min(nx, length(X))
    ti = clamp(t_index, 1, nt)

    xx = X[1:nplot] .+ x_offset

    fig = Figure(size = fig_size)
    
    phi_val = hasproperty(result, :voltages) ? round(result.voltages[ti], digits=2) : "N/A"
    
    ax = Axis(fig[1, 1],
        xlabel = "Distance from electrode [m]",
        ylabel = rich("log", subscript("10"), "(", rich("c", font = :italic), ")"),
        xscale = log10,
        limits = (x_limits, y_limits),
        title  = "t = $(round(result.tsol.t[ti], digits=4)) s | ϕ = $phi_val V",
    )

    for i in 1:nspecies
        conc = tsol[i, 1:nplot, ti] ./ scale
        
        y_log = map(c -> log10(max(c, 1e-25)), conc)
        
        lines!(ax, xx, y_log; color=colors[i], label=string(names[i]))
    end

    axislegend(ax; position=:rt, labelsize=10)
    return fig
end

"""
    cv_conc_gif(pnpresult, bulk, X;
                file="concentrations_cv.gif",
                framerate=10,
                step=5,
                scale=(mol/dm^3),
                x_offset=1e-14,
                x_limits=(1e-12, 1e-3),
                y_limits=(-14, 2),
                legend=true)

Create a GIF of log10(concentration) profiles vs distance (log10 x-scale) over a CV.
Returns the absolute filepath as a string.
"""
function cv_conc_gif(
    pnpresult, m, X;
    file="concentrations_cv.gif",
    framerate=10,
    step=5,
    scale=(mol/dm^3),
    x_offset=1e-14,
    x_limits=(1e-12, 1e-3),
    y_limits=(-14, 2),
    legend=true,
)
    bulk = m.bulk

    raw_data = pnpresult.tsol
    tsol = raw_data ./ scale
    nvar, nx, nt = size(tsol)

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    nspecies = min(nvar, length(names), length(colors))

    nplot = min(nx, length(X))
    xx = X[1:nplot] .+ x_offset

    t_indices = 1:step:nt

    fig = Figure(size=(650, 400))
    ax  = Axis(fig[1, 1],
        xlabel = "Distance from electrode [m]",
        ylabel = "log10(c)",
        xscale = log10,
        limits = (x_limits, y_limits),
    )

    println("Plotting $nspecies species, $nt time steps.")

    ys = [Observable(fill(NaN, nplot)) for _ in 1:nspecies]
    for i in 1:nspecies
        lines!(ax, xx, ys[i]; color=colors[i], label=string(names[i]))
    end
    legend && axislegend(ax)

    record(fig, file, t_indices; framerate=framerate) do ti
        tval = try pnpresult.tsol.t[ti] catch; ti end
        
        phi_str = "N/A"
        if hasproperty(pnpresult, :voltages) && ti <= length(pnpresult.voltages)
            phi_str = "$(round(pnpresult.voltages[ti], digits=2))"
        end
        ax.title = "t = $(round(tval, digits=4)) s | ϕ = $phi_str V"

        for i in 1:nspecies
            conc = tsol[i, 1:nplot, ti]
            # 3. robust log: replace non-positive values with a small floor (1e-20) to avoid breaks in the curve
            ys[i][] = map(c -> (c > 1e-20 ? log10(c) : -20.0), conc)
        end
    end

    return abspath(file)
end



function plot_cv_model_vs_koper_facets(pnpresult; species=iohminus,
                                       koper_csv_relpath="data/Langmuir_CV_data/Figure_1.csv",
                                       fig_size=(1050, 650),
                                       koper_v_shift=-0.4)

    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1],
              ylabel = lab_current,
              xlabel = lab_voltage)

    # --- Gold model ---
    I_model = currents(pnpresult, species) .* (cm^2/mA)
    gold_line = lines!(ax, pnpresult.voltages, I_model;
                       color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))

    # --- Koper data (Figure 1) ---
    root = normpath(joinpath(@__DIR__, ".."))
    csv_path = joinpath(root, splitdir(koper_csv_relpath)...)

    raw_df = CSV.read(csv_path, DataFrame; header=false)
    facet_row = collect(raw_df[1, :])

    numeric_data = [
        parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
        for i in 3:nrow(raw_df)
    ]
    num_df = DataFrame(hcat(numeric_data...)', names(raw_df))

    colors = (:pink, :skyblue, :lightgreen)

    facet_lines = [
        lines!(ax, num_df[!, 1] .+ koper_v_shift, num_df[!, 2], color = colors[1]),
        lines!(ax, num_df[!, 3] .+ koper_v_shift, num_df[!, 4], color = colors[2]),
        lines!(ax, num_df[!, 5] .+ koper_v_shift, num_df[!, 6], color = colors[3]),
    ]
    labels1 = ["CO2RR Gold Model"]
    labels2 = [string(facet_row[1]), string(facet_row[3]), string(facet_row[5])]

    Legend(fig[1, 2],
           [[gold_line], facet_lines],
           [labels1, labels2],
           ["Model", "Koper\nFacets"])

    return fig
end


function plot_pH_varied_sweep(
    pH_recs;
    species = iohminus,
    fig_size = (1600, 900),
    scale = cm^2 / mA,
    legend_title = "Theoretical",
)

    fig = Figure(size = fig_size)
    ax = Axis(
        fig[1, 1],
        xlabel = lab_voltage,
        ylabel = lab_current,
    )

    n = length(pH_recs)
    cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]

    plots  = Any[]
    labels = String[]

    for (j, item) in pairs(pH_recs)

        # ---- support both:
        # 1) old format: (pH, rec)
        # 2) new format: (pH=..., cH=..., cOH=..., c_bulk=..., record=...)
        pH, rec = if item isa NamedTuple
            if haskey(item, :record) && haskey(item, :pH)
                (item.pH, item.record)
            else
                error("NamedTuple input must contain at least :pH and :record")
            end
        elseif item isa Tuple
            if length(item) >= 2
                (item[1], item[2])
            else
                error("Tuple input must have at least 2 elements: (pH, rec)")
            end
        else
            error("Unsupported element type in pH_recs: $(typeof(item))")
        end

        ivres = hasproperty(rec, :ivresult) ? getproperty(rec, :ivresult) : rec
        I = currents(ivres, species) .* scale

        plt = lines!(ax, ivres.voltages, I; color = cols[j], linewidth = 3)
        push!(plots, plt)
        push!(labels, "pH = $(pH)")
    end

    Legend(fig[1, 2], plots, labels, legend_title; framevisible = true)

    return fig
end







"""
plot_iv_with_experiment(pnpresult, ico; csv_path, exp_title="Experimental")

- pnpresult must provide: pnpresult.voltages
- currents(pnpresult, ico) must be defined in the caller environment
- csv is expected to match your Figure_3.csv parsing logic
"""
function plot_iv_with_experiment(pnpresult, ico;
    csv_path::AbstractString = "../data/Langmuir_CV_data/Figure_3.csv",
    exp_title::AbstractString = "Experimental",
    fig_size::Tuple{Int,Int} = (1600, 900),
    markersize::Real = 12,
)
    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1], ylabel = lab_current, xlabel = lab_voltage)

    volts = vec(pnpresult.voltages)
    total_current = currents(pnpresult, ico) .* (cm^2/mA)

    colgrad = RGBf.(range(0, 1, length(volts)), 0.0, 0.0)
    lines!(ax, volts, total_current; color = colgrad)
    scatter!(ax, volts, total_current; markersize = markersize, color = colgrad)

    plot_objs = Any[]
    labels    = String[]

    try
        raw = CSV.read(csv_path, DataFrame; header=false)

        pres = vec(Matrix(raw[1:1, :]))
        sub  = Matrix(raw[4:end, :])

        num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
        num_df = DataFrame(num, :auto)

        npairs = size(num_df, 2) ÷ 2

        pink  = RGB(1.0, 0.7, 0.8)
        pblue = RGB(0.2, 0.5, 1.0)
        cols = [RGB(pink.r + t*(pblue.r-pink.r),
                    pink.g + t*(pblue.g-pink.g),
                    pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]

        for j in 1:npairs
            xcol, ycol = 2j - 1, 2j
            lab = (j == 1) ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
            push!(labels, lab)

            x = num_df[!, xcol]
            y = num_df[!, ycol]
            line = lines!(ax, x, y; color = cols[j])
            push!(plot_objs, line)
        end

        Legend(fig[1, 2], plot_objs, labels, exp_title; framevisible = true)

    catch e
        if e isa UndefVarError
            # skip
        else
            rethrow(e)
        end
    end

    return fig
end

function plot_scanrate_sweeps(
    sweep_vec, scanrates;
    species=iohminus,
    fig_size=(1600, 900),
    scale=cm^2/mA,
    legend_title="Scan Rates (V/s)",
    highlight_index = nothing,
    highlight_color = RGB(1, 0.2, 0.2),
    highlight_lw = 4,
    default_lw = 1,
)
    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1], ylabel = lab_current, xlabel = lab_voltage)

    n = length(sweep_vec)
    cols = [RGB(0.2 + 0.6*(i/n),
                0.3 + 0.5*(1 - i/n),
                0.8 - 0.7*(i/n)) for i in 1:n]

    plot_objs = Any[]
    labels    = String[]

    for (j, rec) in enumerate(sweep_vec)
        push!(labels, "$(scanrates[j])\t\t ")

        color_j = cols[j]
        lw_j    = default_lw

        if highlight_index !== nothing && j == highlight_index
            color_j = highlight_color
            lw_j    = highlight_lw
        end

        I = currents(rec, species) .* scale
        line = lines!(ax, rec.voltages, I; linewidth=lw_j, color=color_j)
        push!(plot_objs, line)
    end

    Legend(fig[1, 2], plot_objs, labels, legend_title; framevisible=true)
    return fig
end

"""
let
	try
	    fig = Figure(size = (1600, 900))
	       ax = Axis(fig[1, 1],
	        xlabel = lab_voltage,
	        ylabel = lab_current,
	        yscale = log10,
	        yminorticksvisible = true,  
	        yminorticks = IntervalsBetween(5),
			#limits = ((-1.3, -0.7),(1e-4, 1e2))
	    )
	
	
	    # Experimental Data Plotting based on M.T.M Koper
	    raw = CSV.read("../data/Langmuir_CV_data/Figure_5.csv", DataFrame; header=false)
	    pres = vec(Matrix(raw[1:1, :]))
	    sub = Matrix(raw[4:end, :])
	    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
	    num_df = DataFrame(num, :auto)
	    npairs = size(num_df, 2) ÷ 2
	    pink, pblue = RGB(0.0, 0.7, 0.8), RGB(0.2, 0.5, 0.0)
	    cols1 = [RGB(pink.r + t*(pblue.r-pink.r),
	                 pink.g + t*(pblue.g-pink.g),
	                 pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]
	
	    plot_objs1 = []
	    labels1 = String[]
	    for j in 1:npairs
	        xcol, ycol = 2j - 1, 2j
	        label = j == 1 ? "\t\t sat" : "\t pCO2(atm)"
	        push!(labels1, label)
	        line = lines!(ax, num_df[!, xcol], (abs.(num_df[!, ycol])); color = cols1[j])
	        #line = lines!(ax, num_df[!, xcol], ((num_df[!, ycol])); color = cols1[j])
	        push!(plot_objs1, line)
	    end
	    Legend(fig[1, 2], plot_objs1, labels1, "Experimental"; framevisible = true)
	
	    # Theoretical Data Plotting based on `LiquidElectrolytes.jl`
	    plot_objs2 = []
	    labels2 = String[]
	    for (j, rec) in enumerate(F5_vec)
	        label2 = j == 1 ? "\t\t sat" : "t pCO2(atm)"
	        push!(labels2, label2)
	        line = lines!(ax, rec.voltages, (abs.(currents(rec, iohminus) .* cm^2/mA)); color = cols1[j])
			#line = lines!(ax, rec.voltages, ((currents(rec, iohminus) .* cm^2/mA)); color = cols2[j])
	        push!(plot_objs2, line)
	    end
	    Legend(fig[1, 3], plot_objs2, labels2, "Theoretical"; framevisible = true)
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
"""


function plot_conc_profile_with_delta(
    pnpresult, m, X;
    ispec::Int = 5,
    frac::Float64 = 0.99,
    t_indices = [1, 10, 50, 100, 140],
    t_index::Int = 140,
    fig_size = (800, 420),
)
    bulk = m.bulk
    tsol = pnpresult.tsol ./ (mol / dm^3)
    nvar, nx, nt = size(tsol)

    species = getproperty.(bulk, :name)
    colors  = getproperty.(bulk, :color)

    nspecies = min(nvar, length(species), length(colors))
    @assert 1 ≤ ispec ≤ nspecies

    nplot = min(nx, length(X))
    xx = X[1:nplot] ./ μm

    t_indices = unique(clamp.(vcat(t_indices, nt), 1, nt))
    t_index   = clamp(t_index, 1, nt)

    c_bulk = tsol[ispec, nplot, 1]
    if !(c_bulk > 0)
        error("c_bulk is not positive (c_bulk=$(c_bulk)). Choose a different ispec or far-field index.")
    end

    function x_at_frac(tsol, ispec, ti, xx, c_bulk, frac)
        c = vec(tsol[ispec, 1:length(xx), ti])
        target = frac * c_bulk
        idx = findfirst(ci -> (ci ≥ target), c)
        return isnothing(idx) ? NaN : xx[idx]
    end

    δs = Float64[]
    ts = Float64[]
    for ti in t_indices
        δ = x_at_frac(tsol, ispec, ti, xx, c_bulk, frac)
        push!(δs, δ)
        push!(ts, pnpresult.tsol.t[ti])
    end

    fig = Figure(size = fig_size)
    ax  = Axis(
        fig[1, 1],
        xlabel = "x / μm",
        ylabel = rich("log", subscript("10"), " ", rich("c", font = :italic), subscript("i"),
                      "  (mol/dm", superscript("3"), ")"),
        title  = "t = $(round(pnpresult.tsol.t[t_index], digits=4)) s  |  δ$(Int(round(frac*100))) for $(species[ispec])",
    )

    conc  = vec(tsol[ispec, 1:nplot, t_index])
    yvals = map(c -> (c > 0 ? log10(c) : NaN), conc)

    lines!(ax, xx, yvals; color = colors[ispec], linewidth = 2, label = species[ispec])

    δ_here = x_at_frac(tsol, ispec, t_index, xx, c_bulk, frac)
    if isfinite(δ_here)
        vlines!(ax, [δ_here]; linestyle = :dash, linewidth = 2)
        ytop = maximum(filter(isfinite, yvals))
        text!(ax, δ_here, ytop;
              text = "  δ$(Int(round(frac*100)))≈$(round(δ_here, digits=4)) μm",
              align = (:left, :top))
    end

    axislegend(ax; position = :rb)

    return (fig = fig, delta = δ_here, ts = ts, deltas = δs)
end

"""
let
	try
	    ic = model.cspecies
	    fig = Figure(size = (1050, 650))
	    ax = Axis(fig[1, 1], 
	  			  #limits = ((-1.25, 0.8),(-0.02, 0.1
										 
										 #)),
	              ylabel = lab_current,
	              xlabel = lab_voltage
	    )
	    colors = [RGB(i/3, 0, 1-(i/3)) for i in 1:3]
	
	    total_current = currents(pnpresult, iohminus) .* (cm^2/mA)
	    gold_line = lines!(ax, pnpresult.voltages, total_current,
	                       color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0))
	    labels1 = ["CO2RR Gold Model"]
	
	    raw_df = CSV.read("../data/Langmuir_CV_data/Figure_5.csv", DataFrame; header=false)
	    pH_row = collect(raw_df[1, :])
	    electrolyte_row = collect(raw_df[2, :])
	
	    numeric_data = [
	        parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
	        for i in 4:nrow(raw_df)
	    ]
	    num_df = DataFrame(hcat(numeric_data...)', names(raw_df))
	
	    conc_lines = [
	        lines!(ax, num_df[!, 1], num_df[!, 2], color = colors[1]),
	        lines!(ax, num_df[!, 3], num_df[!, 4], color = colors[2]),
	        lines!(ax, num_df[!, 5], num_df[!, 6], color = colors[3])
	    ]
	    labels2 = [
	        electrolyte_row[1]*"\t"*pH_row[2]*"pH",
	        electrolyte_row[3]*"\t\t"*pH_row[4]*"pH",
	        electrolyte_row[5]*"\t\t"*pH_row[6]*"pH"
	    ]
	
	    Legend(fig[1, 2],
	        [[gold_line], conc_lines], 
	        [labels1, labels2],         
	        ["Model", "Koper\nElectrolyte"];   
	    )
	
	    fig
	catch e
	   if e isa UndefVarError
			# normal case → skip
	   else
	        println("⚠️ Error occurred: ", e)
	        println(stacktrace(catch_backtrace()))
	    end
	end
end
"""

"""
plot_cv_over_L(results; species=ico, cutoff=-0.4, title="IV vs L")

Plot current vs voltage for each L in `results::Dict{Float64,Any}`.
Works with values that are either IVSweepResult-like or NamedTuple/struct with `ivresult`.
"""
function plot_cv_over_L(results::Dict{Float64,Any};
    species=ico, cutoff=-0.4, title="IV vs L"
)
    vis = GridVisualizer(;
        size   = (800, 500),
        title  = title,
        xlabel = lab_voltage,
        ylabel = lab_current,
        legend = :rt,
        #yscale = :log,
    )

    items = sort(collect(results); by=first)
    n = length(items)
    cols = Makie.resample_cmap(:cool, n)

    for (k, (L, rec)) in enumerate(items)
        ivres = hasproperty(rec, :ivresult) ? getproperty(rec, :ivresult) : rec

        volts = vec(ivres.voltages)
        Iall  = (vec(currents(ivres, species))) .* (cm^2/mA)

        m = min(length(volts), length(Iall))
        volts = volts[1:m]
        Iall  = Iall[1:m]

        mask = volts .< cutoff

        scalarplot!(vis,
            volts[mask],
            Iall[mask];
            clear=false,
            label="L = $(L) μm",
            color=cols[k],
        )
    end

    return reveal(vis)
end

function plot_pressure_varied_sweep(
    P_recs;
    species=iohminus,
    fig_size=(1600, 900),
    scale=cm^2/mA,
    limits=nothing,
    # --- add experimental background (Figure_3.csv) ---
    fig3_csv::Union{Nothing,AbstractString}="../data/Langmuir_CV_data/Figure_3.csv",
    exp_title::AbstractString="Experimental",
    exp_alpha::Real=0.55,
    exp_linewidth::Real=2,
    sim_linewidth::Real=3,
)
    fig = Figure(size = fig_size)

    ax = if limits !== nothing
        Axis(fig[1, 1],
             xlabel = lab_voltage,
             ylabel = lab_current,
             limits = limits)
    else
        Axis(fig[1, 1],
             xlabel = lab_voltage,
             ylabel = lab_current)
    end

    plots  = Any[]
    labels = String[]

    # ------------------------------------------------------------
    # 1) Experimental background (Figure_3.csv)  [optional]
    # ------------------------------------------------------------
    if fig3_csv !== nothing
        try
            raw = CSV.read(fig3_csv, DataFrame; header=false)

            pres = vec(Matrix(raw[1:1, :]))
            sub  = Matrix(raw[4:end, :])

            num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
            num_df = DataFrame(num, :auto)

            npairs = size(num_df, 2) ÷ 2

            pink  = RGB(1.0, 0.7, 0.8)
            pblue = RGB(0.2, 0.5, 1.0)
            cols_exp = [RGB(pink.r + t*(pblue.r-pink.r),
                            pink.g + t*(pblue.g-pink.g),
                            pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]

            for j in 1:npairs
                xcol, ycol = 2j - 1, 2j
                lab = (j == 1) ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"

                x = num_df[!, xcol]
                y = num_df[!, ycol]

                # lighter background curves
                line = lines!(ax, x, y; color = (cols_exp[j], exp_alpha), linewidth = exp_linewidth)
                push!(plots, line)
                push!(labels, "$exp_title | $lab")
            end
        catch e
            # keep behavior: only skip for UndefVarError, otherwise rethrow
            if e isa UndefVarError
                # skip
            else
                rethrow(e)
            end
        end
    end

    # ------------------------------------------------------------
    # 2) Simulation curves (existing logic)
    # ------------------------------------------------------------
    n = length(P_recs)
    cols_sim = [RGB(1 - t, 0, t) for t in LinRange(0, 1, max(n, 1))]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        I = currents(rec, species) .* scale

        line = lines!(ax, rec.voltages, I; color = cols_sim[j], linewidth = sim_linewidth)
        push!(plots, line)
        push!(labels, "Theoretical | $label")
    end

    Legend(fig[1, 2], plots, labels, "Overlay"; framevisible=true)
    return fig
end

function pressure_varied_cvsweep(
    P_recs;
    species=iohminus,
    fig_size=(1600, 900),
    scale=cm^2/mA,
    limits=nothing,
)
    fig = Figure(size = fig_size)
    sim_linewidth = 3
    ax = if limits !== nothing
        Axis(fig[1, 1],
             xlabel = lab_voltage,
             ylabel = lab_current,
             limits = limits)
    else
        Axis(fig[1, 1],
             xlabel = lab_voltage,
             ylabel = lab_current)
    end

    plots  = Any[]
    labels = String[]

    # ------------------------------------------------------------
    # 2) Simulation curves (existing logic)
    # ------------------------------------------------------------
    n = length(P_recs)
    cols_sim = [RGB(1 - t, 0, t) for t in LinRange(0, 1, max(n, 1))]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        I = currents(rec, species) .* scale

        line = lines!(ax, rec.voltages, I; color = cols_sim[j], linewidth = sim_linewidth)
        push!(plots, line)
        push!(labels, "Theoretical | $label")
    end

    Legend(fig[1, 2], plots, labels, "Overlay"; framevisible=true)
    return fig
end


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 1: CV current / scan-rate plotting family.
#
# Species index constants for the Gold CO2RR model, added so that the bare
# `iohminus`/`ico`/`ico2` defaults and bodies (here and in the existing
# functions above) resolve inside the AuCO2RR_plots module.
# =====================================================================
const ikplus   = 1
const ihplus   = 2
const ihco3    = 3
const ico3     = 4
const ico2     = 5
const iohminus = 6
const ico      = 7

# ── moved from cell 7bfe397e (CV_dsp_cap_result) ──
function CV_dsp_cap_result(result, m; scale=cm^2/mA, show_diff=false)
    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(540, 440))
        a = Axis(f[1, 1],
                 ylabel = lab_current,
                 xlabel = lab_voltage)
        return f, a
    end

    j_dsp = result.j_dsp .* scale
    j_cap = result.j_cap .* scale

    lines!(ax, result.voltages, j_dsp; color = :magenta, label = rich(rich("j", font = :italic), subscript("dsp")))
    lines!(ax, result.voltages, j_cap; color = :skyblue, label = rich(rich("j", font = :italic), subscript("cap")))

    if show_diff
        lines!(ax, result.voltages, j_dsp .- j_cap;
               color = :orange, linestyle = :dash, label = rich(rich("j", font = :italic), subscript("dsp"), " − ", rich("j", font = :italic), subscript("cap")))
    end

    axislegend(ax; position = :rt)
    return fig
end

# ── moved from cell ad3a5236 (CV_total_current) ──
function CV_total_current(result, m; species=nothing, scale=cm^2/mA)
    model = m.elydata
    sp = (species === nothing) ? model.cspecies[1] : species

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(540, 440))
        a = Axis(f[1, 1],
                 ylabel = lab_current,
                 xlabel = lab_voltage)
        return f, a
    end

    #i_F   = currents(result, sp) .* scale
    i_cap = result.j_cap         .* scale
    #i_tot = i_F .+ i_cap

   # lines!(ax, result.voltages, i_F;   color = :magenta,  label = L"i_F")
    lines!(ax, result.voltages, i_cap; color = :skyblue,  label = rich(rich("i", font = :italic), subscript("cap")))
   # lines!(ax, result.voltages, i_tot; color = :orange,   label = L"i_{tot}")
    axislegend(ax; position = :rt)
    return fig
end

# ── moved from cell 9949de26 (plot_cv_total_current_tot) ──
function plot_cv_total_current_tot(result, m;
                                   redox_species::Dict{Int,Int},
                                   co_idx,
                                   include_capacitive::Bool = true,
                                   scale = cm^2/mA,
                                   sign::Int = 1,
                                   color_F   = colorant"#F2728A",   # i_F color
                                   color_C   = colorant"#5BA8E8",   # i_C color
                                   mix_mode::Symbol = :mean,         # :sum or :mean
                                   lw = 5)
    electrolyte = m.elydata
    n_t = length(result.voltages)

    # ---- Faradaic current (total, sum over redox species) ----
    I_F = zeros(n_t)
    for (idx, n_e) in redox_species
        I_F .+= n_e .* currents(result, idx)
    end

    # ---- Capacitive current ----
    I_C = zeros(n_t)
    if include_capacitive
        if electrolyte.ircompensation == :ohmicdrop
            icc = electrolyte.icc
            node_we = 1
            I_C = [u[icc, node_we] for u in result.tsol[1:n_t]]
        else
            @warn "Capacitive current only resolved in :ohmicdrop mode " *
                  "(current mode: :$(electrolyte.ircompensation)); j_C set to 0."
        end
    end

    # ---- Scale and sign ----
    I_F_scaled = sign .* I_F          .* scale
    I_C_scaled = sign .* I_C          .* scale
    I_total    = sign .* (I_F .+ I_C) .* scale

    # ---- blended color: RGB average or sum ----
    cF = RGBf(color_F); cC = RGBf(color_C)
    color_tot = mix_mode === :mean ?
        RGBf((cF.r+cC.r)/2, (cF.g+cC.g)/2, (cF.b+cC.b)/2) :
        RGBf(min(cF.r+cC.r, 1f0), min(cF.g+cC.g, 1f0), min(cF.b+cC.b, 1f0))

    # ---- Plot ----
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (900, 800))
        ax_a = Axis(f[1, 2], ylabel = rich(rich("i", font = :italic), subscript("F"), "  (mA cm", superscript("−2"), ")"))
        ax_b = Axis(f[2, 2], ylabel = rich(rich("i", font = :italic), subscript("C"), "  (mA cm", superscript("−2"), ")"))
        ax_c = Axis(f[3, 2], ylabel = rich(rich("i", font = :italic), subscript("tot"), "  (mA cm", superscript("−2"), ")"),
                             xlabel = lab_voltage)

        hidexdecorations!(ax_a; grid = false)
        hidexdecorations!(ax_b; grid = false)
        linkxaxes!(ax_a, ax_b, ax_c)
        return f, (ax_a, ax_b, ax_c)
    end
    f, (ax_a, ax_b, ax_c) = fig

    # ---- panel labels (a), (b), (c) ----
    sub_axes = [ax_a, ax_b, ax_c]
    labels = ["(a)", "(b)", "(c)"]

    for i in 1:3
        Label(f[i, 1], labels[i],
              fontsize = 24,
              font = :bold,
              halign = :left,
              valign = :top,
              padding = (15, 0, 0, 15)
        )
    end

    # ---- draw data lines ----
    lines!(ax_a, result.voltages, I_F_scaled; color = color_F,   linewidth = lw)
    lines!(ax_b, result.voltages, I_C_scaled; color = color_C,   linewidth = lw)
    lines!(ax_c, result.voltages, I_total;    color = color_tot, linewidth = lw)
        rowgap!(f.layout, 15)

        rowsize!(f.layout, 1, Relative(0.30))
        rowsize!(f.layout, 2, Relative(0.30))
        rowsize!(f.layout, 3, Relative(0.40))
    return f
end

# ── moved from cell 02a78c8d (plot_cv_scanrate_grid) ──
function plot_cv_scanrate_grid(result_vec, m;
                               redox_species::Dict{Int,Int},
                               co_idx,
                               scanrates = [0.05, 0.5, 5.0],
                               include_capacitive::Bool = true,
                               scale = cm^2/mA,
                               sign::Int = 1,
                               color_F   = colorant"#F2728A",
                               color_C   = colorant"#5BA8E8",
                               mix_mode::Symbol = :mean,
                               lw = 5.5)
    electrolyte = m.elydata

    unit_i = rich("  (mA cm", superscript("−2"), ")")
    lab_iF   = rich(rich("I", font=:italic), subscript("F"),   unit_i)
    lab_iC   = rich(rich("I", font=:italic), subscript("C"),   unit_i)
    lab_itot = rich(rich("I", font=:italic), subscript("tot"), unit_i)
    lab_U    = rich(rich("U", font=:italic), "  (V vs. SHE)")

    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 700))
        axes_matrix = Matrix{Axis}(undef, 3, 3)
        ylabels = [lab_iF, lab_iC, lab_itot]

        for col in 1:3, row in 1:3
            axes_matrix[row, col] = row == 3 ?
                Axis(f[row, col+1], xlabel = lab_U) :
                Axis(f[row, col+1])

            row < 3 && hidexdecorations!(axes_matrix[row, col]; grid = false)

            if col == 1
                axes_matrix[row, col].ylabel = ylabels[row]
                axes_matrix[row, col].ylabelpadding = 30
            else
                hideydecorations!(axes_matrix[row, col]; ticks = false, grid = false)
            end

            if row == 2
                axes_matrix[row, col].ytickformat = "{:.2f}"
                axes_matrix[row, col].yticks = LinearTicks(3)
            else
                axes_matrix[row, col].yticks = LinearTicks(4)
            end
        end

        for i in 1:3
            linkyaxes!(axes_matrix[i, 1], axes_matrix[i, 2], axes_matrix[i, 3])
            linkxaxes!(axes_matrix[1, i], axes_matrix[2, i], axes_matrix[3, i])
        end

        return f, axes_matrix
    end
    f, axes_matrix = fig

    grid_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"]
    for row in 1:3, col in 1:3
        idx = (row - 1) * 3 + col
        Label(f[row, col+1], grid_labels[idx],
              fontsize = 24, font = :bold, halign = :left, valign = :top,
              padding = (15, 0, 0, 15))
    end

    for col in 1:3
        Label(f[0, col+1],
              rich(string(scanrates[col]), "  V s", superscript("−1"));
              fontsize = 26, font = :bold, halign = :center, valign = :bottom)
    end

    cF, cC = RGBf(color_F), RGBf(color_C)
    color_tot = mix_mode === :mean ?
        RGBf((cF.r + cC.r)/2, (cF.g + cC.g)/2, (cF.b + cC.b)/2) :
        RGBf(min(cF.r + cC.r, 1f0), min(cF.g + cC.g, 1f0), min(cF.b + cC.b, 1f0))

    for col in 1:3
        res = result_vec[col]
        n_t = length(res.voltages)

        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(res, idx)
        end

        I_C = zeros(n_t)
        if include_capacitive && electrolyte.ircompensation == :ohmicdrop
            icc = electrolyte.icc
            I_C = [u[icc, 1] for u in res.tsol[1:n_t]]
        end

        I_F_scaled = sign .* I_F          .* scale
        I_C_scaled = sign .* I_C          .* scale
        I_total    = sign .* (I_F .+ I_C) .* scale

        lines!(axes_matrix[1, col], res.voltages, I_F_scaled; color = color_F,   linewidth = lw)
        lines!(axes_matrix[2, col], res.voltages, I_C_scaled; color = color_C,   linewidth = lw)
        lines!(axes_matrix[3, col], res.voltages, I_total;    color = color_tot, linewidth = lw)
    end

    rowgap!(f.layout, 15)
    colgap!(f.layout, 25)

    rowsize!(f.layout, 1, Relative(0.29))
    rowsize!(f.layout, 2, Relative(0.29))
    rowsize!(f.layout, 3, Relative(0.38))

    colsize!(f.layout, 1, Auto())
    for col in 2:4
        colsize!(f.layout, col, Relative(0.34))
    end

    return f
end

# ── moved from cell 08756476 (plot_cv_scanrate_grid_unc) ──
function plot_cv_scanrate_grid_unc(result_vec, m;
                               redox_species::Dict{Int,Int},
                               co_idx,
                               scanrates = [0.05, 0.5, 5.0],
                               include_capacitive::Bool = true,
                               scale = cm^2/mA,
                               sign::Int = 1,
                               color_tot = "#BAC8FF",
                               lw = 5.5)
    electrolyte = m.elydata

    lab_I = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")")
    lab_U = rich(rich("U", font=:italic), "  (V vs. SHE)")

    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 350))
        axes_vec = Vector{Axis}(undef, 3)

        for col in 1:3
            axes_vec[col] = Axis(f[1, col+1], xlabel = lab_U)

            if col == 1
                axes_vec[col].ylabel = lab_I
                axes_vec[col].ylabelpadding = 30
            else
                hideydecorations!(axes_vec[col]; ticks = false, grid = false)
            end
            axes_vec[col].yticks = LinearTicks(4)
        end

        linkyaxes!(axes_vec[1], axes_vec[2], axes_vec[3])

        return f, axes_vec
    end
    f, axes_vec = fig

    grid_labels = ["(a)", "(b)", "(c)"]
    for col in 1:3
        Label(f[1, col+1], grid_labels[col],
              fontsize = 24, font = :bold, halign = :left, valign = :top,
              padding = (15, 0, 0, 15))
    end

    for col in 1:3
        Label(f[0, col+1],
              rich(string(scanrates[col]), "  V s", superscript("−1"));
              fontsize = 26, font = :bold, halign = :center, valign = :bottom)
    end

    for col in 1:3
        res = result_vec[col]
        n_t = length(res.voltages)

        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(res, idx)
        end

        I_C = zeros(n_t)
        if include_capacitive && electrolyte.ircompensation == :ohmicdrop
            icc = electrolyte.icc
            I_C = [u[icc, 1] for u in res.tsol[1:n_t]]
        end

        I_total = sign .* (I_F .+ I_C) .* scale
        lines!(axes_vec[col], res.voltages, I_total; color = color_tot, linewidth = lw)
    end

    colgap!(f.layout, 25)
    rowsize!(f.layout, 1, Relative(0.98))
    colsize!(f.layout, 1, Auto())
    for col in 2:4
        colsize!(f.layout, col, Relative(0.34))
    end
    return f
end

# ── moved from cell d3b7d864 (plot_scanrate_sweeps_cv_2) ──
function plot_scanrate_sweeps_cv_2(
    sweep_vec, scanrates, m;
    redox_species::Dict{Int,Int},       # all Faradaic species and their electron counts (e.g. Dict(1=>2, 2=>1))
    include_capacitive::Bool = true,
    scale = cm^2/mA,
    sign::Int = 1,
    default_lw = 2
)
    electrolyte = m.elydata
    n = length(sweep_vec)
    pastel2 = cgrad([colorant"#D5F011", colorant"#16D8FF"])
    cols = [pastel2[t] for t in range(0, 1, length=max(n, 1))]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (850, 550))
        a = Axis(f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        n_t = length(rec.voltages)

        # 1. sum Faradaic currents
        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(rec, idx)
        end

        # 2. add capacitive current
        I_C = zeros(n_t)
        if include_capacitive
            if electrolyte.ircompensation == :ohmicdrop
                icc = electrolyte.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in rec.tsol[1:n_t]]
            else
                j == 1 && @warn "Capacitive current only resolved in :ohmicdrop mode; j_C set to 0."
            end
        end

        # 3. total current (apply sign and scale)
        I_total = sign .* (I_F .+ I_C) .* scale

        lines!(ax, rec.voltages, I_total;
            linewidth = default_lw,
            color     = cols[j]
        )

        pos_y = if j == 5
            -2.5
        elseif j == 3
            -1.5
        else
            -0.5 * j
        end

        text!(ax, "$(scanrates[j]) V/s";
            position = (-1.55, pos_y),
            color    = cols[j],
            fontsize = 18,
            font     = :bold
        )
    end

    return fig
end

# ── moved from cell ae50302f (plot_scanrate_sweeps_cv) ──
function plot_scanrate_sweeps_cv(
    sweep_vec, scanrates;
    species=iohminus,
    scale=cm^2/mA,
    default_lw = 4
)
    n = length(sweep_vec)
    pastel2 = cgrad([colorant"#D5F011", colorant"#16D8FF"])
    cols = [pastel2[t] for t in range(0, 1, length=max(n, 1))]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (850, 550))
        a = Axis(f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        I = currents(rec, species) .* scale
        lines!(ax, rec.voltages, I;
            linewidth = default_lw,
            color     = cols[j]
        )
    end

    return fig
end

# ── moved from cell 51064823 (plot_scanrate_sweeps_cap) ──
function plot_scanrate_sweeps_cap(
    sweep_vec, scanrates;
    #species=iohminus,
    scale=cm^2/mA,
    default_lw = 4
)
    n = length(sweep_vec)

    cols = [RGBf(0.1 + 0.6*(i/n),
                 0.2 + 0.6*(1 - i/n),
                 0.7 - 0.3*(i/n)) for i in 1:n]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 550))
        a = Axis(f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        color_j = cols[j]
        lw_j    = default_lw
        J_cap = rec.j_cap .* scale
        lines!(ax, rec.voltages, rec.j_cap; linewidth=lw_j, color=color_j)
    end

    if n > 1
        cgradient = cgrad(cols, categorical = true)

        Colorbar(fig[1, 2],
            limits = (minimum(scanrates), maximum(scanrates)),
            colormap = cgradient,
            label = rich("Scan rate  (V/s)"),
            labelsize = 20,
            ticklabelsize = 18,
            spinewidth = 2.0,
            width = 20
        )
    end

    return fig
end


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 2: contour / multi-panel / overlay family.
# Helpers `powlab` / `lab_time` are moved here (from the panel begin-block).
# Functions that referenced notebook globals `elydata_Gold_unc/odr` were
# parameterized to use the passed `model` / a new `electrolyte` kwarg so they
# resolve inside the module; behavior is preserved at the notebook call sites.
# =====================================================================

# `powlab` and `lab_time` are defined once in plots/struct.jl (shared style helpers).

# ── moved from cell 8e2f8c2c (plot_co2_profiles) ──
function plot_co2_profiles(
    m,
    L_values,
    resL_COMP,
    X_coords_dict;
    L_varied_checkbox::Bool = true
)
    bulk = m.bulk
    if !L_varied_checkbox
        return nothing
    end

    co2_idx = findfirst(s -> s.name == "CO₂", bulk)

    selected_Ls = L_values[1:end]
    num_panels = length(selected_Ls)
    max_L = 1e-2

    fig = Figure(size = (300, 300 * num_panels))

    c_min = log10(0.00001)
    c_max = log10(0.33)
    color_range = (c_min, c_max)

    num_levels = 12
    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    for (idx, L_val) in enumerate(selected_Ls)
        if !haskey(resL_COMP, L_val)
            continue
        end

        res = resL_COMP[L_val]
        tsol = res.tsol
        times = res.times
        X_coords = X_coords_dict[L_val]

        conc_matrix = [(tsol[co2_idx, ix, it] / (mol / dm^3))
                       for ix in 1:length(X_coords), it in 1:length(times)]

        ax = Axis(fig[idx, 1],
            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
            ylabel = "Time [s]",
            title = "Boundary Layer (L) = $(round(L_val/μm)) μm",
            xscale = log10,
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(9)
        )

        hm = heatmap!(ax, X_coords .+ 1e-12, times, conc_matrix;
            colorrange = color_range,
            colormap = discrete_cmap,
            interpolate = false
        )

        vlines!(ax, [L_val], color = :red, linestyle = :dash, linewidth = 2)

        Colorbar(fig[idx, 2], hm, label = rich("log", subscript("10"), "(", rich("c", font = :italic), subscript(rich("CO", subscript("2"))), ")"))
    end

    rowgap!(fig.layout, 35)

    return fig
end

# ── moved from cell d91bcc6a (co2_log_contour) ──
function co2_log_contour(results::Dict, grid_dict, m;
                         scale=mol/dm^3, num_levels=24, L_val=nothing)
    bulk = m.bulk

    if isnothing(L_val)
        L_val = minimum(keys(results))
    end

    result = results[L_val]
    grid = grid_dict[L_val]
    X = grid[Coordinates][1, :]

    co2_idx = findfirst(s -> s.name == "CO₂", bulk)
    times = result.tsol.t

    log_c_bulk = log10(result.tsol[co2_idx, end, 1] / scale)

    M = [log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / scale, 1e-12))
         for ix in 1:length(X), it in 1:length(times)]

    c_min = 0.0
    c_max = log10(result.tsol[co2_idx, end, 1] / scale) - log10(1e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    f = Figure(size = (800, 450))

    ax = Axis(f[1, 1];
        xlabel = lab_time,
        ylabel = rich(rich("x", font = :italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (10.0 .^ (-12:3:-6), [powlab(-12), powlab(-9), powlab(-6)])
    )

    hm = heatmap!(ax, times, X .+ 1e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false)

    cb = Colorbar(f[1, 2], hm;
                  label = rich("log", subscript("10"), "(", rich("c", font = :italic), subscript("bulk"),
                               ") − log", subscript("10"), "(", rich("c", font = :italic), subscript(rich("CO", subscript("2"))), ")"),
                  ticklabelsize = 20,
                  labelsize = 20)

    return f
end

# ── moved from cell e6f43f01 (plot_activity_time_electrode) ──
function plot_activity_time_electrode(result, m;
                                     model_type="DGML_γ", # "DGML_γ" or "Stefan_γ"
                                     nspecies=7,
                                     scale=(mol/dm^3),
                                     ipressure=nothing)
    bulk = m.bulk
    electrolyte = m.elydata
    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    activity_electrode = fill(NaN, nspecies, nt)
    ielectrode = 1

    v0 = electrolyte.v0
    bar_c = 1.0 / v0
    RT = electrolyte.RT
    c_scale = 1.0 / scale

    for t in 1:nt
        c_all = result.tsol[:, ielectrode, t]

        Phi = sum(c_all[ic] * electrolyte.v[ic] for ic in 1:nspecies)
        solvent_frac = max(1.0 - Phi, eps(Float64))

        pnode = (ipressure !== nothing) ? result.tsol[ipressure, ielectrode, t] : 0.0

        for i in 1:nspecies
            c = c_all[i]
            v_i = electrolyte.v[i]
            size_ratio = v_i / v0
            term_conc = c / bar_c

            if model_type == "DGML_γ"
                term_press  = exp((1.0 - size_ratio) * pnode / (bar_c * RT))
                term_steric = solvent_frac^(-size_ratio)
                a_thermo    = term_conc * term_press * term_steric
            elseif model_type == "Stefan_γ"
                a_thermo    = term_conc * (solvent_frac^(-1.0))
            else
                a_thermo    = term_conc
            end

            activity_electrode[i, t] = a_thermo * (bar_c * c_scale)
        end
    end

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (960, 540))

        model_title = model_type == "DGML_γ" ? "Modified DGML Model (Activity)" :
                      model_type == "Stefan_γ" ? "Stefan Model (Activity)" : "Ideal Solution"

        a = Axis(f[1, 1],
            title = model_title,
            xlabel = lab_time,
            ylabel = rich(rich("a", font = :italic), subscript("i, electrode")),
            limits = ((times[1] - (times[end] / 200), times[end] + (times[end] / 100)), (1e-12, 1e4)),
            yscale = log10
        )

        yt_vals = 10.0 .^ (4:-4:-12)
        yt_lbls = [powlab(4), powlab(0), powlab(-4), powlab(-8), powlab(-12)]
        a.yticks = (yt_vals, yt_lbls)

        return f, a
    end

    for i in 1:nspecies
        y = activity_electrode[i, :]
        y_fixed = [isnan(val) || val <= eps(Float64) ? eps(Float64) : val for val in y]

        lines!(ax, times, y_fixed; color=colors[i])
    end

    return fig
end

# ── moved from cell 96c90e9f (CV_overlay_currents) ──
function CV_overlay_currents(results, m;
                             labels   = nothing,
                             species  = nothing,
                             scale    = cm^2/mA,
                             linestyles = [:solid, :dash, :dot])
    model = m.elydata
    sp = (species === nothing) ? model.cspecies[1] : species
    labels = labels === nothing ? ["result $i" for i in 1:length(results)] : labels

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size=(620, 460))
        a = Axis(f[1, 1],
                 ylabel = lab_current,
                 xlabel = lab_voltage)
        return f, a
    end

    col_F   = :magenta
    col_cap = :skyblue
    col_tot = :orange

    for (k, result) in enumerate(results)
        ls = linestyles[mod1(k, length(linestyles))]

        i_F   = currents(result, sp) .* scale
        i_cap = result.j_cap         .* scale
        i_tot = i_F .+ i_cap

        lines!(ax, result.voltages, i_F;   color = col_F,   linestyle = ls,
               label = rich(rich("i", font = :italic), subscript("F"), " (", string(labels[k]), ")"))
        lines!(ax, result.voltages, i_cap; color = col_cap, linestyle = ls,
               label = rich(rich("i", font = :italic), subscript("cap"), " (", string(labels[k]), ")"))
        lines!(ax, result.voltages, i_tot; color = col_tot, linestyle = ls,
               label = rich(rich("i", font = :italic), subscript("tot"), " (", string(labels[k]), ")"))
    end

    axislegend(ax; position = :rt, nbanks = 2)
    return fig
end

# ── moved from cell a2264942 (plot_cv_current_variedL) ──
function plot_cv_current_variedL(results, m;
                                 species = nothing,
                                 scale = cm^2/mA,
                                 color_gradient = true,
                                 linewidth = 3.5,
                                 title = "")
    model = m.elydata
    sp = (species === nothing) ? model.cspecies[1] : species

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 400))
        a = Axis(f[1, 1];
                 title  = title,
                 ylabel = lab_current,
                 xlabel = lab_voltage)
        return f, a
    end

    Lkeys = sort(collect(keys(results)))
    n = length(Lkeys)

	cols = [RGBf(
	    0.6 + 0.3 * (i-1)/max(n-1,1),
	    0.8 - 0.2 * (i-1)/max(n-1,1),
	    0.7 - 0.2 * (i-1)/max(n-1,1)
	) for i in 1:n]


    for (i, L) in enumerate(Lkeys)
        I = currents(results[L], sp) .* scale ./ 2
        lines!(ax, results[L].voltages, I;
               color = cols[i], linewidth = linewidth,
               label = @sprintf("%g μm", L / μm))
    end

    axislegend(ax, "L", position = :rb, framevisible = false, legendtext = 24)
    return fig
end

# ── moved from cell ed4451bd (plot_pressure_varied_sweep_ivc) ──
function plot_pressure_varied_sweep_ivc(
    P_recs;
    species=iohminus,
    scale=cm^2/mA,
    limits=nothing,
    sim_linewidth::Real=5
)
    n = length(P_recs)
	pastel3 = cgrad([colorant"#FFB3BA", colorant"#A3D8FF"])
	cols_sim = [pastel3[t] for t in range(0, 1, length=max(n, 1))]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1050, 500))
        a = if limits !== nothing
            Axis(f[1, 1], xlabel = lab_voltage, ylabel = lab_current, limits = limits)
        else
            Axis(f[1, 1], xlabel = lab_voltage, ylabel = lab_current)
        end
        return f, a
    end

    plot_objs = []
    labels    = String[]

    for j in 1:n
        p, rec = P_recs[j]
        I = currents(rec, species) .* scale

        label_text = "$(p) atm"

        hl = lines!(ax, rec.voltages, I ./ 2;
            color = cols_sim[j],
            linewidth = sim_linewidth)

        push!(plot_objs, hl)
        push!(labels, label_text)
    end

    if n > 0
        leg = Legend(fig[1, 1], plot_objs, labels, rich(rich("p", font = :italic), subscript(rich("CO", subscript("2"))));
            framevisible = false,
            halign = :right, valign = :bottom,
            labelsize = 20, titlesize = 23,
            padding = (0, 0, 0, 0),
            tellwidth = false, tellheight = false)

        translate!(leg.blockscene, -40, 40, 0)
    end

    return fig
end

# ── moved from cell 6f6e779c (plot_cv_current_dict) ──
function plot_cv_current_dict(result_dict, m;
                              species = nothing,
                              scale = 1.0,
                              title = "",
                              linewidth = 3)
    model = m.elydata
    sp = (species === nothing) ? model.cspecies[1] : species

    sorted_keys = sort(collect(keys(result_dict)))

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (400, 300))
        a = Axis(f[1, 1];
                 title  = title,
                 ylabel = lab_current,
                 xlabel = lab_voltage)
        return f, a
    end

    for (i, key) in enumerate(sorted_keys)
        result = result_dict[key]
        I = currents(result, sp) .* scale ./ 2
        lines!(ax, result.voltages, I;
              linewidth = linewidth,
               label = @sprintf("L = %g μm", key / μm))
    end

    axislegend(ax, position = :rb, framevisible = true)
    return fig
end

# ── moved from cell 13a0e5da (panel_conc_time!) ──
function panel_conc_time!(fig, panel_pos, result, m; nspecies=7, scale=mol/dm^3, lw=4,
                          xlabel = lab_time, legend_pos = nothing)
    bulk = m.bulk
    sp_colors = ["#E07B39","#888888","#7B5C3E","#222222",
                 "#C0392B","#27AE60","#2980B9"]
    colors = sp_colors[1:nspecies]
    names  = getproperty.(bulk, :name)
    times  = result.tsol.t
    nt     = length(times)
    conc   = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    ax_c = Axis(panel_pos;
                xlabel = xlabel,
                ylabel = rich(rich("c", font=:italic),
                              subscript(rich("α", font=:italic)),
                              superscript("‡"), "  (M)"),
                yscale = log10,
                yminorticksvisible = true,
                yminorticks = IntervalsBetween(9)
               )

    ax_c.yticks = (10.0 .^ (4:-4:-12),
                   [powlab(4), powlab(0), powlab(-4), powlab(-8), powlab(-12)])

    for i in 1:nspecies
        lines!(ax_c, times, max.(conc[i, :], eps(Float64));
               color = colors[i], linewidth = lw)
    end

    legend_elements = [ [LineElement(color = colors[i], linewidth = lw)] for i in 1:nspecies ]
    legend_labels   = [ rich(string(names[i]), color = colors[i]) for i in 1:nspecies ]

    leg = Legend(legend_pos === nothing ? fig[1, 3] : legend_pos,
                 legend_elements, legend_labels;
                 framevisible = false,
                 labelsize    = 20)

    return ax_c, leg
end

# ── moved from cell 13a0e5da (panel_time_current!)  [ely = elydata_Gold_unc → ely = model] ──
function panel_time_current!(fig, panel_pos, result, m;
                             redox_species = nothing,
                             include_capacitive = true,
                             scale = cm^2/mA,
                             sign = -1,
                             color_gradient = true,
                             lw = 4,
                             xlabel = lab_time)
    model = m.elydata

    ax = Axis(panel_pos; xlabel = xlabel,
              ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"))

    redox = Dict(model.cspecies[ico2] => 2)

    n_t = length(result.times)
    I_F = zeros(n_t)
    for (idx, n_e) in redox
        I_F .+= n_e .* currents(result, idx)
    end

    I_C = zeros(n_t)
    if include_capacitive
        ely = model
        if ely.ircompensation == :ohmicdrop
            icc = ely.icc
            node_we = 1
            I_C = [u[icc, node_we] for u in result.tsol[1:end-1]]
        end
    end

    I = sign .* (I_F .+ I_C) .* scale

    cols = color_gradient ? :skyblue : :black
    lines!(ax, result.times, I; color = cols, linewidth = lw)

    return ax
end

# ── moved from cell 13a0e5da (panel_time_voltage!) ──
# Pass `sawtooth` to plot the APPLIED protocol U_we(t) (the full vmin..vmax ramp).
# Without it, `result.voltages` is plotted, which here is the reaction-plane
# (electrode node) potential — compressed by the gap capacitance, not the sawtooth.
function panel_time_voltage!(fig, panel_pos, result; sawtooth = nothing, lw = 5, xlabel = lab_time)
    ax = Axis(panel_pos; xlabel = xlabel,
              ylabel = rich(rich("U", font = :italic), "  (V vs. SHE)"))
    U = sawtooth === nothing ? result.voltages : sawtooth.(result.times)
    lines!(ax, result.times, U;
           color = parse(Colorant, "#D7C2F0"), linewidth = lw)
    return ax
end

# ── moved from cell 13a0e5da (panel_co2_log_contour!) ──
function panel_co2_log_contour!(fig, panel_pos, cbar_pos, result, X, m;
                                 scale=mol/dm^3, num_levels=24, L_val=nothing)
    bulk = m.bulk
    co2_idx = findfirst(s -> s.name == "CO₂", bulk)
    times = result.tsol.t

    log_c_bulk = log10(result.tsol[co2_idx, end, 1] / scale)

    M = [log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / scale, 1e-12))
         for ix in 1:length(X), it in 1:length(times)]

    c_min = 0.0
    c_max = log10(result.tsol[co2_idx, end, 1] / scale) - log10(1e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    ax = Axis(panel_pos;
        xlabel = lab_time,
        ylabel = rich(rich("x", font=:italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (10.0 .^ (-12:3:-6),
          [powlab(-12), powlab(-9), powlab(-6)])
             )

    hm = heatmap!(ax, times, X .+ 1e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false)

    Colorbar(cbar_pos, hm;
        label = rich("log", subscript("10"), "(", rich("c", font=:italic),
                     subscript("bulk"), ") − log", subscript("10"), "(",
                     rich("c", font=:italic),
                     subscript(rich("CO", subscript("2"))), ")"),
        ticklabelsize = 20, labelsize = 20)
    return ax, hm
end

# ── moved from cell 13a0e5da (panel_time_ph!) ──
function panel_time_ph!(fig, panel_pos, result; ihplus=2, scale=mol/dm^3, lw=5,
                        xlabel = lab_time, color = parse(Colorant, "#7BB661"))
    times = result.tsol.t
    nt    = length(times)
    cH    = [result.tsol[ihplus, 1, t] / scale for t in 1:nt]   # surface (node 1) H⁺ concentration (M)
    pH    = -log10.(max.(cH, eps(Float64)))
    ax = Axis(panel_pos; xlabel = xlabel, ylabel = rich("pH"))
    lines!(ax, times, pH; color = color, linewidth = lw)
    return ax
end

# ── moved from cell ac2bef83 (plot_7species_contours) ──
function plot_7species_contours(result, X, m; scale=mol/dm^3, num_levels=24)
    bulk = m.bulk
    target_species = [
        ("K⁺",    rich("K", superscript("+"))),
        ("H⁺",    rich("H", superscript("+"))),
        ("CO₂",   rich("CO", subscript("2"))),
        ("OH⁻",   rich("OH", superscript("−"))),
        ("HCO₃⁻", rich("HCO", subscript("3"), superscript("−"))),
        ("CO₃²⁻", rich("CO", subscript("3"), superscript("2−"))),
        ("CO",    rich("CO", superscript("−")))
    ]
    fig = Figure(size = (750, 1500))
    times = result.tsol.t
    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    for i in 1:7
        sp_name, sp_label = target_species[i]

        sp_idx = findfirst(s -> s.name == sp_name, bulk)
        if isnothing(sp_idx)
            @warn "Species '$sp_name' not found in bulk; skipping this panel."
            continue
        end

        @views c_matrix = result.tsol[sp_idx, 1:length(X), 1:length(times)] ./ scale

        M = log10.(max.(c_matrix, 1e-12))

        replace!(M, Inf => -12.0, -Inf => -12.0, NaN => -12.0)

        c_min = minimum(M)
        c_max = maximum(M)
        if c_min == c_max
            c_min -= 0.5
            c_max += 0.5
        end

        ax = Axis(fig[i, 1];
            xlabel = (i == 7) ? "Time (s)" : "",
            ylabel = rich(rich("x", font=:italic), "  (m)"),
            yscale = log10,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(9),
            yticks = (10.0 .^ (-12:3:-6), [powlab(-12), powlab(-9), powlab(-6)]),
            title = "$sp_name Concentration Contour"
        )

        hm = heatmap!(ax, times, X .+ 1e-12, M';
            colorrange = (c_min, c_max),
            colormap = discrete_cmap,
            interpolate = false)

        Colorbar(fig[i, 2], hm;
            label = rich("log", subscript("10"), "(", rich("c", font=:italic),
                         subscript(sp_label), " / M)"),
            ticklabelsize = 14, labelsize = 14)
    end
    return fig
end

# ── moved from cell 358b1fba (panel_log_contour!)  [+ electrolyte kwarg] ──
function panel_log_contour!(panel_pos, cbar_pos, result, X, times, sp, m;
                            scale=mol/dm^3, num_levels=24)
    electrolyte = m.elydata
    c_bulk = electrolyte.c_bulk[sp.idx]
    log_c_bulk = log10(c_bulk / scale)

    M = [log_c_bulk - log10(max(result.tsol[sp.idx, ix, it] / scale, 1e-12))
         for ix in 1:length(X), it in 1:length(times)]
    c_min = 0.0
    c_max = log_c_bulk - log10(1e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)
    ax = Axis(panel_pos;
        xlabel = lab_time,
        ylabel = rich(rich("x", font=:italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (10.0 .^ (-12:3:-6),
                  [powlab(-12), powlab(-9), powlab(-6)]))

    hm = heatmap!(ax, times, X .+ 1e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false)

    Colorbar(cbar_pos, hm;
        label = rich("log", subscript("10"), "(", rich("c", font=:italic),
                     subscript("bulk"), ") − log", subscript("10"), "(",
                     rich("c", font=:italic), subscript(sp.label), ")"),
        ticklabelsize = 20, labelsize = 20)

    return ax, hm
end

# ── moved from cell 34857db0 (QoverK)  [+ electrolyte kwarg] ──
function QoverK(fig, panel_pos, result, m; scale=mol/dm^3, lw=4,
                xlabel = lab_time, legend_pos = nothing)
    bulk = m.bulk
    electrolyte = m.elydata

    times = result.tsol.t
    nt    = length(times)

    # Equilibrium Constant
    EqK_hco3 = (electrolyte.c_bulk[ihco3] * electrolyte.c_bulk[iohminus]) /
           		electrolyte.c_bulk[ico3]
	EqK_co3 = (electrolyte.c_bulk[ico2] * electrolyte.c_bulk[iohminus]) /
           electrolyte.c_bulk[ihco3]

    # reaction quotient Q(t)
    cco3  = [result.tsol[ico3,     1, t] for t in 1:nt]
    cohm  = [result.tsol[iohminus, 1, t] for t in 1:nt]
    chco3 = [result.tsol[ihco3,    1, t] for t in 1:nt]
	cco2  = [result.tsol[ico2,     1, t] for t in 1:nt]

    Qt_hco3    = (chco3 .* cohm) ./ cco3
	Qt_co3     = (cco2 .* cohm) ./ chco3

    ax = Axis(panel_pos;
        xlabel = xlabel,
        ylabel = rich(rich("Q", font=:italic), " / ", rich("K", font=:italic)),
        yscale = log10,
        yticks = (10.0 .^ (-6:3:6),
	          [powlab(-6), powlab(-3), powlab(0), powlab(3), powlab(6)]),
				        yminorticksvisible = true,
		        yminorticks = IntervalsBetween(27)
			 )
    hlines!(ax, [1e0]; color = :black, linestyle = :dash, linewidth = 2)

    l1 = lines!(ax, times, max.(Qt_hco3 ./ EqK_hco3, eps(Float64));
                color = "#2980B9", linewidth = lw)
    l2 = lines!(ax, times, max.(Qt_co3  ./ EqK_co3,  eps(Float64));
                color = "#E67E22", linewidth = lw)

	ylims!(ax, low = 1e-7, high = 1e7)

	leg = Legend(legend_pos,
                 [l1, l2],
                 [rich("HCO", subscript("3"), superscript("-"), " ⇌ CO", subscript("3"), superscript("2-")),
                  rich("CO", subscript("2"), " ⇌ HCO", subscript("3"), superscript("-"))],
                 framevisible = false)

    return ax, leg
end

# ── moved from cell f506e83f (plot_cv_total_current)  [ely = elydata_Gold_odr → ely = model] ──
function plot_cv_total_current(result, m;
                               redox_species::Dict{Int,Int},
                               include_capacitive::Bool = true,
                               scale = cm^2/mA,
                               sign::Int = 1,
                               color_gradient::Bool = true)
    model = m.elydata

    n_t = length(result.voltages)

    # ---- Faradaic current: sum over all redox species ----
    I_F = zeros(n_t)
    for (idx, n_e) in redox_species
        I_F .+= n_e .* currents(result, idx)
    end

    # ---- Capacitive current: icc boundary species (only :ohmicdrop) ----
    I_C = zeros(n_t)
    if include_capacitive
        ely = model
        if ely.ircompensation == :ohmicdrop
            icc = ely.icc
            node_we = 1   # working-electrode boundary node (Γ_we = 1)
            I_C = [u[icc, node_we] for u in result.tsol[1:end-1]]
        else
            @warn "Capacitive current only resolved in :ohmicdrop mode " *
                  "(current mode: :$(ely.ircompensation)); j_C set to 0."
        end
    end

    # ---- Total current with sign convention and unit scaling ----
    I_total = sign .* (I_F .+ I_C) .* scale

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (540, 440))
        a = Axis(f[1, 1],
                 ylabel = lab_current,
                 xlabel = lab_voltage)
        return f, a
    end

    if color_gradient
        cols = RGBf.(range(0, 1, length = n_t), 0.0, 0.0)
        lines!(ax, result.voltages, I_total; color = cols)
    else
        lines!(ax, result.voltages, I_total)
    end

    return fig
end

# ── moved from cell 2f0d56cb (plot_combined_exp_sim_ivc)  [+ electrolyte kwarg] ──
function plot_combined_exp_sim_ivc(
    P_recs, m;
    redox_species = Dict(ico => 2),
    include_capacitive = true,
    scale = cm^2/mA,
    sign = 1,
    sim_linewidth::Real = 5
)
    electrolyte = m.elydata
    wanted    = ["0.1", "0.5", "1.0"]
    pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]
    raw_cv = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header=false)
    sub    = Matrix(raw_cv[4:end, :])
    num_cv = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num_cv, :auto)
    npairs = size(num_df, 2) ÷ 2
    keep   = findall(in(wanted), pressures[1:npairs])
    n_exp = length(keep)
    n_sim = length(P_recs)
    pastel2  = cgrad([colorant"#F2728A", colorant"#5BA8E8"])
    cols_exp = [pastel2[t] for t in range(0, 1, length = max(n_exp, 1))]
    pastel3  = cgrad([colorant"#FFB3BA", colorant"#A3D8FF"])
    cols_sim = [pastel3[t] for t in range(0, 1, length = max(n_sim, 1))]
    fig = with_theme(electrochemistry_theme()) do
        return Figure(size = (800, 800))
    end
    ax_exp = Axis(fig[1, 1],
        ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"),
        xticklabelsvisible = false, xticksvisible = false)
    ax_sim = Axis(fig[2, 1],
        xlabel = rich(rich("ϕ", font=:italic), "  (V vs. SHE)"),
        ylabel = rich(rich("I", font=:italic), "  (mA cm", superscript("−2"), ")"))
    linkxaxes!(ax_exp, ax_sim)
    ax_sim.xticks = -1.5:0.3:1.0
    ax_exp.limits = (nothing, (-5.5, 1.8))
    ax_exp.yticks = 1:-2:-5
    # ---- exp ----
    for (k, j) in enumerate(keep)
        xcol, ycol = 2j - 1, 2j
        label_text = "$(pressures[j]) atm"
        lines!(ax_exp, num_df[!, xcol], num_df[!, ycol];
               color = cols_exp[k], linewidth = sim_linewidth)
        text!(ax_exp, label_text;
              position = (-0.79 - 0.03*k, 1.6 - 1.3*k),
              color = cols_exp[k], fontsize = 26, font = :bold)
    end
    # ---- theory ----
    for j in 1:n_sim
        p, rec = P_recs[j]
        n_t = length(rec.voltages)
        I_F = zeros(n_t)
        for (idx, n_e) in redox_species
            I_F .+= n_e .* currents(rec, idx)
        end
        I_C = zeros(n_t)
        if include_capacitive
            ely = electrolyte
            if ely.ircompensation == :ohmicdrop
                icc = ely.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in rec.tsol[1:end-1]]
            end
        end
        I = sign .* (I_F .+ I_C) .* scale .* 2
        label_text = "$(p) atm"
        lines!(ax_sim, rec.voltages, I;
               color = cols_sim[j], linewidth = sim_linewidth)
        text!(ax_sim, label_text;
              position = (-0.99 - 0.05*j, 0.7 - 0.7*j),
              color = cols_sim[j], fontsize = 26, font = :bold)
    end
    text!(ax_sim, "Theory";
          position = (-0.5, -0.8), color = :gray30,
          fontsize = 32, font = :bold)
    text!(ax_exp, "Experiment";
          position = (-0.5, -1.27), color = :gray40,
          fontsize = 32, font = :bold)
    Label(fig[1, 1, TopLeft()], "(a)";
          fontsize = 28, font = :bold, padding = (0, 5, 20, 0))
    Label(fig[2, 1, TopLeft()], "(b)";
          fontsize = 28, font = :bold, padding = (0, 5, 20, 0))
    rowgap!(fig.layout, 1, 15)
    return fig
end
