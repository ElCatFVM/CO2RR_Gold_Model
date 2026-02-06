# plots/cvplot.jl

using CairoMakie
using CSV
using DataFrames


"""
    plot_cv_current(pnpresult, model;
                    species=nothing,
                    fig_size=(650, 400),
                    scale=cm^2/mA,
                    color_gradient=true)

Plot a CV current trace from a `pnpresult` and return only the figure.

- If `species` is provided, plots `currents(pnpresult, species) .* scale`.
- If `species === nothing`, defaults to the first entry of `model.cspecies`.
- `color_gradient=true` applies a red gradient along the voltage axis.
- Returns: `fig`.
"""
function plot_cv_current(pnpresult, model;
                         species=nothing,
                         fig_size=(650, 400),
                         scale=cm^2/mA,
                         color_gradient=true)

    ic = model.cspecies
    sp = (species === nothing) ? ic[6] : species

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
               ylabel = L"I (mA/cm^2)",
               xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
               # limits = ((-0.5, 0.9), (-2e-20, 2e-20)),
    )

    I = currents(pnpresult, sp) .* scale

    if color_gradient
        cols = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0)
        lines!(ax, pnpresult.voltages, I; color=cols)
    else
        lines!(ax, pnpresult.voltages, I)
    end

    return fig
end



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
    plot_conc_time_electrode(result, bulk;
                             nspecies=7,
                             fig_size=(650, 400),
                             scale=(mol/dm^3),
                             log_y=true,
                             legend=true)

Plot electrode-adjacent concentrations vs time for multiple species.

Assumptions
- `result.tsol.t` is the time vector.
- `result.tsol[i, 1, t]` accesses the concentration of species `i` at the electrode node (index 1)
  at time index `t`.
- `bulk` is an array-like object whose entries have `:name` and `:color` fields.

Behavior
- If `log_y=true`, negative/zero values are replaced by `NaN` before plotting.

Returns
- `fig`
"""
function plot_conc_time_electrode(result, bulk;
                                  nspecies=7,
                                  fig_size=(650, 400),
                                  scale=(mol/dm^3),
                                  log_y=true,
                                  legend=true)

    species = getproperty.(bulk, :name)
    colors  = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    conc_at_electrode = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1],
              xlabel = L"time / s",
              ylabel = L"c_{i,\,\mathrm{electrode}} / (\mathrm{mol/dm^3})",
              yscale = (log_y ? log10 : identity))

    for i in 1:nspecies
        y = conc_at_electrode[i, :]
        yplot = log_y ? (map(c -> (c > 0 ? c : NaN), y)) : y
        lines!(ax, times, yplot; color = colors[i], label = string(species[i]))
    end

    if legend
        Legend(fig[1, 2], ax; labelsize=10, backgroundcolor=RGBA(1.0, 1.0, 1.0, 0.5))
    end

    return fig
end


"""
    plot_cv_current(result, model; species=nothing, fig_size=(650, 400), scale=cm^2/mA)

Plot `currents(result, species) .* scale` versus `result.voltages`.
If `species` is `nothing`, use `model.cspecies[1]`.
Returns `fig`.
"""
function plot_cv_current(result, model;
                         species=nothing,
                         fig_size=(650, 400),
                         scale=cm^2/mA,
                         color_gradient=true)

    sp = (species === nothing) ? model.cspecies[1] : species

    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1],
              ylabel = L"I (mA/cm^2)",
              xlabel = L"\phi (V \; \mathrm{vs}\; SHE)")

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
function plot_conc_time_electrode(result, bulk;
                                  nspecies=7,
                                  fig_size=(650, 400),
                                  scale=(mol/dm^3))

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt    = length(times)

    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
               xlabel = L"time / s",
               ylabel = L"c_{i,\,\mathrm{electrode}} / (\mathrm{mol/dm^3})",
               yscale = log10)

    for i in 1:nspecies
        y = conc[i, :]
        y = map(c -> (c > 0 ? c : NaN), y)
        lines!(ax, times, y; color=colors[i], label=string(names[i]))
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
function plot_conc_profile_logx(result, bulk, X, t_index;
                                fig_size=(650, 400),
                                x_limits=(1e-12, 1e-3),
                                y_limits=(-14, 1),
                                x_offset=1e-14,
                                scale=(mol/dm^3))

    tsol = result.tsol ./ scale
    nvar, nx, nt = size(tsol)

    names  = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    nspecies = min(nvar, length(names), length(colors))

    nplot = min(nx, length(X))
    ti = clamp(t_index, 1, nt)

    xx = X[1:nplot] .+ x_offset

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
        xlabel = "Distance from electrode [m]",
        ylabel = "log10(c)",
        xscale = log10,
        limits = (x_limits, y_limits),
        title  = "t = $(round(result.tsol.t[ti], digits=4)) s | ϕ = $(round(result.voltages[ti], digits=2)) V",
    )

    for i in 1:nspecies
        conc = tsol[i, 1:nplot, ti]
        y = map(c -> (c > 0 ? log10(c) : NaN), conc)
        lines!(ax, xx, y; color=colors[i], label=string(names[i]))
    end

    axislegend(ax; position=:rt)
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
    pnpresult, bulk, X;
    file="concentrations_cv.gif",
    framerate=10,
    step=5,
    scale=(mol/dm^3),
    x_offset=1e-14,
    x_limits=(1e-12, 1e-3),
    y_limits=(-14, 2),
    legend=true,
)
    tsol = pnpresult.tsol ./ scale
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

    ys = [Observable(fill(NaN, nplot)) for _ in 1:nspecies]
    for i in 1:nspecies
        lines!(ax, xx, ys[i]; color=colors[i], label=string(names[i]))
    end
    legend && axislegend(ax)

    record(fig, file, t_indices; framerate=framerate) do ti
        tval = pnpresult.tsol.t[ti]
        ϕval = (hasproperty(pnpresult, :voltages) && ti <= length(pnpresult.voltages)) ? pnpresult.voltages[ti] : NaN
        ax.title = "t = $(round(tval, digits=4)) s | ϕ = $(round(ϕval, digits=2)) V"

        for i in 1:nspecies
            conc = tsol[i, 1:nplot, ti]
            ys[i][] = map(c -> (c > 0 ? log10(c) : NaN), conc)
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
              ylabel = L"I (mA/cm^2)",
              xlabel = L"\phi (V \; \mathrm{vs}\; SHE)")

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

function plot_pressure_varied_sweep(P_recs;
                                   species=iohminus,
                                   fig_size=(1600, 900),
                                   scale=cm^2/mA,
                                   limits=nothing
)

    fig = Figure(size = fig_size)
    if limits !== nothing
        ax = Axis(fig[1, 1],
                  xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
                  ylabel = L"I (mA/cm^2)",
                  limits = limits,
                  )
    else
        ax = Axis(fig[1, 1],
                xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
                ylabel = L"I (mA/cm^2)",
                )
    end
    n = length(P_recs)
    cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]

    plots  = Any[]
    labels = String[]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        push!(labels, label)

        I = currents(rec, species) .* scale
        push!(plots, lines!(ax, rec.voltages, I; color = cols[j]))
    end

    Legend(fig[1, 2], plots, labels, "Theoretical"; framevisible=true)
    return fig
end

function plot_pH_varied_sweep(pH_recs;
                              species=iohminus,
                              fig_size=(1600, 900),
                              scale=cm^2/mA)

    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1],
              xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
              ylabel = L"I (mA/cm^2)")

    n = length(pH_recs)
    cols = [RGB(1 - t, 0, t) for t in LinRange(0, 1, n)]

    plots  = Any[]
    labels = String[]

    for j in 1:n
        pH, rec = pH_recs[j]
        label = "pH = $(pH)"
        push!(labels, label)

        ivres = hasproperty(rec, :ivresult) ? getproperty(rec, :ivresult) : rec

        I = currents(ivres, species) .* scale
        push!(plots, lines!(ax, ivres.voltages, I; color = cols[j]))
    end

    Legend(fig[1, 2], plots, labels, "Theoretical"; framevisible=true)
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
    ax  = Axis(fig[1, 1], ylabel = L"I (mA/cm²)", xlabel = L"φ (V vs SHE)")

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
    ax  = Axis(fig[1, 1], ylabel = L"I (mA/cm²)", xlabel = L"φ (V vs SHE)")

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
	        xlabel = L"φ (V vs SHE)",
	        ylabel = L"I (mA/cm²)",
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
    pnpresult, bulk, X;
    ispec::Int = 5,
    frac::Float64 = 0.99,
    t_indices = [1, 10, 50, 100, 140],
    t_index::Int = 140,
    fig_size = (800, 420),
)
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
        ylabel = L"\log_{10} c_i\; (mol/dm^3)",
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
	              ylabel = L"I (mA/cm²)",
	              xlabel = L"φ (V vs SHE)"
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
        xlabel = L"\phi_{we}\;(\mathrm{V\;vs\;SHE})",
        ylabel = L"I\;(\mathrm{mA/cm^2})",
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
             xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
             ylabel = L"I (mA/cm^2)",
             limits = limits)
    else
        Axis(fig[1, 1],
             xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
             ylabel = L"I (mA/cm^2)")
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
