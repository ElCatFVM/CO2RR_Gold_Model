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
function plot_time_voltage_and_dt(
        pnpresult, sawtooth;
        fig_size = (700, 450),
        dt_yscale = log10,
        use_times_field = true,
        t_step = 1.0
    )

    fig = Figure(size = fig_size)

    axV = Axis(
        fig[1, 1],
        xlabel = "Time (s)",
        ylabel = "Voltage"
    )

    axdt = Axis(
        fig[2, 1],
        xlabel = "Time (s)",
        ylabel = "Δt (s)",
        yscale = dt_yscale
    )

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
        dt = T[2:end] .- T[1:(end - 1)]
        lines!(axdt, T[2:end], dt)
    end

    return fig
end


"""
    plot_cv_current(result, model; species=nothing, fig_size=(650, 400), scale=cm^2/mA)

Plot [`cv_current`](@ref) against [`cv_abscissa`](@ref). Returns `fig`.

The old default read the flux of `model.cspecies[1]`, i.e. K⁺ — a spectator that carries
no faradaic current at all. The default is now CO with two electrons.
"""
function plot_cv_current(
        result, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        fig_size = (650, 400),
        scale = cm^2 / mA,
        color_gradient = true
    )
    I = cv_current(result; species, n_e, sgn, include_capacitive) .* scale
    U, xlab = cv_abscissa(result; kind = abscissa)

    fig = Figure(size = fig_size)
    ax = Axis(
        fig[1, 1],
        ylabel = lab_current,
        xlabel = xlab
    )

    if color_gradient
        cols = RGBf.(range(0, 1, length(U)), 0.0, 0.0)
        lines!(ax, U, I; color = cols)
    else
        lines!(ax, U, I)
    end

    return fig
end


"""
    plot_conc_time_electrode(result, bulk; nspecies=7, fig_size=(650, 400), scale=(mol/dm^3))

Plot electrode-adjacent concentrations `result.tsol[i, 1, t] / scale` versus time for `i=1:nspecies`
(using `log10` y-scale; nonpositive values are shown as `NaN`). Returns `fig`.
"""
function plot_conc_time_electrode(
        result, m;
        nspecies = 7,
        fig_size = (800, 500),
        scale = (mol / dm^3)
    ) # if the unit is uncertain, test with 1.0 first
    bulk = m.bulk

    names = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt = length(times)

    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    fig = Figure(size = fig_size)
    ax = Axis(
        fig[1, 1],
        xlabel = lab_time,
        ylabel = rich(
            rich("c", font = :bold_italic), subscript("i, electrode"),
            "  (mol/dm", superscript("3"), ")"
        ),
        limits = ((times[1] - (times[end] / 200), times[end] + (times[end] / 100)), (1.0e-12, 1.0e4)),
        yscale = log10,
        yticks = (
            10.0 .^ (4:-4:-12),
            [powlab(4), powlab(0), powlab(-4), powlab(-8), powlab(-12)],
        )
    )

    for i in 1:nspecies
        y = conc[i, :]

        y_fixed = map(c -> (c > 1.0e-128 ? c : 1.0e-128), y)

        lines!(ax, times, y_fixed; color = colors[i], label = string(names[i]))
    end


    Legend(fig[1, 2], ax; labelsize = 10, backgroundcolor = RGBA(1, 1, 1, 0.5))
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
function plot_conc_profile_logx(
        result, m, X, t_index;
        fig_size = (650, 400),
        x_limits = (1.0e-12, 1.0e-3),
        y_limits = (-14, 4), # range in terms of log10(c)
        x_offset = 1.0e-14,
        scale = 1.0
    )
    bulk = m.bulk

    tsol = result.tsol
    nvar, nx, nt = size(tsol)

    names = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    nspecies = min(nvar, length(names), length(colors))

    nplot = min(nx, length(X))
    ti = clamp(t_index, 1, nt)

    xx = X[1:nplot] .+ x_offset

    fig = Figure(size = fig_size)

    phi_val = hasproperty(result, :voltages) ? round(result.voltages[ti], digits = 2) : "N/A"

    ax = Axis(
        fig[1, 1],
        xlabel = "Distance from electrode [m]",
        ylabel = rich("log", subscript("10"), "(", rich("c", font = :bold_italic), ")"),
        xscale = log10,
        limits = (x_limits, y_limits),
        title = "t = $(round(result.tsol.t[ti], digits = 4)) s | ϕ = $phi_val V",
    )

    for i in 1:nspecies
        conc = tsol[i, 1:nplot, ti] ./ scale

        y_log = map(c -> log10(max(c, 1.0e-25)), conc)

        lines!(ax, xx, y_log; color = colors[i], label = string(names[i]))
    end

    axislegend(ax; position = :rt, labelsize = 10)
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
        file = "concentrations_cv.gif",
        framerate = 10,
        step = 5,
        scale = (mol / dm^3),
        x_offset = 1.0e-14,
        x_limits = (1.0e-12, 1.0e-3),
        y_limits = (-14, 2),
        legend = true,
    )
    bulk = m.bulk

    raw_data = pnpresult.tsol
    tsol = raw_data ./ scale
    nvar, nx, nt = size(tsol)

    names = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)
    nspecies = min(nvar, length(names), length(colors))

    nplot = min(nx, length(X))
    xx = X[1:nplot] .+ x_offset

    t_indices = 1:step:nt

    fig = Figure(size = (650, 400))
    ax = Axis(
        fig[1, 1],
        xlabel = "Distance from electrode [m]",
        ylabel = "log10(c)",
        xscale = log10,
        limits = (x_limits, y_limits),
    )

    println("Plotting $nspecies species, $nt time steps.")

    ys = [Observable(fill(NaN, nplot)) for _ in 1:nspecies]
    for i in 1:nspecies
        lines!(ax, xx, ys[i]; color = colors[i], label = string(names[i]))
    end
    legend && axislegend(ax)

    record(fig, file, t_indices; framerate = framerate) do ti
        tval = try
            pnpresult.tsol.t[ti]
        catch
            ti
        end

        phi_str = "N/A"
        if hasproperty(pnpresult, :voltages) && ti <= length(pnpresult.voltages)
            phi_str = "$(round(pnpresult.voltages[ti], digits = 2))"
        end
        ax.title = "t = $(round(tval, digits = 4)) s | ϕ = $phi_str V"

        for i in 1:nspecies
            conc = tsol[i, 1:nplot, ti]
            # 3. robust log: replace non-positive values with a small floor (1e-20) to avoid breaks in the curve
            ys[i][] = map(c -> (c > 1.0e-20 ? log10(c) : -20.0), conc)
        end
    end

    return abspath(file)
end


function plot_cv_model_vs_koper_facets(
        pnpresult; species = iohminus,
        koper_csv_relpath = "data/Langmuir_CV_data/Figure_1.csv",
        fig_size = (1050, 650),
        koper_v_shift = -0.4
    )

    fig = Figure(size = fig_size)
    ax = Axis(
        fig[1, 1],
        ylabel = lab_current,
        xlabel = lab_voltage
    )

    # --- Gold model ---
    I_model = currents(pnpresult, species) .* (cm^2 / mA)
    gold_line = lines!(
        ax, pnpresult.voltages, I_model;
        color = RGBf.(range(0, 1, length(pnpresult.voltages)), 0.0, 0.0)
    )

    # --- Koper data (Figure 1) ---
    root = normpath(joinpath(@__DIR__, ".."))
    csv_path = joinpath(root, splitdir(koper_csv_relpath)...)

    raw_df = CSV.read(csv_path, DataFrame; header = false)
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

    Legend(
        fig[1, 2],
        [[gold_line], facet_lines],
        [labels1, labels2],
        ["Model", "Koper\nFacets"]
    )

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

    plots = Any[]
    # Not `String[]`: `label_fmt` may return `rich(...)` for a formatted entry.
    labels = Any[]

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

        plt = lines!(ax, ivres.voltages, I; color = cols[j], linewidth = LW_LINE)
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
function plot_iv_with_experiment(
        pnpresult, ico;
        csv_path::AbstractString = "../data/Langmuir_CV_data/Figure_3.csv",
        exp_title::AbstractString = "Experimental",
        fig_size::Tuple{Int, Int} = (1600, 900),
        markersize::Real = 12,
    )
    fig = Figure(size = fig_size)
    ax = Axis(fig[1, 1], ylabel = lab_current, xlabel = lab_voltage)

    volts = vec(pnpresult.voltages)
    total_current = currents(pnpresult, ico) .* (cm^2 / mA)

    colgrad = RGBf.(range(0, 1, length(volts)), 0.0, 0.0)
    lines!(ax, volts, total_current; color = colgrad)
    scatter!(ax, volts, total_current; markersize = markersize, color = colgrad)

    plot_objs = Any[]
    labels = String[]

    try
        raw = CSV.read(csv_path, DataFrame; header = false)

        pres = vec(Matrix(raw[1:1, :]))
        sub = Matrix(raw[4:end, :])

        num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
        num_df = DataFrame(num, :auto)

        npairs = size(num_df, 2) ÷ 2

        pink = RGB(1.0, 0.7, 0.8)
        pblue = RGB(0.2, 0.5, 1.0)
        cols = [
            RGB(
                    pink.r + t * (pblue.r - pink.r),
                    pink.g + t * (pblue.g - pink.g),
                    pink.b + t * (pblue.b - pink.b)
                ) for t in range(0, 1, length = npairs)
        ]

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
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        fig_size = (800, 400),
        scale = cm^2 / mA,
        legend_title = "Scan Rates (V/s)",
        highlight_index = nothing,
        highlight_color = RGB(1, 0.2, 0.2),
        highlight_lw = LW_HIGHLIGHT,
        default_lw = LW_LINE,
    )
    fig = Figure(size = fig_size)

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(sweep_vec) ? lab_voltage : cv_abscissa(sweep_vec[1]; kind = abscissa)[2]
    ax = Axis(fig[1, 1], ylabel = lab_current, xlabel = xlab)

    n = length(sweep_vec)
    cols = [CMAP_SCANRATE[t] for t in range(0, 1, length = max(n, 1))]

    plot_objs = Any[]
    labels = String[]

    for (j, rec) in enumerate(sweep_vec)
        push!(labels, "$(scanrates[j])\t\t ")

        color_j = cols[j]
        lw_j = default_lw

        if highlight_index !== nothing && j == highlight_index
            color_j = highlight_color
            lw_j = highlight_lw
        end

        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)
        line = lines!(ax, U, I; linewidth = lw_j, color = color_j)
        push!(plot_objs, line)
    end

    Legend(fig[1, 2], plot_objs, labels, legend_title; framevisible = true)
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
    colors = getproperty.(bulk, :color)

    nspecies = min(nvar, length(species), length(colors))
    @assert 1 ≤ ispec ≤ nspecies

    nplot = min(nx, length(X))
    xx = X[1:nplot] ./ μm

    t_indices = unique(clamp.(vcat(t_indices, nt), 1, nt))
    t_index = clamp(t_index, 1, nt)

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
    ax = Axis(
        fig[1, 1],
        xlabel = "x / μm",
        ylabel = rich(
            "log", subscript("10"), " ", rich("c", font = :bold_italic), subscript("i"),
            "  (mol/dm", superscript("3"), ")"
        ),
        title = "t = $(round(pnpresult.tsol.t[t_index], digits = 4)) s  |  δ$(Int(round(frac * 100))) for $(species[ispec])",
    )

    conc = vec(tsol[ispec, 1:nplot, t_index])
    yvals = map(c -> (c > 0 ? log10(c) : NaN), conc)

    lines!(ax, xx, yvals; color = colors[ispec], linewidth = LW_LINE, label = species[ispec])

    δ_here = x_at_frac(tsol, ispec, t_index, xx, c_bulk, frac)
    if isfinite(δ_here)
        vlines!(ax, [δ_here]; linestyle = :dash, linewidth = LW_GUIDE)
        ytop = maximum(filter(isfinite, yvals))
        text!(
            ax, δ_here, ytop;
            text = "  δ$(Int(round(frac * 100)))≈$(round(δ_here, digits = 4)) μm",
            align = (:left, :top)
        )
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
function plot_cv_over_L(
        results::Dict{Float64, Any};
        species = ico, cutoff = -0.4, title = "IV vs L"
    )
    vis = GridVisualizer(;
        size = (800, 500),
        title = title,
        xlabel = lab_voltage,
        ylabel = lab_current,
        legend = :rt,
        #yscale = :log,
    )

    items = sort(collect(results); by = first)
    n = length(items)
    cols = Makie.resample_cmap(:cool, n)

    for (k, (L, rec)) in enumerate(items)
        ivres = hasproperty(rec, :ivresult) ? getproperty(rec, :ivresult) : rec

        volts = vec(ivres.voltages)
        Iall = (vec(currents(ivres, species))) .* (cm^2 / mA)

        m = min(length(volts), length(Iall))
        volts = volts[1:m]
        Iall = Iall[1:m]

        mask = volts .< cutoff

        scalarplot!(
            vis,
            volts[mask],
            Iall[mask];
            clear = false,
            label = "L = $(L) μm",
            color = cols[k],
        )
    end

    return reveal(vis)
end

function plot_pressure_varied_sweep(
        P_recs;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        fig_size = (1600, 900),
        scale = cm^2 / mA,
        limits = nothing,
        # --- add experimental background (Figure_3.csv) ---
        fig3_csv::Union{Nothing, AbstractString} = "../data/Langmuir_CV_data/Figure_3.csv",
        exp_title::AbstractString = "Experimental",
        exp_alpha::Real = 0.55,
        exp_linewidth::Real = LW_EXP,
        sim_linewidth::Real = LW_LINE,
    )
    fig = Figure(size = fig_size)

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(P_recs) ? lab_voltage : cv_abscissa(P_recs[1][2]; kind = abscissa)[2]

    ax = if limits !== nothing
        Axis(
            fig[1, 1],
            xlabel = xlab,
            ylabel = lab_current,
            limits = limits
        )
    else
        Axis(
            fig[1, 1],
            xlabel = xlab,
            ylabel = lab_current
        )
    end

    plots = Any[]
    labels = String[]

    # ------------------------------------------------------------
    # 1) Experimental background (Figure_3.csv)  [optional]
    # ------------------------------------------------------------
    if fig3_csv !== nothing
        try
            raw = CSV.read(fig3_csv, DataFrame; header = false)

            pres = vec(Matrix(raw[1:1, :]))
            sub = Matrix(raw[4:end, :])

            num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
            num_df = DataFrame(num, :auto)

            npairs = size(num_df, 2) ÷ 2

            pink = RGB(1.0, 0.7, 0.8)
            pblue = RGB(0.2, 0.5, 1.0)
            cols_exp = [
                RGB(
                        pink.r + t * (pblue.r - pink.r),
                        pink.g + t * (pblue.g - pink.g),
                        pink.b + t * (pblue.b - pink.b)
                    ) for t in range(0, 1, length = npairs)
            ]

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
    cols_sim = [CMAP_PRESSURE[t] for t in LinRange(0, 1, max(n, 1))]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)

        line = lines!(ax, U, I; color = cols_sim[j], linewidth = sim_linewidth)
        push!(plots, line)
        push!(labels, "Theoretical | $label")
    end

    Legend(fig[1, 2], plots, labels, "Overlay"; framevisible = true)
    return fig
end

function pressure_varied_cvsweep(
        P_recs;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        fig_size = (800, 400),
        scale = cm^2 / mA,
        limits = nothing,
    )
    fig = Figure(size = fig_size)
    sim_linewidth = LW_LINE

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(P_recs) ? lab_voltage : cv_abscissa(P_recs[1][2]; kind = abscissa)[2]

    ax = if limits !== nothing
        Axis(
            fig[1, 1],
            xlabel = xlab,
            ylabel = lab_current,
            limits = limits
        )
    else
        Axis(
            fig[1, 1],
            xlabel = xlab,
            ylabel = lab_current
        )
    end

    plots = Any[]
    labels = String[]

    # ------------------------------------------------------------
    # 2) Simulation curves (existing logic)
    # ------------------------------------------------------------
    n = length(P_recs)
    cols_sim = [CMAP_PRESSURE[t] for t in LinRange(0, 1, max(n, 1))]

    for j in 1:n
        p, rec = P_recs[j]
        label = "$(p)\t pCO2(atm)"
        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)

        line = lines!(ax, U, I; color = cols_sim[j], linewidth = sim_linewidth)
        push!(plots, line)
        push!(labels, "Theoretical | $label")
    end

    Legend(fig[1, 2], plots, labels, "Overlay"; framevisible = true)
    return fig
end

"""
    split_axis_pair!(fig, row; ured, uox, kwargs...)

One row of a broken voltage axis: the reductive window `ured` in `fig[row, 1]`, the
oxidative window `uox` in `fig[row, 2]`, butted together with the featureless region
between them omitted and a single black rule at the seam.

Returns `(ax_red, ax_ox)`. The oxidative panel carries its ticks on the right, so the pair
reads as one axis with a piece removed rather than as two plots.

`titles` labels the two windows; pass `nothing` on every row but the first of a stack.
`show_xticklabels = false` likewise suppresses the tick labels on all but the bottom row.
"""
function split_axis_pair!(
        fig, row;
        ured, uox,
        ylabel = lab_current,
        xticks_red = LinearTicks(3),
        xticks_ox = LinearTicks(3),
        titles = nothing,
        show_xticklabels = true,
    )
    ax_red = Axis(
        fig[row, 1];
        ylabel = ylabel, xticks = xticks_red,
        xticklabelsvisible = show_xticklabels, xticksvisible = show_xticklabels,
    )
    ax_ox = Axis(
        fig[row, 2];
        yaxisposition = :right, xticks = xticks_ox,
        xticklabelsvisible = show_xticklabels, xticksvisible = show_xticklabels,
    )
    if titles !== nothing
        ax_red.title = titles[1]
        ax_ox.title = titles[2]
    end
    # One black rule at the seam: keep the left panel's right spine, drop the right
    # panel's left spine. Leaving both on would draw a double-width wall.
    ax_ox.leftspinevisible = false
    xlims!(ax_red, ured...)
    xlims!(ax_ox, uox...)
    return ax_red, ax_ox
end

"""
    split_series!(ax_red, ax_ox, U, I; ured, uox, kwargs...)

Draw one curve into both windows of a [`split_axis_pair!`](@ref) and return
`(plot, anchors)`, where `anchors` is `(; red, ox)`: the current at each window's **vertex**
— the most negative `U` reached inside `ured`, the most positive inside `uox`.

Samples outside a window become `NaN` rather than being dropped, which preserves the
sample order and lets Makie break the line instead of joining across the gap.

A vertex is where a direct label belongs, because that is where the family fans out. The
last in-window sample is the obvious alternative and is wrong: a cycle returns to its
start, so every curve ends on the baseline and all the labels land on one point. `NaN` when
the curve never enters that window.
"""
function split_series!(
        ax_red, ax_ox, U, I;
        ured, uox, anodic_gain = 1.0, color = :black, lw = LW_LINE,
    )
    I_red = [ured[1] <= u <= ured[2] ? x : NaN for (u, x) in zip(U, I)]
    I_ox = [uox[1] <= u <= uox[2] ? x * anodic_gain : NaN for (u, x) in zip(U, I)]

    line = lines!(ax_red, U, I_red; color = color, linewidth = lw)
    lines!(ax_ox, U, I_ox; color = color, linewidth = lw)

    in_red = findall(!isnan, I_red)
    in_ox = findall(!isnan, I_ox)
    anchors = (
        red = isempty(in_red) ? NaN : I_red[in_red[argmin(U[in_red])]],
        ox = isempty(in_ox) ? NaN : I_ox[in_ox[argmax(U[in_ox])]],
    )
    return line, anchors
end

"""
    annotate_split!(ax, anchors, labels, colors; window, at = :end, kwargs...)

Direct labelling for one window of a [`split_axis_pair!`](@ref): each value written at the
end of its own curve, in that curve's colour.

No patch column, no frame sitting on top of the data, and no round trip through a colour
key — the reader never has to match a swatch. The sweep turns around before the axis edge,
so the margin beyond the vertex is empty and the labels cost no data space.

`at` picks the side the labels sit on: `:end` for the oxidative window, whose vertex is at
its high-voltage edge, `:start` for the reductive one, whose vertex is at its low-voltage
edge. Pass the matching field of the [`split_series!`](@ref) anchors — `at = :start` with
the `ox` anchors would put every label at the wrong current.

Bold because a direct label sits on the data rather than in a frame: at this size colour
alone does not carry, and the pale end of a pastel ramp washes out.

`NaN` anchors — a curve that never enters the window — are skipped.
"""
function annotate_split!(
        ax, anchors, labels, colors;
        window, at = :end, annotate_x = nothing, fontsize = 22,
    )
    # 5 % in from the window edge, not 10 %: the labels belong in the margin the sweep
    # leaves beyond its vertex, and the further in they sit the more likely they are to
    # land on the branch coming back. `annotate_x` overrides with an absolute voltage.
    w = window[2] - window[1]
    x = if annotate_x !== nothing
        annotate_x
    elseif at === :start
        window[1] + 0.05w
    else
        window[2] - 0.05w
    end
    for (j, y) in enumerate(anchors)
        isnan(y) && continue
        text!(
            ax, x, y;
            text = labels[j], color = colors[j],
            fontsize = fontsize, font = :bold, align = (:center, :center),
        )
    end
    return nothing
end

"""
    split_widths!(fig, ured, uox; widths = nothing, ncols = 2)

Size the two columns of a broken axis so both windows share the same volts per unit width,
unless `widths` overrides it, and butt them together.
"""
function split_widths!(fig, ured, uox; widths = nothing)
    w = widths === nothing ? (ured[2] - ured[1], uox[2] - uox[1]) : widths
    colsize!(fig.layout, 1, Auto(w[1]))
    colsize!(fig.layout, 2, Auto(w[2]))
    colgap!(fig.layout, 1, 0)
    return nothing
end

"""
    pressure_varied_cvsweep_split(P_recs; ured = (-1.3, -0.6), uox = (0.0, 1.0), kwargs...)

Same data as [`pressure_varied_cvsweep`](@ref), drawn as a **broken voltage axis**: the
reductive window `ured` in `fig[1, 1]` and the oxidative window `uox` in `fig[1, 2]`,
butted together with the flat region between them omitted. A single black rule marks the
seam.

Two things make the anodic feature unreadable on one plain axis: the reduction peak is
one to two orders of magnitude taller, and roughly half the sweep is featureless
baseline. Cutting the baseline out and giving each window its own current axis fixes
both — the left tick labels belong to the reductive panel, the right ones to the
oxidative panel. `anodic_gain` scales the oxidative branch further when auto-scaling is
not enough; any gain other than 1 is written into the panel title, so a magnified branch
cannot be mistaken for a raw one.

By default the panel widths are proportional to the voltage spans they cover, so both
panels share the same volts-per-centimetre and the pair reads as one axis with a piece
removed. Pass `widths` to override.

Samples outside a panel's window become `NaN` rather than being dropped, which preserves
each curve's sample order and lets Makie break the line instead of joining across a gap.

!!! warning
    The two panels share neither a current scale nor a continuous voltage axis, so
    heights and slopes are **not** comparable across the seam. Quote `anodic_gain`, both
    windows and both y-ranges in any caption.

`kwargs` are the usual current/abscissa controls — see [`cv_current`](@ref) and
[`cv_abscissa`](@ref).
"""
function pressure_varied_cvsweep_split(
        P_recs;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        scale = cm^2 / mA,
        ured = (-1.3, -0.6),
        uox = (-0.2, 1.0),
        anodic_gain = 1.0,
        widths = nothing,
        fig_size = (960, 460),
        lw = LW_LINE,
        legend_title = lab_pressure,
        # The unit lives in `legend_title`, so an entry is just its number.
        label_fmt = p -> string(p),
        # Corner of the oxidation panel the legend sits in; `:lt`, `:rt`, `:lb`, `:rb`.
        # Only read when `legend_mode == :axis`.
        legend_position = :lt,
        # `:direct` writes each label in its own curve's colour at the end of the anodic
        # branch, `:axis` draws a boxed legend inside the oxidation panel, `:none` neither.
        legend_mode = :direct,
        # Voltage the direct labels sit at. Default is just right of where the sweep turns
        # around, in the empty margin between the anodic vertex and the axis edge.
        annotate_x = nothing,
        annotate_fontsize = 22,
        colormap = CMAP_PRESSURE,
        xticks_red = LinearTicks(3),
        xticks_ox = LinearTicks(3),
    )
    n = length(P_recs)
    cols_sim = [colormap[t] for t in range(0, 1, length = max(n, 1))]

    xlab = isempty(P_recs) ? lab_voltage : cv_abscissa(P_recs[1][2]; kind = abscissa)[2]
    ox_title = anodic_gain == 1 ? "oxidation" : @sprintf("oxidation  (×%g)", anodic_gain)

    fig = Figure(size = fig_size)
    # Only the left panel carries the current label — the unit is the same on both
    # sides, and the voltage label is a single centred `Label` under the pair, so the
    # broken axis reads as one axis instead of two fully decorated plots.
    ax_red, ax_ox = split_axis_pair!(
        fig, 1; ured, uox, xticks_red, xticks_ox, titles = ("reduction", ox_title),
    )

    plots = Any[]
    labels = String[]
    anchors = Float64[]                 # current at the anodic vertex, one per curve

    for j in 1:n
        p, rec = P_recs[j]
        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)

        line, anchor = split_series!(
            ax_red, ax_ox, U, I; ured, uox, anodic_gain, color = cols_sim[j], lw,
        )
        push!(plots, line)
        push!(labels, label_fmt(p))
        push!(anchors, anchor.ox)
    end

    if legend_mode === :direct
        annotate_split!(
            ax_ox, anchors, labels, cols_sim;
            window = uox, at = :end, annotate_x, fontsize = annotate_fontsize,
        )
    elseif legend_mode === :axis
        # `axislegend`, not `Legend`: it anchors to the axis interior. A `Legend` placed at
        # `fig[1, 2]` would claim layout space and shrink the panel it is meant to sit in.
        axislegend(
            ax_ox, plots, labels, legend_title;
            position = legend_position,
            # One row: the entries are bare numbers of a single swept variable, so reading
            # them left to right matches the colour ramp and costs a fraction of the panel
            # height a five-row column would. `nbanks = 1` pins it to one row — with
            # `:horizontal` alone Makie is free to wrap.
            orientation = :horizontal, nbanks = 1, colgap = 8,
            framevisible = true, backgroundcolor = (:white, 0.85),
            titlesize = FS_LEGEND_TITLE, labelsize = FS_LEGEND,
            patchsize = (18, 10), rowgap = 1, titlegap = 5,
            padding = (8, 8, 5, 5),
        )
    end

    # One voltage label centred under both panels. It is a `Label`, not an axis `xlabel`,
    # so the theme's `xlabelsize` / `xlabelfont` do not reach it — spell both out here or
    # this one string drifts away from every other axis label in the package.
    Label(
        fig[2, 1:2], xlab;
        fontsize = FS_LABEL, font = :regular, padding = (0, 0, 0, 6),
    )

    split_widths!(fig, ured, uox; widths)
    rowgap!(fig.layout, 1, 4)
    return fig
end

"""
    plot_exp_sim_cvsweep_split(P_recs; kwargs...)

Measured against simulated CO₂-pressure CVs, both drawn on the broken voltage axis of
[`pressure_varied_cvsweep_split`](@ref): experiment in row 1, theory in row 2, each split
into the reductive window `ured` and the oxidative window `uox`.

Broken-axis version of [`plot_combined_exp_sim_ivc`](@ref). Same comparison, but the flat
region between the two features is cut out of both rows at once, so the anodic branch —
one to two orders of magnitude below the reduction peak — is legible in the same figure
that shows the peak.

The x axes are linked, the tick labels appear only on the bottom row and the window titles
only on the top, so the four panels read as two rows of one axis. Each row keeps its own
current scale: the measurement and the model do not have to agree in magnitude for their
*shapes* to be compared, which is what this figure is for.

Curves are labelled directly, in their own colour, at the anodic vertex — see
[`annotate_split!`](@ref) — and each row is named by grey text at the top left of the same
panel. Both live in the oxidation window because the sweep turns around before the axis
edge there, leaving a margin the reduction window does not have.

Experiment takes the saturated [`CMAP_PRESSURE_EXP`](@ref) and theory the pale
[`CMAP_PRESSURE`](@ref), so a pressure keeps its hue across the rows while the rows stay
tellable apart.

The experimental CSV is the Langmuir `Figure_3.csv` layout: three header rows, then
`(voltage, current)` column pairs in the order given by `exp_pressures`, of which
`exp_wanted` are drawn. `exp_csv` is resolved relative to the working directory.

!!! warning
    `sim_gain` defaults to 1. [`plot_combined_exp_sim_ivc`](@ref) multiplies its simulated
    current by a hard-coded factor of 2 with no comment; that factor is not carried over
    here, so this figure will not reproduce it unless you pass `sim_gain = 2` and say so in
    the caption.

See [`pressure_varied_cvsweep_split`](@ref) for the current, abscissa and window keywords,
all of which apply here too.
"""
function plot_exp_sim_cvsweep_split(
        P_recs;
        exp_csv = "../data/Langmuir_CV_data/Figure_3.csv",
        exp_pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"],
        exp_wanted = ["0.1", "0.5", "1.0"],
        exp_sgn = 1,
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        scale = cm^2 / mA,
        sim_gain = 1.0,
        ured = (-1.3, -0.6),
        uox = (-0.2, 1.35),
        anodic_gain = 1.0,
        widths = nothing,
        fig_size = (960, 820),
        lw = LW_LINE,
        # The unit rides on every entry here, unlike the swept-family plots, because there
        # is no legend title left to carry it once the labels sit on the curves.
        label_fmt = p -> "$(p) atm",
        annotate_x = nothing,
        annotate_fontsize = 22,
        colormap_exp = CMAP_PRESSURE_EXP,
        colormap_sim = CMAP_PRESSURE,
        xticks_red = LinearTicks(3),
        xticks_ox = LinearTicks(3),
        row_labels = ("Experiment", "Theory"),
        # Relative position inside the oxidation panel, (0,0) bottom left to (1,1) top
        # right. Relative rather than data coordinates so the text does not have to be
        # re-tuned every time the current range changes.
        row_label_pos = (0.04, 0.92),
        row_label_color = :gray35,
        row_label_fontsize = 32,
    )
    # ---- experimental columns ----------------------------------------------------
    raw = CSV.read(exp_csv, DataFrame; header = false)
    # Rows 1-3 are the header block; everything below is numeric with `missing` padding
    # where a trace ends early, which becomes NaN so Makie breaks the line.
    num = map(x -> x === missing ? NaN : parse(Float64, x), Matrix(raw[4:end, :]))
    exp_df = DataFrame(num, :auto)
    npairs = size(exp_df, 2) ÷ 2
    keep = findall(in(exp_wanted), exp_pressures[1:npairs])

    n_exp = length(keep)
    n_sim = length(P_recs)
    cols_exp = [colormap_exp[t] for t in range(0, 1, length = max(n_exp, 1))]
    cols_sim = [colormap_sim[t] for t in range(0, 1, length = max(n_sim, 1))]

    xlab = isempty(P_recs) ? lab_voltage : cv_abscissa(P_recs[1][2]; kind = abscissa)[2]
    ox_title = anodic_gain == 1 ? "oxidation" : @sprintf("oxidation  (×%g)", anodic_gain)

    fig = Figure(size = fig_size)
    ax_exp_red, ax_exp_ox = split_axis_pair!(
        fig, 1; ured, uox, xticks_red, xticks_ox,
        titles = ("reduction", ox_title), show_xticklabels = false,
    )
    ax_sim_red, ax_sim_ox = split_axis_pair!(
        fig, 2; ured, uox, xticks_red, xticks_ox,
    )
    # Linked so the two rows cannot drift apart if a caller sets limits on only one of them.
    linkxaxes!(ax_exp_red, ax_sim_red)
    linkxaxes!(ax_exp_ox, ax_sim_ox)

    # ---- experiment --------------------------------------------------------------
    exp_labels = String[]
    exp_anchors = Float64[]
    for (k, j) in enumerate(keep)
        U = exp_df[!, 2j - 1]
        I = exp_sgn .* exp_df[!, 2j]
        _, anchor = split_series!(
            ax_exp_red, ax_exp_ox, U, I; ured, uox, anodic_gain, color = cols_exp[k], lw,
        )
        push!(exp_labels, label_fmt(exp_pressures[j]))
        push!(exp_anchors, anchor.ox)
    end
    annotate_split!(
        ax_exp_ox, exp_anchors, exp_labels, cols_exp;
        window = uox, at = :end, annotate_x, fontsize = annotate_fontsize,
    )

    # ---- theory ------------------------------------------------------------------
    sim_labels = String[]
    sim_anchors = Float64[]
    for j in 1:n_sim
        p, rec = P_recs[j]
        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale .* sim_gain
        U, _ = cv_abscissa(rec; kind = abscissa)
        _, anchor = split_series!(
            ax_sim_red, ax_sim_ox, U, I; ured, uox, anodic_gain, color = cols_sim[j], lw,
        )
        push!(sim_labels, label_fmt(p))
        push!(sim_anchors, anchor.ox)
    end
    annotate_split!(
        ax_sim_ox, sim_anchors, sim_labels, cols_sim;
        window = uox, at = :end, annotate_x, fontsize = annotate_fontsize,
    )

    if row_labels !== nothing
        # Inside the oxidation panel, not a `Label` in the layout margin: a margin label
        # competes with the window titles above row 1 and pushes the panels apart, and it
        # reads as a figure part number rather than as a name for the data under it. Grey
        # and unnumbered — the row is identified, not enumerated, so it does not fight the
        # coloured pressure labels for attention.
        for (ax, lab) in zip((ax_exp_ox, ax_sim_ox), row_labels)
            text!(
                ax, Point2f(row_label_pos...);
                text = lab, space = :relative,
                color = row_label_color, fontsize = row_label_fontsize, font = :bold,
                align = (:left, :top),
            )
        end
    end

    # One voltage label centred under both columns — see `pressure_varied_cvsweep_split`.
    Label(
        fig[3, 1:2], xlab;
        fontsize = FS_LABEL, font = :regular, padding = (0, 0, 0, 6),
    )

    split_widths!(fig, ured, uox; widths)
    rowgap!(fig.layout, 1, 10)     # between the two data rows
    rowgap!(fig.layout, 2, 4)      # above the shared voltage label
    return fig
end

"""
    plot_scanrate_sweeps_split(sweep_vec, scanrates; kwargs...)

Broken-axis version of [`plot_scanrate_sweeps`](@ref), with the same reduction/oxidation
split as [`pressure_varied_cvsweep_split`](@ref) — which is the shared implementation;
only the legend labelling differs. All of its keywords apply here too.
"""
plot_scanrate_sweeps_split(sweep_vec, scanrates; colormap = CMAP_SCANRATE, kwargs...) =
    pressure_varied_cvsweep_split(
    collect(zip(scanrates, sweep_vec));
    legend_title = lab_scanrate,
    label_fmt = sr -> string(sr),
    colormap = colormap,
    kwargs...
)


# =====================================================================
# Added from scripts/row_interaction_script.jl  (only added, nothing removed)
# Batch 1: CV current / scan-rate plotting family.
#
# Species index constants for the Gold CO2RR model, added so that the bare
# `iohminus`/`ico`/`ico2` defaults and bodies (here and in the existing
# functions above) resolve inside the AuCO2RR_plots module.
# =====================================================================
const ikplus = 1
const ihplus = 2
const ihco3 = 3
const ico3 = 4
const ico2 = 5
const iohminus = 6
const ico = 7

# =====================================================================
# Single source of truth for "the current" and "the voltage" in the CV plots.
#
# These used to be open-coded at every call site, and the sites disagreed on all
# four of: which species the flux is read from, the electron count, the sign, and
# whether the capacitive term is included — while all labelling the result
# identically. Route every CV figure through these so that two figures with the
# same axis label really do show the same quantity.
# =====================================================================

const lab_voltage_rp = rich(
    "Reaction-Plane Potential ", rich("U", font = :bold_italic), subscript("rp"), "\n(V vs. SHE)"
)
const lab_voltage_dl = rich(
    "Double-Layer Voltage ", rich("U", font = :bold_italic), subscript("dl"), "\n(V)"
)

"""
    faradaic_current(result; species = ico, n_e = 2, sgn = 1)

Faradaic current density at the working electrode, from the boundary flux of `species`.

`species` must take part in **no homogeneous reaction**, otherwise its boundary flux also
carries buffer-driven transport, which is not current. In the Gold CO2RR model CO (`ico`)
is the only such species — CO₂, HCO₃⁻, CO₃²⁻, OH⁻ and H⁺ are all buffer-active (see
`buffer_system` in `goldmodel.jl`), and K⁺ is a spectator carrying no faradaic current at
all. Measured on this model `currents(., ico2)` runs roughly 2x the faradaic rate, because
every reduced CO₂ releases 2 OH⁻ which consume further CO₂ through CO₂ + OH⁻ ⇌ HCO₃⁻.

`n_e` is the electron count per molecule of `species` (2 for CO and CO₂, 1 for OH⁻).
`sgn` fixes the display convention: CO and OH⁻ are produced while CO₂ is consumed, so
their fluxes come out with opposite signs.
"""
faradaic_current(result; species = ico, n_e = 2, sgn = 1) =
    sgn .* n_e .* currents(result, species)

"""
    capacitive_current(result)

Capacitive current density at the working electrode, or zeros when the result carries
none.

Reads `result.j_cap`, as the rest of this file does. Do not reach for the `icc` unknown
instead: it only exists when `ircompensation isa OhmicDropEstimation`, so that route
silently evaluates to zero in every other compensation mode and makes currents from
different modes incomparable without any warning.

The field is **present but empty** on a run that accumulated no capacitive term, so testing
`hasproperty` alone is not enough: an `Any[]` passed on to `lines!` fails inside Makie's
argument conversion with `reducing over an empty collection`, several frames away from
anything that names the problem. An empty field is reported and treated as zero here
instead — a figure with a flat capacitive panel is a result, not a crash, but it is never
what the caller expected, so it warns.
"""
function capacitive_current(result)
    n = length(result.times)
    hasproperty(result, :j_cap) || return zeros(n)
    j = result.j_cap
    if isempty(j)
        @warn "`result.j_cap` is empty: this run carried no capacitive current, so the " *
            "capacitive term is taken as zero." maxlog = 1
        return zeros(n)
    end
    length(j) < n && error(
        "`j_cap` has $(length(j)) entries for $(n) time points; it cannot be added to a " *
        "faradaic current of length $(n)"
    )
    # A stored solution often carries one entry more than there are sweep times.
    return length(j) == n ? j : j[1:n]
end

"""
    cv_current(result; include_capacitive = false, kwargs...)

[`faradaic_current`](@ref) plus, optionally, [`capacitive_current`](@ref).
`kwargs` are forwarded to [`faradaic_current`](@ref).
"""
function cv_current(result; include_capacitive = false, kwargs...)
    I = faradaic_current(result; kwargs...)
    include_capacitive || return I
    return I .+ capacitive_current(result)
end

"""
    reaction_plane_potential(result, m)

Electrostatic potential in the electrolyte at the electrode node, ϕ(0), read from the
stored solution.

Use this rather than `result.voltages` whenever compensation modes are compared. That field
is ϕ(0) under `NoIRCompensation` but `ϕ_we + ϕ_DL` — the *compensated applied* potential —
under `OhmicDropEstimation`, because that is what `ohmicdropcompensation` hands to
`potentialbcondition!`. Plotting it against the applied protocol therefore gives a straight
line for the compensated modes and a strongly attenuated curve for the uncompensated one,
which looks like a physical difference and is not.

With the Robin boundary condition, `ϕ_we - ϕ_pzc - ϕ(0)` is the drop across the Helmholtz
gap, so ϕ(0) is the diffuse-layer share of the applied potential. That split is set by
`C_gap` against the diffuse-layer capacitance and is essentially independent of the
compensation factor.
"""
function reaction_plane_potential(result, m)
    iϕ = m.elydata.iϕ
    n = length(result.times)
    return [u[iϕ, 1] for u in result.tsol[1:n]]
end

"""
    electrode_potential(result, m)

Potential of the metal, `φ_M` in the notation of Levey et al. — the value that
`potentialbcondition!` actually imposes at the working electrode.

Under `OhmicDropEstimation` that is `ϕ_we + ϕ_DL`, which is what `result.voltages` holds.
Under `NoIRCompensation` nothing is added, so it is the applied protocol itself; note that
`result.voltages` is *not* usable there, since in that mode it carries ϕ(0) instead.
"""
electrode_potential(result, m) =
    isa(m.elydata.ircompensation, OhmicDropEstimation) ? result.voltages : result.sawtooth

"""
    compensation_potential(result, m)

Potential the ohmic-drop compensation added, `ϕ_DL = factor · Ru · (j_F + j_C)`.

Exactly zero for an uncompensated run and proportional to the factor otherwise. It peaks
where the current peaks and vanishes wherever the current does.
"""
compensation_potential(result, m) = electrode_potential(result, m) .- result.sawtooth

"""
    driving_force(result, m)

`φ_M - φ_PET`, the potential difference actually driving electron transfer, measured from
the potential of zero charge: `electrode_potential - ϕ_pzc - ϕ(0)`.

This is the quantity Levey, Edwards, White & Macpherson plot against the applied potential
(*Phys. Chem. Chem. Phys.* **2023**, 25, 7832, Fig. 3 and 5c). Plotted that way a
loss-free cell gives a straight line of unit slope — their "diffusion model" reference —
and every deviation from it is potential that never reached the reaction plane. Ohmic-drop
compensation pulls the curve back towards that line, which is the whole point of applying
it.

Plot this, not [`reaction_plane_potential`](@ref) on its own: ϕ(0) alone is dominated by
the capacitive split between the Helmholtz gap and the diffuse layer, which compensation
does not touch, so it looks factor-independent and hides the effect.
"""
driving_force(result, m) =
    electrode_potential(result, m) .- m.elydata.ϕ_pzc .- reaction_plane_potential(result, m)

"""
    double_layer_capacitance(result, scanrate; cap_scale = cm^2 / 1.0e-6)

Capacitance from the capacitive current, `C = j_cap / (dU/dt)`.

Divides by the **signed** sweep rate taken from the applied protocol, not by the positive
constant `scanrate`: dU/dt flips at every vertex, so dividing by the constant mirrors one
half of the sweep about zero and the capacitance comes out antisymmetric instead of as a
single curve. `scanrate` is used only to set the threshold below which a sample sits too
close to a vertex to be meaningful; those become `NaN`.

Only interpretable where no faradaic current flows — once the surface reaction turns over,
the coverages change and `j_cap / v` stops being a double-layer capacitance. Mask it with
the faradaic current before reading values off.
"""
function double_layer_capacitance(result, scanrate; cap_scale = cm^2 / 1.0e-6)
    U_we = result.sawtooth
    t = result.times
    dUdt = let d = diff(U_we) ./ diff(t)
        vcat(d, d[end])
    end
    C = capacitive_current(result) ./ dUdt .* cap_scale
    return [abs(v) < 0.5 * abs(scanrate) ? NaN * one(c) : c for (v, c) in zip(dUdt, C)]
end

"""
    cv_abscissa(result; kind = :applied)

Voltage axis for a CV plot, together with a label naming which voltage it is.

Three different voltages live in these results and they are not interchangeable:
`:applied` is the programmed protocol (`result.sawtooth`), `:reaction_plane` is
`result.voltages`, and `:dl` is `result.dlvoltages`, the drop across the double layer.
Compare against experiment on `:applied` — that is the potential a potentiostat controls.
"""
function cv_abscissa(result; kind = :applied)
    if kind === :applied
        hasproperty(result, :sawtooth) ||
            error("result has no `sawtooth` field; use kind = :reaction_plane or :dl")
        return result.sawtooth, lab_voltage
    elseif kind === :reaction_plane
        return result.voltages, lab_voltage_rp
    elseif kind === :dl
        return result.dlvoltages, lab_voltage_dl
    else
        error("unknown abscissa kind: $kind (use :applied, :reaction_plane or :dl)")
    end
end

# ── moved from cell 7bfe397e (CV_dsp_cap_result) ──
function CV_dsp_cap_result(result, m; scale = cm^2 / mA, show_diff = false)
    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (540, 440))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
        return f, a
    end

    j_dsp = result.j_dsp .* scale
    j_cap = result.j_cap .* scale

    lines!(ax, result.voltages, j_dsp; color = :magenta, label = rich(rich("j", font = :bold_italic), subscript("dsp")))
    lines!(ax, result.voltages, j_cap; color = :skyblue, label = rich(rich("j", font = :bold_italic), subscript("cap")))

    if show_diff
        lines!(
            ax, result.voltages, j_dsp .- j_cap;
            color = :orange, linestyle = :dash, label = rich(rich("j", font = :bold_italic), subscript("dsp"), " − ", rich("j", font = :bold_italic), subscript("cap"))
        )
    end

    axislegend(ax; position = :rt)
    return fig
end

# ── moved from cell ad3a5236 (CV_total_current) ──
function CV_total_current(result, m; species = nothing, scale = cm^2 / mA)
    model = m.elydata
    sp = (species === nothing) ? model.cspecies[1] : species

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (540, 440))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
        return f, a
    end

    #i_F   = faradaic_current(result; species = sp) .* scale
    i_cap = capacitive_current(result) .* scale
    #i_tot = i_F .+ i_cap

    # lines!(ax, result.voltages, i_F;   color = :magenta,  label = L"i_F")
    lines!(ax, result.voltages, i_cap; color = :skyblue, label = rich(rich("i", font = :bold_italic), subscript("cap")))
    # lines!(ax, result.voltages, i_tot; color = :orange,   label = L"i_{tot}")
    axislegend(ax; position = :rt)
    return fig
end

# ── moved from cell 9949de26 (plot_cv_total_current_tot) ──
function plot_cv_total_current_tot(
        result, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        redox_species = nothing,      # back-compat: Dict(species => n_e), first entry wins
        co_idx = nothing,
        include_capacitive::Bool = true,
        scale = cm^2 / mA,
        color_F = colorant"#F2728A",   # i_F color
        color_C = colorant"#5BA8E8",   # i_C color
        mix_mode::Symbol = :mean,         # :sum or :mean
        lw = LW_LINE
    )
    sp, ne = if redox_species === nothing
        species, n_e
    else
        first(pairs(redox_species))
    end

    # ---- Faradaic and capacitive current, from the shared definitions ----
    I_F = faradaic_current(result; species = sp, n_e = ne, sgn = sgn)
    I_C = include_capacitive ? capacitive_current(result) : zero(I_F)

    # ---- Scale ----
    I_F_scaled = I_F .* scale
    I_C_scaled = I_C .* scale
    I_total = (I_F .+ I_C) .* scale

    # ---- blended color: RGB average or sum ----
    cF = RGBf(color_F); cC = RGBf(color_C)
    color_tot = mix_mode === :mean ?
        RGBf((cF.r + cC.r) / 2, (cF.g + cC.g) / 2, (cF.b + cC.b) / 2) :
        RGBf(min(cF.r + cC.r, 1.0f0), min(cF.g + cC.g, 1.0f0), min(cF.b + cC.b, 1.0f0))

    # ---- Plot ----
    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (900, 800))
        ax_a = Axis(f[1, 2], ylabel = rich(rich("i", font = :bold_italic), subscript("F"), "  (mA cm", superscript("−2"), ")"))
        ax_b = Axis(f[2, 2], ylabel = rich(rich("i", font = :bold_italic), subscript("C"), "  (mA cm", superscript("−2"), ")"))
        ax_c = Axis(
            f[3, 2], ylabel = rich(rich("i", font = :bold_italic), subscript("tot"), "  (mA cm", superscript("−2"), ")"),
            xlabel = lab_voltage
        )

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
        Label(
            f[i, 1], labels[i],
            fontsize = 24,
            font = :bold,
            halign = :left,
            valign = :top,
            padding = (15, 0, 0, 15)
        )
    end

    # ---- draw data lines ----
    lines!(ax_a, result.voltages, I_F_scaled; color = color_F, linewidth = lw)
    lines!(ax_b, result.voltages, I_C_scaled; color = color_C, linewidth = lw)
    lines!(ax_c, result.voltages, I_total; color = color_tot, linewidth = lw)
    rowgap!(f.layout, 15)

    rowsize!(f.layout, 1, Relative(0.3))
    rowsize!(f.layout, 2, Relative(0.3))
    rowsize!(f.layout, 3, Relative(0.4))
    return f
end

# ── moved from cell 02a78c8d (plot_cv_scanrate_grid) ──
"""
    plot_cv_scanrate_grid(result_vec, m; scanrates, kwargs...)

The current split into its faradaic and capacitive terms, one column per scan rate.

Rows are `I_F`, `I_C` and their sum; there is one column per entry of `result_vec`, headed
by the matching entry of `scanrates`. Rows share a current scale so a term can be read
across scan rates, columns share a voltage scale so the three terms of one sweep line up.

The split is the point of the figure: the faradaic term is nearly scan-rate independent
while the capacitive term grows in proportion to the rate, so the two are separable only by
plotting them apart. At the slowest rate `I_C` is invisible next to `I_F`; by the fastest it
dominates the total.

`I_F` is [`faradaic_current`](@ref) and `I_C` is [`capacitive_current`](@ref) — the same
definitions every other CV figure in this package uses.
"""
function plot_cv_scanrate_grid(
        result_vec, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        abscissa = :applied,
        redox_species = nothing,   # back-compat: Dict(species => n_e), first entry wins
        co_idx = nothing,
        scanrates = [0.05, 0.5, 5.0],
        include_capacitive::Bool = true,
        scale = cm^2 / mA,
        color_F = colorant"#F2728A",
        color_C = colorant"#5BA8E8",
        mix_mode::Symbol = :mean,
        lw = LW_LINE,
        # Below the theme's `FS_LABEL`: a rotated label is bounded by the row height, and a
        # row of a three-high stack is well under half a panel.
        ylabelsize = 20,
        rowgap_px = 6,
        colgap_px = 10,
    )
    sp, ne = redox_species === nothing ? (species, n_e) : first(pairs(redox_species))

    # One column per result. The loops used to be hard-coded to `1:3` while `scanrates`
    # was a keyword, so passing four sweeps silently plotted three of them.
    ncol = length(result_vec)
    ncol == 0 && error("`result_vec` is empty: nothing to plot")
    length(scanrates) < ncol && error(
        "got $(ncol) results but only $(length(scanrates)) scan rates to label them with"
    )

    ylabels = [lab_current_F, lab_current_C, lab_current_tot]
    # The voltage label has to follow `abscissa`; hardcoding it would label a
    # double-layer or reaction-plane potential as if it were the applied one.
    lab_U = cv_abscissa(result_vec[1]; kind = abscissa)[2]

    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (330 * ncol, 700))
        axes_matrix = Matrix{Axis}(undef, 3, ncol)

        for col in 1:ncol, row in 1:3
            axes_matrix[row, col] = row == 3 ?
                Axis(f[row, col + 1], xlabel = lab_U) :
                Axis(f[row, col + 1])

            row < 3 && hidexdecorations!(axes_matrix[row, col]; grid = false)

            if col == 1
                axes_matrix[row, col].ylabel = ylabels[row]
                axes_matrix[row, col].ylabelsize = ylabelsize
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

        # Rows share a current scale so a term can be compared across scan rates; columns
        # share a voltage scale so the three terms of one sweep line up vertically.
        for row in 1:3
            ncol > 1 && linkyaxes!(axes_matrix[row, :]...)
        end
        for col in 1:ncol
            linkxaxes!(axes_matrix[:, col]...)
        end

        return f, axes_matrix
    end
    f, axes_matrix = fig

    for row in 1:3, col in 1:ncol
        Label(
            f[row, col + 1], "($('a' + (row - 1) * ncol + col - 1))",
            fontsize = FS_TITLE, font = :bold, halign = :left, valign = :top,
            padding = (15, 0, 0, 15)
        )
    end

    for col in 1:ncol
        Label(
            f[0, col + 1],
            rich(string(scanrates[col]), " V s", superscript("−1"));
            fontsize = FS_TITLE, font = :bold, halign = :center, valign = :bottom
        )
    end

    cF, cC = RGBf(color_F), RGBf(color_C)
    color_tot = mix_mode === :mean ?
        RGBf((cF.r + cC.r) / 2, (cF.g + cC.g) / 2, (cF.b + cC.b) / 2) :
        RGBf(min(cF.r + cC.r, 1.0f0), min(cF.g + cC.g, 1.0f0), min(cF.b + cC.b, 1.0f0))

    for col in 1:ncol
        res = result_vec[col]

        I_F = faradaic_current(res; species = sp, n_e = ne, sgn = sgn)
        I_C = include_capacitive ? capacitive_current(res) : zero(I_F)
        U, _ = cv_abscissa(res; kind = abscissa)

        lines!(axes_matrix[1, col], U, I_F .* scale; color = color_F, linewidth = lw)
        lines!(axes_matrix[2, col], U, I_C .* scale; color = color_C, linewidth = lw)
        lines!(axes_matrix[3, col], U, (I_F .+ I_C) .* scale;
               color = color_tot, linewidth = lw)
    end

    # Tight: rows 1-2 have their x decorations hidden and columns 2+ their y decorations,
    # so there is nothing between neighbouring panels for a gap to keep apart.
    rowgap!(f.layout, rowgap_px)
    colgap!(f.layout, colgap_px)

    rowsize!(f.layout, 1, Relative(0.29))
    rowsize!(f.layout, 2, Relative(0.29))
    rowsize!(f.layout, 3, Relative(0.38))

    # Column 1 holds only the y labels; the data columns share the rest equally.
    colsize!(f.layout, 1, Auto())
    for col in 2:(ncol + 1)
        colsize!(f.layout, col, Relative(1 / ncol * 0.98))
    end

    return f
end

# ── moved from cell 08756476 (plot_cv_scanrate_grid_unc) ──
function plot_cv_scanrate_grid_unc(
        result_vec, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        abscissa = :applied,
        redox_species = nothing,   # back-compat: Dict(species => n_e), first entry wins
        co_idx = nothing,
        scanrates = [0.05, 0.5, 5.0],
        include_capacitive::Bool = true,
        scale = cm^2 / mA,
        color_tot = "#BAC8FF",
        lw = LW_LINE
    )
    sp, ne = redox_species === nothing ? (species, n_e) : first(pairs(redox_species))

    # shared labels from struct.jl; the voltage one follows `abscissa`
    lab_I = lab_current
    lab_U = isempty(result_vec) ? lab_voltage : cv_abscissa(result_vec[1]; kind = abscissa)[2]

    fig = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 350))
        axes_vec = Vector{Axis}(undef, 3)

        for col in 1:3
            axes_vec[col] = Axis(f[1, col + 1], xlabel = lab_U)

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
        Label(
            f[1, col + 1], grid_labels[col],
            fontsize = 24, font = :bold, halign = :left, valign = :top,
            padding = (15, 0, 0, 15)
        )
    end

    for col in 1:3
        Label(
            f[0, col + 1],
            rich(string(scanrates[col]), "  V s", superscript("−1"));
            fontsize = 26, font = :bold, halign = :center, valign = :bottom
        )
    end

    for col in 1:3
        res = result_vec[col]
        n_t = length(res.voltages)

        I_F = faradaic_current(res; species = sp, n_e = ne, sgn = sgn)
        I_C = include_capacitive ? capacitive_current(res) : zero(I_F)

        I_total = (I_F .+ I_C) .* scale
        U, _ = cv_abscissa(res; kind = abscissa)
        lines!(axes_vec[col], U, I_total; color = color_tot, linewidth = lw)
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
        species = ico,
        n_e = 2,
        sgn = 1,
        abscissa = :applied,
        redox_species = nothing,   # back-compat: Dict(species => n_e), first entry wins
        include_capacitive::Bool = true,
        scale = cm^2 / mA,
        default_lw = LW_LINE
    )
    sp, ne = redox_species === nothing ? (species, n_e) : first(pairs(redox_species))
    n = length(sweep_vec)
    pastel2 = cgrad([colorant"#D5F011", colorant"#16D8FF"])
    cols = [pastel2[t] for t in range(0, 1, length = max(n, 1))]

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(sweep_vec) ? lab_voltage : cv_abscissa(sweep_vec[1]; kind = abscissa)[2]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (850, 550))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = xlab
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        I_F = faradaic_current(rec; species = sp, n_e = ne, sgn = sgn)
        I_C = include_capacitive ? capacitive_current(rec) : zero(I_F)
        n_t = length(I_F)

        # 3. total current
        I_total = (I_F .+ I_C) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)

        lines!(
            ax, U, I_total;
            linewidth = default_lw,
            color = cols[j]
        )

        pos_y = if j == 5
            -2.5
        elseif j == 3
            -1.5
        else
            -0.5 * j
        end

        text!(
            ax, "$(scanrates[j]) V/s";
            position = (-1.55, pos_y),
            color = cols[j],
            fontsize = 18,
            font = :bold
        )
    end

    return fig
end

# ── moved from cell ae50302f (plot_scanrate_sweeps_cv) ──
function plot_scanrate_sweeps_cv(
        sweep_vec, scanrates;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        scale = cm^2 / mA,
        default_lw = LW_LINE
    )
    n = length(sweep_vec)
    pastel2 = cgrad([colorant"#D5F011", colorant"#16D8FF"])
    cols = [pastel2[t] for t in range(0, 1, length = max(n, 1))]

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(sweep_vec) ? lab_voltage : cv_abscissa(sweep_vec[1]; kind = abscissa)[2]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (850, 550))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = xlab
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)
        lines!(
            ax, U, I;
            linewidth = default_lw,
            color = cols[j]
        )
    end

    return fig
end

# ── moved from cell 51064823 (plot_scanrate_sweeps_cap) ──
function plot_scanrate_sweeps_cap(
        sweep_vec, scanrates;
        #species=iohminus,
        scale = cm^2 / mA,
        default_lw = LW_LINE
    )
    n = length(sweep_vec)

    cols = [
        RGBf(
                0.1 + 0.6 * (i / n),
                0.2 + 0.6 * (1 - i / n),
                0.7 - 0.3 * (i / n)
            ) for i in 1:n
    ]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1000, 550))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
        return f, a
    end

    for (j, rec) in enumerate(sweep_vec)
        color_j = cols[j]
        lw_j = default_lw
        J_cap = rec.j_cap .* scale
        lines!(ax, rec.voltages, rec.j_cap; linewidth = lw_j, color = color_j)
    end

    if n > 1
        cgradient = cgrad(cols, categorical = true)

        Colorbar(
            fig[1, 2],
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
    max_L = 1.0e-2

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

        conc_matrix = [
            (tsol[co2_idx, ix, it] / (mol / dm^3))
                for ix in 1:length(X_coords), it in 1:length(times)
        ]

        ax = Axis(
            fig[idx, 1],
            xlabel = idx == num_panels ? "Distance from electrode [m]" : "",
            ylabel = "Time [s]",
            title = "Boundary Layer (L) = $(round(L_val / μm)) μm",
            xscale = log10,
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(9)
        )

        hm = heatmap!(
            ax, X_coords .+ 1.0e-12, times, conc_matrix;
            colorrange = color_range,
            colormap = discrete_cmap,
            interpolate = false
        )

        vlines!(ax, [L_val], color = :red, linestyle = :dash, linewidth = LW_GUIDE)

        Colorbar(fig[idx, 2], hm, label = rich("log", subscript("10"), "(", rich("c", font = :bold_italic), subscript(rich("CO", subscript("2"))), ")"))
    end

    rowgap!(fig.layout, 35)

    return fig
end

# ── moved from cell d91bcc6a (co2_log_contour) ──
function co2_log_contour(
        results::Dict, grid_dict, m;
        scale = mol / dm^3, num_levels = 24, L_val = nothing
    )
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

    M = [
        log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / scale, 1.0e-12))
            for ix in 1:length(X), it in 1:length(times)
    ]

    c_min = 0.0
    c_max = log10(result.tsol[co2_idx, end, 1] / scale) - log10(1.0e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    f = Figure(size = (800, 450))

    ax = Axis(
        f[1, 1];
        xlabel = lab_time,
        ylabel = rich(rich("x", font = :bold_italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (10.0 .^ (-12:3:-6), [powlab(-12), powlab(-9), powlab(-6)])
    )

    hm = heatmap!(
        ax, times, X .+ 1.0e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false
    )

    cb = Colorbar(
        f[1, 2], hm;
        label = rich(
            "log", subscript("10"), "(", rich("c", font = :bold_italic), subscript("bulk"),
            ") − log", subscript("10"), "(", rich("c", font = :bold_italic), subscript(rich("CO", subscript("2"))), ")"
        ),
        ticklabelsize = 20,
        labelsize = 20
    )

    return f
end

# ── moved from cell e6f43f01 (plot_activity_time_electrode) ──
function plot_activity_time_electrode(
        result, m;
        model_type = "DGML_γ", # "DGML_γ" or "Stefan_γ"
        nspecies = 7,
        scale = (mol / dm^3),
        ipressure = nothing
    )
    bulk = m.bulk
    electrolyte = m.elydata
    names = getproperty.(bulk, :name)
    colors = getproperty.(bulk, :color)

    times = result.tsol.t
    nt = length(times)

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
                term_press = exp((1.0 - size_ratio) * pnode / (bar_c * RT))
                term_steric = solvent_frac^(-size_ratio)
                a_thermo = term_conc * term_press * term_steric
            elseif model_type == "Stefan_γ"
                a_thermo = term_conc * (solvent_frac^(-1.0))
            else
                a_thermo = term_conc
            end

            activity_electrode[i, t] = a_thermo * (bar_c * c_scale)
        end
    end

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (960, 540))

        model_title = model_type == "DGML_γ" ? "Modified DGML Model (Activity)" :
            model_type == "Stefan_γ" ? "Stefan Model (Activity)" : "Ideal Solution"

        a = Axis(
            f[1, 1],
            title = model_title,
            xlabel = lab_time,
            ylabel = rich(rich("a", font = :bold_italic), subscript("i, electrode")),
            limits = ((times[1] - (times[end] / 200), times[end] + (times[end] / 100)), (1.0e-12, 1.0e4)),
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

        lines!(ax, times, y_fixed; color = colors[i])
    end

    return fig
end

# ── moved from cell 96c90e9f (CV_overlay_currents) ──
function CV_overlay_currents(
        results, m;
        labels = nothing,
        species = ico,
        n_e = 2,
        sgn = 1,
        abscissa = :applied,
        scale = cm^2 / mA,
        linestyles = [:solid, :dash, :dot]
    )
    labels = labels === nothing ? ["result $i" for i in 1:length(results)] : labels

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(results) ? lab_voltage : cv_abscissa(first(results); kind = abscissa)[2]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (620, 460))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = xlab
        )
        return f, a
    end

    col_F = :magenta
    col_cap = :skyblue
    col_tot = :orange

    for (k, result) in enumerate(results)
        ls = linestyles[mod1(k, length(linestyles))]

        i_F = faradaic_current(result; species, n_e, sgn) .* scale
        i_cap = capacitive_current(result) .* scale
        i_tot = i_F .+ i_cap
        U, _ = cv_abscissa(result; kind = abscissa)

        lines!(
            ax, U, i_F; color = col_F, linestyle = ls,
            label = rich(rich("i", font = :bold_italic), subscript("F"), " (", string(labels[k]), ")")
        )
        lines!(
            ax, U, i_cap; color = col_cap, linestyle = ls,
            label = rich(rich("i", font = :bold_italic), subscript("cap"), " (", string(labels[k]), ")")
        )
        lines!(
            ax, U, i_tot; color = col_tot, linestyle = ls,
            label = rich(rich("i", font = :bold_italic), subscript("tot"), " (", string(labels[k]), ")")
        )
    end

    axislegend(ax; position = :rt, nbanks = 2)
    return fig
end

# ── moved from cell a2264942 (plot_cv_current_variedL) ──
function plot_cv_current_variedL(
        results, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        scale = cm^2 / mA,
        color_gradient = true,
        linewidth = LW_LINE,
        title = ""
    )
    Lkeys = sort(collect(keys(results)))
    n = length(Lkeys)

    # label comes from the abscissa choice, so the axis always names what it shows
    _, xlab = cv_abscissa(results[first(Lkeys)]; kind = abscissa)

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (800, 400))
        a = Axis(
            f[1, 1];
            title = title,
            ylabel = lab_current,
            xlabel = xlab
        )
        return f, a
    end

    cols = if color_gradient
        [
            RGBf(
                    0.6 + 0.3 * (i - 1) / max(n - 1, 1),
                    0.8 - 0.2 * (i - 1) / max(n - 1, 1),
                    0.7 - 0.2 * (i - 1) / max(n - 1, 1)
                ) for i in 1:n
        ]
    else
        fill(RGBf(0.13, 0.13, 0.13), n)
    end

    for (i, L) in enumerate(Lkeys)
        I = cv_current(results[L]; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(results[L]; kind = abscissa)
        lines!(
            ax, U, I;
            color = cols[i], linewidth = linewidth,
            label = @sprintf("%g μm", L / μm)
        )
    end

    axislegend(ax, "L", position = :rb, framevisible = false, legendtext = 24)
    return fig
end

# ── moved from cell ed4451bd (plot_pressure_varied_sweep_ivc) ──
function plot_pressure_varied_sweep_ivc(
        P_recs;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        scale = cm^2 / mA,
        limits = nothing,
        sim_linewidth::Real = LW_LINE
    )
    n = length(P_recs)
    pastel3 = cgrad([colorant"#FFB3BA", colorant"#A3D8FF"])
    cols_sim = [pastel3[t] for t in range(0, 1, length = max(n, 1))]

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(P_recs) ? lab_voltage : cv_abscissa(P_recs[1][2]; kind = abscissa)[2]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (1050, 500))
        a = if limits !== nothing
            Axis(f[1, 1], xlabel = xlab, ylabel = lab_current, limits = limits)
        else
            Axis(f[1, 1], xlabel = xlab, ylabel = lab_current)
        end
        return f, a
    end

    plot_objs = []
    labels = String[]

    for j in 1:n
        p, rec = P_recs[j]
        # the stoichiometric factor lives in `n_e` now; there is no stray /2 here
        I = cv_current(rec; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(rec; kind = abscissa)

        label_text = "$(p) atm"

        hl = lines!(
            ax, U, I;
            color = cols_sim[j],
            linewidth = sim_linewidth
        )

        push!(plot_objs, hl)
        push!(labels, label_text)
    end

    if n > 0
        leg = Legend(
            fig[1, 1], plot_objs, labels, rich(rich("p", font = :bold_italic), subscript(rich("CO", subscript("2"))));
            framevisible = false,
            halign = :right, valign = :bottom,
            labelsize = 20, titlesize = 23,
            padding = (0, 0, 0, 0),
            tellwidth = false, tellheight = false
        )

        translate!(leg.blockscene, -40, 40, 0)
    end

    return fig
end

# ── moved from cell 6f6e779c (plot_cv_current_dict) ──
function plot_cv_current_dict(
        result_dict, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        abscissa = :applied,
        scale = 1.0,
        title = "",
        linewidth = LW_LINE
    )
    sorted_keys = sort(collect(keys(result_dict)))

    # axis label follows the abscissa choice, so it always names what it shows
    xlab = isempty(sorted_keys) ? lab_voltage :
        cv_abscissa(result_dict[first(sorted_keys)]; kind = abscissa)[2]

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (400, 300))
        a = Axis(
            f[1, 1];
            title = title,
            ylabel = lab_current,
            xlabel = xlab
        )
        return f, a
    end

    for (i, key) in enumerate(sorted_keys)
        result = result_dict[key]
        # the stoichiometric factor lives in `n_e` now; there is no stray /2 here
        I = cv_current(result; species, n_e, sgn, include_capacitive) .* scale
        U, _ = cv_abscissa(result; kind = abscissa)
        lines!(
            ax, U, I;
            linewidth = linewidth,
            label = @sprintf("L = %g μm", key / μm)
        )
    end

    axislegend(ax, position = :rb, framevisible = true)
    return fig
end

# ── moved from cell 13a0e5da (panel_conc_time!) ──
function panel_conc_time!(
        fig, panel_pos, result, m; nspecies = 7, scale = mol / dm^3, lw = LW_LINE,
        xlabel = lab_time, legend_pos = nothing
    )
    bulk = m.bulk
    sp_colors = [
        "#E07B39", "#888888", "#7B5C3E", "#222222",
        "#C0392B", "#27AE60", "#2980B9",
    ]
    colors = sp_colors[1:nspecies]
    names = getproperty.(bulk, :name)
    times = result.tsol.t
    nt = length(times)
    conc = [result.tsol[i, 1, t] / scale for i in 1:nspecies, t in 1:nt]

    ax_c = Axis(
        panel_pos;
        xlabel = xlabel,
        ylabel = rich(
            "Surface Concentration ",
            # `α` runs over species, so it is a running index, not a descriptive tag:
            # ISO 80000-1 sets those in italic (unlike the upright `M` of `U_M`). Plain
            # italic rather than bold italic keeps the index subordinate to the `c`.
            rich("c", font = :bold_italic), subscript("α", font = :italic),
            superscript("‡"),
            "\n(M)"
        ),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9)
    )

    ax_c.yticks = (
        10.0 .^ (4:-4:-12),
        [powlab(4), powlab(0), powlab(-4), powlab(-8), powlab(-12)],
    )

    for i in 1:nspecies
        lines!(
            ax_c, times, max.(conc[i, :], eps(Float64));
            color = colors[i], linewidth = lw
        )
    end

    legend_elements = [ [LineElement(color = colors[i], linewidth = lw)] for i in 1:nspecies ]
    legend_labels = [ rich(string(names[i]), color = colors[i]) for i in 1:nspecies ]

    leg = Legend(
        legend_pos === nothing ? fig[1, 3] : legend_pos,
        legend_elements, legend_labels;
        framevisible = false,
        labelsize = 20
    )

    return ax_c, leg
end

# ── moved from cell 13a0e5da (panel_time_current!)  [ely = elydata_Gold_unc → ely = model] ──
function panel_time_current!(
        fig, panel_pos, result, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        redox_species = nothing,       # back-compat: overrides `species` when given
        include_capacitive = true,
        scale = cm^2 / mA,
        color_gradient = true,
        lw = LW_LINE,
        xlabel = lab_time
    )
    ax = Axis(panel_pos; xlabel = xlabel, ylabel = lab_current)

    sp = redox_species === nothing ? species : redox_species
    I = cv_current(result; species = sp, n_e, sgn, include_capacitive) .* scale

    cols = color_gradient ? :skyblue : :black
    lines!(ax, result.times, I; color = cols, linewidth = lw)

    return ax
end

"""
    panel_time_current_diff!(fig, panel_pos, result_1, result_2, m; kwargs...)

Difference between two current traces against time, `I₁(t) − I₂(t)`.

Written for compensated-vs-uncompensated pairs, where the two currents overlap almost
everywhere and the interesting quantity is the residual: the difference resolves a gap
that is invisible when the traces are drawn on top of each other.

The two sweeps do not share a time grid — an adaptive solver puts its steps wherever the
run needed them, so `result_1.times != result_2.times` even for the same protocol.
`result_2` is therefore linearly interpolated onto `result_1.times`. `Line()`
extrapolation covers the ends, so a `result_2` that stops short of `result_1` is
continued along its last slope rather than erroring; keep the two protocols identical if
that matters.

Both currents come from [`cv_current`](@ref), so this panel carries the same definition as
every other CV plot here. `redox_species` overrides `species`, matching
[`panel_time_current!`](@ref).

`xlims` / `ylims` are applied only when given.
"""
function panel_time_current_diff!(
        fig, panel_pos, result_1, result_2, m;
        species = ico,
        n_e = 2,
        sgn = -1,
        redox_species = nothing,       # back-compat: overrides `species` when given
        include_capacitive = true,
        scale = cm^2 / mA,
        lw = LW_LINE,
        color = colorant"#C4844C",
        xlabel = lab_time,
        label_1 = "odr",
        label_2 = "unc",
        title = "$(label_1) − $(label_2)",
        xlims = nothing,
        ylims = nothing,
    )
    ax = Axis(
        panel_pos; xlabel = xlabel, title = title,
        # Δ is an operator, not a quantity, so it stays upright next to the italic symbol.
        ylabel = rich("Current Difference ",rich("Δ", font = :bold), rich("I", font = :bold_italic),
                      "\n(mA cm", superscript("−2"), ")"),
    )

    sp = redox_species === nothing ? species : redox_species
    # `sgn` multiplies the *total*, not just the faradaic part as `cv_current`'s own `sgn`
    # would: it is a display convention (cathodic down), and flipping only one of the two
    # terms would change their sum.
    current(res) = sgn .* cv_current(res; species = sp, n_e, include_capacitive) .* scale

    t1, I_1 = result_1.times, current(result_1)
    t2, I_2 = result_2.times, current(result_2)

    itp2 = Interpolations.linear_interpolation(
        t2, I_2; extrapolation_bc = Interpolations.Line()
    )
    lines!(ax, t1, I_1 .- itp2.(t1); color = color, linewidth = lw)

    xlims === nothing || xlims!(ax, xlims...)
    ylims === nothing || ylims!(ax, ylims...)

    return ax
end

"""
    panel_time_voltage!(fig, panel_pos, result; m = nothing, lw, xlabel)

Electrode potential against time — the same quantity as the `:metal_time` panel of
[`plot_ircomp_compare`](@ref), with the applied protocol drawn dotted behind it.

Given `m`, the curve is [`electrode_potential`](@ref), which resolves to `ϕ_we + ϕ_DL`
under `OhmicDropEstimation` and to the applied protocol itself otherwise. Without `m` it
falls back to `result.voltages`, which means *different things in different compensation
modes* — reaction-plane potential without compensation, compensated applied potential with
it — so pass the model whenever runs are compared.
"""
function panel_time_voltage!(
        fig, panel_pos, result;
        m = nothing, lw = LW_LINE, xlabel = lab_time,
        color = parse(Colorant, "#D7C2F0"),
    )
    ax = Axis(
        panel_pos; xlabel = xlabel,
        ylabel = rich("Electrode Potential ", rich("U", font = :bold_italic),
                      subscript("M"), "\n(V vs. SHE)")
    )
    U = m === nothing ? result.voltages : electrode_potential(result, m)
    # the protocol is the reference the curve is read against; the gap between them is
    # the compensation term, which is invisible unless both are on the axis
    hasproperty(result, :sawtooth) && lines!(
        ax, result.times, result.sawtooth;
        color = (:black, 0.45), linestyle = :dot, linewidth = LW_GUIDE
    )
    lines!(ax, result.times, U; color = color, linewidth = lw)
    return ax
end

# ── moved from cell 13a0e5da (panel_co2_log_contour!) ──
function panel_co2_log_contour!(
        fig, panel_pos, cbar_pos, result, X, m;
        scale = mol / dm^3, num_levels = 24, L_val = nothing
    )
    bulk = m.bulk
    co2_idx = findfirst(s -> s.name == "CO₂", bulk)
    times = result.tsol.t

    log_c_bulk = log10(result.tsol[co2_idx, end, 1] / scale)

    M = [
        log_c_bulk - log10(max(result.tsol[co2_idx, ix, it] / scale, 1.0e-12))
            for ix in 1:length(X), it in 1:length(times)
    ]

    c_min = 0.0
    c_max = log10(result.tsol[co2_idx, end, 1] / scale) - log10(1.0e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    ax = Axis(
        panel_pos;
        xlabel = lab_time,
        ylabel = rich(rich("x", font = :bold_italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (
            10.0 .^ (-12:3:-6),
            [powlab(-12), powlab(-9), powlab(-6)],
        )
    )

    hm = heatmap!(
        ax, times, X .+ 1.0e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false
    )

    Colorbar(
        cbar_pos, hm;
        label = rich(
            "log", subscript("10"), "(", rich("c", font = :bold_italic),
            subscript("bulk"), ") − log", subscript("10"), "(",
            rich("c", font = :bold_italic),
            subscript(rich("CO", subscript("2"))), ")"
        ),
        ticklabelsize = 20, labelsize = 20
    )
    return ax, hm
end

# ── moved from cell 13a0e5da (panel_time_ph!) ──
function panel_time_ph!(
        fig, panel_pos, result; ihplus = 2, scale = mol / dm^3, lw = LW_LINE,
        xlabel = lab_time, color = parse(Colorant, "#7BB661")
    )
    times = result.tsol.t
    nt = length(times)
    cH = [result.tsol[ihplus, 1, t] / scale for t in 1:nt]   # surface (node 1) H⁺ concentration (M)
    pH = -log10.(max.(cH, eps(Float64)))
    # pH is upright by convention — it is not a variable symbol, so it takes no italic
    ax = Axis(panel_pos; xlabel = xlabel, ylabel = rich("Surface pH"))
    lines!(ax, times, pH; color = color, linewidth = lw)
    return ax
end


"""
    plot_ircomp_compare(results, m; reference = nothing, kwargs...)

Compare a family of ohmic-drop compensation factors — the `Dict(factor => result)` that
`cvsweep_over_ircompfactor` returns — in three panels.

`panels` selects which of the three to draw, in the order given:

| `panels` entry | Quantity | Expected behaviour |
|:--|:--|:--|
| `:driving` | `φ_M - φ_PET`, the driving force for electron transfer | approaches the dotted unit-slope line as the factor rises — the reference of Levey et al., Fig. 3 and 5c |
| `:metal` | `φ_M`, the potential imposed at the electrode | rides on the dotted y = x line, displaced by ϕ_DL. This is the "applied and electrode potential track each other" view |
| `:ircomp` | ϕ_DL, the potential the compensation added | **proportional to the factor**, peaking where the current peaks, and identically zero for the uncompensated run |
| `:cv` | the voltammogram itself | the physical consequence: more compensation, more overpotential at the interface, larger and earlier peaks |
| `:metal_time` | `φ_M` against **time** | tracks the applied protocol, drawn as a dotted reference, offset by ϕ_DL |
| `:rp_time` | `ϕ(0)` against **time** | strongly attenuated against the same reference — the capacitive split, not a compensation effect |

The last two put time on the abscissa; the rest use the applied potential. Mixing them in
one call is fine, each axis is labelled for what it carries.

The default is `(:driving, :metal_time, :ircomp)`: what compensation did to the driving
force, whether the electrode still follows the protocol, and how much potential was put
back. Drop to one with `panels = (:driving,)`, or ask for any of the others.

Panels sit in one row unless `layout = (rows, cols)` is given — `layout = (2, 2)` with four
panels wraps them into a square, which keeps the axis labels legible where a single row of
four does not. The figure size follows the grid; override it with `fig_size` and the label
size with `labelsize`.

Pass an uncompensated run as `reference` to overlay it as a dashed black line, **together
with the model it was built from** as `reference_model`. It should sit on top of the
`factor = 0` curve in every panel — at zero factor the two are the same problem, so a
visible gap means the runs were integrated on different time grids.

!!! warning "`reference_model` is not optional in practice"
    Without it the reference is read through `m`, whose electrolyte carries
    `OhmicDropEstimation`, so `electrode_potential` takes `reference.voltages` — which
    under `NoIRCompensation` holds ϕ(0), not φ_M. The driving force then collapses to the
    constant `-ϕ_pzc` and ϕ_DL becomes a spurious straight line.

!!! note "How large should the effect be?"
    ϕ_DL is bounded by `factor · Ru · i`. With `Ru = L/σ` and a boundary layer of a few
    hundred μm this is tens of millivolts against a sweep of volts, so panel (a) looks
    unchanged and panel (b) is where the factor is legible. Raising `L` raises `Ru`
    proportionally — and with it the risk that the positive feedback loop gain
    `factor · Ru · ∂j/∂η` reaches 1 and the solve diverges.
"""
function plot_ircomp_compare(
        results, m;
        reference = nothing,
        reference_model = nothing,
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = false,
        current_scale = cm^2 / mA,
        colormap = CMAP_IRCOMP,
        panels = (:driving, :metal_time, :ircomp),
        layout = nothing,
        fig_size = nothing,
        labelsize = 20,
        lw = LW_LINE,
        legend_title = "f",
    )
    nrow, ncol = layout === nothing ? (1, length(panels)) : layout
    nrow * ncol >= length(panels) ||
        error("layout $(layout) has room for $(nrow * ncol) panels but $(length(panels)) were asked for")
    figsize = fig_size === nothing ? (480 * ncol + 140, 430 * nrow) : fig_size
    ks = sort(collect(keys(results)))
    isempty(ks) && error("plot_ircomp_compare: no results")
    cols = [colormap[t] for t in range(0, 1, length = max(length(ks), 1))]
    ϕ_pzc = m.elydata.ϕ_pzc

    # every run shares the applied protocol; with a fixed grid they share it sample by
    # sample, which is what makes the curves directly comparable
    U_we = results[first(ks)].sawtooth

    # y-value and axis decoration per panel kind, so the layout below is just bookkeeping
    lab_UM = rich("Electrode Potential ", rich("U", font = :bold_italic), subscript("M"), "\n(V vs. SHE)")
    ylabel_of = Dict(
        :driving => rich("Driving Force ",
                         rich("U", font = :bold_italic), subscript("M"), " − ",
                         rich("U", font = :bold_italic), subscript("PET"), "\n(V)"),
        :metal => lab_UM,
        # subscript "DROP", not "DL": U_dl is already the double-layer voltage that
        # `cv_abscissa(:dl)` returns, and the two are unrelated quantities
        :ircomp => rich("Recovered IR Drop ",
                        rich("U", font = :bold_italic), subscript("DROP"), "\n(V)"),
        :cv => lab_current_co,
        :metal_time => lab_UM,
        :rp_time => rich("Reaction-Plane Potential ",
                         rich("U", font = :bold_italic), subscript("PET"), "\n(V vs. SHE)"),
    )
    title_of = Dict(
        :driving => "driving force",
        :metal => "electrode potential",
        :ircomp => "recovered IR drop",
        :cv => "voltammogram",
        :metal_time => "electrode potential vs. time",
        :rp_time => "reaction plane vs. time",
    )
    # panels whose abscissa is time rather than the applied potential
    is_time(kind) = kind in (:metal_time, :rp_time)

    yvalue(kind, r, mm) =
        kind === :driving ? driving_force(r, mm) :
        kind === :metal || kind === :metal_time ? electrode_potential(r, mm) :
        kind === :ircomp ? compensation_potential(r, mm) :
        kind === :rp_time ? reaction_plane_potential(r, mm) :
        kind === :cv ? cv_current(r; species, n_e, sgn, include_capacitive) .* current_scale :
        error("unknown panel :$kind (see the docstring for the list)")

    fig = Figure(size = figsize)
    axes = Dict{Symbol, Axis}()
    for (k, kind) in enumerate(panels)
        r, c = (k - 1) ÷ ncol + 1, (k - 1) % ncol + 1
        tag = length(panels) > 1 ? "($(('a':'z')[k])) " : ""
        axes[kind] = Axis(
            fig[r, c];
            xlabel = is_time(kind) ? lab_time : lab_voltage,
            ylabel = ylabel_of[kind],
            title = tag * title_of[kind],
            xlabelsize = labelsize, ylabelsize = labelsize,
        )
    end

    # Unit slope through the pzc: a cell with no loss between metal and reaction plane —
    # the "diffusion model" reference of Levey et al. Compensation pulls the curves
    # towards it.
    if haskey(axes, :driving)
        ablines!(axes[:driving], -ϕ_pzc, 1; color = (:black, 0.45), linestyle = :dot, linewidth = LW_GUIDE)
        # name the reference on the line itself; a legend entry for a guide line costs a
        # whole row and still leaves the reader matching dash patterns
        ulo, uhi = extrema(U_we)
        x_t = ulo + 0.68 * (uhi - ulo)
        text!(
            axes[:driving], x_t, x_t - ϕ_pzc;
            text = rich(rich("C", font = :bold_italic), subscript("gap"), " → ∞"),
            align = (:left, :bottom), offset = (6, 6), rotation = pi/4,
            fontsize = labelsize, color = (:black, 0.6),
        )
    end
    # For :metal the reference is y = x — the electrode potential with nothing added. The
    # curves ride just above or below it by ϕ_DL, which is why they look glued to it.
    haskey(axes, :metal) &&
        ablines!(axes[:metal], 0, 1; color = (:black, 0.45), linestyle = :dot, linewidth = LW_GUIDE)
    haskey(axes, :ircomp) &&
        hlines!(axes[:ircomp], [0.0]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    # On the time panels the applied protocol is the reference the curves are read
    # against, so draw it once rather than per factor.
    t_ref = results[first(ks)].times
    for kind in (:metal_time, :rp_time)
        haskey(axes, kind) && lines!(axes[kind], t_ref, U_we;
                                     color = (:black, 0.45), linestyle = :dot, linewidth = LW_GUIDE)
    end

    plots = Any[]
    labels = String[]
    for (i, f) in enumerate(ks)
        r = results[f]
        local line
        for kind in panels
            x = is_time(kind) ? r.times : U_we
            line = lines!(axes[kind], x, yvalue(kind, r, m); color = cols[i], linewidth = lw)
        end
        push!(plots, line)
        push!(labels, string(f))
    end

    if reference !== nothing
        # The reference is normally an uncompensated run, i.e. a *different* model, and it
        # has to be read with its own. Passing `m` here would make `electrode_potential`
        # take `reference.voltages`, which under NoIRCompensation holds ϕ(0) rather than
        # φ_M — the driving force then collapses to the constant −ϕ_pzc.
        mref = reference_model === nothing ? m : reference_model
        Uref = reference.sawtooth
        local line
        for kind in panels
            x = is_time(kind) ? reference.times : Uref
            line = lines!(axes[kind], x, yvalue(kind, reference, mref);
                          color = :black, linestyle = :dash, linewidth = LW_GUIDE)
        end
        push!(plots, line)
        push!(labels, "unc")
    end

    Legend(fig[1:nrow, ncol + 1], plots, labels, legend_title; framevisible = false)
    colgap!(fig.layout, 25)
    nrow > 1 && rowgap!(fig.layout, 20)
    return fig
end

"""
    plot_cv_summary(result, m; sawtooth = nothing, kwargs...)

The five-panel summary of a single CV run, stacked on a shared time axis.

| Panel | Content |
|:--|:--|
| (a) | concentrations at the electrode, log scale |
| (b) | current |
| (c) | electrode potential, with the applied protocol dotted behind it — the same quantity as the `:metal_time` panel of [`plot_ircomp_compare`](@ref) |
| (d) | surface pH |
| (e) | reaction quotient over equilibrium constant for the two buffer reactions |

Only the bottom panel carries a time axis; the rest have their x decorations hidden and are
linked to it, so a feature at one time lines up vertically across all five. That is the
point of the figure — the pH excursion in (d), the buffer going off equilibrium in (e) and
the current in (b) are the same event seen three ways.

Panel (c) takes the electrode potential from `m` rather than from `result.voltages`, which
carries a different quantity in each compensation mode. Pass the model the run was built
from — for an uncompensated run that is the uncompensated model, not the odr one.

`conc_row_height` is a fixed pixel height because panel (a) spans some fourteen decades and
needs the room; the other four share the remainder.
"""
function plot_cv_summary(
        result, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        include_capacitive = true,
        conc_limits = (1.0e-14, 1.0e2),
        conc_row_height = 350,
        labelsize = 22,
        fig_size = (800, 1260),
        panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)"],
    )
    fig = Figure(size = fig_size)

    ax1, _ = panel_conc_time!(fig, fig[1, 2], result, m; xlabel = "")
    ax2 = panel_time_current!(fig, fig[2, 2], result, m;
                              species, n_e, sgn, include_capacitive, xlabel = "")
    ax3 = panel_time_voltage!(fig, fig[3, 2], result; m, xlabel = "")
    ax4 = panel_time_ph!(fig, fig[4, 2], result; xlabel = "")
    ax5, _ = QoverK(fig, fig[5, 2], result, m; legend_pos = fig[5, 3])

    for (i, lab) in enumerate(panel_labels)
        Label(fig[i, 1], lab; fontsize = FS_LABEL, font = :bold,
              valign = :top, padding = (0, -20, -10, 0))
    end

    for ax in (ax1, ax2, ax3, ax4, ax5)
        ax.ylabelsize = labelsize
        ax.xlabelsize = labelsize
    end
    for ax in (ax1, ax2, ax3, ax4)
        hidexdecorations!(ax, grid = false)
    end
    ax3.yticks = LinearTicks(3)
    ax4.yticks = LinearTicks(4)
    ylims!(ax1, conc_limits...)

    linkxaxes!(ax1, ax2, ax3, ax4, ax5)
    rowgap!(fig.layout, 15)
    rowsize!(fig.layout, 1, conc_row_height)
    for r in 2:5
        rowsize!(fig.layout, r, Relative(1 / 6))
    end
    return fig
end

"""
    plot_7species_contours(result, X, m; layout = (2, 4), colorrange = nothing, ...)

log₁₀ concentration of every transported species over distance and time, one panel each.

Laid out as `layout = (rows, cols)` with the leftover cell taken by a **shared** colorbar.
Seven species in a 2×4 grid leaves exactly one free slot, which is why that is the default:
a single column of seven panels is far too tall for a page, and splitting into two figures
loses the side-by-side comparison.

`show_potential` fills that free cell with the electrode potential against the same time
axis, which also pushes the colorbar into a column of its own where it can run full height.
Keep it on: seven fields over time with nothing saying where in the sweep a given time
falls makes every feature a bare number rather than "the cathodic vertex". Turn it off and
the colorbar returns to the free cell at half height.

Only the bottom panel of each column carries tick labels, and the time label is written
once under the whole grid rather than per panel.

!!! note "Why this figure is decimated"
    A log-scaled y axis costs CairoMakie its fast path for `heatmap!`: cells of unequal
    height cannot be blitted as one image, so it emits a polygon per cell. Seven panels of
    a run stored at `Δt = 0.05` are millions of vector paths and take minutes to draw.
    `max_time_samples` caps the columns per panel at roughly one per pixel, which is all
    the figure can show; raise it only if a feature narrower than that is being missed, and
    use `x_stride` if the spatial grid is the expensive direction instead.

The colour scale is shared across panels, so a colour means the same concentration
everywhere in the figure and species can be read against one another. Pass `colorrange` to
pin it, e.g. across several figures. Concentrations below `10^floor_exp` are clamped, which
also absorbs the zeros and any negative excursion.
"""
function plot_7species_contours(
        result, X, m;
        scale = mol / dm^3,
        num_levels = 24,
        layout = (2, 4),
        colorrange = nothing,
        floor_exp = -12,
        fig_size = (1250, 520),
        show_potential = true,
        colgap_px = 20,
        rowgap_px = 8,
        max_time_samples = 400,
        x_stride = 1,
    )
    bulk = m.bulk
    target_species = [
        ("K⁺", rich("K", superscript("+"))),
        ("H⁺", rich("H", superscript("+"))),
        ("CO₂", rich("CO", subscript("2"))),
        ("OH⁻", rich("OH", superscript("−"))),
        ("HCO₃⁻", rich("HCO", subscript("3"), superscript("−"))),
        ("CO₃²⁻", rich("CO", subscript("3"), superscript("2−"))),
        ("CO", rich("CO")),
    ]
    times_all = result.tsol.t
    discrete_cmap = cgrad(:jet, num_levels, categorical = true)

    # Decimate before plotting. A log-scaled y axis costs CairoMakie its fast path for
    # `heatmap!`: with unequal cell heights it cannot blit one image and instead emits one
    # polygon per cell. At Δt = 0.05 over 80 s that is 1600 columns per panel, seven
    # panels deep — millions of vector paths for a figure whose panels are ~280 px wide.
    # Nothing below ~1 sample per pixel is visible, so the samples beyond that only cost
    # render time and file size.
    tstride = max(1, cld(length(times_all), max_time_samples))
    tidx = 1:tstride:length(times_all)
    xidx = 1:x_stride:length(X)
    times = times_all[tidx]
    Xp = X[xidx]

    # Build every log field first: a shared colour range cannot be known until all of
    # them exist.
    panels = Any[]
    for (sp_name, sp_label) in target_species
        sp_idx = findfirst(s -> s.name == sp_name, bulk)
        if isnothing(sp_idx)
            @warn "Species '$sp_name' not found in bulk; skipping this panel."
            continue
        end
        @views c_matrix = result.tsol[sp_idx, 1:length(X), 1:length(times_all)] ./ scale
        M = log10.(max.(c_matrix, 10.0^floor_exp))
        replace!(M, Inf => float(floor_exp), -Inf => float(floor_exp), NaN => float(floor_exp))
        # Clamp first, subsample second: the floor has to see every sample, or a spike
        # below it that happens to fall on a dropped column would come back as a hole.
        push!(panels, (sp_label, M[xidx, tidx]))
    end
    isempty(panels) && error("none of the target species were found in `m.bulk`")

    crange = if colorrange === nothing
        lo = minimum(minimum(p[2]) for p in panels)
        hi = maximum(maximum(p[2]) for p in panels)
        lo == hi ? (lo - 0.5, hi + 0.5) : (lo, hi)
    else
        colorrange
    end

    nrow, ncol = layout
    npanel = length(panels)
    # The potential trace claims the first free cell, so it decides which panels have
    # something below them. Without counting it, the last panel of the top row is treated
    # as the bottom of its column and grows a tick row and an axis label, which is what
    # opens the band of white between the two rows.
    has_potential = show_potential && npanel < nrow * ncol
    nfilled = npanel + (has_potential ? 1 : 0)

    fig = Figure(size = fig_size)
    axs = Axis[]
    local hm

    for (k, (sp_label, M)) in enumerate(panels)
        row, col = fldmod1(k, ncol)
        # bottom of its column: nothing sits in the cell below it
        is_bottom = (k + ncol > nfilled)

        ax = Axis(
            fig[row, col];
            # No per-panel `xlabel`: every bottom panel would repeat it, and they all share
            # one linked time axis. A single `Label` under the grid says it once.
            ylabel = col == 1 ? rich(rich("x", font = :bold_italic), "  (m)") : "",
            yscale = log10,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(9),
            yticks = (
                10.0 .^ (floor_exp:3:(floor_exp + 6)),
                [powlab(floor_exp), powlab(floor_exp + 3), powlab(floor_exp + 6)],
            ),
            title = sp_label,
        )
        is_bottom || hidexdecorations!(ax; grid = false)
        col == 1 || hideydecorations!(ax; grid = false, ticks = false)

        hm = heatmap!(
            ax, times, Xp .+ 1.0e-12, M';
            colorrange = crange,
            colormap = discrete_cmap,
            interpolate = false,
        )
        push!(axs, ax)
    end

    length(axs) > 1 && linkaxes!(axs...)

    # The free cell takes the potential protocol. Every contour panel is a field over time,
    # but none of them says where in the sweep a given time falls — with this panel present,
    # a feature at t = 22 s is read as "the cathodic vertex" instead of as a bare number.
    if has_potential
        prow, pcol = fldmod1(npanel + 1, ncol)
        ax_u = Axis(
            fig[prow, pcol];
            ylabel = rich("Electrode Potential ", rich("U", font = :bold_italic),
                          subscript("M"), "\n(V vs. SHE)"),
            ylabelsize = 18,
            yaxisposition = :right,
        )
        lines!(ax_u, result.times, electrode_potential(result, m);
               color = parse(Colorant, "#D7C2F0"), linewidth = LW_LINE)
        # x only: the contour panels share a log-scaled distance axis, and `linkaxes!`
        # would drag this panel's volts onto it.
        isempty(axs) || linkxaxes!(axs[1], ax_u)
    end

    # Shared colorbar. It may size its own column, but not one it shares with a panel: a
    # colorbar reports its 18 px as the column's width requirement, which would collapse
    # the panel above it to a sliver.
    cbar_label = rich(
        "log", subscript("10"), "(", rich("c", font = :bold_italic), " / M)"
    )
    if nfilled < nrow * ncol
        crow, ccol = fldmod1(nfilled + 1, ncol)
        Colorbar(
            fig[crow, ccol], hm;
            label = cbar_label, vertical = true, width = 18,
            tellwidth = false, halign = :left,
        )
    else
        Colorbar(
            fig[1:nrow, ncol + 1], hm;
            label = cbar_label, vertical = true, width = 18,
        )
    end

    # One time label for the whole grid, in the row below it.
    Label(
        fig[nrow + 1, 1:ncol], lab_time;
        fontsize = FS_LABEL, font = :regular, padding = (0, 0, 0, 4),
    )

    colgap!(fig.layout, colgap_px)
    rowgap!(fig.layout, rowgap_px)
    return fig
end

# ── moved from cell 358b1fba (panel_log_contour!)  [+ electrolyte kwarg] ──
function panel_log_contour!(
        panel_pos, cbar_pos, result, X, times, sp, m;
        scale = mol / dm^3, num_levels = 24
    )
    electrolyte = m.elydata
    c_bulk = electrolyte.c_bulk[sp.idx]
    log_c_bulk = log10(c_bulk / scale)

    M = [
        log_c_bulk - log10(max(result.tsol[sp.idx, ix, it] / scale, 1.0e-12))
            for ix in 1:length(X), it in 1:length(times)
    ]
    c_min = 0.0
    c_max = log_c_bulk - log10(1.0e-5)

    discrete_cmap = cgrad(:jet, num_levels, categorical = true)
    ax = Axis(
        panel_pos;
        xlabel = lab_time,
        ylabel = rich(rich("x", font = :bold_italic), "  (m)"),
        yscale = log10,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (
            10.0 .^ (-12:3:-6),
            [powlab(-12), powlab(-9), powlab(-6)],
        )
    )

    hm = heatmap!(
        ax, times, X .+ 1.0e-12, M';
        colorrange = (c_min, c_max),
        colormap = discrete_cmap,
        interpolate = false
    )

    Colorbar(
        cbar_pos, hm;
        label = rich(
            "log", subscript("10"), "(", rich("c", font = :bold_italic),
            subscript("bulk"), ") − log", subscript("10"), "(",
            rich("c", font = :bold_italic), subscript(sp.label), ")"
        ),
        ticklabelsize = 20, labelsize = 20
    )

    return ax, hm
end

# ── moved from cell 34857db0 (QoverK)  [+ electrolyte kwarg] ──
function QoverK(
        fig, panel_pos, result, m; scale = mol / dm^3, lw = LW_LINE,
        xlabel = lab_time, legend_pos = nothing
    )
    bulk = m.bulk
    electrolyte = m.elydata

    times = result.tsol.t
    nt = length(times)

    # Equilibrium Constant
    EqK_hco3 = (electrolyte.c_bulk[ihco3] * electrolyte.c_bulk[iohminus]) /
        electrolyte.c_bulk[ico3]
    EqK_co3 = (electrolyte.c_bulk[ico2] * electrolyte.c_bulk[iohminus]) /
        electrolyte.c_bulk[ihco3]

    # reaction quotient Q(t)
    cco3 = [result.tsol[ico3, 1, t] for t in 1:nt]
    cohm = [result.tsol[iohminus, 1, t] for t in 1:nt]
    chco3 = [result.tsol[ihco3, 1, t] for t in 1:nt]
    cco2 = [result.tsol[ico2, 1, t] for t in 1:nt]

    Qt_hco3 = (chco3 .* cohm) ./ cco3
    Qt_co3 = (cco2 .* cohm) ./ chco3

    ax = Axis(
        panel_pos;
        xlabel = xlabel,
        ylabel = rich(
            "Reaction Quotient\n",
            rich("Q", font = :bold_italic), " / ", rich("K", font = :bold_italic)
        ),
        yscale = log10,
        yticks = (
            10.0 .^ (-6:3:6),
            [powlab(-6), powlab(-3), powlab(0), powlab(3), powlab(6)],
        ),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(27)
    )
    hlines!(ax, [1.0e0]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)

    l1 = lines!(
        ax, times, max.(Qt_hco3 ./ EqK_hco3, eps(Float64));
        color = "#2980B9", linewidth = lw
    )
    l2 = lines!(
        ax, times, max.(Qt_co3 ./ EqK_co3, eps(Float64));
        color = "#E67E22", linewidth = lw
    )

    ylims!(ax, low = 1.0e-7, high = 1.0e7)

    leg = Legend(
        legend_pos,
        [l1, l2],
        [
            rich("HCO", subscript("3"), superscript("-"), " ⇌ CO", subscript("3"), superscript("2-")),
            rich("CO", subscript("2"), " ⇌ HCO", subscript("3"), superscript("-")),
        ],
        framevisible = false
    )

    return ax, leg
end

# ── moved from cell f506e83f (plot_cv_total_current)  [ely = elydata_Gold_odr → ely = model] ──
function plot_cv_total_current(
        result, m;
        species = ico,
        n_e = 2,
        sgn = 1,
        redox_species = nothing,   # back-compat: Dict(species => n_e), first entry wins
        include_capacitive::Bool = true,
        scale = cm^2 / mA,
        color_gradient::Bool = true
    )
    sp, ne = redox_species === nothing ? (species, n_e) : first(pairs(redox_species))

    I_F = faradaic_current(result; species = sp, n_e = ne, sgn = sgn)
    I_C = include_capacitive ? capacitive_current(result) : zero(I_F)
    n_t = length(I_F)

    # ---- Total current with unit scaling ----
    I_total = (I_F .+ I_C) .* scale

    fig, ax = with_theme(electrochemistry_theme()) do
        f = Figure(size = (540, 440))
        a = Axis(
            f[1, 1],
            ylabel = lab_current,
            xlabel = lab_voltage
        )
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
        scale = cm^2 / mA,
        sign = 1,
        sim_linewidth::Real = LW_LINE
    )
    electrolyte = m.elydata
    wanted = ["0.1", "0.5", "1.0"]
    pressures = ["Ar sat", "0.1", "0.2", "0.3", "0.5", "0.6", "1.0"]
    raw_cv = CSV.read("../data/Langmuir_CV_data/Figure_3.csv", DataFrame; header = false)
    sub = Matrix(raw_cv[4:end, :])
    num_cv = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num_cv, :auto)
    npairs = size(num_df, 2) ÷ 2
    keep = findall(in(wanted), pressures[1:npairs])
    n_exp = length(keep)
    n_sim = length(P_recs)
    pastel2 = cgrad([colorant"#F2728A", colorant"#5BA8E8"])
    cols_exp = [pastel2[t] for t in range(0, 1, length = max(n_exp, 1))]
    pastel3 = cgrad([colorant"#FFB3BA", colorant"#A3D8FF"])
    cols_sim = [pastel3[t] for t in range(0, 1, length = max(n_sim, 1))]
    fig = with_theme(electrochemistry_theme()) do
        return Figure(size = (800, 800))
    end
    ax_exp = Axis(
        fig[1, 1],
        ylabel = rich(rich("I", font = :bold_italic), "  (mA cm", superscript("−2"), ")"),
        xticklabelsvisible = false, xticksvisible = false
    )
    ax_sim = Axis(
        fig[2, 1],
        xlabel = rich(rich("ϕ", font = :bold_italic), "  (V vs. SHE)"),
        ylabel = rich(rich("I", font = :bold_italic), "  (mA cm", superscript("−2"), ")")
    )
    linkxaxes!(ax_exp, ax_sim)
    ax_sim.xticks = -1.5:0.3:1.0
    ax_exp.limits = (nothing, (-5.5, 1.8))
    ax_exp.yticks = 1:-2:-5
    # ---- exp ----
    for (k, j) in enumerate(keep)
        xcol, ycol = 2j - 1, 2j
        label_text = "$(pressures[j]) atm"
        lines!(
            ax_exp, num_df[!, xcol], num_df[!, ycol];
            color = cols_exp[k], linewidth = sim_linewidth
        )
        text!(
            ax_exp, label_text;
            position = (-0.79 - 0.03 * k, 1.6 - 1.3 * k),
            color = cols_exp[k], fontsize = 26, font = :bold
        )
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
            if isa(ely.ircompensation , OhmicDropEstimation)
                icc = ely.icc
                node_we = 1
                I_C = [u[icc, node_we] for u in rec.tsol[1:(end - 1)]]
            end
        end
        I = sign .* (I_F .+ I_C) .* scale .* 2
        label_text = "$(p) atm"
        lines!(
            ax_sim, rec.voltages, I;
            color = cols_sim[j], linewidth = sim_linewidth
        )
        text!(
            ax_sim, label_text;
            position = (-0.99 - 0.05 * j, 0.7 - 0.7 * j),
            color = cols_sim[j], fontsize = 26, font = :bold
        )
    end
    text!(
        ax_sim, "Theory";
        position = (-0.5, -0.8), color = :gray30,
        fontsize = 32, font = :bold
    )
    text!(
        ax_exp, "Experiment";
        position = (-0.5, -1.27), color = :gray40,
        fontsize = 32, font = :bold
    )
    Label(
        fig[1, 1, TopLeft()], "(a)";
        fontsize = 28, font = :bold, padding = (0, 5, 20, 0)
    )
    Label(
        fig[2, 1, TopLeft()], "(b)";
        fontsize = 28, font = :bold, padding = (0, 5, 20, 0)
    )
    rowgap!(fig.layout, 1, 15)
    return fig
end


function plottsol(
        grd, celldata, tsol; xsplit = 5ufac"nm",
        species = celldata.iϕ,
        limits = nothing,
        figscale = 1,
        stride = 1,
        levels = 5
    )
    label = "C/(mol/dm^3)"

    X = grd[XCoordinates]
    wl = 100
    hy = 250
    wr = log10(X[end] * figscale / ufac"nm") * 40
    figsize = (wl + wr + 200, hy)
    fig = Figure(size = figsize, fontsize = 5)
    axl = Axis(
        fig[1, 1],
        xlabel = L"x/nm",
        ylabel = L"t/s",
    )
    xlims!(axl, 0, 10)
    axr = Axis(
        fig[1, 2],
        yticksvisible = false,
        yticklabelsvisible = false,
        xscale = log10,
        xlabel = L"x/nm",
    )


    xlims!(axr, 10, X[end] / ufac"nm")
    #	hidedecorations!(axl)
    #	hidedecorations!(axr)
    colorscale = identity
    colormap = :terrain
    scale = 1 / ufac"mol/dm^3"
    if species == celldata.iϕ
        colorscale = identity
        colormap = :seismic
        scale = 1
    end
    T = tsol.t
    nx = length(X)
    nt = length(T)
    u = [tsol[species, ix, it] for ix in 1:stride:nx, it in 1:stride:nt ] * scale
    x = [X[ix] / ufac"nm" for ix in 1:stride:nx, it in 1:stride:nt ]
    y = [T[it] for ix in 1:stride:nx, it in 1:stride:nt ]
    if isnothing(limits)
        limits = extrema(u)
    end
    lvs = range(limits..., length = levels)

    contourf!(axl, x, y, u; levels = lvs, colormap, colorscale)
    cr = contourf!(axr, x, y, u; levels = lvs, colormap, colorscale)
    Colorbar(fig[1, 3], cr; label)

    colgap!(fig.layout, 2)
    colsize!(fig.layout, 1, Fixed(wl))
    colsize!(fig.layout, 2, Fixed(wr))
    colsize!(fig.layout, 3, Fixed(40))
    rowsize!(fig.layout, 1, Fixed(hy))

    return fig
end
