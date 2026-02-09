using CairoMakie
using VoronoiFVM
using LiquidElectrolytes

"""
plot_caps_comparison(; is_Landstorfer, model_key, result_pb, result_pnp,
                      molarities, ϕ0_pzc=0.972,
                      ref_low=nothing, ref_high=nothing)

- Landstorfer mode: reference curves + sweep results
- Comparison mode: PB vs PNP via capsplot_v
- ref_low/ref_high: reference data as (voltages, dlcaps)
"""
function plot_caps_comparison(;
    is_Landstorfer::Bool = false,
    model_key::AbstractString,
    result_pb,
    result_pnp=nothing,
    molarities=nothing,
    ϕ0_pzc::Real=0.972,
    ref_low=nothing,
    ref_high=nothing,
)
    if is_Landstorfer
        f = Figure(size=(900, 900))
        ax = Axis(
            f[1, 1],
            xlabel = "φ / (V vs φ_pzc)",
            ylabel = "dlcaps / (μF / cm²)",
            title  = "CSV Plot",
            limits = ((-1.2, 1.2), (0, 120))
        )

        handles = Any[]
        labels  = String[]

        if ref_low !== nothing
            vref, cref = ref_low
            h = lines!(ax, vref .+ ϕ0_pzc, cref, linestyle=:dash)
            push!(handles, h)
            push!(labels,  model_key * "\t 5mM")
        end

        if ref_high !== nothing
            vref, cref = ref_high
            h = lines!(ax, vref .+ ϕ0_pzc, cref, linestyle=:dash)
            push!(handles, h)
            push!(labels,  model_key * "\t 100mM")
        end

        l = 1 / max(length(result_pb), 1)
        ks = Any[]

        for i in 1:length(result_pb)
            c = RGB(i*l, 0.0, 1 - i*l)

            v   = result_pb[i].voltage_range
            cdl = result_pb[i].dlcaps / (μF / cm^2)

            n = min(length(v), length(cdl))
            hi = lines!(ax, v[1:n], cdl[1:n], color=c, label="Result $i")
            push!(ks, hi)

            if molarities !== nothing
                m = molarities[i]
                push!(labels, "LiquidElectrolyte $m M")
            else
                push!(labels, "LiquidElectrolyte (case $i)")
            end
        end

        Legend(
            f[1, 1],
            [handles... , ks...],
            labels,
            halign = :left,
            valign = :top,
            tellheight = false,
            tellwidth  = false,
            framevisible = false
        )

        return f
    else
        results = [
            ("Poisson-Boltzmann",     result_pb),
            ("Poisson-Nernst-Planck", result_pnp),
        ]

        plots = Tuple{String,Any}[]
        for (name, res) in results
            try
                if res !== nothing && !isempty(res)
                    push!(plots, (name, res))
                end
            catch
                continue
            end
        end

        vis = GridVisualizer(Plotter=CairoMakie, legend=:lt, title="Compare PB & PNP")
        capsplot_v(vis, plots)

        return reveal(vis)
    end
end

function capsplot_v(vis, result_named)
    color = [:magenta, :blue]

    for (i, (name, res)) in enumerate(result_named)
        v = res[1].voltage_range
        c = res[1].dlcaps

        v = vec(v)               
        c = vec(c)               

        n = min(length(v), length(c))
        v = v[1:n]
        c = c[1:n] ./ (μF / cm^2)

        scalarplot!(
            vis, v, c;
            limits = (0, 400),    
            xlimits = (-1.1, 1.1),
            color = color[i],
            clear = (i == 1),
            label = name,
            markershape = :none,
            yscale = 10,
            xlabel = "φ / (V vs φ_pzc)",
            ylabel = "dlcaps / (μF / cm²)"
        )
    end

    return vis
end


function capsplot_κ(vis, sys; κ_values = [0, 1, 5, 10, 20, 30, 40])
    n = length(κ_values)
    colors = [RGB(0, 0, i/n) for i in 1:n]
    dls = LiquidElectrolytes.DLCapSweepResult[]

    sys_local = deepcopy(sys)
    data = electrolytedata(sys_local)

    κ_original = copy(data.κ)
    cbulk_original = copy(data.c_bulk)

    try
        for (j, κval) in enumerate(κ_values)
            data.κ .= κval
            data.c_bulk .= [0.05, 0.05] * ufac"mol/dm^3"

            result = caps(sys_local)
            push!(dls, result)

            scalarplot!(
                vis,
                result.voltages,
                result.dlcaps / (μF / cm^2);
                linestyle = :solid,
                color = colors[j],
                clear = false,
                label = "κ = $(κval)"
            )
        end
    catch e
        @warn "capsplot_κ failed" exception=e
        rethrow()
    finally
        data.κ .= κ_original
        data.c_bulk .= cbulk_original
    end

    return dls
end

function capsplot_κ(vis, sys; κ_values = [0, 1, 5, 10, 20, 30, 40], c_bulk = [0.05, 0.05])
    n = length(κ_values)
    colors = [RGB(0, 0, i/n) for i in 1:n]
    dls = LiquidElectrolytes.DLCapSweepResult[]

    sys_local = deepcopy(sys)
    data = electrolytedata(sys_local)

    κ_original = copy(data.κ)
    cbulk_original = copy(data.c_bulk)

    try
        for (j, κval) in enumerate(κ_values)
            data.κ .= κval
            data.c_bulk .= c_bulk * ufac"mol/dm^3"

            result = caps(sys_local)
            push!(dls, result)

            scalarplot!(
                vis,
                result.voltages,
                result.dlcaps / (μF / cm^2);
                linestyle = :solid,
                color = colors[j],
                clear = false,
                label = "κ = $(κval)"
            )
        end
    finally
        data.κ .= κ_original
        data.c_bulk .= cbulk_original
    end

    return dls
end

function plot_κ_sweep_with_refs(sys;
    ref_low,
    ref_high,
    κ_values = [0, 1, 5, 10, 20, 30, 40],
    c_bulk = [0.05, 0.05],
    size = (650, 650),
    xlimits = (-0.5, 0.5),
    xlabel = "φ / (V vs φ_pzc)",
    ylabel = "dlcaps / (μF / cm²)",
    legend = :lt,
)
    vis = GridVisualizer(Plotter=CairoMakie, legend=legend, size=size)

    scalarplot!(
        vis,
        ref_low.voltages,
        ref_low.dlcaps;
        color = :black,
        linestyle = :dot,
        label = "κ = 0",
        xlimits = xlimits,
        xlabel = xlabel,
        ylabel = ylabel,
        clear = true,
    )

    scalarplot!(
        vis,
        ref_high.voltages,
        ref_high.dlcaps;
        color = :blue,
        linestyle = :dot,
        label = "κ = 40",
        clear = false,
    )

    dls = capsplot_κ(vis, sys; κ_values=κ_values, c_bulk=c_bulk)

    reveal(vis)
    return (vis = vis, results = dls)
end

function capsplot(vis, ::Nothing, title)
    scalarplot!(
        vis, [NaN], [NaN];
        clear = false,
        markershape = :none,
        label = "(no data)",
        title = title,
        xlabel = L"φ / (V vs φ_{pzc})",
        ylabel = L"dlcaps / (μF / cm²)",
        xlimits = (-1.1, 1.1),
        ylimits = (-1, 250)
    )
    return vis
end

function capsplot_fixed(
    vis, result, title;
    is_Landstorfer::Bool = false, # molarity 
    nshow::Int = 201,
    xlimits_L=(-1.0, 1.0),
    ylimits_L=(0, 100),
    show_cdl0::Bool=true,
)
    nres = length(result)
    hmol = 1 / max(nres, 1)

    for i in 1:nres
        c = RGB(i * hmol, 0, 1 - i * hmol)

        v   = result[i].voltage_range
        cap = vec(result[i].dlcaps)

        n = min(nshow, length(v), length(cap))
        v   = v[1:n]
        cap = cap[1:n] / (μF / cm^2)

        # --- Legend(라벨)---
        lbl = if hasproperty(result[i], :molarity)
            "$(result[i].molarity) M"
        elseif hasproperty(result[i], :comb)
            "×$(result[i].comb)"
        else
            "Run $i"
        end
        # ----------------------------

        scalarplot!(
            vis, v, cap;
            clear=(i == 1),
            color=c,
            label=lbl, 
            title=title,
            markershape=:none,
            xlabel=L"\phi / (V vs \phi_{pzc})",
            ylabel=L"C_{dl} / (\mu F / cm^2)",
            xlimits=xlimits_L,
            limits=(ylimits_L[1], ylimits_L[2]),
        )

        # PZC 
        if show_cdl0
            scalarplot!(
                vis, [0.0], [result[i].cdl0] / (μF / cm^2);
                clear=false,
                color=c, 
                markershape=:circle,
                markersize=5,
                label=""
            )
        end
    end

    return vis
end


"""
overlay_csv!(vis, dfs; v=:Voltage, y=:Cdl, lab=nothing, ls=:dash, lw=2)

Overlay one or more CSV DataFrames onto an existing `vis`.
- `dfs` is a Vector of DataFrames (each from CSV.read).
- Default columns are `:Voltage` and `:Cdl`. Change with `v=` and `y=`.
- `lab` can be a Vector of labels (same length as dfs).
"""
function overlay_csv!(
    vis,
    dfs;
    ϕ_pzc::Real,
    v::Symbol = :Voltage,
    y::Symbol = :Cdl,
    lab = nothing,
    ls = :dashdot,
    lw = 2,
)
    n = length(dfs)
    labels = lab === nothing ? ["CSV $i" for i in 1:n] : lab
    length(labels) == n || error("lab length must match dfs length")

    for i in 1:n
        colors = [RGB(0.5, 0.3, i / n) for i in 1:n]
        df = dfs[i]
        scalarplot!(
            vis,
            Float64.(df[!, v]) .+ ϕ_pzc,
            Float64.(df[!, y]);
            color = colors[i],
            label = labels[i],
            linestyle = ls,
            linewidth = lw,
            markershape = :none,
        )
    end
    return vis
end



"""
capsplot_with_csv!(vis, result, title, dfs; ...)

Draw simulation (capsplot_fixed) and overlay multiple CSV DataFrames on top.
Returns `vis`.
"""
function capsplot_with_csv!(
    vis,
    result,
    title,
    ϕ_pzc,
    dfs;
    nshow::Int = 201,
    xlimits = (-1.0, 1.0),
    ylimits = (0, 100),
    show_cdl0::Bool = true,
    v::Symbol = :Voltage,
    y::Symbol = :Cdl,
    lab = nothing,
    ls = :dash,
    lw = 2,
)
    capsplot_fixed(
        vis, result, title;
        nshow = nshow,
        xlimits_L = xlimits,
        ylimits_L = ylimits,
        show_cdl0 = show_cdl0,
    )

    overlay_csv!(
        vis, dfs;
        ϕ_pzc = ϕ_pzc,
        v = v,
        y = y,
        lab = lab,
        ls = ls,
        lw = lw,
    )

    return vis
end



