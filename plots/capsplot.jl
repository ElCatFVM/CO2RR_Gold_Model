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
    result, title;
    is_Landstorfer::Bool = false,
    nshow::Int = 201,
    xlimits_L=(-1.0, 1.0),
    ylimits_L=(0, 100),
    show_cdl0::Bool=true,
)

    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1];
        xlabel = L"\phi~(\mathrm{V~vs~}\phi_{pzc})",
        ylabel = L"C_{dl}~(\mu \mathrm{F\,cm^{-2}})",
        title = title,
        titlesize = 25,
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 25,
        yticklabelsize = 25,
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 4.5,
        xtickwidth = 4.5,
        ytickwidth = 4.5,
        limits = (xlimits_L, ylimits_L)
    )

    nres = length(result)
    hmol = 1 / max(nres, 1)

    plot_objs = Any[]
    labels = String[]

    for i in 1:nres
        c = RGB(i * hmol, 0, 1 - i * hmol)

        v   = result[i].voltage_range
        cap = vec(result[i].dlcaps)

        n = min(nshow, length(v), length(cap))
        v   = v[1:n]
        cap = cap[1:n] / (μF / cm^2)

        lbl = if hasproperty(result[i], :molarity)
            "$(result[i].molarity) M"
        elseif hasproperty(result[i], :comb)
            "×$(result[i].comb)"
        else
            "Run $i"
        end

        line = lines!(ax, v, cap; color=c, linewidth=4)
        push!(plot_objs, line)
        push!(labels, lbl)

        if show_cdl0
            scatter!(ax, [0.0], [result[i].cdl0] / (μF / cm^2); color=c, markersize=10)
        end
    end

    leg = Legend(fig[1, 1], plot_objs, labels;
        framevisible = false,
        halign = :right,
        valign = :bottom,
        labelsize = 20,
        titlesize = 23,
        tellwidth = false,
        tellheight = false
    )
    translate!(leg.blockscene, -40, 40, 0)

    return fig, ax
end



"""
overlay_csv!(vis, dfs; v=:Voltage, y=:Cdl, lab=nothing, ls=:dash, lw=2)

Overlay one or more CSV DataFrames onto an existing `vis`.
- `dfs` is a Vector of DataFrames (each from CSV.read).
- Default columns are `:Voltage` and `:Cdl`. Change with `v=` and `y=`.
- `lab` can be a Vector of labels (same length as dfs).
"""
function overlay_csv_on_axis!(
    ax,
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

    colors = [RGB(0.5, 0.3, i / max(n, 1)) for i in 1:n]

    for i in 1:n
        df = dfs[i]

        # x축이 "V vs φ_pzc" 라면 보통 -가 맞습니다 (필요하면 +로 바꾸세요)
        x = Float64.(df[!, v]) .+ ϕ_pzc
        yv = Float64.(df[!, y])

        lines!(ax, x, yv; color=colors[i], linestyle=ls, linewidth=lw, label=labels[i])
    end

    return ax
end



"""
capsplot_with_csv!(vis, result, title, dfs; ...)

Draw simulation (capsplot_fixed) and overlay multiple CSV DataFrames on top.
Returns `vis`.
"""
function capsplot_with_csv(
    result, title,
    ϕ_pzc, dfs;
    nshow::Int = 201,
    xlimits = (-1.0, 1.0),
    ylimits = (0, 100),
    show_cdl0::Bool = true,
    v::Symbol = :Voltage,
    y::Symbol = :Cdl,
    lab = nothing,
    ls = :dashdot,
    lw = 2,
)
    fig, ax = capsplot_fixed(result, title;
        nshow=nshow, xlimits_L=xlimits, ylimits_L=ylimits, show_cdl0=show_cdl0
    )

    overlay_csv_on_axis!(ax, dfs; ϕ_pzc=ϕ_pzc, v=v, y=y, lab=lab, ls=ls, lw=lw)

    # CSV legend도 보고 싶으면 (기존 legend랑 겹치면 position만 조정)
    axislegend(ax; position=:lt, framevisible=false)

    return fig
end




