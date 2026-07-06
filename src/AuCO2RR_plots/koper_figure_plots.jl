"""
    plot_koper_fig1(; csv_relpath="data/Langmuir_CV_data/Figure_1.csv",
                     fig_size=(650, 400),
                     pairs=[(1,2), (3,4), (5,6)],
                     colors=(:pink, :skyblue, :lightgreen),
                     legend_pos=:rt,
                     log_y=false)

Load a Koper (Langmuir) Figure 1–style CSV and plot facet CV curves.
If `log_y=true`, plot `log10(abs(y))`. Returns `fig`.
"""
function plot_koper_fig1(; csv_relpath="data/Langmuir_CV_data/Figure_1.csv",
                         fig_size=(650, 400),
                         pairs=[(1,2), (3,4), (5,6)],
                         colors=(:pink, :skyblue, :lightgreen),
                         legend_pos=:rt,
                         log_y=false)

    root = normpath(joinpath(@__DIR__, ".."))
    csv_path = joinpath(root, splitdir(csv_relpath)...)

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
               ylabel = "Current Density (mA/cm²)",
               xlabel = "Voltage (ϕ-ϕₚ)")

    raw_df = CSV.read(csv_path, DataFrame; header=false)

    facet_row = collect(raw_df[1, :])

    numeric_rows = [
        parse.(Float64, coalesce.(collect(raw_df[i, :]), "NaN"))
        for i in 3:nrow(raw_df)
    ]
    num_df = DataFrame(hcat(numeric_rows...)', names(raw_df))

    nplot = min(length(pairs), length(colors))

    for i in 1:nplot
        xcol, ycol = pairs[i]
        y = num_df[!, ycol]
        ydata = log_y ? log10.(abs.(y)) : y
        lines!(ax, num_df[!, xcol], ydata; color=colors[i], label=string(facet_row[xcol]))
    end

    axislegend(ax; position=legend_pos)
    return fig
end

"""
    plot_koper_fig3(; csv_relpath="data/Langmuir_CV_data/Figure_3.csv", fig_size=(1600, 900))

Load a Koper (Langmuir) Figure 3–style CSV file and plot CV curves.
- `csv_relpath` is a path relative to the project root (assumed to be the parent of this file's directory).
- Returns: `(fig, ax)`.
"""
function plot_koper_fig3(; csv_relpath="data/Langmuir_CV_data/Figure_3.csv", fig_size=(1600, 900), log_y = false)
    # Assume the project root is the parent directory of this file's folder (i.e., <root>/plots/cvplot.jl).
    root = normpath(joinpath(@__DIR__, ".."))

    # Build an absolute path from the project root and the provided relative path.
    csv_path = joinpath(root, splitdir(csv_relpath)...)

    # Create figure and axis (no computation; only plotting loaded data).
    fig = Figure(size = fig_size)
    ax  = Axis(
        fig[1, 1],
        ylabel = L"I (mA/cm^2)",
        xlabel = L"\phi (V \; \mathrm{vs}\; SHE)",
    )

    # Read CSV as raw table (no header in the original file format).
    raw = CSV.read(csv_path, DataFrame; header=false)

    # The first row contains labels/pressure metadata.
    pres = vec(Matrix(raw[1:1, :]))

    # Data starts at row 4 in the original convention.
    sub = Matrix(raw[4:end, :])

    # Convert each cell to Float64; treat missing values as NaN.
    num = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)

    # Interpret columns as (x, y) pairs: (1,2), (3,4), ...
    npairs = size(num_df, 2) ÷ 2

    # Construct a simple color gradient from pink to blue.
    pink  = RGB(1.0, 0.7, 0.8)
    pblue = RGB(0.2, 0.5, 1.0)
    cols  = [RGB(pink.r + t*(pblue.r-pink.r),
                 pink.g + t*(pblue.g-pink.g),
                 pink.b + t*(pblue.b-pink.b)) for t in range(0, 1, length=npairs)]

    plot_objs = Any[]
    labels    = String[]

    for j in 1:npairs
        xcol, ycol = 2j - 1, 2j

        # Keep the original label logic from your notebook code.
        label = (j == 1) ? "$(pres[1])\t\t sat" : "$(pres[2j])\t pCO2(atm)"
        push!(labels, label)
        if log_y
            ydata = log10.(abs.(num_df[!, ycol]))
        else
            ydata = num_df[!, ycol]
        end
        line = lines!(ax, num_df[!, xcol], ydata; color = cols[j])
        push!(plot_objs, line)
    end

    # Place the legend to the right.
    Legend(fig[1, 2], plot_objs, labels, "Experimental"; framevisible=true)

    return fig
end

"""
    plot_koper_fig5(; csv_relpath="data/Langmuir_CV_data/Figure_5.csv",
                     fig_size=(900, 600),
                     legend_title="Experimental",
                     log_y=true)

Load a Koper (Langmuir) Figure 5–style CSV file and plot CV curves.
If `log_y=true`, plot `log10(abs(y))`.

Returns: `(fig, ax)`.
"""
function plot_koper_fig5(; csv_relpath="data/Langmuir_CV_data/Figure_5.csv",
                         fig_size=(900, 600),
                         legend_title="Experimental",
                         log_y=true)

    root = normpath(joinpath(@__DIR__, ".."))
    csv_path = joinpath(root, splitdir(csv_relpath)...)

    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1],
        xlabel = L"\phi \, (\mathrm{V \; vs \; SHE})",
        ylabel = L"I \; (\mathrm{mA/cm^2})",
    )

    raw = CSV.read(csv_path, DataFrame; header=false)

    pres_labels = vec(Matrix(raw[1:1, :]))
    sub  = Matrix(raw[4:end, :])
    num  = map(x -> x === missing ? NaN : parse(Float64, x), sub)
    num_df = DataFrame(num, :auto)

    npairs = size(num_df, 2) ÷ 2
    colors = [RGB(0.9 - (0.05 * i/npairs), 0.15 * (1 - i/npairs), 0.3 * i/npairs) for i in 1:npairs]

    plots  = Any[]
    labels = String[]

    for j in 1:npairs
        xcol, ycol = 2j - 1, 2j

        label_idx = min(2j, length(pres_labels))
        label = "$(pres_labels[label_idx]) pH"
        push!(labels, label)

        y = num_df[!, ycol]
        ydata = log_y ? log10.(abs.(y)) : y

        h = lines!(ax, num_df[!, xcol], ydata; color = colors[j])
        push!(plots, h)
    end

    Legend(fig[1, 2], plots, labels, legend_title; framevisible=true)
    return fig
end

