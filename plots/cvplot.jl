# plots/cvplot.jl

using CairoMakie
using CSV
using DataFrames

"""
    plot_koper_fig3(; csv_relpath="data/Langmuir_CV_data/Figure_3.csv", fig_size=(1600, 900))

Load a Koper (Langmuir) Figure 3–style CSV file and plot CV curves.
- `csv_relpath` is a path relative to the project root (assumed to be the parent of this file's directory).
- Returns: `(fig, ax)`.
"""
function plot_koper_fig3(; csv_relpath="data/Langmuir_CV_data/Figure_3.csv", fig_size=(1600, 900))
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

        line = lines!(ax, num_df[!, xcol], num_df[!, ycol]; color = cols[j])
        push!(plot_objs, line)
    end

    # Place the legend to the right.
    Legend(fig[1, 2], plot_objs, labels, "Experimental"; framevisible=true)

    return fig, ax
end
