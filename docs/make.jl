push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "examples"))

using Documenter
using Capacitance_Code
using ExampleJuggler
using CairoMakie
using LiquidElectrolytes
using PlutoStaticHTML
using VoronoiFVM

ExampleJuggler.verbose!(true)
const thisdir = @__DIR__

function make(; with_notebooks = true)
    pages = [
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Code Structure" => "code_structure.md",
        "Model" => "model.md",
        "Usage" => "usage.md",
        "API" => "api.md",
    ]

    cleanexamples()

    exampledir  = joinpath(@__DIR__, "..", "examples")
    notebookdir = joinpath(@__DIR__, "..", "notebooks")

    size_threshold_ignore = []

    if with_notebooks
        notebooks = [
            "EquilibriumCheck.jl",
            "plot.jl",
        ]

        notebook_examples = @docplutonotebooks(
            notebookdir,
            notebooks,
            iframe = false,
            append_build_context = false,
        )

        size_threshold_ignore = last.(notebook_examples)
        push!(pages, "Notebooks" => notebook_examples)
    end

    DocMeta.setdocmeta!(
        Capacitance_Code,
        :DocTestSetup,
        :(using Capacitance_Code, LiquidElectrolytes, Unitful, LessUnitful),
        recursive = true,
    )

    cd(thisdir)

    makedocs(
        sitename = "Capacitance_Code.jl",
        modules = [Capacitance_Code],
        format = Documenter.HTML(
            mathengine = MathJax3(),
            size_threshold_ignore = size_threshold_ignore,
        ),
        clean = false,
        doctest = true,
        warnonly = true,
        draft = false,
        authors = "S.M. Choi",
        repo = "https://github.com/ElCatFVM/Capacitance_Code.jl",
        pages = pages,
    )

    cleanexamples()

    if !isinteractive()
        deploydocs(
            repo = "github.com/ElCatFVM/Capacitance_Code.jl.git",
            devbranch = "main",
        )
    end

    return nothing
end

if isinteractive() || (haskey(ENV, "DOCSONLY") && ENV["DOCSONLY"] == "true")
    make(; with_notebooks = false)
else
    make()
end