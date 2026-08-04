using Documenter
using DocumenterMermaid          # mermaid diagrams inside markdown
using AuCO2RR
# using PlutoStaticHTML          # uncomment to render the Pluto notebooks into the docs

function mkdocs()
    DocMeta.setdocmeta!(
        AuCO2RR, :DocTestSetup, :(using AuCO2RR); recursive = true,
    )

    makedocs(
        sitename = "AuCO2RR.jl",
        modules = [AuCO2RR, AuCO2RR.GoldModel, AuCO2RR.AuCO2RR_plots],
        authors = "Sumin Choi",
        # This repository has no GitHub remote configured, so source links are disabled.
        # To enable them, drop `remotes` and set
        #     repo = Remotes.GitHub("USER", "Capacitance_Code")
        remotes = nothing,
        format = Documenter.HTML(
            size_threshold = nothing,          # allow large pages
            mathengine = MathJax3(
                Dict(
                    :tex => Dict(
                        "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                        "tags" => "ams",
                        # mhchem => \ce{CO2 + OH^- <=> HCO3^-}
                        "packages" => ["base", "ams", "autoload", "mhchem"],
                    ),
                ),
            ),
        ),
        clean = false,
        doctest = false,   # set true to actually run the docstring examples
        draft = false,     # set true to skip code blocks for a fast preview
        # Most functions in this package are still undocumented; without this the
        # build fails on the first missing docstring instead of producing a site.
        warnonly = [:missing_docs, :cross_references],
        pages = [
            "Home" => "index.md",
            "Guide" => "guide.md",
            "Public API" => "public.md",
        ],
    )
end

mkdocs()

# Skip deployment from the REPL; only CI publishes to gh-pages.
# Enable once a GitHub remote exists.
# if !isinteractive()
#     deploydocs(
#         repo = "github.com/USER/Capacitance_Code.git",
#         devbranch = "main",
#     )
# end
