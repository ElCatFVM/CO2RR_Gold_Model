using Documenter
using DocumenterMermaid          # mermaid diagrams inside markdown
using AuCO2RR
using Pkg
using TOML
# using PlutoStaticHTML          # uncomment to render the Pluto notebooks into the docs

"Version string of the package itself, read from the root Project.toml."
const PKG_VERSION = VersionNumber(
    TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"]
)

"Label shown in the sidebar version selector, e.g. `dev (0.3.1)`."
const DOC_VERSION_LABEL = "dev ($(PKG_VERSION))"

function mkdocs()
    DocMeta.setdocmeta!(
        AuCO2RR, :DocTestSetup, :(using AuCO2RR); recursive = true,
    )

    makedocs(
        sitename = "AuCO2RR.jl",
        modules = [AuCO2RR, AuCO2RR.GoldModel, AuCO2RR.AuCO2RR_plots],
        authors = "Sumin Choi",
        version = string(PKG_VERSION),
        # This repository has no GitHub remote configured, so source links are disabled.
        # To enable them, drop `remotes` and set
        #     repo = Remotes.GitHub("USER", "Capacitance_Code")
        remotes = nothing,
        format = Documenter.HTML(
            size_threshold = nothing,          # allow large pages
            # No git remote, so there is nothing to link "edit this page" or the
            # navbar icon at. Both must be nothing explicitly or Documenter guesses
            # a branch and warns on every build.
            edit_link = nothing,
            repolink = nothing,
            footer = "AuCO2RR.jl $(DOC_VERSION_LABEL) — built with " *
                "[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
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
            "AuCO2RR.jl" => "index.md",
            "Notations" => "notations.md",
            "Standard calculations" => "calculations.md",
            "Plots" => "plots.md",
            "API" => "api.md",
            "Internal API" => "internal.md",
        ],
    )
    return nothing
end

"""
    write_version_selector(builddir)

Populate the sidebar version selector for a **local** build.

`deploydocs` normally writes `versions.js` next to the version directories on gh-pages and
a `siteinfo.js` inside each one; the selector at the bottom left of the sidebar is filled
in from those two files by `documenter.js`. A local build has neither, so the selector
stays empty. Writing them by hand gives the same result without deploying.

Delete this call once `deploydocs` is enabled — it would otherwise overwrite what the real
deployment produced.
"""
#function write_version_selector(builddir)
#    label = DOC_VERSION_LABEL
#    write(
#        joinpath(builddir, "siteinfo.js"),
#        """
#        var DOCUMENTER_CURRENT_VERSION = "$(label)";
#        """
#    )
#    write(
#        joinpath(builddir, "versions.js"),
#        """
#        var DOC_VERSIONS = ["$(label)"];
#        var DOCUMENTER_NEWEST = "$(label)";
#        var DOCUMENTER_STABLE = "$(label)";
#        """
#    )
#    return nothing
#end

mkdocs()
write_version_selector(joinpath(@__DIR__, "build"))

# Skip deployment from the REPL; only CI publishes to gh-pages.
# Enable once a GitHub remote exists — and drop `write_version_selector` above,
# since deploydocs writes versions.js / siteinfo.js itself.
 if !isinteractive()
     deploydocs(
         repo = "github.com/USER/Capacitance_Code.git",
         devbranch = "main",
         versions = ["stable" => "v^", "dev" => "dev", "v#.#"],
     )
 end
