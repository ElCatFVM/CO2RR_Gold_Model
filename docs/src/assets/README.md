# Figure assets

PNGs referenced by the documentation pages. Documenter copies this directory into the built
site as-is, so a file dropped here is reachable from a page as `assets/<name>.png`.

[Cyclic voltammetry walkthrough](../cv_notebook.md) expects these, in the order it uses
them:

| File | Produced by |
|:--|:--|
| `cv_summary.png` | `plot_cv_summary(cv_odr, elystruct_odr)` |
| `conc_time_electrode.png` | `plot_conc_time_electrode(cv_odr, elystruct_odr)` |
| `contours_qk_odr.png` | `plot_7species_contours_qk(cv_odr, X, elystruct_odr)` |
| `contours_qk_unc_slow.png` | `plot_7species_contours_qk(cv_unc_low_sr, X, elystruct_unc)` |
| `ircomp_compare.png` | `plot_ircomp_compare(cv_odr, m; reference, reference_model)` |
| `summary_compare_ircomp.png` | `plot_cv_summary_compare(cv_odr, …, cv_unc, …)` |
| `ircomp_diff.png` | `panel_time_current_diff!(f, f[1, 1], cv_odr, cv_unc, m)` |
| `ircomp_factors.png` | `cvsweep_over_ircompfactor(...)` + `plot_ircomp_compare` |
| `summary_compare_scanrate.png` | `plot_cv_summary_compare(cv_exp, …, cv_odr, …)` |
| `scanrate_split.png` | `plot_scanrate_sweeps_split(SR_vec, scanrates)` |
| `pressure_split.png` | `pressure_varied_cvsweep_split(P_recs)` |
| `exp_sim_split.png` | `plot_exp_sim_cvsweep_split(P_recs)` |

Export from the notebook with `save(joinpath(@__DIR__, "..", "docs", "src", "assets",
"cv_summary.png"), fig; px_per_unit = 2)`. `px_per_unit = 2` keeps the text sharp on a
high-DPI display; the figure sizes in `struct.jl` are in points, not pixels.

The PNGs here are committed (`.gitignore` excludes `*.png` everywhere else), so the
published site shows them. Figures exported from the Pluto HTML export were rendered from its
SVG outputs at 2× scale.

A page whose image is missing still builds — it renders as a broken image rather than
failing the build, so pages can be written before their figures exist.
