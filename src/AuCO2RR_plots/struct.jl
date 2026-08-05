# Central plotting style for every AuCO2RR figure.
#
# `electrochemistry_theme()` returns the shared Makie Theme. It is also applied
# globally in `__init__` (via `set_theme!`), so figures that do NOT wrap
# themselves in `with_theme(...)` still inherit the same style and font. This
# keeps all plots consistent without repeating styling code in each function
# (efficiency), and lets us drop per-function `with_theme` wrappers over time.
#
# Typography follows the row_interaction_script convention: a single sans-serif
# family (DejaVu Sans, the script's `sans_font`) for all text, with `rich(...)`
# labels that render physical quantities in italic and units in roman. Build
# labels with the `lab_*` / `powlab` helpers below to stay consistent.
#
# To reset to Makie defaults in a session: `CairoMakie.set_theme!()`.

"Sans-serif family shared by every figure (matches row_interaction_script's `sans_font`)."
const PLOT_FONT_FAMILY = "DejaVu Sans"

# --------------------------------------------------------------------------
# Line widths.
#
# These used to be open-coded at every plotting function, which left seven
# different values (1, 2, 3, 3.5, 4, 5, 5.5) in circulation under six different
# keyword names — two panels of the *same* multi-panel figure could disagree.
# Default every width keyword to one of these instead of to a bare number, and
# reach for a new constant rather than a literal when a new role appears.
# --------------------------------------------------------------------------

"Width of a data series. Also the `Lines` default of [`electrochemistry_theme`](@ref)."
const LW_LINE = 4

"Width of an emphasised series in an overlay; must stay visibly above `LW_LINE`."
const LW_HIGHLIGHT = 6

"Width of experimental reference data drawn behind the simulation."
const LW_EXP = 2

"Width of dashed guides — `vlines!`, `hlines!`, boundary-layer markers."
const LW_GUIDE = 2

# --------------------------------------------------------------------------
# Colour sequences for a swept family of curves (scan rate, L, pressure, …).
#
# Index them as `CMAP[t]` with `t in range(0, 1, length = n)` so a family always
# spans the full sequence regardless of how many members it has.
# --------------------------------------------------------------------------

"""
Teal → orange sequence for scan-rate families.

Equivalent to the `RGB(0.2 + 0.6t, 0.8 - 0.5t, 0.8 - 0.7t)` ramp that
[`plot_scanrate_sweeps`](@ref) used to build inline — that expression is linear in `t`,
so a two-stop gradient reproduces it exactly.
"""
const CMAP_SCANRATE = cgrad([RGB(0.2, 0.8, 0.8), RGB(0.8, 0.3, 0.1)])

"Muted red → blue sequence for CO₂ partial-pressure families."
const CMAP_PRESSURE = cgrad([colorant"#F2A6A6", colorant"#A9C6EE"])

function electrochemistry_theme()
    Theme(
        # Global font family: every textual element (ticks, titles, legends and
        # `rich(...)` labels) renders in DejaVu Sans. `:regular/:bold/:italic`
        # used inside `rich(...)` resolve to the matching variant here.
        fonts = (
            regular     = PLOT_FONT_FAMILY,
            bold        = "DejaVu Sans Bold",
            italic      = "DejaVu Sans Oblique",
            bold_italic = "DejaVu Sans Bold Oblique",
        ),
        size = (960, 540),
        Axis = (
            spinewidth = 5.5,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            xticksize = 8,
            yticksize = 8,
            xlabelsize = 25,
            ylabelsize = 25,
            xticklabelsize = 25,
            yticklabelsize = 25,
            xgridvisible = false,
            ygridvisible = false,
            xlabelpadding = 10,
            ylabelpadding = 10,
            xlabelfont = :bold,
            ylabelfont = :bold,
            # NOTE: do NOT force LinearTicks here. As a global theme this would
            # also hit log-scale axes, crowding linear tick values (e.g. 2000,
            # 4000, … 10000) at the top of a log axis. Let each axis pick its
            # own locator (log axes get log ticks); set ticks per-figure when a
            # specific count is needed.
        ),
        Lines = (
            linewidth = LW_LINE,
        ),
        Legend = (
            framevisible = false,
        ),
        Colorbar = (
            labelsize = 20,
            ticklabelsize = 18,
        ),
    )
end

# --------------------------------------------------------------------------
# rich-text label helpers (italic quantity + roman unit), matching the
# row_interaction_script typography. Prefer these over `L"..."` so the labels
# share the global font instead of the TeX math font.
# --------------------------------------------------------------------------

"`10^n` as rich text, e.g. for log-scaled tick labels."
powlab(n) = rich("10", superscript(string(n)))

# Pre-built rich labels (values, matching the row_interaction_script usage
# `xlabel = lab_time`). Built once and reused.
"Time axis label: t (s)."
const lab_time = rich(rich("t", font = :italic), "  (s)")

"Voltage axis label: U (V vs. SHE)."
const lab_voltage = rich(rich("U", font = :italic), "  (V vs. SHE)")

"Current-density axis label: I (mA cm^-2)."
const lab_current = rich(rich("I", font = :italic), "  (mA cm", superscript("−2"), ")")

"Apply the electrochemistry style globally. Called from `__init__`; reset with `set_theme!()`."
apply_electrochemistry_style!() = set_theme!(electrochemistry_theme())

# Apply the shared style as soon as the module is loaded (runtime, not precompile).
function __init__()
    apply_electrochemistry_style!()
    return nothing
end
