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
            linewidth = 3,
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
