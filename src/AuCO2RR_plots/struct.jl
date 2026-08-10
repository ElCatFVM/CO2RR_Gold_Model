# Central plotting style for every AuCO2RR figure.
#
# `electrochemistry_theme()` returns the shared Makie Theme. It is also applied
# globally in `__init__` (via `set_theme!`), so figures that do NOT wrap
# themselves in `with_theme(...)` still inherit the same style and font. This
# keeps all plots consistent without repeating styling code in each function
# (efficiency), and lets us drop per-function `with_theme` wrappers over time.
#
# Typography: a single sans-serif family (`PLOT_FONT_FAMILY`) for all text, and
# one label pattern — `<Name> <Symbol>\n(<unit>)`, e.g. `Voltage U (V vs. SHE)`.
# Weight and slant carry meaning and nothing else does: the physical symbol is
# bold italic, the name and the unit are regular upright. Build labels with the
# `lab_*` / `powlab` helpers below rather than passing a bare string, so a new
# panel cannot drift from the pattern.
#
# Upright even though they look like symbols: `pH` (not a variable), species
# subscripts (`c_α`), and unit letters.
#
# `set_theme!` runs once, from `__init__`. Revise does NOT re-run `__init__`, so
# after editing this file call `AuCO2RR_plots.apply_electrochemistry_style!()` or
# the session keeps the old theme.
#
# To reset to Makie defaults in a session: `CairoMakie.set_theme!()`.

"""
Sans-serif family shared by every figure.

Was `"DejaVu Sans"` (row_interaction_script's `sans_font`), but a stock Debian/Ubuntu
installs only `DejaVuSans.ttf` and `DejaVuSans-Bold.ttf` — the oblique faces live in the
separate `fonts-dejavu-extra` package. Without them `:italic` / `:bold_italic` fail to
resolve and Makie silently falls back to the upright face, so a symbol set in
`:bold_italic` came out **bold but not slanted**.

TeX Gyre Heros ships *inside* Makie (`Makie.assetpath("fonts")`) with all four faces, so
it resolves on any machine with no system font install. It is a Helvetica clone, which is
what most journals ask for anyway.

Glyphs missing from Heros (`⇌`, some superscripts) are filled in per-glyph by Makie's
fallback chain, so the reaction-arrow legends still render.
"""
const PLOT_FONT_FAMILY = "TeX Gyre Heros Makie"

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

"""
Dusty rose → dusty blue sequence for CO₂ partial-pressure families.

Deliberately dark for a pastel: the same curve colour also sets the direct label written on
top of the data, and a pale pastel that reads fine as a 4 pt line disappears as text. Both
stops sit near 50 % lightness so either end carries at both sizes.
"""
const CMAP_PRESSURE = cgrad([colorant"#C97B8A", colorant"#7B9BC4"])

"""
Darker counterpart of [`CMAP_PRESSURE`](@ref), for measured data shown beside simulated.

Same two hues, so a pressure keeps its identity across the panels, but a step deeper — the
rows are told apart by weight of ink rather than by hue, which leaves hue free to mean
pressure and nothing else.
"""
const CMAP_PRESSURE_EXP = cgrad([colorant"#A34E60", colorant"#4E7196"])

"""
Deep blue → purple → red sequence for IR-compensation factor families.

Every stop stays dark enough to read as a line on white, unlike the pale end of the
perceptual maps, and the hue path is far enough from [`CMAP_SCANRATE`](@ref) and
[`CMAP_PRESSURE`](@ref) that a compensation figure is not mistaken for either.
"""
const CMAP_IRCOMP = cgrad([colorant"#2C3E8C", colorant"#7B3FA0", colorant"#C0392B"])

# --------------------------------------------------------------------------
# Font sizes. Axis labels and ticks were already set here; panel titles and
# legends were not, so they fell back to Makie's defaults and came out far
# smaller than everything around them.
# --------------------------------------------------------------------------

"Axis label and tick label size."
const FS_LABEL = 25

"Panel title size — a step below the axis labels so it reads as a caption."
const FS_TITLE = 22

"Legend title size."
const FS_LEGEND_TITLE = 22

"Legend entry size."
const FS_LEGEND = 20

function electrochemistry_theme()
    Theme(
        # Global font family: every textual element (ticks, titles, legends and
        # `rich(...)` labels) renders in one family. The `:regular/:bold/:italic/
        # :bold_italic` symbols used inside `rich(...)` resolve through this table.
        #
        # All four faces must exist, otherwise Makie falls back to the nearest
        # available one *without warning* and the style silently degrades — a
        # missing italic face is invisible in code review and only shows up as
        # upright text in the rendered figure. Check a face with
        #     Makie.to_font("TeX Gyre Heros Makie Italic")
        # which errors if the name does not resolve.
        fonts = (
            regular     = PLOT_FONT_FAMILY,
            bold        = "$(PLOT_FONT_FAMILY) Bold",
            italic      = "$(PLOT_FONT_FAMILY) Italic",
            bold_italic = "$(PLOT_FONT_FAMILY) Bold Italic",
        ),
        size = (960, 540),
        Axis = (
            spinewidth = 5.5,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            xticksize = 8,
            yticksize = 8,
            titlesize = FS_TITLE,
            titlefont = :bold,
            xlabelsize = FS_LABEL,
            ylabelsize = FS_LABEL,
            xticklabelsize = FS_LABEL,
            yticklabelsize = FS_LABEL,
            xgridvisible = false,
            ygridvisible = false,
            # Gap between the axis label and the tick labels. The label is measured
            # from the tick-label block, not the spine, so a long tick label (`10⁻¹²`,
            # `−10`) pushes the label outwards but leaves the *gap* unchanged — 10 px
            # under a 25 pt label reads as touching. Keep these two paddings equal to
            # the tick padding or larger.
            xlabelpadding = 16,
            ylabelpadding = 20,
            xticklabelpad = 6,
            yticklabelpad = 6,
            # Regular, not bold: weight is reserved for the physical symbol inside the
            # label, which is set to `:bold_italic` in the `lab_*` constants. A bold
            # label makes the whole string heavy and the symbol stops standing out.
            xlabelfont = :regular,
            ylabelfont = :regular,
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
            titlesize = FS_LEGEND_TITLE,
            # A legend title is an axis label for the swept variable, so it follows the
            # same rule: the name and the unit regular, only the symbol bold italic. Makie
            # defaults this to :bold, which overrode the `rich` parts and set the whole
            # string heavy — the symbol then stopped standing out from its own label.
            titlefont = :regular,
            labelsize = FS_LEGEND,
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
"""
Axis labels follow one pattern: **name, italic symbol, unit in brackets**, e.g.
`Voltage U (V vs. SHE)`. Spelling the quantity out keeps a figure readable on its own,
and every potential in this package is written `U`, never `φ`, so that panels can be put
side by side without the reader translating between two notations.
"""
const lab_time = rich("Time ", rich("t", font = :bold_italic), "\n(s)")

"Applied potential: Voltage U (V vs. SHE)."
const lab_voltage = rich("Voltage ", rich("U", font = :bold_italic), "\n(V vs. SHE)")

"Generic current density: Current I (mA cm⁻²)."
const lab_current = rich("Current ", rich("I", font = :bold_italic), "\n(mA cm", superscript("−2"), ")")

"Faradaic current read from the CO flux, the definition [`cv_current`](@ref) defaults to."
const lab_current_co = rich(
    "CO Partial Current ", rich("I", font = :bold_italic), "\n(mA cm", superscript("−2"), ")"
)

"""
The three terms a voltammetric current splits into, for figures that show them separately.

`I_F` is [`faradaic_current`](@ref), `I_C` is [`capacitive_current`](@ref), `I_tot` their
sum. The subscripts name which term is meant rather than standing for a quantity, so they
are upright — unlike the `α` of a species index, which is a running index and is italic.

Built from the same pieces as [`lab_current`](@ref) so a split figure and a plain one cannot
drift apart.

These three break after the **name** rather than after the symbol, unlike every other label
here. A rotated y label competes with the panel *height*, and these only ever appear on the
rows of a stacked grid, where a row is a couple of hundred pixels tall while
`Faradaic Current I_F` on one line needs more than that and spills into the row above. The
information and its order are unchanged; only the wrap point moves.
"""
const lab_current_F = rich(
    "Faradaic Current\n", rich("I", font = :bold_italic), subscript("F"),
    " (mA cm", superscript("−2"), ")"
)

"Capacitive term of the current. See [`lab_current_F`](@ref)."
const lab_current_C = rich(
    "Capacitive Current\n", rich("I", font = :bold_italic), subscript("C"),
    " (mA cm", superscript("−2"), ")"
)

"Faradaic plus capacitive. See [`lab_current_F`](@ref)."
const lab_current_tot = rich(
    "Total Current\n", rich("I", font = :bold_italic), subscript("tot"),
    " (mA cm", superscript("−2"), ")"
)

"""
CO₂ partial pressure: `CO₂ Partial Pressure p (atm)`.

Legend *titles* follow the same `<Name> <Symbol> (<unit>)` pattern as axis labels, so the
entries underneath can be bare numbers. Repeating the unit on every entry — `0.2 pCO2(atm)`,
`0.4 pCO2(atm)`, … — is what the swept-family legends used to do, and it both wastes the
legend width and drops the subscript.

The `2` of CO₂ is a stoichiometric index, not a running number, so it stays upright.
"""
const lab_pressure = rich(
    "CO", subscript("2"), " Partial Pressure\n",
    rich("p", font = :bold_italic), " (atm)"
)

"Scan rate: `Scan Rate v (V s⁻¹)`. Unit as a product of powers, not `V/s`."
const lab_scanrate = rich(
    "Scan Rate\n", rich("v", font = :bold_italic), " (V s", superscript("−1"), ")"
)

# --------------------------------------------------------------------------
# Symbol-only forms, for panels too short to fit a spelled-out label.
#
# A rotated y-label competes with the panel *height*, not its width. In a
# stacked figure a row is a couple of hundred pixels tall while a full label
# needs three or four hundred, so the text is clipped. Use these there and keep
# the spelled-out forms for single wide panels.
#
# Every label puts the unit on a second line: it halves the length a rotated
# label needs and keeps the bold-italic symbol on a line of its own.
# --------------------------------------------------------------------------

"Time axis, symbol only."
const lab_time_short = rich(rich("t", font = :bold_italic), "\n(s)")

"Voltage axis, symbol only."
const lab_voltage_short = rich(rich("U", font = :bold_italic), "\n(V vs. SHE)")

"Current-density axis, symbol only."
const lab_current_short = rich(rich("I", font = :bold_italic), "\n(mA cm", superscript("−2"), ")")

"Apply the electrochemistry style globally. Called from `__init__`; reset with `set_theme!()`."
apply_electrochemistry_style!() = set_theme!(electrochemistry_theme())

# Apply the shared style as soon as the module is loaded (runtime, not precompile).
function __init__()
    apply_electrochemistry_style!()
    return nothing
end
