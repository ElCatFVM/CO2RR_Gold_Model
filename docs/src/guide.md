```@meta
CurrentModule = AuCO2RR
```

# Guide

## Species layout

`GoldModel.SpeciesLayout` fixes the index of every unknown. Several parts of the code
slice the unknown vector by range rather than by name, so the ordering is load-bearing.

| Index | Species | Transported | Buffer reactions | Surface mechanism |
|:--|:--|:--|:--|:--|
| 1 | K⁺ | yes | — | — |
| 2 | H⁺ | yes | yes | — |
| 3 | HCO₃⁻ | yes | yes | — |
| 4 | CO₃²⁻ | yes | yes | — |
| 5 | CO₂ | yes | yes | yes |
| 6 | OH⁻ | yes | yes | yes |
| 7 | CO | yes | — | yes |
| 8–10 | CO\*, COOH\*, CO₂\* | boundary only (coverages) | — | yes |

Two contiguous ranges matter:

- `ibufferstart:ibufferend` = `2:6` — the slice handed to the buffer rate law
- `isurfacestart:isurfaceend` = `5:10` — the slice handed to the microkinetic rate law

!!! warning "Mixed units inside `5:10`"
    Indices 5–7 are aqueous concentrations, 8–10 are dimensionless coverages. The
    microkinetic function writes turnover frequencies into all six, and only the aqueous
    ones are afterwards multiplied by the site density. Changing `SpeciesLayout` without
    updating that conversion will silently mis-scale the fluxes.

## Buffer network

Five reversible reactions act in the volume (`GoldModel.buffer_system`):

```math
\begin{aligned}
\ce{CO2 + OH^- &<=> HCO3^-} \\
\ce{HCO3^- + OH^- &<=> CO3^{2-} + H2O} \\
\ce{CO2 + H2O &<=> HCO3^- + H^+} \\
\ce{HCO3^- &<=> CO3^{2-} + H^+} \\
\ce{H2O &<=> H^+ + OH^-}
\end{aligned}
```

The alkaline and acidic routes are thermodynamically consistent — the equilibrium
constants satisfy $K_\mathrm{b} = K_\mathrm{a}/K_\mathrm{w}$ — but they are *kinetically*
independent, so the network can sit far from the water equilibrium locally.

!!! note "Water is not locally equilibrated at the electrode"
    Under a driven CV the ion product at the electrode departs from $K_\mathrm{w}$ by
    several decades ($Q_\mathrm{w}/K_\mathrm{w} \approx 4\times10^{-6}$ has been measured
    on the anodic branch). Anything that infers one of H⁺/OH⁻ from the other through
    $K_\mathrm{w}$ is therefore invalid near the electrode, even though it is fine in the
    bulk.

Reaction 2 is the stiff one: with a surface bicarbonate concentration of order 1 M its
pseudo-first-order rate constant for OH⁻ reaches $\sim 10^{8}\,\mathrm{s^{-1}}$, i.e. a
time constant of nanoseconds against time steps of order 0.1 s. It is also the main
*supplier* of OH⁻ near the electrode — slowing it down makes the OH⁻ balance worse, not
better.

## From turnover frequency to boundary flux

`GoldModel.surface_reaction` returns a rate law whose output is a **turnover frequency**,
i.e. per catalytic site per second. Coverages want exactly that. Aqueous species want a
flux per unit electrode area, so they are multiplied by the site density $S$:

```math
\underbrace{r}_{\mathrm{s^{-1}}} \times \underbrace{S}_{\mathrm{mol\,sites\,m^{-2}}}
= \underbrace{j}_{\mathrm{mol\,m^{-2}\,s^{-1}}}
```

With `ReactionData().S ≈ 1.6\times10^{-8}\ \mathrm{mol\,m^{-2}}` — about $10^{16}$ sites
per m², roughly 0.07 % of the gold surface atom density, so $S$ encodes an *active site
fraction* rather than a geometric count.

!!! important "$S$ scales the whole faradaic current"
    Steady-state coverages do not depend on $S$, so $I_\mathrm{F} = n_e F \, r \, S$ is
    linear in $S$. Matching an absolute measured current is therefore, to a large extent,
    fitting $S$.

## Booking the proton stoichiometry

The mechanism in `data/catmap_CO2R_data/catmap_CO2R_template.mkm` writes the two
electron transfers with water as the proton donor:

```math
\ce{CO2* + H2O + e^- <=> COOH* + OH^-}, \qquad
\ce{COOH* + e^- <=> CO* + OH^-}
```

Read literally, the reverse (anodic) direction consumes free OH⁻. That is not something
the electrolyte can support:

| Reservoir | Concentration | Max diffusive supply | vs. demand $\sim 10^{-4}\,\mathrm{mol\,m^{-2}\,s^{-1}}$ |
|:--|--:|--:|--:|
| OH⁻ | $6\times10^{-5}$ mol m⁻³ | $3\times10^{-9}$ | 1/37 000 |
| H⁺ | $1.6\times10^{-4}$ mol m⁻³ | $1.5\times10^{-8}$ | 1/8 200 |
| H₂O | $5.6\times10^{4}$ mol m⁻³ | $5.5\times10^{-1}$ | 4 600× surplus |

A boundary sink that does not vanish as its species is depleted has no non-negative
solution once it exceeds the maximum local supply; the Newton solve then converges to a
negative concentration. This is a property of the discrete system, not of the solver
settings — refining the time step, tightening tolerances, refining the grid and slowing
the buffer all leave it unchanged.

`we_breactions` therefore books the stoichiometry on whichever species is being
*produced*, branch-free so the sparsity tracer can walk both paths:

```julia
r_oh  = f[iohminus] * S
f[iohminus] = min(r_oh, zero(r_oh))    # cathodic: water reduced   → OH⁻ produced
f[ihplus]  -= max(r_oh, zero(r_oh))    # anodic:   proton released → H⁺ produced
```

Charge balance is unaffected (one ion per electron either way), and the split between H⁺
and OH⁻ is left to the autoprotolysis reaction that is already in the buffer network.
Where water *is* locally equilibrated the two formulations are identical, so the cathodic
branch is bit-for-bit unchanged.

!!! warning "This is a mechanism choice, not bookkeeping"
    Microscopic reversibility ties the reverse of a water-donor step to OH⁻ consumption.
    Routing the anodic direction to H⁺ production is the *proton route*, which is not what
    the `.mkm` text says. Either rewrite steps 2–3 in the `.mkm` as proton transfers, or
    feed the transported OH⁻ activity into the rate law and accept a smaller anodic
    current. The code and the mechanism file should not disagree.

## Defining "the current"

Boundary flux is only a faradaic current for a species with **no homogeneous reaction**.
In this model that is CO alone — CO₂, HCO₃⁻, CO₃²⁻, OH⁻ and H⁺ are all buffer-active, and
K⁺ is a spectator.

Measured on this model, `currents(result, ico2)` runs about twice the faradaic rate,
because each reduced CO₂ releases two OH⁻ which consume further CO₂ through
$\ce{CO2 + OH^- <=> HCO3^-}$.

All CV figures go through one definition:

```julia
faradaic_current(result; species = ico, n_e = 2, sgn = 1)
capacitive_current(result)                      # result.j_cap, valid in every IR mode
cv_current(result; include_capacitive = false, kwargs...)
cv_abscissa(result; kind = :applied)            # returns (values, axis label)
```

`cv_abscissa` exists because three different voltages live in a result and are not
interchangeable:

| `kind` | Quantity | Use for |
|:--|:--|:--|
| `:applied` | `result.sawtooth`, the programmed protocol | comparing against experiment |
| `:reaction_plane` | `result.voltages` | reaction-plane potential |
| `:dl` | `result.dlvoltages` | drop across the double layer |

Each returns its own axis label, so a figure always names the voltage it plots.

!!! tip "Do not read the `icc` unknown for the capacitive current"
    `icc` only exists when `ircompensation isa OhmicDropEstimation`. Reading it makes the
    capacitive term evaluate to zero in every other mode, without any warning, and
    currents from different compensation modes stop being comparable. Use
    `capacitive_current`, which reads `result.j_cap`.

## Validity limits

- **Potential window.** The CatMAP input declares
  `descriptor_ranges = [[-1.5, 0.0], [298, 298]]`. Sweeping to $+0.8$ V vs. SHE
  extrapolates the surface energetics roughly 0.8 V beyond where they were parameterised.
  The anodic branch should be read as extrapolation.
- **Positivity.** Check it rather than assume it. Scan rate, boundary-layer thickness $L$
  and the anodic vertex all move the boundary between well-posed and not; a negative
  concentration anywhere invalidates the run.
- **Absolute currents.** These scale with $S$ and with the choice of `species`/`n_e`;
  curve *shapes* are far more robust than magnitudes.
