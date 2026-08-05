```@meta
CurrentModule = AuCO2RR
```

# Notations

## Species indexing

`GoldModel.SpeciesLayout` fixes the index of every unknown. Several routines slice the
unknown vector by **range** rather than by name, so the ordering is load-bearing.

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

Two contiguous ranges carry the reaction terms:

| Range | Field | Contents |
|:--|:--|:--|
| `2:6` | `ibufferstart:ibufferend` | the slice handed to the buffer rate law |
| `5:10` | `isurfacestart:isurfaceend` | the slice handed to the microkinetic rate law |

!!! warning "Mixed units inside `5:10`"
    Indices 5–7 are aqueous concentrations, 8–10 are dimensionless coverages. The
    microkinetic rate law writes turnover frequencies into all six, and only the aqueous
    ones are afterwards multiplied by the site density. Changing `SpeciesLayout` without
    updating that conversion mis-scales the fluxes silently.

The plotting layer re-declares the same indices as `ikplus`, `ihplus`, `ihco3`, `ico3`,
`ico2`, `iohminus`, `ico`.

## Symbols

| Symbol | Code | Meaning |
|:--|:--|:--|
| $c_i$ | `u[i]` | concentration of species $i$, mol m⁻³ |
| $c_i^\mathrm{bulk}$ | `c_bulk[i]` | bulk concentration |
| $z_i$ | `z[i]` | charge number |
| $D_i$ | `D[i]` | diffusion coefficient, m² s⁻¹ |
| $\gamma_i$ | `γ[i]` | activity coefficient |
| $a_i$ | — | activity, $\gamma_i c_i$ referenced to mol dm⁻³ |
| $v_i$ | `v[i]` | molar volume |
| $\kappa_i$ | `κ[i]` | solvation number |
| $\phi$ | `u[iϕ]` | electrostatic potential |
| $p$ | `u[ip]` | pressure |
| $\phi_\mathrm{we}$ | `ϕ_we`, `u[iϕ_we]` | working-electrode potential |
| $\phi_\mathrm{pzc}$ | `ϕ_pzc` | potential of zero charge |
| $\sigma$ | `σ` | electrode surface charge |
| $C_\mathrm{gap}$ | `C_gap` | Helmholtz gap capacitance |
| $S$ | `S` | catalytic site density, mol m⁻² |
| $\theta_j$ | `u[8:10]` | surface coverage |
| $\bar c$ | `c̄` | total site density of the solvent lattice |

## Rate constants

Named `k` + route + direction, following the `.mkm` and `ReactionData`:

| Prefix | Route |
|:--|:--|
| `kb…` | alkaline (OH⁻ as reactant) |
| `ka…` | acidic (H⁺ as product) |
| `kw…` | water autoprotolysis |

with suffix `e` for the equilibrium constant, `f` forward, `r` reverse — so `kbf2` is the
forward rate constant of the second alkaline reaction and `kbe2` its equilibrium constant.

## The three voltages

These are **not** interchangeable, and a result carries all three. `cv_abscissa` returns
each with a label that names it, so a figure cannot present one as another.

| `kind` | Field | Meaning |
|:--|:--|:--|
| `:applied` | `result.sawtooth` | the programmed protocol — what a potentiostat controls. Use this to compare against experiment. |
| `:reaction_plane` | `result.voltages` | potential at the reaction plane |
| `:dl` | `result.dlvoltages` | drop across the double layer |

Under an uncompensated sweep the difference is large: a $\pm 1.2$ V applied ramp can
appear as a $\pm 0.2$ V swing at the reaction plane, with the remainder lost to the
electrolyte.

## Units

Physical quantities are SI internally, with `LessUnitful`'s `@unitfactors` providing the
conversion factors. Two conventions to keep in mind:

- **Concentrations are mol m⁻³ internally**, while the literature and the plots use
  mol dm⁻³ (= M). The factor `ufac"mol/dm^3"` is 1000.
- **Currents are A m⁻² internally**; the plots multiply by `cm^2/mA` to reach the
  conventional mA cm⁻².
