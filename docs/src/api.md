```@meta
CurrentModule = AuCO2RR
```

# API

```@contents
Pages = ["api.md"]
Depth = 3
```

## Current and voltage

Every CV figure is built from these four. [Plots](@ref) explains why `species` defaults to
CO and why the abscissa carries its own label.

```@docs
AuCO2RR_plots.faradaic_current
AuCO2RR_plots.capacitive_current
AuCO2RR_plots.cv_current
AuCO2RR_plots.cv_abscissa
```

## Model construction

```@autodocs
Modules = [AuCO2RR.GoldModel]
Order = [:type, :function, :constant]
```

## Sweep drivers

```@autodocs
Modules = [AuCO2RR]
Order = [:type, :function, :constant]
```

## Plotting

Line widths, colour sequences and axis labels are shared constants rather than literals at
each call site — reach for one of these instead of a number when adding a figure.

```@autodocs
Modules = [AuCO2RR.AuCO2RR_plots]
Order = [:type, :function, :constant]
Filter = x -> !(x isa Function && nameof(x) in (:faradaic_current, :capacitive_current, :cv_current, :cv_abscissa))
```

## Index

```@index
```
