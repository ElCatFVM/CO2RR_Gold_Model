```@meta
CurrentModule = AuCO2RR
```

# Public API

```@contents
Pages = ["public.md"]
Depth = 3
```

## Current and voltage definitions

Every CV figure is built from these. See
[Defining "the current"](@ref) for why `species` defaults to CO.

```@docs
AuCO2RR_plots.faradaic_current
AuCO2RR_plots.capacitive_current
AuCO2RR_plots.cv_current
AuCO2RR_plots.cv_abscissa
```

## Model construction

```@autodocs
Modules = [AuCO2RR.GoldModel]
Order = [:type, :function]
```

## Sweep drivers

```@autodocs
Modules = [AuCO2RR]
Order = [:function]
```

## Plotting

```@autodocs
Modules = [AuCO2RR.AuCO2RR_plots]
Order = [:function]
Filter = f -> !(nameof(f) in (:faradaic_current, :capacitive_current, :cv_current, :cv_abscissa))
```

## Index

```@index
```
