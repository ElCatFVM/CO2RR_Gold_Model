# CO2RR_Gold (Au) Model

This repository provides a **CO₂ electroreduction reaction model on gold (CO2RR on Au)**, implemented in **Julia**.

The model was originally developed in **COMSOL**, and has been **reconstructed to Julia**. The current Julia implementation is validated against the COMSOL version, showing good agreement in key outputs (e.g., CV/IV trends and equilibrium behavior).

Our primary development goal is to extend the model beyond the COMSOL baseline by implementing a more advanced **activity-coefficient** framework, enabling more realistic electrolyte thermodynamics and improved predictive capability.

It supports running standard electrochemical simulations such as:

- CV (Cyclic Voltammetry)
- IV / polarization curves (I–V curve)
- DLCap (Double Layer Capacitance curve)
- (and related analysis/plot utilities)

The main workflow is: use the package code in `src/`, then open a runnable notebook in `notebooks/` (or a plain script in `scripts/`) that performs the simulation and generates outputs.

---

## 1. How to start

- The package entry point `src/AuCO2RR.jl` **automatically loads** the model, including all relevant modules from `src/` and the plotting helpers in `src/AuCO2RR_plots/`. Each notebook in `notebooks/` loads this package at the top, so running a notebook executes the full simulation workflow without requiring manual includes.
- Each simulation now lives in its **own notebook** (e.g., `Capacitance_notebook.jl` for DLCap, `CyclicVoltammetry_notebook.jl` for CV, `polarization_notebook.jl` for IV/polarization). Simulations are **triggered interactively** inside each notebook.

### 1) Clone

```bash
git clone https://github.com/ElCatFVM/Capacitance_Code
cd CO2RR_Gold
```

### 2) Instantiate Julia environment

Open Julia in the repository root, then enter Pkg mode (`]`) and run:

```julia
] activate .
```

```julia
(AuCO2RR)> instantiate
```

(Optional but recommended)

```julia
(AuCO2RR)> precompile
```

### 3) Add CatmapInterface.jl (required)

In the same Julia environment (still in Pkg mode):

```julia
(AuCO2RR)> add https://github.com/ElCatFVM/CatmapInterface.jl
```

### 4) Run a notebook

Open Pluto and run the notebook for the simulation you want, e.g.:

```julia
using Pluto; Pluto.run()
```

then open `notebooks/CyclicVoltammetry_notebook.jl` (CV), `notebooks/Capacitance_notebook.jl` (DLCap), or `notebooks/polarization_notebook.jl` (IV).

---

## 2. Repository structure

### data/

- Contains the **experimental datasets** used for model validation and comparison (e.g., measured CV/IV/capacitance data).

### notebooks/

- Holds the **Pluto notebooks**, which is where **most of the runnable execution scripts**.
- Each notebook loads the `AuCO2RR` package and runs a specific simulation interactively:
  - `Capacitance_notebook.jl` — double-layer capacitance (DLCap)
  - `CyclicVoltammetry_notebook.jl` — cyclic voltammetry (CV)
  - `polarization_notebook.jl` — IV / polarization curves
  - `plot_publication_notebook.jl` — publication-quality figures
  - `row_interaction_script.jl`, `sweepcompare_plots.jl` — row_intraction_script(Equilibriumcheck.jl)

### scripts/

- Plain runnable Julia scripts (non-notebook), e.g. `SweepCompare1.jl`, `SweepCompare2.jl`, `Sweeps0.jl`, used for sweep comparisons and batch runs.

### src/

- `AuCO2RR.jl` — the package entry point for the $\text{CO}_2\text{RR}$ Au model. It includes (`include(".../-.jl")`) the modules/functions from the sub-files below and the plotting helpers in `AuCO2RR_plots/`. Intended to be loaded by the notebooks/scripts.
- `cv.jl`, `iv.jl`, `dlcap.jl` — CV / IV / double-layer capacitance solvers.
- `goldmodel.jl` — the gold (Au) electrode / interface model definition.
- `sweeps_csv.jl` — sweep setup and CSV export helpers.
- `overwritten_functions.jl` — local overrides of upstream functions.
- `AuCO2RR_plots/` — plotting helpers used across the simulations.

---

## 3. Typical usage flow

1. Edit/confirm model definitions in `src/AuCO2RR.jl` (or its included files).

2. Open the notebook for the target simulation (e.g., `notebooks/CyclicVoltammetry_notebook.jl`) in Pluto. After, you can checks:

   - simulation results (data structs / exported files)
   - figures (CV/IV/DLCap plots)
   - logs (timing, convergence, parameter summary)

---

## 4. E-acta Link

Electrochemica acta paper overleaf:
https://www.overleaf.com/project/6952a316a06a94de96ff001d
