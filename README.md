# AuCO2RR.jl

**Fully integrated, time-dependent multiscale model of CO₂ electroreduction on gold**

Code accompanying *Maaß, Choi, Fuhrmann & Ringe, Electrochimica Acta (in press)* — see [Publication](#4-publication).

This repository provides a **CO₂ electroreduction reaction model on gold (CO2RR on Au)**, implemented in **Julia**.

The model was originally developed in **COMSOL** and has been **reconstructed in Julia**. The Julia implementation is validated against the COMSOL version, showing good agreement in key outputs (e.g. CV/IV trends and equilibrium behaviour).

Our primary development goal is to extend the model beyond the COMSOL baseline by implementing a more advanced **activity-coefficient** framework, enabling more realistic electrolyte thermodynamics and improved predictive capability.

It supports running standard electrochemical simulations:

- **CV** — cyclic voltammetry
- **IV** — polarization curves (I–V)
- **DLCap** — double-layer capacitance curves
- related analysis and plotting utilities

The main workflow is: load the package code in `src/`, then open a runnable notebook in `notebooks/` (or a plain script in `scripts/`) that performs the simulation and generates outputs.

---

## 1. Getting started

The package entry point `src/AuCO2RR.jl` **automatically loads** the model, including all relevant modules from `src/` and the plotting helpers in `src/AuCO2RR_plots/`. Each notebook in `notebooks/` loads this package at the top, so running a notebook executes the full simulation workflow without requiring manual includes.

Each simulation lives in its **own notebook** (e.g. `CyclicVoltammetry_notebook.jl` for CV). Simulations are **triggered interactively** inside each notebook.

### Requirements

| | |
|---|---|
| Julia | 1.9 or later |
| Required external package | [`CatmapInterface.jl`](https://github.com/ElCatFVM/CatmapInterface.jl) (see step 3) |

### 1) Clone

```bash
git clone https://github.com/ElCatFVM/Capacitance_Code
cd Capacitance_Code
```

### 2) Instantiate the Julia environment

Open Julia in the repository root, enter Pkg mode with `]`, then:

```julia
] activate .
```

```julia
(AuCO2RR)> instantiate
```

Optional but recommended:

```julia
(AuCO2RR)> precompile
```

### 3) Add `CatmapInterface.jl` (required)

`CatmapInterface.jl` is not registered in the Julia General registry, so it must be added by URL. In the same environment (still in Pkg mode):

```julia
(AuCO2RR)> add [https://github.com/ElCatFVM/CatmapInterface.jl](https://github.com/ElCatFVM/AuCO2RR.jl)
```

### 4) Run a notebook

```julia
using Pluto; Pluto.run()
```

Then open the notebook for the simulation you want:

| Notebook | Simulation |
|---|---|
| `notebooks/CyclicVoltammetry_notebook.jl` | CV |
| `notebooks/Capacitance_notebook.jl` | DLCap |
| `notebooks/polarization_notebook.jl` | IV / polarization |

---

## 2. Repository structure

### `data/`

**Experimental datasets** used for model validation and comparison (measured CV / IV / capacitance data).

### `notebooks/`

**Pluto notebooks** — this is where most of the runnable execution lives. Each notebook loads the `AuCO2RR` package and runs a specific simulation interactively.

| File | Purpose |
|---|---|
| `Capacitance_notebook.jl` | Double-layer capacitance (DLCap) |
| `CyclicVoltammetry_notebook.jl` | Cyclic voltammetry (CV) |
| `polarization_notebook.jl` | IV / polarization curves |
| `plot_publication_notebook.jl` | Publication-quality figures |
| `row_interaction_script.jl` | Row-interaction study (see also `Equilibriumcheck.jl`) |
| `sweepcompare_plots.jl` | Sweep comparison plots |

### `scripts/`

Plain runnable Julia scripts (non-notebook) for sweep comparisons and batch runs: `SweepCompare1.jl`, `SweepCompare2.jl`, `Sweeps0.jl`.

### `src/`

| File | Purpose |
|---|---|
| `AuCO2RR.jl` | Package entry point. Includes the modules below and the plotting helpers in `AuCO2RR_plots/`. Loaded by the notebooks and scripts |
| `cv.jl`, `iv.jl`, `dlcap.jl` | CV / IV / double-layer capacitance solvers |
| `goldmodel.jl` | Gold (Au) electrode and interface model definition |
| `sweeps_csv.jl` | Sweep setup and CSV export helpers |
| `overwritten_functions.jl` | Local overrides of upstream functions |
| `AuCO2RR_plots/` | Plotting helpers shared across simulations |

---

## 3. Typical usage flow

1. Edit or confirm the model definitions in `src/AuCO2RR.jl` (or its included files).
2. Open the notebook for the target simulation in Pluto — e.g. `notebooks/CyclicVoltammetry_notebook.jl`.
3. Run the notebook and inspect:
   - simulation results (data structs / exported files)
   - figures (CV / IV / DLCap plots)
   - logs (timing, convergence, parameter summary)

---

## 4. Publication

The methodology implemented in this repository is described in:

> S. Maaß\*, S. Choi\*, J. Fuhrmann, S. Ringe,
> *"Fully integrated, finite-volume-based, time-dependent multi-scale modeling of electrochemical CO₂ reduction"*,
> **Electrochimica Acta** (in press).
> \*These authors contributed equally.



---

## 5. Authors

- **Sumin Choi** — Korea University — model implementation, simulations, validation
- **Jürgen Fuhrmann** — Weierstrass Institute (WIAS), Berlin — numerical methods and solver framework

This package builds on [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl), [LiquidElectrolytes.jl](https://github.com/j-fu/LiquidElectrolytes.jl) and [CatmapInterface.jl](https://github.com/ElCatFVM/CatmapInterface.jl).
