# CO2RR_Gold (Au) Model 

This repository provides a **CO₂ electroreduction reaction model on gold (CO2RR on Au)**, implemented in **Julia**.
The model was originally developed in **COMSOL**, and has been **reconstructed to Julia**. The current Julia implementation is validated against the COMSOL version, showing good agreement in key outputs (e.g., CV/IV trends and equilibrium behavior).

Our primary development goal is to extend the model beyond the COMSOL baseline by implementing a more advanced **activity-coefficient** framework, enabling more realistic electrolyte thermodynamics and improved predictive capability.

It supports running standard electrochemical simulations such as:

- CV (Cyclic Voltammetry)
- IV / polarization curves (I–V curve)
- DLCap (Double Layer Capacitance curve) 
- (and related analysis/plot utilities)

The main workflow is: use the package code in src/, then execute a runnable script in script/ that performs the simulation and generates outputs.

---
## 1. How to start

- In the `script/EquilibriumCheck.jl` script, the **AuCO2RR code is automatically loaded**, which includes all relevant modules from `src/` and `plots/`. Running this script executes the full simulation workflow without requiring manual includes.
- All simulations (e.g., **double-layer capacitance (DLCap)**, **IV**, **CV**, and various parametric studies) are **triggered interactively**

### 1) Clone
```bash
git clone https://github.com/ElCatFVM/Capacitance_Code
cd Capacitance_Code
```
### 2) Instantiate Julia environment

Open Julia in the repository root, then enter Pkg mode (]) and run:
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
### 4) Run the main script
running script using `Pluto script/Equilibriumcheck.jl`

---
## 2. Repository structure

### src/CO2RR_Au.jl

- Acts as the package entry point for the $\text{CO}_2\text{RR}$ Au model.
- Includes (“-.jl”) modules/functions from sub-files in the downward directory, e.g. CV/IV solvers, plotting helpers(`../plots/`), constants/parameters.    
- Intended to be imported/loaded by runnable scripts.

### script/EquilibriumCheck.jl

- The primary execution script (driver) to run the model.
    
- Typical responsibilities:
    
    - load the package (src/CO2RR_Au.jl)
    - define simulation settings (electrolyte/parameters, voltage range, solver options)
    - run equilibrium checks and sweeps (e.g., CV/IV)
    - save results / generate plots
        
    
- In addition, the initial execution performs **no computations by default**. The script only **loads the model and initializes the interface**.
- All simulations (double-layer capacitance (DLCap), IV, CV, and various parametric plots) are triggered interactively:  
when a user selects the corresponding **checkbox** above a plot, the associated calculation is executed, and the relevant plots are generated automatically after the computation completes.

This design enables efficient, on-demand simulations without unnecessary recomputation.

---
## 3. Typical usage flow

1. Edit/confirm model definitions in src/CO2RR_Au.jl (or its included files)
    
2. Run the main driver:
    
    - julia script/EquilibriumCheck.jl
    
3. Script produces:
    
    - simulation results (data structs / exported files)
    - figures (CV/IV plots)
    - logs (timing, convergence, parameter summary)
  
---
## 4. E-acta Link

Electrochemica acta paper overleaf URL :
 https://www.overleaf.com/project/6952a316a06a94de96ff001d
    
