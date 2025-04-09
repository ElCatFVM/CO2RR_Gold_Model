<<<<<<< HEAD
### Explain
 - For DLCap_03.jl, you may encounter an error saying that v_0 is not defined if you run the code directly. To avoid this, please make sure to run it within the liquidelectrolyte.jl package environment.

 - All reactions are currently turned off, and mass transport has not been defined yet.

 - As for equilibrium_pluto, I’ve attached a slightly modified version of the code used to solve dlcapsweep_equi.
=======
[![linux-macos-windows](https://github.com/j-fu/LiquidElectrolytes.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/j-fu/LiquidElectrolytes.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/LiquidElectrolytes.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/LiquidElectrolytes.jl/stable)


LiquidElectrolytes.jl
=====================

This package is a solver for generalized Poisson-Nernst-Planck models taking into account finite ion size, solvation and electroosmotic pressure based on formulations [derived from first principles of nonequilibrium thermodymamics](https://doi.org/10.1016/j.elecom.2014.03.015).  It utilizes a [thermodynamically consistent finite volume space discretization approach](https://doi.org/10.1007/s00211-022-01279-y) which uses the sum of the electrotstatic potential and the excess chemical potential as convective terms. It is realized on top of the solver kernel [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl) for coupled nonlinear PDE systems which takes advantage of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to generate the full linearization of the coupled nonlinear system as a basis for a robust Newton solver for the discrete nonlinear systems.

Parts of the package are still work in progress.

>>>>>>> Edited Based on LiquidElectrolyte.jl/notebook/EquilibriumCheck
