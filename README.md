# 6-Stroke Internal Combustion Engine Simulation

<img src="animation/porsche_6stroke_engine.gif" width="600" alt="6-Stroke Engine Animation">

A MATLAB implementation for steady-state simulation of a novel 6-stroke internal combustion engine architecture. This codebase provides a complete simulation framework, from geometric analysis to thermodynamic modeling and dynamic visualization.

## Overview

This repository contains a simulation framework for analyzing a 6-stroke engine cycle inspired by [US Patent US20240301817A1](https://patents.google.com/patent/US20240301817A1/en). The engine employs a planetary gear mechanism to enable variable stroke lengths, resulting in a unique 1080° crank-angle cycle comprising six distinct strokes (two long, two short, and two long strokes in alternating pattern).

### Technical Approach

The simulation employs a **single-zone thermodynamics model** with three state variables per cylinder:
- **Cylinder pressure** (p)
- **Trapped gas mass** (m)
- **Oxygen mass** (mO₂)

The thermodynamic modeling combines:
- **Wiebe-based heat release** functions for two separate combustion events
- **Woschni correlation** for convective heat transfer
- **Compressible nozzle flow** models for intake, exhaust, and scavenge ports

The thermodynamics are coupled to an **8-DOF multibody dynamics model** representing a flat-six boxer engine configuration:
- 6 rotational DOFs for individual planet gears
- 1 DOF for the carrier/crankshaft
- 1 DOF for the output shaft

A dual-mass flywheel (DMF) model dampens torque fluctuations between the engine and gearbox. The coupled system's equations of motion are integrated in the crank-angle domain (0-1080°), with a cycle-to-cycle fixed-point iteration used to find periodic steady-state solutions.

## Code Structure

The codebase is organized into three main directories:

### `preprocessing/`
Functions for parameter loading, geometric analysis, and initialization:
- `parameters.m` - Engine parameters and operating conditions
- `derived_parameters.m` - Computed geometric and kinematic parameters
- `extrema_6stroke.m` - TDC/BDC computation for all six strokes
- `manifold_geometry.m` - Intake/exhaust/scavenge port geometry setup
- `initialize_particles.m` - Particle initialization for visualization

### `simulation/`
Core simulation engine including thermodynamics and dynamics:
- `engine_cycle_steady_6stroke.m` - Single-cylinder cycle simulation with fixed-point iteration
- `full_engine_simulation_6stroke.m` - 8-DOF flat-six engine dynamics
- `piston_kinematics.m` - Kinematic relationships for planetary mechanism
- `generalized_gas_force.m` - Gas pressure → generalized force conversion
- `generalized_mesh_force.m` - Gear mesh stiffness modeling
- `generalized_DMF_force.m` - Dual-mass flywheel model

### `animation/`
Visualization and post-processing:
- `plot_results_6stroke.m` - Comprehensive plotting of simulation results
- `animate_results_6stroke.m` - Animated visualization of engine cycle
- `plot_engine.m` - 2D cross-section engine visualization
- `particle_simulation.m` - Particle-based flow visualization
- Various `Draw*.m` functions for rendering engine components

## Getting Started

### Requirements

- MATLAB R2016b or later
- Optimization Toolbox (for `fzero` and ODE solvers)

### Running the Simulation

1. **Clone or download this repository**

2. **Open MATLAB and navigate to the repository directory**

3. **Run the main script:**
   ```matlab
   Main_engine
   ```

The script automatically adds the required subdirectories to MATLAB's path and executes the complete simulation workflow:
1. Parameter loading and geometric preprocessing
2. Periodic steady-state calculation (single cylinder)
3. Full 6-cylinder engine simulation with dynamics
4. Result visualization and animation

### Customization

Engine parameters can be modified in `parameters.m`, including:
- Operating conditions (RPM, number of cylinders)
- Geometric parameters (bore, stroke lengths, gear ratios)
- Combustion timing and fuel split
- Valve and port timing
- Heat transfer coefficients

## Model Limitations and Notes

This codebase represents a streamlined version intended for public release. While it captures the essential physics and dynamics of the 6-stroke cycle, certain advanced features and refinements present in more detailed implementations may not be included. The modeling approach prioritizes computational efficiency while maintaining physical accuracy for steady-state analysis.

Users seeking to extend or modify the code may find that some simplifying assumptions have been made to keep the codebase accessible and maintainable. The core algorithms remain intact and can serve as a solid foundation for further development.

## Converting to Other Languages

While this codebase is written in MATLAB, the algorithms and mathematical formulations are language-agnostic. Users comfortable with Python, C++, Julia, or other programming languages may find it straightforward to port this codebase. Modern AI-assisted coding tools can significantly streamline the translation process, particularly for the vectorized numerical operations and ODE integration routines that form the core of this simulation.

The modular structure of the codebase—with clear separation between preprocessing, simulation, and visualization—also facilitates translation into object-oriented frameworks where appropriate.

## Contact

Bart Blockmans  
Email: bart@blockmans.net

---

**Note:** This code is provided as-is for educational and research purposes. Users should verify model assumptions and results against experimental data or more detailed simulations before making engineering decisions based on these simulations.

