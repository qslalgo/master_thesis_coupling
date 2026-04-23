# Multiphysics Coupling: OpenMC – Thermal-Hydraulics Solver

## Overview

This repository contains the coupling framework developed for my master’s thesis.
The objective is to study the impact of thermal-hydraulic feedback on neutron transport in a research reactor configuration.

The workflow couples:

* **OpenMC** (neutronics)
* **Custom 2.5D thermal-hydraulics solver** (natural circulation, FEM-based)

The coupling is iterative and performed on an **assembly-wise and axial level**.

---

## Coupling Methodology

1. OpenMC computes **kappa-fission distributions** (power profiles) per pin and axial slice
2. Power is converted to **linear heat rate (q')**
3. Thermal-hydraulics solver computes:

   * temperature fields (fuel, cladding, coolant)
   * coolant density and flow velocities (natural circulation)
4. Updated temperatures and densities are fed back into OpenMC
5. Steps are repeated until convergence

---

## Repository Structure

```
coupling/
├── src/            # Core solver and coupling logic
├── meshes/         # 2D FEM meshes for subchannel / assembly modeling
├── scripts/        # Pre- and post-processing scripts
├── configs/        # Case definitions and parameters
├── results/        # Processed results (plots, small data)
```

---

## Requirements

* Python 3.x
* NumPy, SciPy
* OpenMC (with nuclear data libraries)
* meshio (for mesh handling)
* IAPWS (water properties)

---

## How to Run

1. Prepare OpenMC input and run neutronics simulation
2. Extract axial power profiles (kappa-fission tallies)
3. Run thermal-hydraulics solver:

   ```
   python run_coupling.py
   ```
4. Iterate coupling until convergence criteria are met

---

## Notes

* Large simulation outputs (`.h5`, statepoints) are not included in this repository
* Only processed data and essential inputs are stored for reproducibility

---

## Author

Enida Fatic
MSc Nuclear Engineering (Erasmus Mundus SARENA)

---

## Status

Work in progress – part of ongoing master’s thesis research.
