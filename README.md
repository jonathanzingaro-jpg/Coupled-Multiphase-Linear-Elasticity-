# Coupled Vertically Averaged Multiphase Flow – Linear Poroelasticity (FEM)

MATLAB finite element code for simulating CO₂ injection into a deep saline reservoir and the resulting surface uplift. A vertically averaged (vertical equilibrium) two-phase flow model computes the pressure build-up and CO₂ plume, and is one-way coupled to an axisymmetric linear poroelastic model of the reservoir and surrounding rock to compute ground deformation.

This is the numerical model used to generate the synthetic dataset in:

> G. Zingaro, S. Hatefi Ardakani, R. Gracie, Y. Leonenko (2025). **Deep learning assisted monitoring framework for geological carbon sequestration.** *International Journal of Greenhouse Gas Control*, 144, 104372. https://doi.org/10.1016/j.ijggc.2025.104372

---

## Model Overview

The problem is a single vertical injection well in a horizontal, homogeneous reservoir bounded above and below by low-permeability seals (overburden and underburden). Because of radial symmetry:

- the vertically averaged flow model reduces to a **1D radial** problem, and
- the poroelastic model is solved as a **2D axisymmetric** (r–z) problem.

Both are discretized in space with the **Galerkin finite element method**. The flow model is integrated in time with finite differences.

### Solution strategy

```
Initialize h_B, compute initial p
for each time step:
    compute pressure build-up p(r, t)        ┐ sequentially coupled
    compute brine/CO₂ interface h_B(r, t)    ┘  (vertically averaged flow)
    compute displacement d(r, z, t)            one-way coupled (poroelasticity)
```

### Vertically Averaged Two-Phase Flow

Mass conservation for the brine (B) and CO₂ (C) phases:

```math
\phi\,(1 - S_B)\,\frac{\partial h_B}{\partial t} + \nabla\cdot\hat{\mathbf{q}}_B = 0
```

```math
\phi\,(1 - S_B)\,\frac{\partial h_C}{\partial t} + \nabla\cdot\hat{\mathbf{q}}_C = Q_{inj}
```

The vertically averaged fluxes follow the multiphase extension of Darcy's law:

```math
\hat{\mathbf{q}}_B = -h_B\,\frac{\mathbf{k}}{\mu_B}\,\nabla p
```

```math
\hat{\mathbf{q}}_C = -h_C\,\frac{\mathbf{k}\,k_C^{rel}}{\mu_C}\left(\nabla p + \Delta\rho\, g\,\nabla h_C\right)
```

where

| Symbol | Description |
|---|---|
| $`r`$ | Radial distance from the wellbore |
| $`\phi`$ | Porosity of the reservoir |
| $`S_B`$ | Residual saturation of the brine phase |
| $`h_B(r,t)`$, $`h_C(r,t)`$ | Height of the brine and CO₂ phases |
| $`\hat{\mathbf{q}}_B(r,t)`$, $`\hat{\mathbf{q}}_C(r,t)`$ | Vertically averaged brine and CO₂ fluxes |
| $`Q_{inj}`$ | CO₂ injection source term |
| $`\mathbf{k}`$ | Permeability tensor of the rock mass |
| $`k_C^{rel}`$ | Relative permeability of CO₂ |
| $`\mu_B`$, $`\mu_C`$ | Viscosity of brine and CO₂ |
| $`p(r,t)`$ | Pressure build-up in the reservoir |
| $`\Delta\rho = \rho_B - \rho_C`$ | Density difference between brine and CO₂ |
| $`g`$ | Gravitational acceleration |
| $`H = h_B + h_C`$ | Reservoir thickness |

Because the phases are immiscible with a sharp interface, the CO₂ saturation field is expressed in terms of the plume height $`h_C`$. Streamline Upwind (SU) stabilization is used for numerical stability (Ladubec et al., 2015; Ladubec & Gracie, 2017).

**Assumptions:** infinite, fully saturated reservoir; sharp, immiscible CO₂–brine interface; incompressible brine, CO₂ and rock; constant fluid viscosities and densities; no capillary pressure, viscous fingering, thermal or chemical effects; homogeneous, linear isotropic permeability.

**Initial and boundary conditions:** reservoir initially fully brine-saturated, $`s_B(r, t=0) = 1`$; impermeable top and bottom; infinite-acting far-field boundary with $`s_B(r_\infty, t) = 1`$ and $`p(r_\infty, t) = 0`$.

### Linear Poroelasticity

Quasi-static balance of linear momentum for the solid skeleton:

```math
\nabla\cdot\boldsymbol{\sigma} + \mathbf{b} = \mathbf{0}
```

Principle of effective stress:

```math
\boldsymbol{\sigma} = \boldsymbol{\sigma}^{\prime} - \alpha\,\mathbf{m}\,p
```

Linear elastic constitutive law:

```math
\boldsymbol{\sigma}^{\prime} = \mathbf{C} : \boldsymbol{\epsilon}
```

Small-strain kinematics:

```math
\boldsymbol{\epsilon} = \frac{1}{2}\left(\nabla\mathbf{d} + \nabla\mathbf{d}^{\top}\right)
```

where

| Symbol | Description |
|---|---|
| $`\boldsymbol{\sigma}`$ | Cauchy total stress tensor |
| $`\boldsymbol{\sigma}^{\prime}`$ | Effective stress tensor |
| $`\mathbf{b}`$ | Body force |
| $`\alpha`$ | Biot's coefficient |
| $`\mathbf{m}`$ | Vectorized Kronecker delta |
| $`p`$ | Pressure build-up from the flow model |
| $`\mathbf{C}`$ | Fourth-order stiffness tensor of the solid skeleton |
| $`\boldsymbol{\epsilon}`$ | Strain tensor |
| $`\mathbf{d}`$ | Displacement vector |

The pressure build-up from the flow model is applied as a load on the poroelastic mesh (reservoir, overburden, underburden), giving the displacement field and surface uplift at the measurement sensor locations.

---

## Repository structure

```
├── Vertically Averaged Two-Phase Flow/
│   ├── RunSim_SW_INJ.m              # Run script: CO₂ injection simulation
│   ├── Functions/
│   │   ├── Main/                    # Solver driver, stabilization selection
│   │   ├── Assembly/                # Element/global matrices and force vectors
│   │   ├── Adaptive Time Step/      # Time-step control
│   │   ├── Linear Solver/           # Dirichlet BCs and linear solve
│   │   ├── Mesh Generation/         # Structured, imported, Excel and HyperMesh meshes
│   │   ├── Shape Functions and Quadrature/
│   │   ├── Error Computation/       # Mass balance and manufactured-solution errors
│   │   └── Post Processing/         # VTK/text output, pore-pressure fitting
│   ├── Convergence Study/           # Method of manufactured solutions (MMS)
│   └── Nodal Coordinates/
│
└── Linear Poroelasticity/
    ├── RunSim.m                     # Run script: axisymmetric poroelastic model
    ├── Functions/
    │   ├── main.m
    │   ├── Assembly/                # Stiffness matrix, external force, strain
    │   ├── Constitutive Law/        # Elastic stiffness matrix D
    │   ├── Fitting/                 # Maps 1D pressure field onto the 2D mesh
    │   ├── Linear Solver/
    │   ├── Mesh Generation/
    │   ├── Shape Functions and Quadrature/
    │   ├── Error Computation/       # Manufactured-solution check
    │   └── Post Processing/         # VTK and displacement output
    └── Mesh/                        # 25 km × 100 m subsurface meshes (Q4, Q8), sensor locations
```

---

## Citation

If you use this code, please cite:

```bibtex
@article{zingaro2025deep,
  title   = {Deep learning assisted monitoring framework for geological carbon sequestration},
  author  = {Zingaro, Giovanni and Hatefi Ardakani, Saeed and Gracie, Robert and Leonenko, Yuri},
  journal = {International Journal of Greenhouse Gas Control},
  volume  = {144},
  pages   = {104372},
  year    = {2025},
  doi     = {10.1016/j.ijggc.2025.104372}
}
```

## References

- Celia, M.A., Nordbotten, J.M. (2009). Practical modeling approaches for geological storage of carbon dioxide. *Groundwater*, 47(5), 627–638.
- Ladubec, C., Gracie, R., Craig, J. (2015). An extended finite element method model for carbon sequestration. *IJNME*, 102(3–4), 316–331.
- Ladubec, C., Gracie, R. (2017). Stabilized finite element methods for vertically averaged multiphase flow for carbon sequestration. *IJNME*, 111(8), 724–752.
- Nordbotten, J.M., Celia, M.A. (2011). *Geological Storage of CO₂: Modeling Approaches for Large-Scale Simulation*. Wiley.

## Acknowledgments

Funded by an NSERC Alliance Missions Grant supporting collaboration with Stelco Inc. Department of Civil and Environmental Engineering, University of Waterloo.

## Contact

Giovanni Zingaro, gzingaro@uwaterloo.ca
