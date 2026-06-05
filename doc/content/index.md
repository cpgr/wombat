!config navigation breadcrumbs=False scrollspy=False

# wombat — WOrkflow for Modelling Buoyant Aquifer Transport

**A vertical-equilibrium (VE) solver for geological CO2 storage, built on the
[MOOSE](https://mooseframework.inl.gov) framework and its PorousFlow module.**

`wombat` depth-integrates the multiphase porous-flow equations under the assumption of
vertical pressure equilibrium, reducing a 3D CO2-storage problem to a 2D problem on the
top surface of the storage formation.  This is appropriate for CO2 plumes that are thin
relative to their lateral extent — the common case at field scale and for long-timescale,
post-injection migration studies — and is dramatically cheaper than a full 3D simulation
while retaining the dominant physics: buoyant up-dip migration, structural trapping,
residual trapping, and dissolution.

## Governing equations

Per fluid component $c$, depth-integrated mass conservation reads:

\begin{equation}
\frac{\partial}{\partial t}\!\left(H\,\bar{\phi}\,\rho_c\,\bar{S}_c\right)
  + \nabla \cdot F_c - q_c = 0
\end{equation}

\begin{equation}
F_c = -H\,\mathbf{K}_{up}\,\frac{k_{r,c}^{up}\,\rho_c}{\mu_c}
      \left(\nabla p_{top} + \rho_c\,g\,\nabla z_T\right)
\end{equation}

The $\rho_c\,g\,\nabla z_T$ term — the **sloped-top-surface buoyancy drive** — is the
dominant up-dip migration mechanism at field scale.  Primary variables are $p_{top}$
(pore pressure at the formation top) and $\bar{S}_n$ (depth-averaged CO2 saturation).

## Documentation

- [Object Reference](source/index.md) — all registered kernels, materials, user objects,
  and physics actions grouped by type with full parameter listings
- [Cross-code benchmark](benchmarks/cross_code.md) — comparison against the MRST co2lab
  VE reference implementation

## Capabilities

| Area | Objects |
|------|---------|
| **Physics actions** | `VEFlowFE` (continuous Galerkin), `VEFlowFV` (cell-centred FV) |
| **Mass storage / flux (FE)** | `VEMassTimeDerivative`, `VEAdvectiveFluxS`, `VEAdvectiveFluxP` |
| **Mass storage / flux (FV)** | `VEFVMassTimeDerivative`, `VEFVAdvectiveFlux` |
| **Geometry / petrophysics** | `VEGeometry`, `VEPorosity`, `VEPermeability` |
| **Plume reconstruction** | `VEPlumeReconstruction` (sharp-interface and capillary-fringe) |
| **Relative permeability** | `VERelPermSharpUO` (analytical), `VERelPermTableUO` (tabulated), `VERelPermHysteresisUO` (Land/Killough) |
| **Capillary pressure** | `VEUpscaledCapPressure`, `VEFVCapPressure`, `VECapillaryPressureTableUO` |
| **Fluid properties** | `VEFluidProperties`, `VEFVFluidProperties`, `ConstantFluidProperties`, `SimpleBrineFluidProperties` |
| **Open boundary (FV)** | `VEFVAdvectiveOutflowBC` |
| **Trapping and dissolution** | `VERelPermHysteresisUO`, `VEDissolution`, `VEDissolutionSink`, `VEFVDissolutionSink` |
| **Diagnostics** | `VEPlumeHeightAux`, `VEGravityNumberAux`, `VESaturationMaxAux`, `VETrappedSaturationAux`, `VEDissolvedCO2Aux` |
| **Postprocessors** | `VEMobileCO2MassPostprocessor`, `VETrappedCO2MassPostprocessor`, `VEDissolvedCO2MassPostprocessor`, `VEPlumeFootprintPostprocessor` |

## Building

```bash
# Activate the MOOSE conda environment first
source ~/mambaforge3/etc/profile.d/conda.sh && conda activate moose

make -j8               # build wombat-opt
./run_tests -j8        # run the full test suite
```

## References

- Nordbotten & Celia, *Geological Storage of CO2: Modeling Approaches for Large-Scale Simulation* — the canonical VE reference
- Nordbotten & Dahle (2011) — capillary-fringe VE formulation
- Pau et al.; Hesse et al. — convective dissolution correlations
- MRST co2lab — reference VE implementation used for cross-code comparison
- [MOOSE PorousFlow documentation](https://mooseframework.inl.gov/modules/porous_flow/)
