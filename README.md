# wombat

### WOrkflow for Modelling Buoyant Aquifer Transport

**A vertical-equilibrium (VE) solver for geological CO2 storage, built on the
[MOOSE](https://mooseframework.inl.gov) framework and its PorousFlow module.**

`wombat` depth-integrates the multiphase porous-flow equations under the assumption
of vertical pressure equilibrium, reducing a 3D CO2-storage problem to a 2D problem
on the top surface of the storage formation. This is appropriate for CO2 plumes that
are thin relative to their lateral extent — the common case at field scale and for
long-timescale, post-injection migration studies — and is dramatically cheaper than a
full 3D simulation while retaining the dominant physics (buoyant up-dip migration,
structural trapping, residual trapping, dissolution).

It is a research application under active development.

---

## Why vertical equilibrium?

At field scale a CO2 plume is often metres-to-tens-of-metres thick but kilometres
wide. Resolving the thin vertical dimension in 3D is expensive and largely wasted: the
fluids segregate gravitationally almost instantly compared to the lateral migration
timescale, so the vertical pressure distribution is hydrostatic to good approximation.
Integrating the governing equations over the formation thickness `H(x,y)` collapses
the problem onto the 2D top surface, with the vertical structure (the CO2-brine
interface, the capillary fringe) reconstructed analytically or from upscaled tables.

Per fluid component `c`, depth-integrated mass conservation reads:

```
d/dt ( H * phi_bar * rho_c * S_c )  +  div(F_c)  -  q_c  =  0

F_c = - H * K_up * kr_c^up * rho_c / mu_c * ( grad(pp_top) + rho_c * g * grad(z_T) )
```

The `rho_c * g * grad(z_T)` term — the **sloped-top-surface buoyancy drive** — is the
dominant up-dip migration mechanism at field scale and is treated as a first-class
verification target throughout the code.

## Design

`wombat` is built **on top of** PorousFlow rather than as a standalone tool, so it
inherits battle-tested CO2-brine equations of state, fluid-property objects, the PETSc
solver stack, parallel partitioning, and Exodus I/O. VE-specific physics lives in
custom MOOSE kernels, materials, and user objects that wrap or extend PorousFlow
primitives — PorousFlow itself is not forked.

- **Primary variables:** pore pressure at the top surface `pp_top`, and depth-averaged
  non-wetting (CO2) saturation `sat_n`. The (P, S) formulation keeps saturation bounded
  in [0, 1] and plays well with PorousFlow's residual/Jacobian patterns.
- **Two solver paths:** a continuous finite-element (FE) path and a finite-volume (FV)
  path. The FV path (with optional TVD flux limiting) is the robust choice for the
  advection-dominated migration that VE problems present; the FE path is convenient for
  smooth/diffusive cases and some verification work.
- **Geometry from the mesh:** formation geometry (`H`, `z_T`, `z_B`, `grad(z_T)`) is
  supplied by the `VEGeometry` material (FE) or by coupled `z_top`/`z_bottom` variables
  (FV), sourced from the upscaled Exodus mesh. The VE "flavour" (`sharp_interface` vs.
  `capillary_fringe`) is selected by which relperm / capillary-pressure objects are wired
  in, not by a central configuration object.
- **Field-scale input:** petrophysical fields (`H`, `z_T`, `phi_bar`, `K_up`, upscaled
  rel-perm / capillary tables) are read directly from an Exodus mesh produced by an
  external Petrel-to-2D upscaling workflow.

## Capabilities

| Area | Objects |
| --- | --- |
| **Mass storage / flux (FE)** | `VEMassTimeDerivative`, `VEAdvectiveFluxS` (CO2), `VEAdvectiveFluxP` (brine) |
| **Mass storage / flux (FV)** | `VEFVMassTimeDerivative`, `VEFVAdvectiveFlux` (unified; `upwind` + `vanLeer`/`min_mod`/`sou`/`quick`/`venkatakrishnan` limiters) |
| **Geometry / petrophysics** | `VEGeometry` (H, z_T, grad z_T, grad H), `VEPorosity`, `VEPermeability` (anisotropic `K_up`) |
| **Plume reconstruction** | `VEPlumeReconstruction` (sharp-interface closed form; capillary-fringe Newton inversion, Nordbotten & Dahle 2011) |
| **Relative permeability** | swappable `VERelativePermeability` user objects — `VERelPermSharpUO` (analytical), `VERelPermTableUO` (tabulated/upscaled), `VERelPermHysteresisUO` (Land/Killough); FE/FV adapters `VERelPerm`/`VEFVRelPerm` |
| **Capillary pressure** | `VEUpscaledCapPressure` / `VEFVCapPressure` (two-pressure VE), `VECapillaryPressureTable` |
| **Fluids** | `VEFluidProperties` (density/viscosity from any `SinglePhaseFluidProperties` UO via `fp_nw`/`fp_w`) driven by `ConstantFluidProperties` (verification) or `CO2FluidProperties`/`BrineFluidProperties` (full P/T EOS) |
| **Wells & boundaries** | `ConstantPointSource` Dirac wells (FE + FV), `VEFVAdvectiveOutflowBC` (open up-dip outflow), standard Neumann/Dirichlet flux BCs |
| **Trapping & dissolution** | residual (hysteretic) trapping via `VESaturationMaxAux` + `VERelPermHysteresisUO`; convective dissolution via `VEDissolution` + `VEDissolutionSink`/`VEFVDissolutionSink` |
| **Diagnostics** | `VEPlumeHeightAux`, `VEGravityNumberAux` (VE-validity indicator), `VEMobileCO2MassPostprocessor`, `VETrappedCO2MassPostprocessor`, `VEDissolvedCO2MassPostprocessor`, `VEPlumeFootprintPostprocessor` |

## Verification

Every kernel and material has a corresponding test in `test/tests/`, and the core
physics is anchored by an analytical / cross-code **verification ladder**:

- **1D Buckley-Leverett with gravity** (`buckley_leverett/`) — advective flux kernel.
- **Nordbotten-Celia sloped aquifer** (`nordbotten_celia/`) — the sloped-surface gravity
  drive (analytical plume-nose velocity).
- **McWhorter-Sunada** (`mcwhorter_sunada/`) — counter-current capillary imbibition
  (sign and magnitude of the upscaled capillary term).
- **Capillary equilibrium** (`capillary_equilibrium/`) — the `grad(H)` chain-rule term.
- **Real-field input pipeline** (`exodus/`), **wells** (`dirac_source/`), **outflow**
  (`outflow/`), **hysteresis/trapping** (`hysteresis/`), **dissolution**
  (`dissolution/`), and more.
- **Cross-code benchmark** (`cross_code/`, in progress) — comparison against the MRST
  co2lab VE reference (see `doc/content/benchmarks/cross_code.md`).

## Building

`wombat` is a standard MOOSE application. It requires a built MOOSE with the
**PorousFlow** module enabled. MOOSE is expected as a sibling checkout (`../moose`) or
as a `moose` submodule.

```sh
# from the wombat root, with a MOOSE conda environment activated
make -j8                 # builds wombat-opt (and wombat-devel for debugging)
./run_tests -j8          # runs the regression + verification test suite
```

See the MOOSE [Getting Started](https://mooseframework.inl.gov/getting_started/)
guide for installing the environment and building the framework.

## Repository layout

```
src/ , include/      C++ source and headers for the VE objects
  kernels/ fvkernels/ fvbcs/ materials/ userobjects/ auxkernels/ postprocessors/
test/tests/          regression + verification inputs and gold files
doc/content/         design notes (incl. benchmarks/cross_code.md)
Makefile , run_tests standard MOOSE app build/test entry points
```

## Status and scope

This is a v1 research code focused on the CO2-brine system with constant fluid
properties, sharp-interface and capillary-fringe VE, residual trapping, and
constant-flux convective dissolution. Roadmap items include a pressure-dependent
`PorousFlowBrineCO2` EOS wrap, transported dissolved CO2, the full cross-code
benchmark, and coupling to TensorMechanics for caprock-integrity studies.
Standalone-style optimisations (front tracking, GPU kernels) are out of scope.

## References

- Nordbotten & Celia, *Geological Storage of CO2: Modeling Approaches for Large-Scale
  Simulation* — the canonical VE reference.
- Nordbotten & Dahle (2011) — capillary-fringe VE formulation.
- Pau et al.; Hesse et al. — convective dissolution correlations.
- MRST co2lab — reference VE implementation (cross-code comparison).
- [MOOSE PorousFlow documentation](https://mooseframework.inl.gov/modules/porous_flow/)
  — the framework being extended.

---
