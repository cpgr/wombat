# VEFlowFE

!syntax description /Physics/VEFlow/FiniteElement/VEFlowFE

## Overview

`VEFlowFE` is the continuous-Galerkin (Galerkin finite-element) Physics action for the
vertical-equilibrium two-phase CO2–brine flow equations.  It collapses the standard VE
input boilerplate by automatically creating:

- Two `LAGRANGE` primary variables: `pp_top` (pore pressure at top surface) and `sat_n`
  (depth-averaged CO2 saturation).
- Two geometry `AuxVariables` (`z_top`, `z_bottom`) unless `define_geometry_variables =
  false`.
- Four FE flow kernels: `VEMassTimeDerivative` and `VEAdvectiveFluxS` on `sat_n`, and
  `VEMassTimeDerivative` and `VEAdvectiveFluxP` on `pp_top`.
- The full FE material chain: `VEGeometry`, `VEPorosity`, `VEPermeability`, `VESaturation`,
  `VEFluidProperties`, `VERelPerm`.
- Optionally, capillary materials (`VEPlumeReconstruction`, `VEUpscaledCapPressure`) when
  `capillary` is set.
- Optionally, the dissolution chain (`VEDissolution` material + `VEDissolvedCO2Aux`)
  when `dissolution = true`.

Syntax: `[Physics/VEFlow/FiniteElement/<name>]`.

### Geometry variables

By default (`define_geometry_variables = true`) `VEFlowFE` declares `z_top` and
`z_bottom` as `LAGRANGE FIRST ORDER` AuxVariables and applies the
`define_geometry_variables = false` path when the user has already declared them in
`[AuxVariables]`.  In either case the variables must be `LAGRANGE FIRST ORDER` (a
`CONSTANT MONOMIAL` declaration would zero the buoyancy drive — see `VEGeometry`).

### Capillary option

Set `capillary` to `sharp_interface` or `capillary_fringe` to enable the capillary
diffusion term on the CO2 equation (`VEAdvectiveFluxS capillary=true`) and add the
supporting materials.  In `capillary_fringe` mode `pc_uo` is required.

## Example Input File Syntax

Minimal FE setup (verification — flat formation, constant properties):

```text
[FluidProperties]
  [co2_fp]   type = ConstantFluidProperties   density = 700.0   viscosity = 1.0e-3 []
  [brine_fp] type = ConstantFluidProperties   density = 1000.0  viscosity = 1.0e-3 []
[]
[UserObjects]
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.2   krn_max = 1.0   krw_max = 1.0
  []
[]
[Physics/VEFlow/FiniteElement]
  [ve]
    z_top    = z_top
    z_bottom = z_bottom
    phi_bar  = 0.2
    K_up_xx  = 1.0e-10
    K_up_yy  = 1.0e-10
    fp_nw    = co2_fp
    fp_w     = brine_fp
    relperm_uo = relperm_uo
  []
[]
```

With capillary diffusion (sharp-interface) and user-declared geometry variables:

```text
[Physics/VEFlow/FiniteElement]
  [ve]
    z_top    = z_top
    z_bottom = z_bottom
    phi_bar  = phi_bar
    K_up_xx  = K_up_xx
    K_up_yy  = K_up_yy
    fp_nw    = co2_fp
    fp_w     = brine_fp
    relperm_uo = relperm_uo
    capillary  = sharp_interface
    S_wr       = 0.2
    define_geometry_variables = false
  []
[]
```

!syntax parameters /Physics/VEFlow/FiniteElement/VEFlowFE

!syntax inputs /Physics/VEFlow/FiniteElement/VEFlowFE

!syntax children /Physics/VEFlow/FiniteElement/VEFlowFE
