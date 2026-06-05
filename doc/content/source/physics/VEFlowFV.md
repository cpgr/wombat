# VEFlowFV

!syntax description /Physics/VEFlow/FiniteVolume/VEFlowFV

## Overview

`VEFlowFV` is the cell-centered finite-volume Physics action for the
vertical-equilibrium two-phase CO2–brine flow equations.  It collapses the FV input
boilerplate by automatically creating:

- Two `MooseVariableFVReal` primary variables: `pp_top` and `sat_n`.
- Two FV geometry variables (`z_top`, `z_bottom`, as `MooseVariableFVReal`) unless
  `define_geometry_variables = false`.
- Four FV flow kernels: `VEFVMassTimeDerivative` and `VEFVAdvectiveFlux` on `sat_n`,
  and `VEFVMassTimeDerivative` and `VEFVAdvectiveFlux` on `pp_top`.
- Elemental materials: `VEPorosity`, `VEPermeability`, `VESaturation`, `VEFluidProperties`.
- Functor materials: `VEFVRelPerm`, `VEFVFluidProperties`.
- Optionally, `VEFVCapPressure` (when `capillary` is set).
- Optionally, the dissolution chain when `dissolution = true`.

Syntax: `[Physics/VEFlow/FiniteVolume/<name>]`.

### No VEGeometry

The FV flux and storage kernels build $H$ and $\nabla z_T$ inline from the FV-variable
cell values rather than from `VEGeometry`'s `ve_H`, because FE-coupled materials are
not reinitialised on FV faces.  `VEGeometry` is therefore **not** added by `VEFlowFV`.

### Petrophysics and geometry variable types

For FV models, `phi_bar`, `K_up_xx`, etc. must be `MooseVariableFVReal` when supplied
as coupled field variables (not constants).  A regular FE AuxVariable is not
reinitialised on FV faces and would cause the harmonic permeability average to collapse
to zero, silently killing the advective flux.  `VEFlowFV` checks this and emits an
actionable error if a mismatch is detected.

Similarly, `z_top` / `z_bottom` must be `MooseVariableFVReal`; a mismatched type is
caught by MOOSE's variable-type consistency check.

### Preconditioner

The coupled FV solve requires `SMP full = true` in `[Preconditioning]`.

## Example Input File Syntax

Minimal FV setup (verification — flat formation, constant properties):

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
[Physics/VEFlow/FiniteVolume]
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
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
```

Field model with spatially varying petrophysics and capillary diffusion:

```text
[Physics/VEFlow/FiniteVolume]
  [ve]
    z_top    = z_top
    z_bottom = z_bottom
    phi_bar  = phi_bar      # MooseVariableFVReal from Exodus
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

!syntax parameters /Physics/VEFlow/FiniteVolume/VEFlowFV

!syntax inputs /Physics/VEFlow/FiniteVolume/VEFlowFV

!syntax children /Physics/VEFlow/FiniteVolume/VEFlowFV
