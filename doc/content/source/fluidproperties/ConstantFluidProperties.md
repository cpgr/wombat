# ConstantFluidProperties

!syntax description /FluidProperties/ConstantFluidProperties

## Overview

`ConstantFluidProperties` is a `SinglePhaseFluidProperties` UserObject that returns a
constant density and viscosity independent of pressure and temperature (zero
derivatives).

Its purpose is to provide a verification-case EOS that plugs into `VEFluidProperties`
(and `VEFVFluidProperties`) using the same "plug and play" interface as the physical
EOS classes (`CO2FluidProperties`, `SimpleBrineFluidProperties`).  This means:

- All analytical verification golds (Buckley–Leverett, Nordbotten–Celia, McWhorter–
  Sunada, …) are run with `ConstantFluidProperties` and reproduce the pre-AD numbers
  bit-for-bit.
- Switching to a full pressure/temperature-dependent EOS for field runs is a UserObject
  swap in the `[FluidProperties]` block — no kernel or material changes.

Only `rho_from_p_T` and `mu_from_p_T` are implemented; all other properties inherit
the base-class "not implemented" error.

## Example Input File Syntax

```text
[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density   = 700.0    # kg/m^3
    viscosity = 5.0e-5   # Pa.s
  []
  [brine_fp]
    type = ConstantFluidProperties
    density   = 1020.0
    viscosity = 1.0e-3
  []
[]
[Materials]
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w  = brine_fp
    pp_top = pp_top
    temperature = 330.0
  []
[]
```

!syntax parameters /FluidProperties/ConstantFluidProperties

!syntax inputs /FluidProperties/ConstantFluidProperties

!syntax children /FluidProperties/ConstantFluidProperties
