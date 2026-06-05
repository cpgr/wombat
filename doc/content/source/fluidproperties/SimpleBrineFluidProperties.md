# SimpleBrineFluidProperties

!syntax description /FluidProperties/SimpleBrineFluidProperties

## Overview

`SimpleBrineFluidProperties` exposes the FluidProperties-module `BrineFluidProperties`
through the two-variable $(p, T)$ `SinglePhaseFluidProperties` interface consumed by
`VEFluidProperties`, at a fixed (constant) salt (NaCl) mass fraction.

`BrineFluidProperties` is a three-variable EOS $(p, T, X_{NaCl})$; it cannot be
handed directly to `VEFluidProperties`, which carries only pressure and temperature.
`SimpleBrineFluidProperties` bridges the gap by holding a constant `xnacl` and
forwarding $\rho$ and $\mu$ evaluations to the underlying three-variable EOS at the
fixed salinity.  Derivatives with respect to salinity are discarded; only $\partial/\partial p$
and $\partial/\partial T$ are exposed — exactly what a constant-salinity VE model needs.

Variable salinity is intentionally out of scope; if ever needed, brine salinity would
have to become a solution/coupled field.

Only `rho_from_p_T` and `mu_from_p_T` are implemented.

## Example Input File Syntax

```text
[FluidProperties]
  [brine_fp]
    type = SimpleBrineFluidProperties
    xnacl = 0.05    # 5 wt% NaCl
  []
  [co2_fp]
    type = CO2FluidProperties
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

!syntax parameters /FluidProperties/SimpleBrineFluidProperties

!syntax inputs /FluidProperties/SimpleBrineFluidProperties

!syntax children /FluidProperties/SimpleBrineFluidProperties
