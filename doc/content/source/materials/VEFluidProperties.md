# VEFluidProperties

!syntax description /Materials/VEFluidProperties

## Overview

`VEFluidProperties` evaluates the phase density and viscosity for both fluid phases
(CO2 and brine) from two `SinglePhaseFluidProperties` UserObjects and declares the
results as `ADMaterialProperty<std::vector<Real>>`:

- `ve_density[0]` / `ve_density[1]` — CO2 and brine densities (kg/m$^3$)
- `ve_viscosity[0]` / `ve_viscosity[1]` — CO2 and brine dynamic viscosities (Pa·s)

Phase 0 is non-wetting (CO2) and phase 1 is wetting (brine), matching the convention
used by the relperm and capillary-pressure objects.

The "plug and play" EOS pattern mirrors PorousFlow: supply `ConstantFluidProperties`
for analytical verification cases, `CO2FluidProperties` + `SimpleBrineFluidProperties`
for field runs using real EOS — without touching any kernel or material.

Temperature is a constant scalar parameter (isothermal assumption); it is not a
solution variable.

`initQpStatefulProperties()` seeds the old-value storage read by
`VEMassTimeDerivative` and `VEFVMassTimeDerivative`, so the storage term is correct
at $t = 0$.

### EOS reference depth

The `eos_reference_depth` parameter selects where the EOS is evaluated:

- `top_surface` (default) — evaluate at the top-surface pressure `pp_top`.  Reproduces
  all existing verification gold files bit-for-bit.
- `interface` — evaluate at the CO2–brine contact $z_T + h$, using the hydrostatic
  pressure $p = p_{top} + \rho_n(p_{top}) \cdot |g| \cdot h$ with the
  sharp-interface thickness $h = S_n H / (1 - S_{wr})$.  Provides a more physically
  consistent EOS pressure for buoyancy-driven CO2 columns.  Requires `sat_n`,
  `z_top`, `z_bottom`, `S_wr`, and `gravity` to be specified.

!alert note title=Required in both FE and FV inputs
`VEFVMassTimeDerivative` reads the stateful `ve_density` old-value from
`VEFluidProperties` even in FV models.  The FV flux kernel (`VEFVAdvectiveFlux`)
additionally needs `VEFVFluidProperties` for face-correct density evaluation.

## Example Input File Syntax

Constant properties (verification):

```text
[FluidProperties]
  [co2_fp]   type = ConstantFluidProperties   density = 700.0   viscosity = 1.0e-3 []
  [brine_fp] type = ConstantFluidProperties   density = 1000.0  viscosity = 1.0e-3 []
[]
[Materials]
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w  = brine_fp
    pp_top = pp_top
    temperature = 330.0     # K
  []
[]
```

Real EOS with interface reference depth:

```text
[FluidProperties]
  [co2_fp]   type = CO2FluidProperties []
  [brine_fp] type = SimpleBrineFluidProperties  xnacl = 0.05 []
[]
[Materials]
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w  = brine_fp
    pp_top = pp_top
    temperature = 330.0
    eos_reference_depth = interface
    sat_n    = sat_n
    z_top    = z_top
    z_bottom = z_bottom
    S_wr     = 0.2
  []
[]
```

!syntax parameters /Materials/VEFluidProperties

!syntax inputs /Materials/VEFluidProperties

!syntax children /Materials/VEFluidProperties
