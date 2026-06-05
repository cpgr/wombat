# VEFVFluidProperties

!syntax description /Materials/VEFVFluidProperties

## Overview

`VEFVFluidProperties` is a `FunctorMaterial` that publishes phase density and viscosity
as the functors

- `ve_density_n` / `ve_density_w` (kg/m$^3$)
- `ve_viscosity_n` / `ve_viscosity_w` (Pa·s)

consumed by `VEFVAdvectiveFlux` for face-correct EOS evaluation, delegating to the same
`SinglePhaseFluidProperties` UserObjects as `VEFluidProperties`.

Because `VEFVAdvectiveFlux` evaluates these functors at a central-difference
(face-averaged) spatial argument, the density and viscosity depend on both adjacent cell
pressures and carry AD derivatives with respect to `pp_top` on both sides of the face.
For constant fluid properties the face average reduces to the constant, so this is exact.

The `eos_reference_depth` parameter mirrors `VEFluidProperties`:

- `top_surface` (default) — evaluate at the `pp_top` functor value.
- `interface` — evaluate at $p = p_{top} + \rho_n(p_{top}) \cdot |g| \cdot h$,
  where $h = S_n (z_T - z_B) / (1 - S_{wr})$ is the sharp-interface plume thickness
  evaluated at the same face argument.  Requires `sat_n`, `z_top`, `z_bottom`, and
  `S_wr` to be specified.

!alert note title=Used alongside VEFluidProperties
`VEFVFluidProperties` feeds only the FV flux kernel.  The elemental `VEFluidProperties`
material is still required in every FV input because `VEFVMassTimeDerivative` reads
the stateful `ve_density` old-value from it.

## Example Input File Syntax

```text
[FunctorMaterials]
  [fv_fluid_props]
    type = VEFVFluidProperties
    fp_nw  = co2_fp
    fp_w   = brine_fp
    pp_top = pp_top
    temperature = 330.0
  []
[]
```

!syntax parameters /Materials/VEFVFluidProperties

!syntax inputs /Materials/VEFVFluidProperties

!syntax children /Materials/VEFVFluidProperties
