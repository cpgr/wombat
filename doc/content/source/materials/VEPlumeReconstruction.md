# VEPlumeReconstruction

!syntax description /Materials/VEPlumeReconstruction

## Overview

`VEPlumeReconstruction` reconstructs the local CO2 plume thickness `ve_h` (m) from
the depth-averaged saturation `sat_n`.  Two modes are supported:

### sharp_interface (default)

Nordbotten–Celia analytical formula:

$$
h = \frac{S_n \, H}{1 - S_{wr}}
$$

Closed form, cheap, and used in all standard verification cases.

### capillary_fringe (Nordbotten & Dahle 2011)

Newton inversion of the depth-averaged saturation integral:

$$
\bar{S}_n(h) = \frac{1}{H} \int_0^h \left[1 - S_w\!\left(P_c(\zeta)\right)\right]
\mathrm{d}\zeta
$$

where $P_c(\zeta) = (\rho_w - \rho_{CO_2}) \, g \, \zeta$ is the buoyancy-driven
capillary pressure at height $\zeta$ above the CO2–brine contact, and $S_w(P_c)$ is
provided by a `PorousFlowCapillaryPressure` UserObject (`pc_uo`).

The Newton derivative $dF/dh = S_n(h)/H$ follows from Leibniz rule.  The integral
uses the trapezoidal rule (64 intervals).  The sharp-interface estimate is the warm
start.

The density difference $\Delta\rho = \rho_w - \rho_{CO_2}$ is read per quadrature
point from `ve_density` (from `VEFluidProperties`), so spatial density variation is
handled correctly.

### AD seeding

`ve_h` is declared as `ADMaterialProperty<Real>` so that $\partial h / \partial S_n$
propagates into the Jacobian of `VEUpscaledCapPressure` and `VEFVCapPressure` without
hand-coded derivatives.  In sharp-interface mode the seed is $H/(1-S_{wr})$; in
capillary-fringe mode it is $H/S_n(h)$ (inverse-function-theorem result on the
saturation integral).

The result is clamped to $[0, H]$ in both modes.

!alert note title=Pc UserObject must represent the upscaled column curve
For real multi-facies formations the `PorousFlowCapillaryPressure` UO fed to
`VEPlumeReconstruction` must be the column-averaged $S_w(P_c)$ table produced by the
upscaling workflow — not a parametric fit to a single rock type.  Use
`VECapillaryPressureTableUO` for this purpose.  Parametric UOs (`VG`, `BC`) are
correct for verification only.

## Example Input File Syntax

Sharp-interface (default):

```text
[Materials]
  [plume_recon]
    type = VEPlumeReconstruction
    sat_n = sat_n
    S_wr  = 0.2
  []
[]
```

Capillary-fringe with tabulated Pc:

```text
[UserObjects]
  [pc_uo]
    type = VECapillaryPressureTableUO
    sat_lr   = 0.2
    pc_points = '0  5000  20000  60000'
    sw_points = '1.0  0.8  0.5  0.2'
  []
[]
[Materials]
  [plume_recon]
    type = VEPlumeReconstruction
    sat_n = sat_n
    S_wr  = 0.2
    mode  = capillary_fringe
    pc_uo = pc_uo
  []
[]
```

!syntax parameters /Materials/VEPlumeReconstruction

!syntax inputs /Materials/VEPlumeReconstruction

!syntax children /Materials/VEPlumeReconstruction
