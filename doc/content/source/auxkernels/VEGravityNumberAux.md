# VEGravityNumberAux

!syntax description /AuxKernels/VEGravityNumberAux

## Overview

`VEGravityNumberAux` computes the gravity number $\Gamma$, a dimensionless
VE-validity indicator that quantifies the ratio of vertical gravity segregation to
lateral viscous flow:

$$
\Gamma = \frac{k_v \, \Delta\rho \, g \, H^2}{\mu_n \, \bar{\phi} \, Q \, L}
$$

where:

| Symbol | Description | Source |
|--------|-------------|--------|
| $k_v$ | Vertical permeability (m$^2$) | `k_v` (coupled value or AuxVariable) |
| $\Delta\rho = \rho_w - \rho_n$ | Density difference (kg/m$^3$) | parameter |
| $g$ | Gravitational acceleration (m/s$^2$) | `gravity` parameter |
| $H$ | Formation thickness (m) | `ve_H` from `VEGeometry` |
| $\mu_n$ | CO2 viscosity (Pa·s) | parameter |
| $\bar{\phi}$ | Depth-averaged porosity (-) | `ve_phi_bar` from `VEPorosity` |
| $Q$ | Characteristic injection rate (m$^3$/s) | parameter |
| $L$ | Characteristic lateral length (m) | parameter |

Large $\Gamma$ (>> 1) means fast vertical segregation — the VE assumption is valid.
Small $\Gamma$ (< ~10) flags regions where vertical equilibrium may not hold; those
elements should be treated with caution and reported to stakeholders.

$H$ and $\bar{\phi}$ vary spatially, so $\Gamma$ is computed per quadrature point and
averaged to the element.

## Example Input File Syntax

```text
[AuxVariables]
  [k_v] ... []    # vertical permeability [m^2] (constant or from Exodus)
  [gamma]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [gravity_number]
    type = VEGravityNumberAux
    variable    = gamma
    k_v         = k_v
    delta_rho   = 300.0   # kg/m^3
    mu_n        = 5.0e-5  # Pa.s
    Q           = 1.0     # m^3/s (characteristic injection rate)
    L           = 5000.0  # m (characteristic length)
    gravity     = 9.81    # m/s^2
    execute_on  = TIMESTEP_END
  []
[]
```

!syntax parameters /AuxKernels/VEGravityNumberAux

!syntax inputs /AuxKernels/VEGravityNumberAux

!syntax children /AuxKernels/VEGravityNumberAux
