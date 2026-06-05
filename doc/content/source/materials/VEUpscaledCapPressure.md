# VEUpscaledCapPressure

!syntax description /Materials/VEUpscaledCapPressure

## Overview

`VEUpscaledCapPressure` computes the macroscopic (column-averaged) capillary pressure
and its partial derivatives, consumed by `VEAdvectiveFluxS` when `capillary != none`.

The upscaled capillary pressure is:

\begin{equation}
P_c^{up} = (\rho_w - \rho_n) \, |g| \, h + P_{c,entry}
\end{equation}

where $h$ is the plume thickness from `VEPlumeReconstruction` (`ve_h`) and
$P_{c,entry}$ is an optional constant entry pressure (default 0).

`VEAdvectiveFluxS` assembles the full chain-rule gradient:

\begin{equation}
\nabla P_c^{up} = \underbrace{\frac{\partial P_c^{up}}{\partial S_n}}_{\texttt{ve\_dPcup\_dsatn}} \nabla S_n
                 + \underbrace{\frac{\partial P_c^{up}}{\partial H}}_{\texttt{ve\_dPcup\_dH}} \nabla H
\end{equation}

with, writing $S_n(h) = 1 - S_w(\Delta\rho \, g \, h)$ for the CO2 saturation at the
plume top:

|Mode|$\partial h/\partial S_n$|$\partial h/\partial H$|
|----|---|---|
|`sharp_interface`|$H/(1-S_{wr})$|$S_n/(1-S_{wr})$|
|`capillary_fringe`|$H/S_n(h)$|$S_n(h_0)/S_n(h)$ (via the UO)|

Both coefficient properties are `ADMaterialProperty<Real>`.  In sharp-interface mode
the AD derivative is zero (constant density → zero sat_n Jacobian on the $\nabla S_n$
term, matching the original behaviour).  In capillary-fringe mode the AD chain from
`ve_h` makes the $\nabla S_n$ Jacobian exact.

### Declares

- `ve_pc_up` — upscaled capillary pressure (Pa) (diagnostic)
- `ve_dPcup_dsatn` — $\partial P_c^{up}/\partial S_n$ (Pa)
- `ve_dPcup_dH` — $\partial P_c^{up}/\partial H$ (Pa/m)

### Reads

- `ve_h` from `VEPlumeReconstruction`
- `ve_density` from `VEFluidProperties`
- `ve_H` from `VEGeometry`

## Example Input File Syntax

```text
[Materials]
  [cap_pressure]
    type = VEUpscaledCapPressure
    sat_n = sat_n
    S_wr  = 0.2
    mode  = sharp_interface   # or capillary_fringe
    gravity = '0 0 -9.81'
  []
[]
```

!syntax parameters /Materials/VEUpscaledCapPressure

!syntax inputs /Materials/VEUpscaledCapPressure

!syntax children /Materials/VEUpscaledCapPressure
