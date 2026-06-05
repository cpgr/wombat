# VEFVCapPressure

!syntax description /Materials/VEFVCapPressure

## Overview

`VEFVCapPressure` is the FV counterpart of `VEUpscaledCapPressure`: a `FunctorMaterial`
that publishes the chain-rule gradient coefficients of the upscaled capillary pressure
as functors consumed by `VEFVAdvectiveFlux` when `capillary != none`.

The gradient of the upscaled capillary pressure is decomposed as:

$$
\nabla P_c^{up} \cdot \hat{n}
   = \underbrace{\frac{\partial P_c^{up}}{\partial S_n}}_{\texttt{ve\_dPcup\_dsatn}}
     \nabla S_n \cdot \hat{n}
   + \underbrace{\frac{\partial P_c^{up}}{\partial H}}_{\texttt{ve\_dPcup\_dH}}
     \nabla H \cdot \hat{n}
$$

The kernel uses the `sat_n` variable gradient (boundary-aware: picks up Dirichlet BCs)
rather than a material-functor Green–Gauss gradient, so `ve_pc_up` is published for
diagnostics only.

The partial-derivative coefficients are computed from the chain-rule partials of the
plume thickness $h(S_n, H)$, identical to `VEUpscaledCapPressure`.  The density
difference $\Delta\rho = \rho_w - \rho_n$ is read per evaluation from the density
functors, so spatially varying EOS is handled correctly.

### Declares

- `ve_pc_up` — upscaled capillary pressure (Pa) (diagnostic functor)
- `ve_dPcup_dsatn` — $\partial P_c^{up}/\partial S_n$ (Pa) (functor)
- `ve_dPcup_dH` — $\partial P_c^{up}/\partial H$ (Pa/m) (functor)

## Example Input File Syntax

```text
[FunctorMaterials]
  [cap_pressure]
    type = VEFVCapPressure
    sat_n      = sat_n
    z_top      = z_top
    z_bottom   = z_bottom
    mode       = sharp_interface
    S_wr       = 0.2
    gravity    = '0 0 -9.81'
  []
[]
```

!syntax parameters /Materials/VEFVCapPressure

!syntax inputs /Materials/VEFVCapPressure

!syntax children /Materials/VEFVCapPressure
