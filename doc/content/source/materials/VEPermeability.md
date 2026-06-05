# VEPermeability

!syntax description /Materials/VEPermeability

## Overview

`VEPermeability` provides the upscaled permeability tensor `ve_K_up` (m$^3$)
(units: m$^2$ × m because the formation thickness $H$ is already folded in by the
upscaling workflow) consumed by `VEAdvectiveFluxP`, `VEAdvectiveFluxS`, and
`VEFVAdvectiveFlux`.

The symmetric 2×2 horizontal tensor is assembled from up to three independent
components:

| Parameter | Description | Required |
|-----------|-------------|---------|
| `K_up_xx` | $K_{xx}$ component (m$^3$) | yes |
| `K_up_yy` | $K_{yy}$ component (m$^3$) | yes |
| `K_up_xy` | Off-diagonal $K_{xy}$ (m$^3$) | no (default 0) |

Each component accepts either a constant scalar or a spatially varying
`AuxVariable` from the Exodus mesh.  The $K_{zz}$ component is zero (no vertical
flow in a VE model).

`ve_K_up` is a plain `Real` tensor (not AD): permeability is time-invariant.

!alert note title=FV heterogeneous-K face averaging
`VEFVAdvectiveFlux` uses a distance-weighted harmonic mean of the element-side and
neighbour-side `ve_K_up` values for interior faces, giving the correct TPFA
two-point transmissibility across permeability contrasts.

## Example Input File Syntax

Isotropic homogeneous case (verification):

```text
[Materials]
  [permeability]
    type = VEPermeability
    K_up_xx = 1.0e-10
    K_up_yy = 1.0e-10
  []
[]
```

Anisotropic spatially varying case (field model):

```text
[AuxVariables]
  [K_up_xx] [K_up_yy] [K_up_xy] ...
[]
[Materials]
  [permeability]
    type = VEPermeability
    K_up_xx = K_up_xx
    K_up_yy = K_up_yy
    K_up_xy = K_up_xy
  []
[]
```

!syntax parameters /Materials/VEPermeability

!syntax inputs /Materials/VEPermeability

!syntax children /Materials/VEPermeability
