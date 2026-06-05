# VERelPermTableUO

!syntax description /UserObjects/VERelPermTableUO

## Overview

`VERelPermTableUO` provides piecewise-linear tabulated upscaled relative permeabilities
$k_{r,n}(S_n)$ and $k_{r,w}(S_n)$ from tables produced by the fine-scale upscaling
workflow.  It is the field-case counterpart of the analytical `VERelPermSharpUO`.

Both curves share the same saturation axis `sat_n_points`.  Values outside the tabulated
range are clamped to the endpoint values (zero derivative), the correct behaviour for a
bounded upscaled curve.

Being a `VERelativePermeability` subclass, `VERelPermTableUO` is consumed identically
by `VERelPerm` (FE) and `VEFVRelPerm` (FV) — swapping the UserObject changes the kr
model without touching any kernel.

## Example Input File Syntax

```text
[UserObjects]
  [relperm_uo]
    type = VERelPermTableUO
    sat_n_points = '0.0  0.1  0.3  0.6  0.8  1.0'
    krn_points   = '0.0  0.01 0.08 0.35 0.65 0.80'
    krw_points   = '1.0  0.85 0.60 0.25 0.05 0.00'
  []
[]
[Materials]
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
[]
```

!syntax parameters /UserObjects/VERelPermTableUO

!syntax inputs /UserObjects/VERelPermTableUO

!syntax children /UserObjects/VERelPermTableUO
