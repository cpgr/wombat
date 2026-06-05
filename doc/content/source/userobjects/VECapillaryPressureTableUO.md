# VECapillaryPressureTableUO

!syntax description /UserObjects/VECapillaryPressureTableUO

## Overview

`VECapillaryPressureTableUO` is a `PorousFlowCapillaryPressure` UserObject backed by a
piecewise-linear $(P_c, S_w)$ lookup table.  It is the production companion to
`VEPlumeReconstruction` (capillary_fringe mode) and `VEUpscaledCapPressure` /
`VEFVCapPressure` for real field models.

### Motivation

Parametric Pc UOs (`PorousFlowCapillaryPressureVG`, `BC`) assume a single rock type
throughout the entire vertical column.  Real storage formations have multiple stacked
facies with different Pc curves.  The fine-scale upscaling workflow already integrates
through all vertical facies and produces a column-averaged $S_w(P_c)$ table as output.
`VECapillaryPressureTableUO` feeds this table directly into `VEPlumeReconstruction`
without forcing the data through a parametric fit that may not represent the multi-facies
column accurately.

Different material blocks can reference different `VECapillaryPressureTableUO` instances
to handle lateral variation in the facies sequence across the domain.

### Derivatives

First derivatives (`dEffectiveSaturation`, `dCapillaryPressureCurve`) are the
piecewise-constant slopes of the linear segments — exact for the interpolant.  Second
derivatives are returned as 0.0 (acceptable for Newton convergence in practice with a
piecewise-linear interpolant).

Set `log_extension = false` (the base-class logarithmic low-saturation extension is not
meaningful for a table with a finite lower $S_w$ bound).

### Input

- `pc_points` — capillary pressure values (Pa), strictly increasing from 0
- `sw_points` — water saturation values (-), strictly decreasing, one-to-one with `pc_points`

## Example Input File Syntax

```text
[UserObjects]
  [pc_uo]
    type = VECapillaryPressureTableUO
    sat_lr   = 0.2
    pc_points = '0      5000   20000  60000  150000'
    sw_points = '1.0    0.85   0.65   0.40   0.20'
  []
[]
```

Different lateral zones:

```text
[UserObjects]
  [pc_zone_A]
    type = VECapillaryPressureTableUO
    sat_lr   = 0.15
    pc_points = '0  8000  30000'
    sw_points = '1.0  0.75  0.20'
  []
  [pc_zone_B]
    type = VECapillaryPressureTableUO
    sat_lr   = 0.25
    pc_points = '0  3000  12000'
    sw_points = '1.0  0.80  0.30'
  []
[]
```

!syntax parameters /UserObjects/VECapillaryPressureTableUO

!syntax inputs /UserObjects/VECapillaryPressureTableUO

!syntax children /UserObjects/VECapillaryPressureTableUO
