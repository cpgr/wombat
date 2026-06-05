# VEDissolution

!syntax description /Materials/VEDissolution

## Overview

`VEDissolution` computes the convective-dissolution rate of free-phase CO2 into the
underlying brine as an areal mass flux `ve_dissolution_rate` (kg/m$^2$/s):

\begin{equation}
\dot{m}_{diss} = q_0 \cdot \text{gate}(S_n) \cdot \text{capacity}(c_{diss})
\end{equation}

**$q_0$** (`dissolution_flux`) is the constant convective dissolution flux.  In the
density-driven convective regime (Pau et al. 2010; Hesse et al.) the flux per unit
plume footprint becomes approximately constant ($q_0 \sim \alpha k \Delta\rho_c g
C_{sat} / \mu$).  For v1 it is supplied directly; deriving it from the correlation is
a later refinement.

**gate**$(S_n) = \min(1,\, S_n / s_{ref})$ rate-limits dissolution as the plume thins:
full rate when $S_n \geq s_{ref}$, tapering linearly to zero as $S_n \to 0$ (no plume
→ no dissolution).

**capacity** $= \max(0,\ 1 - c_{diss} / c_{cap})$ optionally halts dissolution once the
column has absorbed its CO2 capacity.  Omit `dissolved_co2` / `c_cap` to disable the
cap (unbounded dissolution).

This is a **sink-only** formulation (Option A): dissolved CO2 is tracked as an immobile
areal inventory `c_diss` (written by `VEDissolvedCO2Aux`) and is not transported.  The
freed pore volume is supplied by brine inflow through the coupled pressure solve (open
system).

The rate is consumed by `VEDissolutionSink` (FE) and `VEFVDissolutionSink` (FV) on
the CO2 saturation equation, and is also read by `VEDissolvedCO2Aux` to advance the
dissolved inventory.

## Example Input File Syntax

```text
[AuxVariables]
  [c_diss]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  []
[]
[Materials]
  [dissolution]
    type = VEDissolution
    dissolution_flux = 1.4    # kg/m^2/s
    s_ref = 0.05
    # Optional capacity cap:
    # dissolved_co2 = c_diss
    # c_cap = 500.0           # kg/m^2
  []
[]
```

!syntax parameters /Materials/VEDissolution

!syntax inputs /Materials/VEDissolution

!syntax children /Materials/VEDissolution
