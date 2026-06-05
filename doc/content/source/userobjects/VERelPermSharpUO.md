# VERelPermSharpUO

!syntax description /UserObjects/VERelPermSharpUO

## Overview

`VERelPermSharpUO` implements the Nordbotten–Celia sharp-interface upscaled relative
permeability curves:

\begin{equation}
s_{eff} = \text{clamp}\!\left(\frac{S_n}{1 - S_{wr}},\, 0,\, 1\right)
\end{equation}

\begin{equation}
k_{r,n}(S_n) = k_{r,n}^{max} \cdot s_{eff} \qquad (\text{CO}_2,\; \text{phase 0})
\end{equation}

\begin{equation}
k_{r,w}(S_n) = k_{r,w}^{max} \cdot (1 - s_{eff}) \qquad (\text{brine},\; \text{phase 1})
\end{equation}

These are linear in effective saturation — the exact result for a vertically uniform
column with a sharp CO2–brine contact.  This is the standard analytical verification
curve.  For tabulated kr from the fine-scale upscaling workflow, use `VERelPermTableUO`.
For hysteretic trapping, use `VERelPermHysteresisUO`.

All three UserObjects are `VERelativePermeability` subclasses and are consumed
identically by `VERelPerm` (FE) and `VEFVRelPerm` (FV).

## Example Input File Syntax

```text
[UserObjects]
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr    = 0.2     # residual water saturation in CO2 zone
    krn_max = 0.8     # max CO2 relative permeability
    krw_max = 1.0     # max brine relative permeability
  []
[]
[Materials]
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
[]
```

!syntax parameters /UserObjects/VERelPermSharpUO

!syntax inputs /UserObjects/VERelPermSharpUO

!syntax children /UserObjects/VERelPermSharpUO
