# VERelPermHysteresisUO

!syntax description /UserObjects/VERelPermHysteresisUO

## Overview

`VERelPermHysteresisUO` implements hysteretic upscaled relative permeability with
residual (capillary) trapping via the Land–Killough scanning-curve model.

In effective CO2 saturation $s = S_n / (1 - S_{wr})$, with Land coefficient
$C = (1 - S_{wr}) / S_{gr,max} - 1$:

**Drainage** ($S_n \geq S_{n,max}$):

\begin{equation}
k_{r,n} = k_{r,n}^{max} \cdot s
\end{equation}

**Imbibition** ($S_n < S_{n,max}$):

\begin{equation}
s_{gr} = \frac{s_{max}}{1 + C \, s_{max}}, \qquad
k_{r,n} = k_{r,n}^{max} \cdot s_{max} \cdot \frac{s - s_{gr}}{s_{max} - s_{gr}}
\end{equation}

The imbibition curve is continuous with the drainage line at the turning point and
vanishes at $s = s_{gr}$, immobilising the residual CO2.  A monotonically increasing
saturation history (pure drainage) never enters the imbibition branch, so this UO
reduces **exactly** to `VERelPermSharpUO` during injection.

**Brine** (wetting phase): $k_{r,w} = k_{r,w}^{max} \cdot (1 - s)$ on both branches
(trapped CO2 still occupies pore space and blocks brine).

The trapped saturation is:

\begin{equation}
S_{n,trap} = (1 - S_{wr}) \cdot s_{gr} \cdot \frac{s_{max} - s}{s_{max} - s_{gr}}
\end{equation}

reported by `trappedSaturation()` and written by `VETrappedSaturationAux`.

The turning point `sat_n_max` is supplied by `VESaturationMaxAux` (a stateful AuxVariable
frozen within each solve) and enters all history-aware calls as a plain `Real` — the AD
seed remains purely $d(k_r)/d(S_n)$.

## Example Input File Syntax

```text
[UserObjects]
  [relperm_uo]
    type = VERelPermHysteresisUO
    S_wr     = 0.2
    krn_max  = 0.8
    krw_max  = 1.0
    S_gr_max = 0.25   # max trapped CO2 saturation after full drainage
  []
[]
```

Wire with `VESaturationMaxAux` and `VETrappedSaturationAux` as described in the
hysteresis example tests.

!syntax parameters /UserObjects/VERelPermHysteresisUO

!syntax inputs /UserObjects/VERelPermHysteresisUO

!syntax children /UserObjects/VERelPermHysteresisUO
