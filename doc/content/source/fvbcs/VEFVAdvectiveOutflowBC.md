# VEFVAdvectiveOutflowBC

!syntax description /FVBCs/VEFVAdvectiveOutflowBC

## Overview

`VEFVAdvectiveOutflowBC` is an open / advective-outflow boundary condition for the
FV VE solver.  It lets a phase leave the domain at the Darcy rate rather than piling
up against the default no-flow wall.

Assign it to `variable = sat_n` (`fluid_phase = 0`) on the updip or downstream
boundary to allow migrating CO2 to exit; the brine equation is typically handled by a
Dirichlet `pp_top` on the same boundary.

The flux is evaluated using the same functors as `VEFVAdvectiveFlux`:

\begin{equation}
F_c \cdot \hat{n} = -H \, K_{nn} \, \frac{k_{r,c} \, \rho_c}{\mu_c}
  \left( \nabla p_{top} \cdot \hat{n} + \rho_c \, |g| \, \nabla z_T \cdot \hat{n} \right)
\end{equation}

The boundary pressure drive picks up the `pp_top` Dirichlet BC (so the Darcy velocity
is set by prescribed pressure + buoyancy).  The relative permeability is taken from the
interior cell (the upwind side for outflow).

**Strictly outflow only**: when the boundary Darcy velocity points inward, the flux is
set to zero — the BC never draws the phase back into the domain.

Capillary and TVD options are intentionally omitted: an open far boundary does not need
front sharpening.

## Example Input File Syntax

```text
[FVBCs]
  [co2_outflow]
    type = VEFVAdvectiveOutflowBC
    variable   = sat_n
    boundary   = right
    fluid_phase = 0
    pp_top     = pp_top
    z_top      = z_top
    z_bottom   = z_bottom
  []
[]
```

!syntax parameters /FVBCs/VEFVAdvectiveOutflowBC

!syntax inputs /FVBCs/VEFVAdvectiveOutflowBC

!syntax children /FVBCs/VEFVAdvectiveOutflowBC
