# VETrappedCO2MassPostprocessor

!syntax description /Postprocessors/VETrappedCO2MassPostprocessor

## Overview

`VETrappedCO2MassPostprocessor` integrates the depth-averaged residually trapped CO2
mass over the 2-D domain:

\begin{equation}
M_{trap} = \int_\Omega \rho_{CO_2} \, \bar{\phi} \, H \, S_{n,trap} \, dA \quad (\text{kg})
\end{equation}

`sat_n_trap` is coupled from the `VETrappedSaturationAux` AuxVariable.  When not
coupled it defaults to zero, making this postprocessor a no-op until hysteresis is
active.

Pair with `VEMobileCO2MassPostprocessor` and `VEPlumeFootprintPostprocessor` for a
complete CO2 inventory summary.

## Example Input File Syntax

```text
[Postprocessors]
  [trapped_co2_mass]
    type = VETrappedCO2MassPostprocessor
    sat_n_trap = sat_n_trap
    rho_co2    = 700.0
    execute_on = TIMESTEP_END
  []
[]
```

!syntax parameters /Postprocessors/VETrappedCO2MassPostprocessor

!syntax inputs /Postprocessors/VETrappedCO2MassPostprocessor

!syntax children /Postprocessors/VETrappedCO2MassPostprocessor
