# VEDissolvedCO2MassPostprocessor

!syntax description /Postprocessors/VEDissolvedCO2MassPostprocessor

## Overview

`VEDissolvedCO2MassPostprocessor` integrates the CO2 mass dissolved in the brine
column:

$$
M_{diss} = \int_\Omega \rho_{brine} \, \bar{\phi} \, H \, (1 - S_n) \, X_{CO_2} \, dA
  \quad (\text{kg})
$$

where $X_{CO_2}$ is the dissolved CO2 mass fraction in brine (kg CO2 / kg brine),
supplied as a coupled `AuxVariable`.  It defaults to zero when not coupled, making
this postprocessor a no-op until dissolution transport is active.

!alert note title=Sink-only vs. full dissolution transport
In the current v1 sink-only formulation, the dissolved inventory is tracked by
`VEDissolvedCO2Aux` as an elemental areal mass variable `c_diss` (kg/m$^2$).  For
that case, use `ElementIntegralVariablePostprocessor` on `c_diss` directly.
`VEDissolvedCO2MassPostprocessor` is intended for future full dissolution-transport
solves where $X_{CO_2}$ is a field variable.

## Example Input File Syntax

```text
[Postprocessors]
  [dissolved_co2_mass]
    type = VEDissolvedCO2MassPostprocessor
    sat_n    = sat_n
    X_co2    = X_co2    # dissolved CO2 mass fraction AuxVariable
    rho_brine = 1020.0
    execute_on = TIMESTEP_END
  []
[]
```

!syntax parameters /Postprocessors/VEDissolvedCO2MassPostprocessor

!syntax inputs /Postprocessors/VEDissolvedCO2MassPostprocessor

!syntax children /Postprocessors/VEDissolvedCO2MassPostprocessor
