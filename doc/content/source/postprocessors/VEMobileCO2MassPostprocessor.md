# VEMobileCO2MassPostprocessor

!syntax description /Postprocessors/VEMobileCO2MassPostprocessor

## Overview

`VEMobileCO2MassPostprocessor` integrates the depth-averaged mobile CO2 mass over the
2-D domain:

$$
M_{mobile} = \int_\Omega \rho_{CO_2} \, \bar{\phi} \, H \, S_{n,mobile} \, dA \quad (\text{kg})
$$

where $S_{n,mobile} = S_n - S_{n,trap}$ (free CO2 saturation).

The optional `sat_n_trap` AuxVariable (written by `VETrappedSaturationAux`) provides
the trapped fraction.  When not coupled, it defaults to zero and the full `sat_n` is
treated as mobile — correct for non-hysteretic verification cases.

`rho_co2` is a uniform scalar parameter, appropriate when `ConstantFluidProperties` is
used.  Reads `ve_H` from `VEGeometry` and `ve_phi_bar` from `VEPorosity`.

## Example Input File Syntax

```text
[Postprocessors]
  [mobile_co2_mass]
    type = VEMobileCO2MassPostprocessor
    sat_n    = sat_n
    rho_co2  = 700.0
    execute_on = TIMESTEP_END
  []
[]
```

With hysteresis:

```text
[Postprocessors]
  [mobile_co2_mass]
    type = VEMobileCO2MassPostprocessor
    sat_n      = sat_n
    sat_n_trap = sat_n_trap
    rho_co2    = 700.0
    execute_on = TIMESTEP_END
  []
[]
```

!syntax parameters /Postprocessors/VEMobileCO2MassPostprocessor

!syntax inputs /Postprocessors/VEMobileCO2MassPostprocessor

!syntax children /Postprocessors/VEMobileCO2MassPostprocessor
