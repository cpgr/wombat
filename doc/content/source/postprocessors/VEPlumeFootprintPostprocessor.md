# VEPlumeFootprintPostprocessor

!syntax description /Postprocessors/VEPlumeFootprintPostprocessor

## Overview

`VEPlumeFootprintPostprocessor` integrates the area of elements where the
depth-averaged CO2 saturation `sat_n` exceeds a threshold, reporting the 2-D
footprint of the CO2 plume:

$$
A_{plume} = \int_\Omega \mathbf{1}(S_n > \varepsilon) \, dA \quad (\text{m}^2)
$$

The threshold $\varepsilon$ defaults to $10^{-6}$ to exclude numerical noise at the
plume front.

Combine with `VEMobileCO2MassPostprocessor` and `VETrappedCO2MassPostprocessor` for
a complete CO2 inventory summary output at each timestep.

## Example Input File Syntax

```text
[Postprocessors]
  [plume_footprint]
    type = VEPlumeFootprintPostprocessor
    sat_n     = sat_n
    threshold = 1.0e-5
    execute_on = TIMESTEP_END
  []
[]
```

!syntax parameters /Postprocessors/VEPlumeFootprintPostprocessor

!syntax inputs /Postprocessors/VEPlumeFootprintPostprocessor

!syntax children /Postprocessors/VEPlumeFootprintPostprocessor
