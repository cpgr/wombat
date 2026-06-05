# VEPlumeHeightAux

!syntax description /AuxKernels/VEPlumeHeightAux

## Overview

`VEPlumeHeightAux` exposes the CO2 plume thickness `ve_h` (m) computed by
`VEPlumeReconstruction` as an elemental `AuxVariable` for Exodus output and
postprocessing.

The plume height is not a primary solution variable; it is a derived material property.
This kernel bridges the gap between the material layer and the output layer.

## Example Input File Syntax

```text
[AuxVariables]
  [plume_height]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [plume_height_aux]
    type = VEPlumeHeightAux
    variable = plume_height
    execute_on = 'TIMESTEP_END'
  []
[]
```

!syntax parameters /AuxKernels/VEPlumeHeightAux

!syntax inputs /AuxKernels/VEPlumeHeightAux

!syntax children /AuxKernels/VEPlumeHeightAux
