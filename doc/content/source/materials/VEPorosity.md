# VEPorosity

!syntax description /Materials/VEPorosity

## Overview

`VEPorosity` provides the depth-averaged porosity `ve_phi_bar` (-) consumed by
`VEMassTimeDerivative`, `VEFVMassTimeDerivative`, and the diagnostic postprocessors.

Supply `phi_bar` as a constant scalar for homogeneous verification cases, or couple it
to a spatially varying `AuxVariable` read from the upscaled Exodus mesh for field
models.  MOOSE resolves both through the same `coupledValue` interface.

`ve_phi_bar` is a plain `Real` (not AD): porosity is time-invariant and carries no
Jacobian contribution with respect to the primary flow variables.

## Example Input File Syntax

Constant porosity (verification):

```text
[Materials]
  [porosity]
    type = VEPorosity
    phi_bar = 0.20
  []
[]
```

Spatially varying porosity from an Exodus AuxVariable:

```text
[AuxVariables]
  [phi_bar]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[Materials]
  [porosity]
    type = VEPorosity
    phi_bar = phi_bar
  []
[]
```

!syntax parameters /Materials/VEPorosity

!syntax inputs /Materials/VEPorosity

!syntax children /Materials/VEPorosity
