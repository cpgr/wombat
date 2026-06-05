# VESaturation

!syntax description /Materials/VESaturation

## Overview

`VESaturation` wraps the primary variable `sat_n` (depth-averaged CO2 saturation) into
the two-element `ADMaterialProperty<std::vector<Real>>` named `ve_saturation` expected
by the VE kernels and `VERelPerm`:

$$
\text{ve\_saturation}(0) = S_n \quad (\text{CO}_2,\;\text{non-wetting})
$$
$$
\text{ve\_saturation}(1) = 1 - S_n \quad (\text{brine},\;\text{wetting})
$$

Because the property is declared as `ADMaterialProperty`, derivatives of `sat_n`
propagate automatically into `VEMassTimeDerivative` and `VERelPerm` without any
hand-coded chain rule.

`initQpStatefulProperties()` seeds the old-value storage from the initial condition
so that `VEMassTimeDerivative` computes a correct storage term at the first timestep.

This material is always required in an FE input.  In FV inputs, `VEFVMassTimeDerivative`
reads `ve_saturation` from the elemental `VESaturation` material (the FV functor
materials `VEFVRelPerm` and `VEFVFluidProperties` evaluate `sat_n` directly as a
functor and do not consume this property).

## Example Input File Syntax

```text
[Variables]
  [sat_n] []
[]
[Materials]
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
[]
```

!syntax parameters /Materials/VESaturation

!syntax inputs /Materials/VESaturation

!syntax children /Materials/VESaturation
