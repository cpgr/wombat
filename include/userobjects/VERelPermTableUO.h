#pragma once

#include "VERelativePermeability.h"
#include "LinearInterpolation.h"

/**
 * VERelPermTableUO
 *
 * Tabulated upscaled relative permeability: piecewise-linear lookup of
 * kr_n(sat_n) and kr_w(sat_n) from (sat_n, kr) tables supplied by the
 * fine-scale upscaling workflow. This is the field-case counterpart of the
 * analytical VERelPermSharpUO; both are VERelativePermeability UserObjects, so
 * the FE (VERelPerm) and FV (VEFVRelPerm) adapters consume either with no
 * change.
 *
 * The two curves share the saturation axis (sat_n_points). Outside the tabulated
 * range LinearInterpolation clamps to the end values (zero derivative), which is
 * the desired behaviour for a bounded upscaled curve.
 */
class VERelPermTableUO : public VERelativePermeability
{
public:
  static InputParameters validParams();
  VERelPermTableUO(const InputParameters & parameters);

  Real relativePermeability(Real sat_n, unsigned int phase) const override;
  Real dRelativePermeability(Real sat_n, unsigned int phase) const override;

protected:
  /// kr_n(sat_n): CO2 (non-wetting) relative permeability table.
  LinearInterpolation _krn;
  /// kr_w(sat_n): brine (wetting) relative permeability table.
  LinearInterpolation _krw;
};
