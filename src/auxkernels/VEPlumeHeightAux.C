#include "VEPlumeHeightAux.h"

registerMooseObject("wombatApp", VEPlumeHeightAux);

InputParameters
VEPlumeHeightAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Exposes the ve_h material property (CO2 plume thickness (m)) computed by "
      "VEPlumeReconstruction as an elemental AuxVariable for output and postprocessing.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPlumeHeightAux::VEPlumeHeightAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _h(getADMaterialProperty<Real>("ve_h"))
{
}

Real
VEPlumeHeightAux::computeValue()
{
  return MetaPhysicL::raw_value(_h[_qp]);
}
