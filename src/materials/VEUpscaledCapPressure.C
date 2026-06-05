#include "VEUpscaledCapPressure.h"
#include "PorousFlowCapillaryPressure.h"

registerMooseObject("wombatApp", VEUpscaledCapPressure);

InputParameters
VEUpscaledCapPressure::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Upscaled column capillary pressure ve_pc_up = (rho_w - rho_n) * |g| * h + pc_entry "
      "and its partial derivatives ve_dPcup_dsatn = delta_rho*|g|*dh/dsat_n and "
      "ve_dPcup_dH = delta_rho*|g|*dh/dH. Both partials are AD; in sharp_interface mode "
      "dh/dsat_n = H/(1-S_wr) and dh/dH = sat_n/(1-S_wr), in capillary_fringe mode "
      "dh/dsat_n = H/S_n(h) and dh/dH = sat_n/S_n(h) with S_n(h) = 1 - Sw(delta_rho*|g|*h). "
      "Consumed by VEAdvectiveFluxS when capillary != none to assemble "
      "grad(Pc^up) = ve_dPcup_dsatn*grad(sat_n) + ve_dPcup_dH*grad(H).");

  MooseEnum modes("sharp_interface capillary_fringe", "sharp_interface");
  params.addParam<MooseEnum>(
      "mode",
      modes,
      "Reconstruction mode the coefficients are consistent with. Must match the mode of "
      "VEPlumeReconstruction (it supplies ve_h).");

  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector (m/s2). Only the magnitude enters the buoyancy term.");
  params.addParam<Real>(
      "pc_entry", 0.0, "Capillary entry/fringe pressure added at the plume top (Pa).");
  params.addRangeCheckedParam<Real>(
      "S_wr",
      0.0,
      "S_wr >= 0 & S_wr < 1",
      "Residual water saturation (-). Used in sharp_interface mode; must match the value "
      "used in VEPlumeReconstruction so that the partials are consistent with ve_h.");

  params.addParam<UserObjectName>(
      "pc_uo",
      "PorousFlowCapillaryPressure UserObject defining the Sw(Pc) curve. Required in "
      "capillary_fringe mode (supplies S_n(h) = 1 - Sw(delta_rho*|g|*h)); must be the same "
      "UO used by VEPlumeReconstruction.");
  params.addCoupledVar(
      "sat_n",
      0.0,
      "Depth-averaged CO2 saturation (primary variable). Required in capillary_fringe mode "
      "for ve_dPcup_dH = delta_rho*|g|*sat_n/S_n(h); unused in sharp_interface mode (the 0.0 "
      "default is then never read).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEUpscaledCapPressure::VEUpscaledCapPressure(const InputParameters & parameters)
  : Material(parameters),
    _mode(getParam<MooseEnum>("mode")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pc_entry(getParam<Real>("pc_entry")),
    _S_wr(getParam<Real>("S_wr")),
    _h(getADMaterialProperty<Real>("ve_h")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _H(getMaterialProperty<Real>("ve_H")),
    _pc_uo(nullptr),
    _sat_n(adCoupledValue("sat_n")),
    _pc_up(declareADProperty<Real>("ve_pc_up")),
    _dPcup_dsatn(declareADProperty<Real>("ve_dPcup_dsatn")),
    _dPcup_dH(declareADProperty<Real>("ve_dPcup_dH"))
{
  if (_mode == "capillary_fringe")
  {
    if (!isParamValid("pc_uo"))
      paramError("pc_uo", "Must be supplied for capillary_fringe mode.");
    if (!isCoupled("sat_n"))
      paramError("sat_n", "Must be coupled for capillary_fringe mode.");

    _pc_uo = &getUserObject<PorousFlowCapillaryPressure>("pc_uo");
  }
}

ADReal
VEUpscaledCapPressure::saturationAtPlumeTop(const ADReal & delta_rho) const
{
  // S_n(h) = 1 - Sw(Pc(h)), Pc(h) = delta_rho*|g|*h. The AD saturation overload
  // propagates ve_h's d/d(sat_n) through dSw/dPc, so S_n carries the correct
  // Jacobian. Floor near zero (plume nose, where Sw -> 1 and dh/dsat_n -> inf) to
  // keep the H/S_n divisions finite, mirroring VEPlumeReconstruction's guard.
  const ADReal pc_h = delta_rho * _gravity_magnitude * _h[_qp];
  const ADReal Sn = 1.0 - _pc_uo->saturation(pc_h);
  if (MetaPhysicL::raw_value(Sn) < 1.0e-12)
    return ADReal(1.0e-12);
  return Sn;
}

void
VEUpscaledCapPressure::computeQpProperties()
{
  const ADReal delta_rho = _density[_qp][1] - _density[_qp][0];
  _pc_up[_qp] = delta_rho * _gravity_magnitude * _h[_qp] + _pc_entry;

  if (_mode == "sharp_interface")
  {
    // dh/dsat_n = H/(1-S_wr) (constant in sat_n; the grad(sat_n) Jacobian rides
    // on the AD variable gradient _grad_u in the kernel).
    _dPcup_dsatn[_qp] = delta_rho * _gravity_magnitude * _H[_qp] / (1.0 - _S_wr);
    // dh/dH = sat_n/(1-S_wr) = h/H (since sharp h = sat_n*H/(1-S_wr)). Kept AD so
    // the grad(H) term's d/d(sat_n) Jacobian rides on this coefficient.
    _dPcup_dH[_qp] = delta_rho * _gravity_magnitude * _h[_qp] / _H[_qp];
  }
  else
  {
    // Capillary-fringe partials: dh/dsat_n = H/S_n(h), dh/dH = sat_n/S_n(h).
    const ADReal Sn_top = saturationAtPlumeTop(delta_rho);
    _dPcup_dsatn[_qp] = delta_rho * _gravity_magnitude * _H[_qp] / Sn_top;
    _dPcup_dH[_qp] = delta_rho * _gravity_magnitude * _sat_n[_qp] / Sn_top;
  }
}
