#include "VEPlumeReconstruction.h"
#include "PorousFlowCapillaryPressure.h"

registerMooseObject("wombatApp", VEPlumeReconstruction);

InputParameters
VEPlumeReconstruction::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Reconstructs CO2 plume thickness ve_h from depth-averaged saturation. "
      "sharp_interface mode: h = sat_n * H / (1 - S_wr). "
      "capillary_fringe mode: Newton inversion of the vertically-integrated "
      "capillary-gravity equilibrium profile using a PorousFlowCapillaryPressure "
      "UserObject (Nordbotten & Dahle 2011). Result is clamped to (0, H).");

  MooseEnum modes("sharp_interface capillary_fringe", "sharp_interface");
  params.addParam<MooseEnum>(
      "mode", modes, "Reconstruction mode: sharp_interface or capillary_fringe.");

  params.addRequiredCoupledVar("sat_n", "Depth-averaged CO2 saturation (primary variable).");
  params.addRangeCheckedParam<Real>(
      "S_wr", "S_wr >= 0 & S_wr < 1", "Residual water saturation in the CO2 zone (-).");

  params.addParam<UserObjectName>(
      "pc_uo",
      "PorousFlowCapillaryPressure UserObject defining the Pc(Sw) curve "
      "(capillary_fringe mode only).");
  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector (m/s2) (capillary_fringe mode only). Defaults to (0, 0, -9.81).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPlumeReconstruction::VEPlumeReconstruction(const InputParameters & parameters)
  : Material(parameters),
    _mode(getParam<MooseEnum>("mode")),
    _sat_n(adCoupledValue("sat_n")),
    _H(getMaterialProperty<Real>("ve_H")),
    _S_wr(getParam<Real>("S_wr")),
    _pc_uo(nullptr),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _h(declareADProperty<Real>("ve_h"))
{
  if (_mode == "capillary_fringe")
  {
    if (!isParamValid("pc_uo"))
      paramError("pc_uo", "Must be supplied for capillary_fringe mode.");

    _pc_uo = &getUserObject<PorousFlowCapillaryPressure>("pc_uo");
  }
}

Real
VEPlumeReconstruction::satNBar(Real h, Real H, Real delta_rho) const
{
  // Trapezoidal integration of S_n(zeta) = 1 - Sw(Pc(zeta)) over [0, h].
  // Pc(zeta) = delta_rho * g * zeta is the buoyancy pressure at height zeta
  // above the CO2-brine contact.
  const unsigned int N = 64;
  const Real dz = h / N;
  Real sum = 0.0;
  for (unsigned int i = 0; i <= N; ++i)
  {
    const Real pc = delta_rho * _gravity_magnitude * i * dz;
    const Real Sn = 1.0 - _pc_uo->saturation(pc);
    sum += (i == 0 || i == N) ? 0.5 * Sn : Sn;
  }
  return sum * dz / H;
}

void
VEPlumeReconstruction::computeQpProperties()
{
  const Real H = _H[_qp];

  if (_mode == "sharp_interface")
  {
    // Natural AD propagation: d(h)/d(sat_n_dof) = H/(1-S_wr) * d(sat_n)/d(dof).
    const ADReal h = _sat_n[_qp] * H / (1.0 - _S_wr);
    if (MetaPhysicL::raw_value(h) <= 0.0)
      _h[_qp] = ADReal(0.0);
    else if (MetaPhysicL::raw_value(h) >= H)
      _h[_qp] = ADReal(H);
    else
      _h[_qp] = h;
  }
  else
  {
    // Newton inversion on raw values. Derivative seeded manually via the
    // inverse function theorem: d(h)/d(sat_n) = H / S_n(h).
    const Real raw_sat = MetaPhysicL::raw_value(_sat_n[_qp]);
    const Real delta_rho = MetaPhysicL::raw_value(_density[_qp][1]) -
                           MetaPhysicL::raw_value(_density[_qp][0]);

    Real h_val = std::max(0.0, std::min(raw_sat * H / (1.0 - _S_wr), H)); // warm start

    for (unsigned int iter = 0; iter < 20; ++iter)
    {
      const Real F = satNBar(h_val, H, delta_rho) - raw_sat;
      if (std::abs(F) < 1.0e-12)
        break;
      const Real Sn_top = 1.0 - _pc_uo->saturation(delta_rho * _gravity_magnitude * h_val);
      if (Sn_top < 1.0e-14)
        break;
      h_val -= F * H / Sn_top;
      h_val = std::max(0.0, std::min(h_val, H));
    }

    if (h_val <= 0.0 || h_val >= H)
    {
      _h[_qp] = ADReal(h_val);
    }
    else
    {
      // Linearise: h(sat_n) ~= h_val + dh/dsatn * (sat_n - raw_sat).
      // Value part = h_val, derivative = (H/Sn_h) * d(sat_n)/d(dof).
      const Real Sn_h = 1.0 - _pc_uo->saturation(delta_rho * _gravity_magnitude * h_val);
      const Real dh_dsatn = (Sn_h > 1.0e-14) ? H / Sn_h : H / 1.0e-14;
      _h[_qp] = h_val + dh_dsatn * (_sat_n[_qp] - raw_sat);
    }
  }
}
