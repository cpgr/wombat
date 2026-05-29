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
      "UserObject (Nordbotten & Dahle 2011). Result is clamped to [0, H].");

  MooseEnum modes("sharp_interface capillary_fringe", "sharp_interface");
  params.addParam<MooseEnum>(
      "mode", modes, "Reconstruction mode: sharp_interface or capillary_fringe.");

  params.addRequiredCoupledVar("sat_n", "Depth-averaged CO2 saturation (primary variable).");
  params.addRangeCheckedParam<Real>(
      "S_wr", "S_wr >= 0 & S_wr < 1", "Residual water saturation in the CO2 zone [-].");

  params.addParam<UserObjectName>(
      "pc_uo",
      "PorousFlowCapillaryPressure UserObject defining the Pc(Sw) curve "
      "(capillary_fringe mode only).");
  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2] (capillary_fringe mode only). Defaults to (0, 0, -9.81).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPlumeReconstruction::VEPlumeReconstruction(const InputParameters & parameters)
  : Material(parameters),
    _mode(getParam<MooseEnum>("mode")),
    _sat_n(coupledValue("sat_n")),
    _H(getMaterialProperty<Real>("ve_H")),
    _S_wr(getParam<Real>("S_wr")),
    _pc_uo(nullptr),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _h(declareProperty<Real>("ve_h"))
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
  Real h;

  if (_mode == "sharp_interface")
  {
    h = _sat_n[_qp] * _H[_qp] / (1.0 - _S_wr);
  }
  else
  {
    // Newton inversion: find h such that satNBar(h, H, delta_rho) = sat_n_qp.
    //   F(h)  = satNBar(h, H, delta_rho) - sat_n_qp
    //   F'(h) = S_n(h) / H                          (Leibniz rule on the integral)
    const Real H = _H[_qp];
    const Real target = _sat_n[_qp];
    const Real delta_rho = MetaPhysicL::raw_value(_density[_qp][1]) -
                           MetaPhysicL::raw_value(_density[_qp][0]);

    h = std::max(0.0, std::min(target * H / (1.0 - _S_wr), H)); // warm start

    for (unsigned int iter = 0; iter < 20; ++iter)
    {
      const Real F = satNBar(h, H, delta_rho) - target;
      if (std::abs(F) < 1.0e-12)
        break;
      const Real Sn_top = 1.0 - _pc_uo->saturation(delta_rho * _gravity_magnitude * h);
      if (Sn_top < 1.0e-14)
        break; // plume top is at residual CO2; can't improve further
      h -= F * H / Sn_top;
      h = std::max(0.0, std::min(h, H));
    }
  }

  _h[_qp] = std::max(0.0, std::min(h, _H[_qp]));
}
