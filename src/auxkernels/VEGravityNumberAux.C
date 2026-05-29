#include "VEGravityNumberAux.h"

registerMooseObject("wombatApp", VEGravityNumberAux);

InputParameters
VEGravityNumberAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes the VE gravity number Gamma = k_v * delta_rho * g * H^2 "
      "/ (mu_n * phi_bar * Q * L).  Large Gamma (>> 1) supports the VE "
      "assumption; small Gamma (< ~10) should be flagged for review.");

  params.addRequiredParam<Real>("k_v", "Vertical permeability [m^2].");
  params.addRequiredParam<Real>("delta_rho",
                                "Density difference rho_brine - rho_CO2 [kg/m^3].");
  params.addRequiredParam<Real>("mu_n", "CO2 (non-wetting phase) viscosity [Pa*s].");
  params.addRequiredParam<Real>(
      "Q",
      "Characteristic volumetric injection rate [m^3/s].  Used as the "
      "reference flow scale; typically the total well rate for the domain.");
  params.addRequiredParam<Real>("L",
                                "Characteristic horizontal length of the domain [m].");
  RealVectorValue g_default(0.0, 0.0, -9.81);
  params.addParam<RealVectorValue>(
      "gravity", g_default, "Gravity vector [m/s^2]. Defaults to (0, 0, -9.81).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEGravityNumberAux::VEGravityNumberAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _H(getMaterialProperty<Real>("ve_H")),
    _phi_bar(getMaterialProperty<Real>("ve_phi_bar")),
    _k_v(getParam<Real>("k_v")),
    _delta_rho(getParam<Real>("delta_rho")),
    _mu_n(getParam<Real>("mu_n")),
    _Q(getParam<Real>("Q")),
    _L(getParam<Real>("L")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm())
{
  if (_k_v <= 0.0)
    paramError("k_v", "Vertical permeability must be positive.");
  if (_delta_rho <= 0.0)
    paramError("delta_rho", "Density difference must be positive.");
  if (_mu_n <= 0.0)
    paramError("mu_n", "CO2 viscosity must be positive.");
  if (_Q <= 0.0)
    paramError("Q", "Reference injection rate must be positive.");
  if (_L <= 0.0)
    paramError("L", "Characteristic length must be positive.");
}

Real
VEGravityNumberAux::computeValue()
{
  const Real denom = _mu_n * _phi_bar[_qp] * _Q * _L;
  return _k_v * _delta_rho * _gravity_magnitude * _H[_qp] * _H[_qp] / denom;
}
