#include "VEFVAdvectiveOutflowBC.h"

registerMooseObject("wombatApp", VEFVAdvectiveOutflowBC);

InputParameters
VEFVAdvectiveOutflowBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addClassDescription(
      "Open advective-outflow BC for the depth-integrated VE flow (boundary "
      "counterpart of VEFVAdvectiveFlux). Imposes the interior Darcy advective flux "
      "F_c.n = -H*K_nn*(kr_c*rho_c/mu_c)*(grad(pp_top).n + rho_c*|g|*grad(z_T).n) on "
      "the boundary face, with kr_c taken from the interior (upwind-for-outflow) "
      "cell, so a phase leaves the domain instead of accumulating against a no-flow "
      "wall. Strictly outflow: the flux is zero where the boundary Darcy velocity "
      "points inward, so it never draws the phase back in. Assign to variable=sat_n "
      "(fluid_phase=0) on the updip boundary.");

  params.addParam<unsigned int>("fluid_phase", 0,
      "Fluid phase index (0 = CO2, 1 = brine).");
  RealVectorValue g_default(0.0, 0.0, -9.81);
  params.addParam<RealVectorValue>("gravity", g_default,
      "Gravity vector [m/s2]. Only the magnitude enters the buoyancy term. Must "
      "match the VEFVAdvectiveFlux setting.");

  params.addRequiredParam<MooseFunctorName>("pp_top",
      "Pore pressure at the top surface [Pa]. Its boundary-face-normal gradient "
      "(which picks up the pp_top Dirichlet BC) is the Darcy pressure drive.");
  params.addRequiredParam<MooseFunctorName>("z_top",
      "Top-surface elevation z_T [m] as an FV variable. Supplies H and grad(z_T).");
  params.addRequiredParam<MooseFunctorName>("z_bottom",
      "Bottom-surface elevation z_B [m] as an FV variable. H = z_top - z_bottom.");

  // Two ghost layers so functor gradient reconstruction is available across the
  // boundary face, including in parallel (matches VEFVAdvectiveFlux).
  params.set<unsigned short>("ghost_layers") = 2;
  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

VEFVAdvectiveOutflowBC::VEFVAdvectiveOutflowBC(const InputParameters & parameters)
  : FVFluxBC(parameters),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pp_top(getFunctor<ADReal>("pp_top")),
    _z_top(getFunctor<ADReal>("z_top")),
    _z_bottom(getFunctor<ADReal>("z_bottom")),
    _relperm(getFunctor<ADReal>(
        getParam<unsigned int>("fluid_phase") == 0 ? "ve_relperm_n" : "ve_relperm_w")),
    _K_up(getMaterialProperty<RealTensorValue>("ve_K_up")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _viscosity(getADMaterialProperty<std::vector<Real>>("ve_viscosity"))
{
  // VE is always a two-phase CO2-brine system (phase 0 = CO2, 1 = brine).
  if (_fluid_phase > 1)
    paramError("fluid_phase", "fluid_phase=", _fluid_phase,
               " but a VE simulation has only 2 fluid phases (0 = CO2, 1 = brine).");
}

ADReal
VEFVAdvectiveOutflowBC::computeQpResidual()
{
  const auto state = determineState();
  const auto face = singleSidedFaceArg();

  // --- Geometry and pressure gradient via boundary-aware functors. The pp_top
  //     gradient incorporates the pp_top Dirichlet BC on this boundary. ---
  const ADReal H_face = _z_top(face, state) - _z_bottom(face, state);
  const ADReal grad_zt_n = _z_top.gradient(face, state) * _normal;
  const ADReal grad_pp_n = _pp_top.gradient(face, state) * _normal;

  // --- Permeability projected on the (outward) face normal; constant -> elem-side ---
  const ADReal K_nn = _normal * (_K_up[_qp] * _normal);

  // --- Darcy potential gradient . n and the boundary Darcy velocity ---
  const ADReal rho_c = _density[_qp][_fluid_phase];
  const ADReal mu_c = _viscosity[_qp][_fluid_phase];
  const ADReal dphi_dn = grad_pp_n + rho_c * _gravity_magnitude * grad_zt_n;
  const ADReal darcy_velocity_n = -K_nn * dphi_dn;

  // --- Outflow only. If the boundary Darcy velocity points inward (or is zero) the
  //     flux is zero: the exterior holds none of this phase, so nothing flows in.
  //     This is what makes it an OUTFLOW BC rather than a transparent face. There is a
  //     derivative kink at darcy_velocity_n = 0, but the flux magnitude vanishes there
  //     so it does not impede the Newton solve in practice (cf. the FV upwind-switch
  //     kink, HANDOVER Gotcha C). ---
  if (MetaPhysicL::raw_value(darcy_velocity_n) <= 0.0)
    return 0.0;

  // --- Mobility from the interior cell (the upwind side for outflow). On this open
  //     boundary sat_n carries no Dirichlet BC, so the relperm functor extrapolates
  //     the interior cell value: the flux uses the mobility of the phase actually
  //     present inside, which the outward velocity carries out. ---
  const ADReal kr = _relperm(face, state);
  const ADReal mobility = H_face * kr * rho_c / mu_c;

  // Outward-normal flux, same +div(F) sign convention as VEFVAdvectiveFlux.
  return mobility * darcy_velocity_n;
}
