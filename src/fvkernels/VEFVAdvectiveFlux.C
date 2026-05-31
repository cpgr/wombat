#include "VEFVAdvectiveFlux.h"

registerMooseObject("wombatApp", VEFVAdvectiveFlux);

InputParameters
VEFVAdvectiveFlux::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription(
      "FV depth-integrated Darcy advective flux for one fluid phase. Assign to "
      "variable=sat_n (fluid_phase=0) for CO2 and variable=pp_top (fluid_phase=1) "
      "for brine. Geometry and pressure gradients use boundary-aware functors; the "
      "relative permeability is upwinded at the face via a functor (ve_relperm_n/w "
      "from VEFVRelPerm), so Dirichlet saturation inlets drive inflow.");

  params.addRequiredParam<UserObjectName>("VEDictator", "The VEDictator UserObject.");
  params.addParam<unsigned int>("fluid_phase", 0,
      "Fluid phase index (0 = CO2, 1 = brine).");
  RealVectorValue g_default(0.0, 0.0, -9.81);
  params.addParam<RealVectorValue>("gravity", g_default,
      "Gravity vector [m/s2]. Only the magnitude enters the buoyancy term.");

  params.addRequiredParam<MooseFunctorName>("pp_top",
      "Pore pressure at the top surface [Pa]. Its face-normal gradient is the "
      "Darcy pressure drive (carries AD derivatives wrt pp_top).");
  params.addRequiredParam<MooseFunctorName>("z_top",
      "Top-surface elevation z_T [m] as an FV variable. Supplies H and the "
      "buoyancy slope grad(z_T).");
  params.addRequiredParam<MooseFunctorName>("z_bottom",
      "Bottom-surface elevation z_B [m] as an FV variable. H = z_top - z_bottom.");

  MooseEnum advected_interp_method("average upwind vanLeer min_mod sou quick venkatakrishnan",
                                   "upwind");
  params.addParam<MooseEnum>("advected_interp_method", advected_interp_method,
      "Face interpolation for the advected relperm kr_c. 'upwind' (default, first "
      "order) takes the upstream value; 'average' is unstabilised central "
      "(verification only); 'vanLeer'/'min_mod'/'sou'/'quick'/'venkatakrishnan' are "
      "second-order TVD limiters that sharpen the plume front (transient only).");

  // Two ghost layers so functor gradient reconstruction and upwind face
  // evaluation are available across the face, including in parallel.
  params.set<unsigned short>("ghost_layers") = 2;
  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

VEFVAdvectiveFlux::VEFVAdvectiveFlux(const InputParameters & parameters)
  : FVFluxKernel(parameters),
    _dictator(getUserObject<VEDictator>("VEDictator")),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _advected_interp_method(
        Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("advected_interp_method"))),
    _pp_top(getFunctor<ADReal>("pp_top")),
    _z_top(getFunctor<ADReal>("z_top")),
    _z_bottom(getFunctor<ADReal>("z_bottom")),
    _relperm(getFunctor<ADReal>(
        getParam<unsigned int>("fluid_phase") == 0 ? "ve_relperm_n" : "ve_relperm_w")),
    _K_up(getMaterialProperty<RealTensorValue>("ve_K_up")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _viscosity(getADMaterialProperty<std::vector<Real>>("ve_viscosity"))
{
  if (_fluid_phase >= _dictator.numPhases())
    paramError("fluid_phase", "fluid_phase=", _fluid_phase,
               " but VEDictator reports only ", _dictator.numPhases(), " phases.");
}

ADReal
VEFVAdvectiveFlux::computeQpResidual()
{
  using namespace Moose::FV;

  const auto state = determineState();
  const auto face = singleSidedFaceArg();

  // --- Geometry and pressure gradient via boundary-aware functors ---
  const ADReal H_face = _z_top(face, state) - _z_bottom(face, state);
  const ADReal grad_zt_n = _z_top.gradient(face, state) * _normal;
  const ADReal grad_pp_n = _pp_top.gradient(face, state) * _normal;

  // --- Permeability projected on the face normal (constant -> elem-side) ---
  const Real K_nn = _normal * (_K_up[_qp] * _normal);

  // --- Darcy potential gradient . normal, and advecting flux ---
  const ADReal rho_c = _density[_qp][_fluid_phase];
  const ADReal mu_c = _viscosity[_qp][_fluid_phase];
  const ADReal dphi_dn = grad_pp_n + rho_c * _gravity_magnitude * grad_zt_n;
  const ADReal darcy_velocity_n = -K_nn * dphi_dn;

  // --- Upwind relperm via the functor at the upwind face ---
  // elem_is_upwind: flux outward from elem (darcy_velocity_n > 0) => elem is
  // upstream. On an inflow boundary (darcy_velocity_n < 0) the face/BC side is
  // upstream, so a Dirichlet sat_n inlet supplies the mobility.
  const bool elem_is_upwind = MetaPhysicL::raw_value(darcy_velocity_n) > 0;
  const auto upwind_face =
      makeFace(*_face_info, limiterType(_advected_interp_method), elem_is_upwind, false);
  const ADReal kr = _relperm(upwind_face, state);

  const ADReal mobility = H_face * kr * rho_c / mu_c;

  // Outward-normal flux. The leading minus (in darcy_velocity_n) matches the FE
  // kernel sign: the FE ADKernelGrad form represents -div(F), the FV form +div(F).
  return mobility * darcy_velocity_n;
}
