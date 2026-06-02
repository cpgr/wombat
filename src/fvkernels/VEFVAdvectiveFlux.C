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

  params.addParam<bool>(
      "capillary",
      false,
      "If true, adds grad(Pc^up).n to the CO2 (fluid_phase=0) Darcy potential "
      "gradient as ve_dPcup_dsatn(face)*grad(sat_n).n + ve_dPcup_dH(face)*grad(H).n. "
      "Has no effect on the brine equation. Requires VEFVCapPressure in "
      "[FunctorMaterials] (declares ve_dPcup_dsatn, ve_dPcup_dH). The grad(H).n term "
      "vanishes for constant thickness, so flat/constant-H golds are unchanged. "
      "Default OFF so existing inputs are unchanged.");

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
    _capillary(getParam<bool>("capillary")),
    _advected_interp_method(
        Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("advected_interp_method"))),
    _pp_top(getFunctor<ADReal>("pp_top")),
    _z_top(getFunctor<ADReal>("z_top")),
    _z_bottom(getFunctor<ADReal>("z_bottom")),
    _relperm(getFunctor<ADReal>(
        getParam<unsigned int>("fluid_phase") == 0 ? "ve_relperm_n" : "ve_relperm_w")),
    _density(getFunctor<ADReal>(
        getParam<unsigned int>("fluid_phase") == 0 ? "ve_density_n" : "ve_density_w")),
    _viscosity(getFunctor<ADReal>(
        getParam<unsigned int>("fluid_phase") == 0 ? "ve_viscosity_n" : "ve_viscosity_w")),
    _dPcup_dsatn(_capillary && _fluid_phase == 0 ? &getFunctor<ADReal>("ve_dPcup_dsatn")
                                                 : nullptr),
    _dPcup_dH(_capillary && _fluid_phase == 0 ? &getFunctor<ADReal>("ve_dPcup_dH")
                                              : nullptr),
    _K_up(getMaterialProperty<RealTensorValue>("ve_K_up"))
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
  // Density and viscosity are smooth in pressure, so they are taken at the same
  // boundary-aware face argument as the geometry/pressure above. singleSidedFaceArg
  // uses CentralDifference, so on an interior face this is the two-sided face average
  // (face-correct for a real EOS, AD wrt pp_top on both sides) and on a boundary face
  // it is the BC-extrapolated value. For constant fluid properties it is exactly that
  // constant. Being flow-direction independent, it creates no circular dependence with
  // the upwind selection below.
  const ADReal rho_c = _density(face, state);
  const ADReal mu_c = _viscosity(face, state);
  ADReal dphi_dn = grad_pp_n + rho_c * _gravity_magnitude * grad_zt_n;
  if (_capillary && _fluid_phase == 0)
  {
    // grad(Pc^up).n = d(Pc^up)/d(sat_n)*grad(sat_n).n + d(Pc^up)/d(H)*grad(H).n.
    // gradUDotNormal() is the sat_n variable's boundary-aware face-normal gradient
    // (the FVDiffusion idiom), so a Dirichlet sat_n inlet drives capillary
    // imbibition. The coefficients are material-functor VALUEs (no material-functor
    // gradient is taken). grad(H).n uses the z_top/z_bottom functor gradients --
    // fixed geometry, so taking their Green-Gauss gradient is safe (unlike sat_n).
    const ADReal grad_H_n = grad_zt_n - _z_bottom.gradient(face, state) * _normal;
    dphi_dn += (*_dPcup_dsatn)(face, state) * gradUDotNormal(state, false) +
               (*_dPcup_dH)(face, state) * grad_H_n;
  }
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
