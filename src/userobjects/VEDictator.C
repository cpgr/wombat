#include "VEDictator.h"

registerMooseObject("wombatApp", VEDictator);

InputParameters
VEDictator::validParams()
{
  // Inherit all PorousFlowDictator parameters (porous_flow_vars,
  // number_fluid_phases, number_fluid_components, ...) so input files look
  // familiar to PorousFlow users.
  InputParameters params = PorousFlowDictator::validParams();

  params.addClassDescription(
      "Dictator for a vertical-equilibrium (VE) CO2-brine simulation. "
      "Always uses 2 fluid phases (CO2=0, brine=1) and 2 fluid components "
      "(CO2=0, water=1). Carries the VE physics flavour flag used by all VE "
      "Materials and Kernels.");

  MooseEnum flavour_enum("sharp_interface capillary_fringe with_hysteresis",
                         "sharp_interface");
  params.addParam<MooseEnum>(
      "ve_flavour",
      flavour_enum,
      "VE physics variant. 'sharp_interface' uses a classical sharp CO2-brine interface with "
      "analytical upscaled relperm/Pc (default, enables the Nordbotten-Celia analytical tests). "
      "'capillary_fringe' adds the Nordbotten-Dahle 2011 fringe via a precomputed lookup table. "
      "'with_hysteresis' extends capillary_fringe with residual trapping scanning curves.");

  // Force the fixed phase/component counts so that the parent-class
  // error-checking logic rejects mis-configured input files early.
  params.setParameters<unsigned int>("number_fluid_phases", 2);
  params.suppressParameter<unsigned int>("number_fluid_phases");

  params.setParameters<unsigned int>("number_fluid_components", 2);
  params.suppressParameter<unsigned int>("number_fluid_components");

  return params;
}

VEDictator::VEDictator(const InputParameters & parameters)
  : PorousFlowDictator(parameters),
    _flavour([&]() {
      const std::string s = getParam<MooseEnum>("ve_flavour").operator std::string();
      if (s == "sharp_interface")
        return VEFlavour::SHARP_INTERFACE;
      if (s == "capillary_fringe")
        return VEFlavour::CAPILLARY_FRINGE;
      return VEFlavour::WITH_HYSTERESIS;
    }())
{
  if (numVariables() != 2)
    paramError("porous_flow_vars",
               "VEDictator requires exactly two primary variables: pp_top and sat_n. "
               "Provide them in that order.");
}
