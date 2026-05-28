#include "wombatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
wombatApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

wombatApp::wombatApp(const InputParameters & parameters) : MooseApp(parameters)
{
  wombatApp::registerAll(_factory, _action_factory, _syntax);
}

wombatApp::~wombatApp() {}

void
wombatApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<wombatApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"wombatApp"});
  Registry::registerActionsTo(af, {"wombatApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
wombatApp::registerApps()
{
  registerApp(wombatApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
wombatApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  wombatApp::registerAll(f, af, s);
}
extern "C" void
wombatApp__registerApps()
{
  wombatApp::registerApps();
}
