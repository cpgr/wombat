//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "wombatTestApp.h"
#include "wombatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
wombatTestApp::validParams()
{
  InputParameters params = wombatApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

wombatTestApp::wombatTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  wombatTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

wombatTestApp::~wombatTestApp() {}

void
wombatTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  wombatApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"wombatTestApp"});
    Registry::registerActionsTo(af, {"wombatTestApp"});
  }
}

void
wombatTestApp::registerApps()
{
  registerApp(wombatApp);
  registerApp(wombatTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
wombatTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  wombatTestApp::registerAll(f, af, s);
}
extern "C" void
wombatTestApp__registerApps()
{
  wombatTestApp::registerApps();
}
