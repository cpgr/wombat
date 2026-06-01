# Generates a small 2-D Exodus mesh carrying spatially-varying ELEMENTAL petrophysical
# fields, standing in for the output of the Petrel-to-2D upscaling workflow. Read back by
# exodus_fields.i to test the real-field input pipeline (FileMesh + initial_from_file_var ->
# VEGeometry / VEPorosity / VEPermeability).
#
# Fields are linear (so the CONSTANT MONOMIAL element value equals the centroid value
# exactly), giving analytic gold values for the reader. On the 4x2 mesh over [0,400]x[0,200]
# the element containing (50,50) has centroid (50,50):
#   H        = 100 + 0.05*x        -> 102.5
#   phi_bar  = 0.2 + 0.0002*x      -> 0.21
#   K_up_xx  = 1e-12*(1 + 0.001*x) -> 1.05e-12
#   K_up_yy  = 5e-13               -> 5e-13
#   K_up_xy  = 1e-13*(y/200)       -> 2.5e-14
#   z_top    = -1000 + 0.02*x      -> structural surface (not verified here)

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 2
  xmin = 0
  xmax = 400
  ymin = 0
  ymax = 200
[]

[Problem]
  solve = false
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  [H]
    family = MONOMIAL
    order = CONSTANT
  []
  [phi_bar]
    family = MONOMIAL
    order = CONSTANT
  []
  [K_up_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [K_up_yy]
    family = MONOMIAL
    order = CONSTANT
  []
  [K_up_xy]
    family = MONOMIAL
    order = CONSTANT
  []
  [z_top]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [H_fn]
    type = ParsedFunction
    expression = '100 + 0.05*x'
  []
  [phi_fn]
    type = ParsedFunction
    expression = '0.2 + 0.0002*x'
  []
  [kxx_fn]
    type = ParsedFunction
    expression = '1.0e-12*(1 + 0.001*x)'
  []
  [kyy_fn]
    type = ParsedFunction
    expression = '5.0e-13'
  []
  [kxy_fn]
    type = ParsedFunction
    expression = '1.0e-13*(y/200)'
  []
  [ztop_fn]
    type = ParsedFunction
    expression = '-1000 + 0.02*x'
  []
[]

[AuxKernels]
  [H_aux]
    type = FunctionAux
    variable = H
    function = H_fn
    execute_on = 'INITIAL'
  []
  [phi_aux]
    type = FunctionAux
    variable = phi_bar
    function = phi_fn
    execute_on = 'INITIAL'
  []
  [kxx_aux]
    type = FunctionAux
    variable = K_up_xx
    function = kxx_fn
    execute_on = 'INITIAL'
  []
  [kyy_aux]
    type = FunctionAux
    variable = K_up_yy
    function = kyy_fn
    execute_on = 'INITIAL'
  []
  [kxy_aux]
    type = FunctionAux
    variable = K_up_xy
    function = kxy_fn
    execute_on = 'INITIAL'
  []
  [ztop_aux]
    type = FunctionAux
    variable = z_top
    function = ztop_fn
    execute_on = 'INITIAL'
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 1
[]

[Outputs]
  execute_on = 'INITIAL'
  exodus = true
[]
