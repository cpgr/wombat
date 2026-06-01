# Real-field input pipeline test: read spatially-varying ELEMENTAL petrophysical fields from
# an Exodus mesh (produced by exodus_fields_generator.i, standing in for the Petrel-to-2D
# upscaling workflow) and verify they reach the VE materials correctly.
#
# Reads H, phi_bar, K_up_xx/yy/xy, z_top via initial_from_file_var, feeds them to VEGeometry
# (thickness mode), VEPorosity and VEPermeability, and outputs the resulting material
# properties (ve_H, ve_phi_bar, ve_K_up) via MaterialReal[Tensor]Aux. No flow solve is needed
# (no VEDictator) -- this isolates the field-reading + material-plumbing path.
#
# At the element containing (50,50) the analytic fields give (see generator):
#   ve_H = 102.5, ve_phi_bar = 0.21, ve_K_up = [[1.05e-12, 2.5e-14],[2.5e-14, 5e-13]].

[Mesh]
  type = FileMesh
  file = exodus_fields_generator_out.e
[]

[Problem]
  solve = false
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  # --- Fields read from the Exodus mesh (elemental) ---
  [H]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = H
    initial_from_file_timestep = LATEST
  []
  [phi_bar]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = phi_bar
    initial_from_file_timestep = LATEST
  []
  [K_up_xx]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = K_up_xx
    initial_from_file_timestep = LATEST
  []
  [K_up_yy]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = K_up_yy
    initial_from_file_timestep = LATEST
  []
  [K_up_xy]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = K_up_xy
    initial_from_file_timestep = LATEST
  []
  [z_top]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = z_top
    initial_from_file_timestep = LATEST
  []

  # --- Material-property readouts ---
  [ve_H_out]
    family = MONOMIAL
    order = CONSTANT
  []
  [ve_phi_out]
    family = MONOMIAL
    order = CONSTANT
  []
  [ve_Kxx_out]
    family = MONOMIAL
    order = CONSTANT
  []
  [ve_Kxy_out]
    family = MONOMIAL
    order = CONSTANT
  []
  [ve_Kyy_out]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [ve_H_aux]
    type = MaterialRealAux
    property = ve_H
    variable = ve_H_out
    execute_on = 'INITIAL'
  []
  [ve_phi_aux]
    type = MaterialRealAux
    property = ve_phi_bar
    variable = ve_phi_out
    execute_on = 'INITIAL'
  []
  [ve_Kxx_aux]
    type = MaterialRealTensorValueAux
    property = ve_K_up
    row = 0
    column = 0
    variable = ve_Kxx_out
    execute_on = 'INITIAL'
  []
  [ve_Kxy_aux]
    type = MaterialRealTensorValueAux
    property = ve_K_up
    row = 0
    column = 1
    variable = ve_Kxy_out
    execute_on = 'INITIAL'
  []
  [ve_Kyy_aux]
    type = MaterialRealTensorValueAux
    property = ve_K_up
    row = 1
    column = 1
    variable = ve_Kyy_out
    execute_on = 'INITIAL'
  []
[]

[Materials]
  [geometry]
    type = VEGeometry
    z_top = z_top
    H = H
  []
  [porosity]
    type = VEPorosity
    phi_bar = phi_bar
  []
  [permeability]
    type = VEPermeability
    K_up_xx = K_up_xx
    K_up_yy = K_up_yy
    K_up_xy = K_up_xy
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 1
[]

[Postprocessors]
  [H_probe]
    type = PointValue
    variable = ve_H_out
    point = '50 50 0'
    execute_on = 'INITIAL'
  []
  [phi_probe]
    type = PointValue
    variable = ve_phi_out
    point = '50 50 0'
    execute_on = 'INITIAL'
  []
  [Kxx_probe]
    type = PointValue
    variable = ve_Kxx_out
    point = '50 50 0'
    execute_on = 'INITIAL'
  []
  [Kxy_probe]
    type = PointValue
    variable = ve_Kxy_out
    point = '50 50 0'
    execute_on = 'INITIAL'
  []
  [Kyy_probe]
    type = PointValue
    variable = ve_Kyy_out
    point = '50 50 0'
    execute_on = 'INITIAL'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL'
  []
[]
