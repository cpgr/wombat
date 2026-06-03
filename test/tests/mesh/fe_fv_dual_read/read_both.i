# Regression guard for the FE+FV dual-read contract that wombat relies on when
# consuming an ex2ve mesh written with --z-discretization both.
#
# dual_field.e (from make_dual_field_mesh.py) carries z_top and z_bottom in BOTH
# the nodal and the elemental Exodus namespaces under the same names, with
#   z_top    = -900 + 0.5x + 0.3y   (flat mesh; topography only in the field)
#   z_bottom = z_top - 100
#
# MOOSE must resolve which namespace to read from the receiving variable's TYPE:
#   LAGRANGE  variable  -> nodal field
#   MooseVariableFVReal -> elemental field
#
# The nodal and elemental forms differ (node values vs centroid means), so the
# min over each is a fingerprint of which namespace was read:
#   z_top    nodal min = -900 (corner 0,0)   ; elemental min = -896 (centroid 5,5)
#   z_bottom nodal min = -1000               ; elemental min = -996
# And VEGeometry (consuming the FE/nodal z_top, z_bottom) must recover the exact
# horizontal slope (0.5, 0.3) on the flat mesh. If a future MOOSE change breaks
# the namespace disambiguation, these collide or move and the CSVDiff fails.

[Mesh]
  type = FileMesh
  file = dual_field.e
[]

[Problem]
  solve = false
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  # FE: nodal LAGRANGE -> must read the NODAL z_top/z_bottom.
  [z_top_fe]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = z_top
    initial_from_file_timestep = LATEST
  []
  [z_bottom_fe]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = z_bottom
    initial_from_file_timestep = LATEST
  []
  # FV: cell-centred -> must read the ELEMENTAL z_top/z_bottom.
  [z_top_fv]
    type = MooseVariableFVReal
    initial_from_file_var = z_top
    initial_from_file_timestep = LATEST
  []
  [z_bottom_fv]
    type = MooseVariableFVReal
    initial_from_file_var = z_bottom
    initial_from_file_timestep = LATEST
  []
  [grad_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [grad_y]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [gx]
    type = MaterialRealVectorValueAux
    variable = grad_x
    property = ve_grad_z_top
    component = 0
    execute_on = 'INITIAL'
  []
  [gy]
    type = MaterialRealVectorValueAux
    variable = grad_y
    property = ve_grad_z_top
    component = 1
    execute_on = 'INITIAL'
  []
[]

[Materials]
  # top_bottom mode exercises BOTH nodal fields: ve_H = z_top_fe - z_bottom_fe
  # and ve_grad_z_top = grad(z_top_fe).
  [geometry]
    type = VEGeometry
    z_top = z_top_fe
    z_bottom = z_bottom_fe
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 1
[]

[Postprocessors]
  [zt_fe_min]                 # expect -900 (nodal read)
    type = NodalExtremeValue
    variable = z_top_fe
    value_type = min
    execute_on = 'INITIAL'
  []
  [zt_fv_min]                 # expect -896 (elemental read)
    type = ElementExtremeValue
    variable = z_top_fv
    value_type = min
    execute_on = 'INITIAL'
  []
  [zb_fe_min]                 # expect -1000 (nodal read)
    type = NodalExtremeValue
    variable = z_bottom_fe
    value_type = min
    execute_on = 'INITIAL'
  []
  [zb_fv_min]                 # expect -996 (elemental read)
    type = ElementExtremeValue
    variable = z_bottom_fv
    value_type = min
    execute_on = 'INITIAL'
  []
  [grad_x_avg]               # expect 0.5
    type = ElementAverageValue
    variable = grad_x
    execute_on = 'INITIAL'
  []
  [grad_y_avg]               # expect 0.3
    type = ElementAverageValue
    variable = grad_y
    execute_on = 'INITIAL'
  []
[]

[Outputs]
  execute_on = 'INITIAL'
  [csv]
    type = CSV
  []
[]
