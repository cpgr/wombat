# Tests VEGravityNumberAux with a SPATIALLY-VARYING coupled k_v (instead of a
# constant param). k_v is supplied as an AuxVariable that varies linearly in x,
# standing in for an upscaled PERMZ field read from an ex2ve mesh, so the
# VE-validity map honours the real vertical permeability.
#
# Flat formation: z_top = 0, z_bottom = -20 -> H = 20 m, phi_bar = 0.2.
#   gamma(x) = k_v(x) * delta_rho * g * H^2 / (mu_n * phi_bar * Q * L)
# with k_v(x) = 1e-12 * (1 + x/1000) over x in [0, 1000] (so k_v doubles across
# the domain). The constant-k_v base value (k_v = 1e-12) is gamma0 = 0.005886
# (see ve_auxkernels.i), so gamma(x) = gamma0 * (1 + x/1000). On the 10-element
# mesh (centroids 50..950):
#   gamma_min = gamma0 * (1 + 50/1000)  = 0.0061803   (x_c = 50)
#   gamma_max = gamma0 * (1 + 950/1000) = 0.0114777    (x_c = 950)
#   gamma_avg = gamma0 * (1 + 500/1000) = 0.008829     (mean centroid x = 500)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 1000
[]

[Problem]
  solve = false
[]

[Variables]
  [dummy]
  []
[]

[Functions]
  [kv_func]
    type = ParsedFunction
    expression = '1.0e-12 * (1 + x/1000)'
  []
[]

[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
    initial_condition = -20.0
  []
  [kv_field]
    order = CONSTANT
    family = MONOMIAL
  []
  [gamma_ve]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [kv_aux]
    type = FunctionAux
    variable = kv_field
    function = kv_func
    execute_on = 'INITIAL'
  []
  [gravity_number]
    type = VEGravityNumberAux
    variable = gamma_ve
    k_v = kv_field            # spatially-varying coupled k_v
    delta_rho = 300.0
    gravity = '0 0 -9.81'
    mu_n = 1.0e-3
    Q = 1.0e-3
    L = 1000.0
    execute_on = 'INITIAL'
  []
[]

[Materials]
  [geometry]
    type = VEGeometry
    z_top = z_top
    z_bottom = z_bottom
  []
  [porosity]
    type = VEPorosity
    phi_bar = 0.2
  []
[]

[Executioner]
  type = Transient
  end_time = 1.0
  dt = 1.0
[]

[Postprocessors]
  [gamma_avg]
    type = ElementAverageValue
    variable = gamma_ve
    execute_on = 'INITIAL'
  []
  [gamma_min]
    type = ElementExtremeValue
    variable = gamma_ve
    value_type = min
    execute_on = 'INITIAL'
  []
  [gamma_max]
    type = ElementExtremeValue
    variable = gamma_ve
    value_type = max
    execute_on = 'INITIAL'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL'
  []
[]
