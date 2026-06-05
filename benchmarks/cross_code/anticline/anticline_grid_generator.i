# Cross-code benchmark, Tier A -- conceptual anticline grid generator.
#
# Writes a 2D top-surface Exodus mesh standing in for the shared VE grid exchanged
# between MRST co2lab and wombat (verification-ladder step 3; see
# doc/content/benchmarks/cross_code.md). The SAME analytic geometry is reproduced
# in MRST so the two codes run on an identical grid -- the non-negotiable rule for a
# cross-code comparison.
#
# Geometry: a regionally tilted plane (updip = +x) carrying a Gaussian dome (a
# structural culmination near (7 km, 5 km)) that forms the trap the plume migrates
# into:
#   z_top(x,y) = -800 + tan(1.5 deg)*x + 40*exp(-((x-7000)^2+(y-5000)^2)/(2*1200^2))
#                tan(1.5 deg) = 0.0261859
#   H          = 50 m (constant)               -> grad(H) = 0
#   phi_bar    = 0.20
#   K_up       = 200 mD isotropic = 1.9738e-13 m^2   (200 * 9.869e-16)
#
# FIELD LAYOUT (matters -- cf. test/tests/exodus and HANDOVER):
#   z_top, H : NODAL  (LAGRANGE FIRST) so the reader gets non-zero grad(z_T);
#              the buoyancy drive rho_c*g*grad(z_T) is the whole point of the test.
#   phi_bar, K_up_xx/yy/xy : ELEMENTAL (CONSTANT MONOMIAL), as upscaled cell data.
#
# Mesh is 50x50 here for a fast generator; the full benchmark uses 100x100 (+ a
# 200x200 refinement run) -- bump nx/ny and regenerate. Read back by anticline_a1.i.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = 0
  xmax = 10000
  ymin = 0
  ymax = 10000
[]

[Problem]
  solve = false
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  # --- Nodal geometry (LAGRANGE FIRST -> non-zero gradient in the reader) ---
  [z_top]
    order = FIRST
    family = LAGRANGE
  []
  [H]
    order = FIRST
    family = LAGRANGE
  []

  # --- Elemental upscaled petrophysics (CONSTANT MONOMIAL) ---
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
[]

[Functions]
  [ztop_fn]
    type = ParsedFunction
    expression = '-800 + 0.0261859*x + 40*exp(-((x-7000)^2 + (y-5000)^2)/(2*1200^2))'
  []
  [H_fn]
    type = ParsedFunction
    expression = '50'
  []
  [phi_fn]
    type = ParsedFunction
    expression = '0.20'
  []
  [kxx_fn]
    type = ParsedFunction
    expression = '1.9738e-13'
  []
  [kyy_fn]
    type = ParsedFunction
    expression = '1.9738e-13'
  []
  [kxy_fn]
    type = ParsedFunction
    expression = '0'
  []
[]

[AuxKernels]
  [ztop_aux]
    type = FunctionAux
    variable = z_top
    function = ztop_fn
    execute_on = 'INITIAL'
  []
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
