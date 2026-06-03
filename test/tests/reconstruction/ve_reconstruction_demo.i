# Vertical-reconstruction DEMO -- 2D VE injection that carries the fields the
# scripts/ve_reconstruct_vertical.py post-processor needs to rebuild the 3D plume.
#
# Physics is the known-good 2D updip injection from mesh/dip/dip_injection.i: CO2
# injected through the lower-left (downdip) edge of a 100 x 50 m aquifer that dips in
# x and y, migrating updip as a buoyant tongue. The only additions here are the
# diagnostic AuxVariables the reconstruction reads:
#
#   z_top, z_bottom  (already present) -- column geometry; H = z_top - z_bottom = 20 m
#   h_plume          -- mobile-plume thickness ve_h (VEPlumeReconstruction + VEPlumeHeightAux)
#   rho_n, rho_w     -- the two phase densities (ADMaterialStdVectorAux on ve_density,
#                       index 0 = non-wetting/CO2, 1 = wetting/brine). Output BOTH
#                       densities, not a precomputed delta_rho, so the file stays
#                       self-documenting; the script forms delta_rho = rho_w - rho_n.
#
# Exodus writes every 25 steps -> ~20 frames of a growing plume.
#
# ----------------------------------------------------------------------------------
# Reconstruct the 3D vertical saturation AFTER the run (needs SEACAS exodus.py):
#
#   # sharp interface (box column; mass-consistent with this sharp run):
#   scripts/ve_reconstruct_vertical.py ve_reconstruction_demo_out.e plume3d_sharp.e \
#       --mode sharp --s-wr 0.0 --n-layers 40 --vertical up --sat-min 1e-2 \
#       --z-top-var z_top --z-bottom-var z_bottom --h-var h_plume
#
#   # capillary fringe (graded column; supply an Sw(Pc) curve -- here analytic VG):
#   scripts/ve_reconstruct_vertical.py ve_reconstruction_demo_out.e plume3d_fringe.e \
#       --mode fringe --pc vg --vg-alpha 1.0e-3 --vg-m 0.5 --s-wr 0.0 \
#       --n-layers 40 --vertical up --sat-min 1e-2 \
#       --z-top-var z_top --z-bottom-var z_bottom --h-var h_plume \
#       --rho-n-var rho_n --rho-w-var rho_w
#
#   # drape onto the ORIGINAL 3D geology instead of an extruded column mesh -- exact
#   # faults/layering/topography, with the mesh's porosity/permeability carried through:
#   scripts/ve_reconstruct_vertical.py ve_reconstruction_demo_out.e plume3d_orig.e \
#       --original-mesh field_model.e --mode sharp --s-wr 0.0 --vertical up --sat-min 1e-2 \
#       --z-top-var z_top --z-bottom-var z_bottom --h-var h_plume
#
#   field_model.e is the em2ex 3D grid the 2D run was upscaled from (via ex2ve). This
#   synthetic GeneratedMesh demo has no such original grid, so the command above is shown
#   for the field workflow; on a real case each original node is mapped to its nearest VE
#   column and sat_co2 is evaluated at the node's true depth below z_top (add --no-carry
#   to write saturation only, without the porosity/permeability fields).
#
#   --sat-min clips the plume FOOTPRINT: continuous-Galerkin FE advection leaves a thin
#   diffuse halo of tiny sat_n ahead of the front, which would otherwise reconstruct as
#   a faint sheet over the whole top surface (worst in sharp mode, where any h > 0 paints
#   the top layer at 1 - S_wr). --sat-min 1e-2 (= ~5% of the peak sat_n here) renders only
#   columns with a real plume. Default --sat-min 0 keeps everything (faithful, but haloed).
#
#   The plume here is thin (h ~ 2 m in a 20 m column), so the buoyancy pressure
#   Pc = delta_rho*g*zeta is small; vg-alpha = 1e-3 (low entry pressure) keeps the
#   fringe legible (max sat_co2 ~ 0.84). A realistic, higher entry pressure
#   (vg-alpha ~ 1e-5) gives a fainter grade (max ~ 0.06) -- correct VE physics for a
#   thin plume, just less to look at.
#
# Open plume3d_*.e in ParaView and colour by sat_co2 (co2_region flags 0 brine /
# 1 mobile / 2 trapped). The sharp box and the fringe grade share the same plume
# SHAPE h(x,y); the fringe additionally shows the smooth vertical saturation profile.
#
# NOTE: this run uses sharp-interface dynamics, so the sharp reconstruction is the
# mass-consistent one. The VG fringe command overlays an equilibrium fringe of the
# same thickness for visualisation; for an end-to-end fringe-consistent result, switch
# VEPlumeReconstruction to mode = capillary_fringe with a matching pc_uo and reconstruct
# with the SAME curve.
# ----------------------------------------------------------------------------------

[Mesh]
  [base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 50
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 50
  []
  [injection_zone]
    type = ParsedGenerateSideset
    input = base_mesh
    combinatorial_geometry = 'x < 1e-6 & y <= 5.0'
    normal = '-1 0 0'
    new_sideset_name = 'injection_zone'
  []
[]

[GlobalParams]
  VEDictator = ve_dictator
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 5.0e-5
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1020.0
    viscosity = 8.0e-4
  []
[]

[UserObjects]
  [ve_dictator]
    type = VEDictator
    porous_flow_vars = 'pp_top sat_n'
    ve_flavour = sharp_interface
  []
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '0.034899 * x + 0.017450 * y'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * x + 0.017450 * y - 20'
  []
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - x)'
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
  []
[]

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_init
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.0
  []
  [z_top_ic]
    type = FunctionIC
    variable = z_top
    function = z_top_func
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = z_bottom_func
  []
[]

[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
  []
  # Reconstruction inputs (elemental: they read qp material properties).
  [h_plume]
    order = CONSTANT
    family = MONOMIAL
  []
  [rho_n]
    order = CONSTANT
    family = MONOMIAL
  []
  [rho_w]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [co2_storage]
    type = VEMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
  []
  [co2_flux]
    type = VEAdvectiveFluxS
    variable = sat_n
    fluid_phase = 0
    pp_top = pp_top
  []
  [brine_storage]
    type = VEMassTimeDerivative
    variable = pp_top
    fluid_phase = 1
  []
  [brine_flux]
    type = VEAdvectiveFluxP
    variable = pp_top
    fluid_phase = 1
  []
[]

[AuxKernels]
  [plume_height]
    type = VEPlumeHeightAux
    variable = h_plume
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rho_n_aux]
    type = ADMaterialStdVectorAux
    variable = rho_n
    property = ve_density
    index = 0 # non-wetting (CO2)
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rho_w_aux]
    type = ADMaterialStdVectorAux
    variable = rho_w
    property = ve_density
    index = 1 # wetting (brine)
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[BCs]
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [co2_injection]
    type = NeumannBC
    variable = sat_n
    boundary = injection_zone
    value = 2.8e-3
  []
  [s_right]
    type = DirichletBC
    variable = sat_n
    boundary = right
    value = 0.0
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
  [permeability]
    type = VEPermeability
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
  []
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
  [plume_reconstruction]
    type = VEPlumeReconstruction
    sat_n = sat_n
    S_wr = 0.0
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
[]

[Executioner]
  type = Transient
  end_time = 5.0e6
  dt = 1.0e4

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

[Postprocessors]
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_plume_height]
    type = ElementExtremeValue
    variable = h_plume
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  # ~20 frames of the growing plume.
  exodus = true
  time_step_interval = 25
[]
