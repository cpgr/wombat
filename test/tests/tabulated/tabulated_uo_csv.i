# data_file variant of tabulated_uo.i: the table UOs read the ve_upscale_curves.py
# CSV output (curves_relperm.csv, curves_pc.csv) DIRECTLY via data_file, instead of
# !include'ing the generated [UserObjects] vector block. Same sampling and the same
# independent numpy gold (gold/tabulated_uo_csv_out.csv == gold/tabulated_uo_out.csv),
# so this proves the data_file path reads/interpolates identically to the inline path.
# The script-regression test (prereq) writes the two CSVs.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  solve = false
[]

[UserObjects]
  [relperm_uo]
    type = VERelPermTableUO
    data_file = curves_relperm.csv # columns sat_n, kr_n_up, kr_w_up
  []
  [pc_uo]
    type = VECapillaryPressureTableUO
    data_file = curves_pc.csv # columns pc, sw
    sat_lr = 0.2
  []
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 1.0e-3
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1000.0 # delta_rho = 300 kg/m3 (matches tabulated_gold.py)
    viscosity = 1.0e-3
  []
[]

[Functions]
  [sat_n_drive]
    type = ParsedFunction
    expression = '0.1 * t' # t = 1..6 -> sat_n = 0.1..0.6
  []
[]

[AuxVariables]
  [sat_n]
    family = MONOMIAL
    order = CONSTANT
  []
  [pp_top]
    family = MONOMIAL
    order = CONSTANT
  []
  [z_top]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
    initial_condition = -20.0 # H = 20 m
  []
  [krn]
    family = MONOMIAL
    order = CONSTANT
  []
  [krw]
    family = MONOMIAL
    order = CONSTANT
  []
  [ve_h]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e6
  []
[]

[AuxKernels]
  [sat_n_drive_aux]
    type = FunctionAux
    variable = sat_n
    function = sat_n_drive
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [krn_aux]
    type = ADMaterialStdVectorAux
    variable = krn
    property = ve_relperm
    index = 0
    execute_on = 'TIMESTEP_END'
  []
  [krw_aux]
    type = ADMaterialStdVectorAux
    variable = krw
    property = ve_relperm
    index = 1
    execute_on = 'TIMESTEP_END'
  []
  [ve_h_aux]
    type = VEPlumeHeightAux
    variable = ve_h
    execute_on = 'TIMESTEP_END'
  []
[]

[Materials]
  [geometry]
    type = VEGeometry
    z_top = z_top
    z_bottom = z_bottom
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
    mode = capillary_fringe
    sat_n = sat_n
    S_wr = 0.2
    pc_uo = pc_uo
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 6.0
  dt = 1.0
[]

[Postprocessors]
  [krn]
    type = ElementAverageValue
    variable = krn
    execute_on = 'TIMESTEP_END'
  []
  [krw]
    type = ElementAverageValue
    variable = krw
    execute_on = 'TIMESTEP_END'
  []
  [ve_h]
    type = ElementAverageValue
    variable = ve_h
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  execute_on = 'TIMESTEP_END' # one CSV row per step (t = 1..6); no INITIAL row
  csv = true
[]
