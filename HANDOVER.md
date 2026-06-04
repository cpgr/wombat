# Session handover — VE CCS app

Working note to bootstrap the next coding session. Read this plus `CLAUDE.md`
(project memory). Everything below is committed on `main`; the tree is clean
except the intentionally-untracked `test/tests/mesh/topography_demo.i` and this
note.

## Where things stand (all on `main`, full suite passing)

**Solver paths**
- FE: `VEMassTimeDerivative` + `VEAdvectiveFluxS/P`.
- FV: `VEFVMassTimeDerivative` + unified `VEFVAdvectiveFlux` (boundary-aware
  functor geometry/pressure; functor-upwinded mobility; `advected_interp_method`
  supports `upwind` (default) and the TVD limiters `vanLeer`/`min_mod`/`sou`/
  `quick`/`venkatakrishnan`).

**Relative permeability — swappable UserObject** (the refactor that pays off)
- `VERelativePermeability` (base UO): `relativePermeability(sat_n, phase)`,
  `dRelativePermeability(...)`, and `relativePermeabilityAD()` (seeds the AD
  sat_n derivative once).
- `VERelPermSharpUO` (analytical sharp), `VERelPermTableUO` (tabulated,
  `LinearInterpolation`).
- Adapters: `VERelPerm` (FE qp material) and `VEFVRelPerm` (FV functor material).
  A new relperm model = one new UO subclass, no kernel/adapter changes.

**Capillary pressure — now in the solve (two-pressure VE)**

This session wired `Pc^up` into the flux as an opt-in capillary-diffusion term:

- `ve_h` in `VEPlumeReconstruction` is now an `ADMaterialProperty<Real>`:
  - sharp_interface mode: natural AD propagation (`h = sat_n * H / (1-S_wr)`).
  - capillary_fringe mode: Newton inversion on raw values, derivative seeded
    manually via inverse function theorem: `_h[qp] = h_val + (H/Sn_h) * (sat_n - raw_sat)`.
  - `VEPlumeHeightAux` updated to read `ADMaterialProperty<Real>` + `raw_value()`.

- `VEUpscaledCapPressure` upgraded:
  - `ve_pc_up` is now an `ADMaterialProperty<Real>` (Jacobian flows through
    density, ready for future pressure-dependent EOS).
  - New `ve_dPcup_dsatn` (`MaterialProperty<Real>`) =
    `(rho_w - rho_n) * g * H / (1 - S_wr)` — the FE kernel gradient coefficient.
  - New `S_wr` parameter (optional, default 0.0; must match VEPlumeReconstruction).
  - `ve_auxkernels.i` updated: `MaterialRealAux` → `ADMaterialRealAux` for `ve_pc_up`.

- **New `VEFVCapPressure`** (`FunctorMaterial`): declares the `ve_pc_up` functor
  for FV kernels using `sat_n`, `z_top`, `z_bottom` functors and a `delta_rho`
  parameter. `VEFVAdvectiveFlux` calls `.gradient(face, state) * _normal` on it.

- `VEFVAdvectiveFlux`: new `capillary = false` param (bool). When true and
  `fluid_phase == 0`: `dphi_dn += _pc_up->gradient(face, state) * _normal`.
  Default OFF — all existing sharp-interface golds are byte-identical.

- `VEAdvectiveFluxS`: new `capillary = false` param (bool). When true:
  `potential_gradient += (*_dPcup_dsatn)[qp] * _grad_u[qp]`
  (= `d(Pc^up)/dS * grad(sat_n)`, with `_grad_u[qp]` carrying the diffusion
  Jacobian automatically via AD). Default OFF.

**Verification tests (CSVDiff, shortened via cli_args)**
- `buckley_leverett`: `bl_flat`, `bl_flat_fv`, `bl_flat_fv_vanleer`,
  `bl_flat_table`.
- `nordbotten_celia`: `nc_sloped`, `nc_sloped_fv`, `nc_sloped_fv_vanleer`.
- `mesh/dip`: `dip` (gradient), `dip_injection{,_x,_y}` (FE) and `_fv` variants,
  `dip_injection_fv_vanleer`.
- `auxkernels/ve_auxkernels`: diagnostics — plume height, gravity number, and
  `Pc^up` (= 44145 Pa for the uniform field).
- `plume_reconstruction`, `postprocessors`.
- `buckley_leverett/bl_flat_capillary_fe_jacobian`,
  `buckley_leverett/bl_flat_capillary_fv_jacobian` — `PetscJacobianTester` tests
  for the capillary-term Jacobian (ratio_tol = 1e-7).
- `mcwhorter_sunada/mcwhorter_fe`, `mcwhorter_sunada/mcwhorter_fv` — McWhorter-
  Sunada semi-analytical CSVDiff (capillary sign + magnitude). **All pass.**

**Recent commits**
```
226bd68 Add McWhorter-Sunada capillary verification; fix FV capillary drive
d7ef756 Replace Unicode arrows with ASCII in VECapillaryPressureTable comments
d54973f Wire upscaled capillary pressure into the VE flux (two-pressure VE)
7cac172 Add VEUpscaledCapPressure diagnostic material
a7225e8 Add tabulated relative permeability UserObject
ac45dd6 Add FV dip tests and TVD-limiter regression tests
```

## Capillary term — COMPLETE (committed, all tests green)

The two-pressure VE capillary term is wired, AD-consistent (Jacobian tests), and
physically verified (McWhorter-Sunada). Two gotchas resolved this session, both
worth remembering:

**Gotcha A (relperm clamp kink — Jacobian tests):** the `PetscJacobianTester`s *must*
be evaluated at a strictly-interior uniform `sat_n` (the spec overrides IC and BC to 0.5 via
cli_args). The sharp relperm derivative `dRelativePermeability` returns 0 at the
clamp boundaries (sat_n = 0 and the s_eff = 1 point); at sat_n = 0 (the physical
IC) every interior node sits on that kink, AD reports slope 0, the one-sided FD
steps into the active region and reports 1/(1-S_wr) -> a ~21% Frobenius mismatch
that has NOTHING to do with the capillary term. Evaluating at sat_n = 0.5 makes
s_eff interior, AD = FD. The capillary term (C*grad(sat_n), C = ve_dPcup_dsatn
constant) is linear and fully AD-tracked, so its Jacobian is exact; the brine
flux and mass storage kernels are pure-AD too. Do not "fix" the capillary code
in response to a kink-state Jacobian failure -- fix the test state.

What a Jacobian test does NOT prove: sign/magnitude of the residual (a wrong-sign
Pc^up is still self-consistent). That is what McWhorter-Sunada checks — and it did:

**Gotcha B (FV capillary drive silently zero — found by McWhorter):** the FV kernel
originally took the gradient of the `ve_pc_up` MATERIAL functor. A material functor's
Green-Gauss gradient does NOT pick up the `sat_n` Dirichlet BC at a boundary face, so
the inlet capillary drive was zero, CO2 never imbibed, and `mcwhorter_fv` measured
`||sat_n - S_ref||` with `sat_n` stuck at the IC. The FV Jacobian test missed it
entirely (it runs at uniform `sat_n`, where the term vanishes either way) — a sharp
reminder to pair every Jacobian test with a physics test. Fix: `VEFVCapPressure`
publishes the coefficient `ve_dPcup_dsatn`, and the kernel forms
`grad(Pc^up).n = ve_dPcup_dsatn(face) * gradUDotNormal()` using the sat_n VARIABLE
gradient (boundary-aware, the FVDiffusion idiom). Mirrors the FE kernel exactly.

### McWhorter-Sunada test design (for reference / extension)

`test/tests/mcwhorter_sunada/` holds the verification. The mapping that makes it
exact: constant per-phase densities + flat formation + closed far boundary => total
volume flux q_t = 0 => pure counter-current imbibition => the CO2 equation collapses
to `phi dS/dt = d/dx(D(S) dS/dx)`, `D(S) = (K*C/mu)*S*(1-S)`, `C = (rho_w-rho_n)*g*H`.
`mcwhorter_reference.py` solves the Boltzmann similarity ODE (`scipy.solve_bvp`) and
writes `mcwhorter_reference.csv` (committed input data; read via `data_file` with
`axis = x` — PiecewiseLinear defaults to a function of TIME, a trap that cost a
debug cycle). `ElementL2Error` (`co2_l2_error`, gold 0, `abs_zero = 0.03`) checks the
match; observed ~2e-3, cross-checked by an independent method-of-lines solve. A wrong
capillary sign gives anti-diffusion, L2 ~ O(0.3). Regenerate the reference after any
parameter change: `python3 mcwhorter_reference.py`.

## grad(H) capillary term -- COMPLETE (this session)

The `sat_n*grad(H)` chain-rule part of the capillary drive is now in both solver
paths. Previously the kernels used only `C*grad(sat_n)` (`C = d(Pc^up)/d(sat_n)`),
exact only for constant H. The full gradient is

    grad(Pc^up) = d(Pc^up)/d(sat_n) * grad(sat_n) + d(Pc^up)/d(H) * grad(H)
                = delta_rho*g*grad(h)    (sharp interface, h = sat_n*H/(1-S_wr))

- `VEGeometry` now declares **`ve_grad_H`** (`MaterialProperty<RealGradient>`)
  `= grad(z_T) - grad(z_B)`. z_bottom must be LAGRANGE FIRST for a non-zero gradient
  (same caveat as z_top).
- `VEUpscaledCapPressure` (FE) declares **`ve_dPcup_dH`** (`ADMaterialProperty<Real>`)
  `= delta_rho*g*ve_h/H` (= sharp `delta_rho*g*sat_n/(1-S_wr)`). It is **AD** (unlike the
  plain-Real `ve_dPcup_dsatn`) because grad(H) is fixed geometry and carries no
  Jacobian -- so the term's d/d(sat_n) must ride on this coefficient.
- `VEFVCapPressure` (FV) declares the **`ve_dPcup_dH`** functor `= delta_rho*g*sat_n(face)/(1-S_wr)`.
- `VEAdvectiveFluxS`: `potential_gradient += ve_dPcup_dH * ve_grad_H`.
- `VEFVAdvectiveFlux`: `dphi_dn += ve_dPcup_dH(face) * grad(H).n`, where
  `grad(H).n = (z_top.gradient - z_bottom.gradient).n` (z_top/z_bottom functor
  gradients -- fixed geometry, BC-safe; NOT a material-functor gradient, so Gotcha
  B does not apply). The boundary-aware sat_n enters via the `ve_dPcup_dH(face)` VALUE.
- **No new input parameter.** Folded into the existing `capillary=true`. grad(H)=0 for
  constant thickness, so all flat/constant-H golds (bl_flat_capillary, McWhorter)
  are byte-identical.

**Verification (`test/tests/capillary_equilibrium/`):** no-flow capillary
equilibrium on a laterally varying thickness. Flat top (kills buoyancy), H(x)=10+0.2x,
IC set to the analytic equilibrium `sat_n*H = const` (sat_n in [0.2,0.6], interior).
`grad(Pc^up)=delta_rho*g*grad(const)=0` is a fixed point ONLY with the grad(H) term;
without it the kernel drives sat_n toward uniform and the L2 drift balloons -- so the
test checks the SIGN and MAGNITUDE the Jacobian test cannot. `cap_equilibrium_{fe,fv}`
(CSVDiff, `satn_drift` gold 0, `abs_zero=1e-3`) + `cap_equilibrium_{fe,fv}_jacobian`
(PetscJacobianTester at the non-uniform interior state, both capillary terms active,
no relperm kink). **All four pass.**

Observed drift: FE 1.3e-4 (under tol), FV 3.6e-16 (machine zero). The FV scheme
preserves this equilibrium EXACTLY: with arithmetic-mean face interpolation the
discrete product rule  b_bar*(a_N-a_P) + a_bar*(b_N-b_P) == a_N*b_N - a_P*b_P  makes
the two chain-rule halves telescope, H_face*grad(sat_n).n + sat_n_face*grad(H).n =
grad(sat_n*H).n = grad(6).n = 0. FE interpolates sat_n=C/H linearly inside each
element (true profile is 1/H), so sat_n*H at interior qps isn't exactly C -- a small
quadratic interpolation residual (the 1.3e-4, shrinks ~h^2 on refinement). Both
correct; FV is exact for this case by construction.

**Gotcha C (FV upwind-switch kink at zero velocity -- Jacobian tests):** the FV
Jacobian test FAILED at first (||J-Jfd||/||J|| ~ 9.7e-6, > ratio_tol 1e-7) while the
FE Jacobian passed at the IDENTICAL state. Cause: at the no-flow fixed point the face
Darcy velocity is ~0 everywhere, sitting exactly on the raw_value(velocity) > 0 branch
that picks the upwind cell for the advected relperm. An FD perturbation flips that
branch -- a finite jump in the relperm source -- which AD does not (and should not)
track. Nothing to do with the grad(H) term (FE has no upwind switch -> clean at the
same state). Fix: run the FV Jacobian test with advected_interp_method=average
(central, no upwind branch, fully differentiable); the capillary term is still active
and tested. Same lesson as Gotcha A: evaluate the Jacobian away from a kink; don't
"fix" correct code in response to a kink-state mismatch.

Housekeeping: the `*_jacobian` tests share the input file with the physics tests, so
they carry `Outputs/file_base=cap_equilibrium_{fe,fv}_jac` in cli_args -- their short
end_time=2000 run writes to a separate base instead of clobbering the physics
`*_csv.csv` (and it avoids a parallel-run I/O race on the shared output).

## Capillary-fringe flux path -- DONE (this session)

`capillary = capillary_fringe` now drives the two-pressure VE flux in both FE and FV
(previously sharp_interface only; the fringe option `paramError`ed). The unifying idea:
in BOTH modes `grad(Pc^up) = delta_rho*g*grad(h)`, so the flux coefficients are just the
chain-rule partials of `h(sat_n, H)`:

    ve_dPcup_dsatn = delta_rho*g * dh/dsat_n
    ve_dPcup_dH    = delta_rho*g * dh/dH

with `dh/dsat_n = H/S_n(h)`, `dh/dH = sat_n/S_n(h)`, and `S_n(h) = 1 - Sw(delta_rho*g*h)`
the CO2 saturation at the plume top. Sharp interface is the special case `S_n = 1 - S_wr`
(constant), which reproduces the old `H/(1-S_wr)` / `sat_n/(1-S_wr)` formulae exactly.

What changed:
- `VEUpscaledCapPressure` (FE) and `VEFVCapPressure` (FV) gained a `mode` enum + optional
  `pc_uo`. In fringe mode they evaluate `S_n(h)` from the Pc UO. `ve_h` is the AD plume
  thickness from `VEPlumeReconstruction` (FE) or an inline Newton inversion (FV, mirroring
  VEPlumeReconstruction's loop + inverse-function-theorem AD seed).
- **`ve_dPcup_dsatn` is now an `ADMaterialProperty<Real>`** (was plain `Real`), and
  `VEAdvectiveFluxS` reads it with `getADMaterialProperty`. In sharp mode the AD coefficient
  carries zero sat_n-derivative (constant density), so the `ADcoeff * _grad_u` product is
  byte-identical to the old `Real * _grad_u` -- all sharp golds unchanged. In fringe mode the
  coefficient depends on sat_n through `S_n(h)`, and its AD derivative (propagated from AD
  `ve_h` via `PorousFlowCapillaryPressure::saturation(const ADReal&)`, which chains through
  `dSw/dPc`) makes the grad(sat_n) Jacobian exact. FV is AD-consistent for free (the
  coefficient functor VALUE carries d/dsat_n and multiplies the AD `gradUDotNormal`).
- `VEFlowPhysicsBase`: new `pc_uo` param (required when fringe), `assignCapillaryMode()`
  helper threads `mode`+`pc_uo` to the downstream materials; the FE wiring also couples
  `sat_n` into `VEUpscaledCapPressure` (needed for `dh/dH = sat_n/S_n`).

Tests (`test/tests/physics/`): `ve_flow_fe_capillary_fringe` (PetscJacobianTester at a
NON-uniform cap_equilibrium-style interior state -- both grad(sat_n) and grad(H) fringe terms
active; FE has no upwind kink so the no-flow geometry is safe) and `ve_flow_fv_capillary_fringe`
(flowing uniform state: a pressure gradient keeps the Darcy velocity non-zero so the FV upwind
branch is differentiable -- the Physics does not expose `advected_interp_method=average`, so a
no-flow FV state would trip Gotcha C). Negative `err_capillary_fringe_requires_pc_uo`.

Both fringe Jacobian tests use a **LINEAR** `VECapillaryPressureTableUO` (the
`ve_plume_recon_cf` table): `dSw/dPc` is piecewise-constant and exact, so `S_n(h)` and its
derivative are mutually consistent and the strict `ratio_tol = 1e-7` holds. (The FV fringe
test keeps `ratio_tol = 1e-7` but uses `difference_tol = 1e-3` -- the per-face Newton
inversion makes the absolute J-vs-Jfd entries slightly noisier than the FE qp path.) A *nonlinear* table would
amplify the value-vs-analytic-derivative inconsistency (Gotcha-D family) -- prefer a fine
piecewise-linear table or relax `ratio_tol` if a curved table is ever Jacobian-tested. The
near-zero floor on `S_n` (plume nose, `dh/dsat_n -> inf`) is a kink; keep test states off it
(thin but non-vanishing plumes, `S_n(h) ~ 0.3` in these inputs).

FOLLOW-UPS (both DONE this session; verified by the user's `./run_tests`):
- **Fringe IN-SOLVE physics verification (sign/magnitude) -- DONE.** McWhorter-fringe FE+FV in
  `test/tests/mcwhorter_sunada/` (`mcwhorter_fringe_{fe,fv}.i`). Same counter-current
  imbibition as the sharp McWhorter but the upscaled-Pc slope is the fringe one, so the
  effective diffusion is `D(S) = D_sharp(S)/S_n(h(S))` (column imbibes FASTER). The ONLY change
  from the sharp derivation is `dPc/dS`, so the proven sharp framework carries over. Reference
  `mcwhorter_fringe_reference.py` (Boltzmann similarity ODE with the fringe `D`, the table
  inversion for `h(S)`) -> `mcwhorter_fringe_reference.csv`; **cross-checked against an
  independent method-of-lines solve to L2 ~ 1e-7**. `H = 20 m` keeps `S_n(h) ~ 0.6-0.7` (away
  from the `D ~ sqrt(S)` degeneracy at `S -> 0`, which otherwise destabilises both the MOL and
  the solve). `co2_l2_error` gold 0 with **`abs_zero = 0.012`**: the fringe-vs-sharp-coefficient
  separation is L2 ~ 0.022, so the test discriminates the `S_n(h)` factor (a wrong sign blows
  up; the sharp coefficient fails). Keep the `.i` params in sync with the reference `.py`.
- **FV grad(H) fringe Jacobian -- DONE.** `VEFlowFV` now exposes `advected_interp_method`
  (forwarded to both `VEFVAdvectiveFlux` kernels; default `upwind`, so existing golds are
  unchanged). `test/tests/physics/ve_flow_fv_capillary_fringe_gradh.i` uses
  `advected_interp_method = average` on the no-flow varying-H equilibrium so the
  zero-velocity upwind kink (Gotcha C) is avoided and the FV grad(H) fringe term
  `ve_dPcup_dH = delta_rho*g*sat_n/S_n(h)` is Jacobian-tested.
- STILL OPEN: FV fringe runs a per-face Newton inversion inside the functor (CLAUDE.md
  key-subtlety 8); profile on a real mesh.

## Hysteresis / residual trapping -- v1 ESSENTIALLY DONE (Phases 0-4 committed)

Saturation-space Land/Killough, kr+trapping only (NO Pc hysteresis in v1). Most
important remaining physics for the post-injection migration regime VE targets.
All 33 tests pass. Commits: d6e901c, ed5e3e5, fcb5269, 25586d0, 3b683b0, 6637b87.

DONE (committed):
- `VERelativePermeability` history-aware overloads `relativePermeability(sat_n,
  sat_n_max, phase)`, `dRelativePermeability(...)`, `trappedSaturation(...)`, 3-arg
  `relativePermeabilityAD`. Base defaults IGNORE history -> existing UOs unchanged.
- `VESaturationMaxAux`: stateful `sat_n_max <- max(sat_n, uOld)`, `execute_on=
  TIMESTEP_END` (lagged/explicit). The sat_n_max VARIABLE must carry an
  `initial_condition` = initial sat_n (set it on the variable block, NOT a separate
  [ICs] ConstantIC -- that did not reliably apply to an aux-with-kernel and left
  sat_n_max=0 -> drainage branch -> no trapping). Family matches discretisation.
- `VERelPermHysteresisUO` (param `S_gr_max`, Land `C=(1-S_wr)/S_gr_max-1`):
  drainage `kr_n=krn_max*s` while `sat_n>=sat_n_max`; imbibition linear scanning
  `kr_n=krn_max*s_max*(s-s_gr)/(s_max-s_gr)`, `s_gr=s_max/(1+C*s_max)`.
  `trappedSaturation=(1-S_wr)*s_gr*(s_max-s)/(s_max-s_gr)` -- chosen so
  `sat_n - sat_n_trap` == kr-implied mobile sat (continuous at reversal, -> Land
  residual at full imbibition). Reduces EXACTLY to VERelPermSharpUO under monotonic
  drainage. All three scanning methods GUARD s_max->0 (denom<=1e-12) -> drainage /
  zero-trapping, else 0/0 NaN when sat_n_max lags at its zero IC during drainage.
- Adapters `VERelPerm`/`VEFVRelPerm`: OPTIONAL `sat_n_max` coupling, passed FROZEN
  (AD seed stays wrt sat_n). Absent -> 2-arg drainage path -> golds byte-identical.
- `VETrappedSaturationAux` writes `sat_n_trap`; mass PPs already consume it (no
  change): `VEMobileCO2MassPostprocessor` (sat_n - sat_n_trap),
  `VETrappedCO2MassPostprocessor` (sat_n_trap).
- Kernels UNCHANGED: storage carries total sat_n (incl. trapped); flux immobilises
  trapped CO2 via scanning kr_n -> 0.
- Tests (`test/tests/hysteresis/`): `sat_n_max_tracking` (running-max),
  `column_cycle` (analytic Land cycle, solve=false), `hysteresis_jacobian`
  (PetscJacobianTester through UO->adapter->kernel on the imbibition branch, away
  from the turning-point/clamp/trapped kinks -- cf. Gotcha A, Gotcha C),
  `hysteresis_drainage_fv` (FV adapter VEFVRelPerm+sat_n_max in a real transient
  solve; monotonic drainage reduces to sharp BL, co2_volume=5e-4*t).

TODO (Phase 4b -- real-flow IMBIBITION-displacement demo, deferred; HARD):
- Goal: a transient flow solve where hysteretic kr immobilises CO2 during imbibition
  (the one thing not yet shown in a solve -- column_cycle covers the algebra, the FV
  drainage test covers the adapter under drainage only).
- Attempted a 1D throughflow drainage->imbibition sweep; it kept hitting GENUINE
  multiphase pathologies, not code bugs:
  (1) Starting fully CO2-saturated (sat_n=1, S_wr=0) -> kr_w=0 everywhere ->
      DEGENERATE pressure equation -> blow-up. Need both phases mobile at t=0.
  (2) FE continuous-Galerkin advection of the nonlinear hysteretic fractional flow
      is UNSTABLE (sat_n oscillates to +/-large). Use the FV path.
  (3) FV upwind is bounded, but injecting brine to displace mobile CO2 forms a CO2
      BANK that reaches max saturation where kr_w=0 -> brine blocked -> displacement
      LOCKS UP (sat_n froze ~0.6-1.0, inlet 0.25 never propagated). Plus starting the
      whole domain exactly at the turning point sat_n=sat_n_max is a derivative kink
      (first dt cut to ~4 s).
  Path forward for a working version: FV + TVD limiter, S_wr > 0 (connate water keeps
  kr_w > 0 so no lock-up), start partially saturated and OFF the turning point, gentle
  ramp; verify trapped vs swept-volume integral / MacMinn-Juanes migration distance.
  Needs live iteration (the user runs ./run_tests; not tunable blind).
- `VEDictator::WITH_HYSTERESIS` flavour enum exists; not yet used for gating.

## Test-coverage audit (done) -- remaining gaps

Audited every registered object vs test inputs. All objects instantiated; PPs, both
plume-reconstruction modes, relperm/capillary/hysteresis verified analytically.
CLOSED earlier (commit dc33c87): VEGeometry thickness mode
(`geometry/thickness_mode`, since REMOVED -- see open-items note) and
anisotropic/off-diagonal K_up
(`permeability/anisotropic_permeability`, manufactured solution, QUAD9/2nd-order).
Remaining gaps, lower priority:
- **capillary_fringe in-solve**: `ve_plume_recon_cf` verifies the Newton-inversion
  VALUE only. The fringe `ve_h` AD seeding (inverse-function-theorem, subtlety 7) is
  never Jacobian-tested (all PetscJacobianTesters run sharp_interface) and fringe is
  never wired into a flux/relperm in any test. Untested the moment it enters a residual.
- **pc_entry** (entry/fringe pressure offset, VEUpscaledCapPressure / VEFVCapPressure):
  always 0 -> untested.
- **Buoyancy term Jacobian**: `rho_c*g*grad(z_T)` contributes 0 to the Jacobian with
  constant density. Now COVERED for the pressure-dependent EOS: `fv_eos_jacobian`
  (sloped top, CO2+brine) and the new `eos_interface_fe_jacobian` (sloped top, interface
  mode) both exercise it.

## Convective dissolution -- v1 DONE (Option A, sink-only; commits 949cc17, 8d4cb79)

Constant-flux convective dissolution, FE + FV. All tests pass.
- `VEDissolution` (material): `ve_dissolution_rate = q0*gate(sat_n)*capacity(c_diss)`
  [kg/m^2/s], AD in sat_n via the gate. q0 = constant convective areal flux (direct
  param for v1); gate = min(1, sat_n/s_ref) prevents over-depletion; optional capacity
  factor (needs dissolved_co2 + c_cap) stops at column saturation.
- `VEDissolutionSink` (FE, ADKernelValue) / `VEFVDissolutionSink` (FV, FVElementalKernel):
  +ve_dissolution_rate on the sat_n equation. Mass units, no rho scaling (matches storage
  residual). Same elemental sign FE/FV. Off by default (separate kernels).
- `VEDissolvedCO2Aux`: stateful areal dissolved mass c_diss += rate*dt (TIMESTEP_END,
  converged rate => mass-conservative with the implicit sink). MUST be elemental (CONSTANT
  MONOMIAL) -- reads the QP rate material. Integrate for total dissolved mass.
- Tests (`test/tests/dissolution/`): constant_rate (FE, single elem + pp Dirichlet both
  nodes -> pure ODE), constant_rate_fv (FV, krn_max=0 -> CO2 immobile -> pure ODE, brine
  fills via open BCs), + Jacobian variants (gate AD-active at s_ref>sat_n). Analytic:
  sat_n->0.4, dissolved->q0*t=140, free CO2 lost = dissolved (conservation).
- Sink-only: dissolved CO2 is an immobile inventory (c_diss), not transported; freed pore
  volume supplied by brine inflow through the pressure solve. The X_co2-based
  VEDissolvedCO2MassPostprocessor is left for a future transported/equilibrium model (B).
- FOLLOW-UPS: (1) compute q0 from the convective correlation alpha*k*delta_rho_c*g/mu*C_sat
  (Pau/Neufeld) instead of a direct input; (2) the capacity-cap path (dissolved_co2 + c_cap)
  is implemented but NOT yet exercised by a test; (3) Option B transported dissolved CO2.

## Exodus real-field input pipeline -- DONE (commit da8d780)

The production input path (read upscaled petrophysics from an Exodus mesh) is now
exercised; previously ALL tests used GeneratedMesh + uniform initial_condition scalars.
`test/tests/exodus/`: a generator writes a 2D Exodus with spatially-varying ELEMENTAL
fields (z_top, z_bottom, phi_bar, K_up_xx/yy/xy; linear so element value == centroid value);
the reader (prereq on generator) reads them via `[Mesh] type=FileMesh` +
`initial_from_file_var` into VEGeometry / VEPorosity / VEPermeability and
checks ve_H (= z_top - z_bottom), ve_phi_bar, ve_K_up at a sample element. Notes for extending:
- No flow solve / VEDictator needed -- `[Problem] solve=false` + MaterialReal[Tensor]Aux
  at INITIAL computes the materials and reads them out (this WORKS -- "Solve Skipped").
- Two-test pattern: generator (RunApp) writes `<gen>_out.e`; reader has `prereq=<gen>`.
  The .e is a regenerated run artifact (gitignored), not committed.
- Confirmed field LAYOUT used: elemental (CONSTANT MONOMIAL). z_top here is elemental too
  (value only); in production z_top/H want nodal LAGRANGE for non-zero grad(z_T)/grad(H).
- relperm/Pc tables were OUT OF SCOPE (kept analytic). Next: read VERelPermTableUO /
  VECapillaryPressureTableUO tables (CSV today) and/or upscaled curve fields from the mesh.
- NATURAL FOLLOW-ON: cross-code benchmark (MRST co2lab Johansen / Sleipner-like), verif.
  ladder step 3 -- builds directly on this read pipeline + a real sloped/varying mesh.

## EOS reference depth z_T + h -- DONE (key-subtlety 2)

Switchable `eos_reference_depth` on `VEFluidProperties` (FE) and `VEFVFluidProperties`
(FV): `top_surface` (default, evaluate at `pp_top`; all golds byte-identical) or
`interface` (evaluate at the CO2-brine contact via the hydrostatic pressure
`p = pp_top + rho_n(pp_top)*|g|*h`, sharp-interface `h = sat_n*H/(1-S_wr)` clamped to
`[0,H]`). Design points:

- BOTH materials source `H = z_top - z_bottom` from COUPLED geometry variables (params
  `z_top`/`z_bottom`), NOT from `VEGeometry`'s `ve_H`. This was a deliberate revision:
  no FV input runs VEGeometry (FV kernels build H inline), so an ve_H dependency would
  have made interface mode unusable in FV -- forcing the inconsistent "flux at interface,
  storage at top_surface" split. Coupling z_top/z_bottom (VALUES only; no gradient
  needed) works identically in FE (LAGRANGE aux) and FV (FV vars), so the FV STORAGE can
  be interface too. The FE flux still uses VEGeometry for grad(z_T); the fluid material
  no longer does.
- `h` is built LOCALLY from `sat_n`, NOT from `VEPlumeReconstruction`'s `ve_h` --
  `VEPlumeReconstruction` consumes `ve_density`, so pulling its `ve_h` into the fluid
  material would be a material-dependency CYCLE. The sharp `h` is adequate for the EOS
  pressure (the fringe difference is 2nd order).
- `initQpStatefulProperties` seeds the old-value density at the SAME interface pressure as
  the regular compute. This is fine because H/sat_n now come from coupled VARIABLE values,
  which are populated from the ICs before stateful init (unlike material properties such
  as the former ve_H, which are NOT computed during initStatefulProps -- that earlier
  blocker is gone with the ve_H dependency).

`p` carries AD wrt both `pp_top` and `sat_n` (the latter through `h`). FV mirrors the FE
math via a templated `refPressure(r,t)` helper inlined in the header (header now
`#include`s SinglePhaseFluidProperties.h). Tests: `fluidproperties/eos_interface_fe_jacobian`
(full FE interface path, sloped top) and `eos_interface_fv_jacobian` (fully consistent FV
model: BOTH flux and storage at interface, H from coupled FV z_top/z_bottom). Both
PetscJacobianTester, ratio_tol 1e-7, SimpleFluidProperties (see Gotcha D), PASS.

**Convergence with a real (CO2) EOS in interface mode.** The Gotcha-D Jacobian
inconsistency does NOT stop a real solve. The residual is exact (true depth-integrated
mass balance with the real CO2 density), so Newton's fixed point -- and the converged
answer -- is correct regardless of Jacobian quality; an inexact Jacobian only costs
convergence RATE (quadratic -> superlinear), typically an extra iteration or two per step.
The amplification factor is `|g|*h`, so it grows with plume thickness: a very thick plume
(h ~ hundreds of m) could reach ~1e-2 Jacobian error and want a smaller dt, but still
converges. ESCAPE HATCH if it ever bites a real case: freeze `rho_n_top` within the Newton
loop (lag the increment density `rho_n*g*h`). That drops the nested `|g|*h*drho/dp`
coefficient, so `d(rho)/d(pp_top)` collapses to plain `drho/dp` (top_surface-quality
consistency), at the cost of the increment's pp_top-derivative lagging one iteration. Not
implemented (1e-3 is fine for now); noted for the field runs.

**Gotcha D (interface mode amplifies EOS derivative inconsistency -- Jacobian tests):**
the FE interface Jacobian test first FAILED at ||J-Jfd||/||J|| ~ 1.4e-3 with
CO2FluidProperties + SimpleBrineFluidProperties, then PASSED at <1e-7 with
SimpleFluidProperties (same geometry/state). Root cause is NOT a code bug. The interface
density derivative is a composition with a nested coefficient
`d(rho)/d(pp_top) = drho/dp|_p * (1 + |g|*h*drho/dp|_{p_top})`, so the small inconsistency
`d` between an EOS's value and its analytic `drho/dp` is amplified by `|g|*h` (~100). For
CO2 (Newton-inverted Span-Wagner) `d` is tiny -- invisible wherever the derivative is just
`drho/dp` (all prior top_surface tests, incl. heavily-weighted FV mass storage in
fv_eos_jacobian, pass <1e-7) -- but `|g|*h*d` reaches ~1e-3 once interface mode puts that
density in the heavily-weighted mass-STORAGE term. Two conditions, both needed: the |g|*h
amplification (interface-only) and the storage weight. (NB the earlier guess that "FV
storage doesn't include density" was wrong -- FV storage DOES carry density; the original
FV interface test merely happened to leave VEFluidProperties in top_surface, so its storage
never saw the amplified term. Both interface Jacobian tests now run interface storage with
SimpleFluidProperties.) LESSON: use SimpleFluidProperties (exact analytic derivatives) for
storage-weighted / interface-mode Jacobian TESTS; CO2FluidProperties is fine for real SOLVES
(residual exact; only Newton rate is affected -- see the convergence note above) and for
value/physics checks. Same family as Gotchas A/C: don't "fix" correct AD code in response to
an EOS- or kink-induced FD mismatch -- change the test config to remove the confounder.

## [Physics/VEFlow] system -- option test-coverage plan (proposed)

The `[Physics/VEFlow/FiniteElement|FiniteVolume]` action system landed this session
(commits 7f15b1c, 463569b, + the variable-ownership expansion): `VEFlowPhysicsBase` +
`VEFlowFE`/`VEFlowFV` generate the primary variables, the 4 flow kernels, the material
chain, AND the auxiliary variables the action can default.

**Action now owns its variables.** It declares the primary `pp_top`/`sat_n` (names are
the user-facing handles for BCs/wells/ICs/PPs -- documented in the class descriptions;
override with pressure_variable/saturation_variable), the geometry `z_top`/`z_bottom`
(plain `VariableName` params, defaults "z_top"/"z_bottom"; FE LAGRANGE / FV FVReal), and
`c_diss` when dissolution is on. The user supplies VALUES via ICs (or overrides a name to
a mesh field). `enable_dissolution = true` now fully wires VEDissolution + c_diss +
VEDissolvedCO2Aux + the sink kernel; `dissolution_flux` is required when enabled.

Geometry ownership is an explicit opt-out: `define_geometry_variables` (default true). True
= action declares z_top/z_bottom (FE LAGRANGE / FV FVReal), user sets values via IC. False =
the user provides them (e.g. Exodus elemental fields via initial_from_file_var) and the action
validates existence + type. This replaces relying on MOOSE's opaque same-type "merge" / the
nondeterministic action-vs-AddAuxVariableAction ordering (`shouldCreateVariable` does not
reliably skip a user-declared aux). Type validation is `checkGeometryVariableType` (called from
`checkRequiredFields`): FV requires `MooseVariableFVReal`, FE rejects FV vars. (`isNodalDefined()`
is a runtime per-node query -- false during setup -- so it's NOT usable for a static nodal
check; use `isVariableFV`.) c_diss is always action-owned regardless of the flag.

Double-declaration guard: with define=true, `checkGeometryNotUserDeclared()` errors if the user
also declares z_top/z_bottom in [AuxVariables] (detected via the ActionWarehouse,
`_awh.getActions<AddAuxVariableAction>()` -- reliable at parse time, unlike `_problem->hasVariable`).
So an input must not both set define=true and declare the geometry; the pre-declaring verification
inputs (`ve_flow_*_capillary`, `ve_flow_anisotropic`) use define=false. Negative test
`err_geometry_double_declared`.

Field validation: `checkRequiredFields()` (at add_material) errors with an actionable
"rebuild the Exodus file with fields z_top/z_bottom/phi_bar/K_up_*" message when a coupled
field names a variable that is not present -- e.g. an upscaling field the workflow did not
emit, or a misspelled coupling. Constants and the unset optional K_up_xy are skipped.
(A field missing from the Exodus *file* itself, read via initial_from_file_var, is already
caught by MOOSE's CopyNodalVars.) Negative test `err_missing_field`.

**Committed coverage so far (full suite 63):** `physics/ve_flow_{fe,fv}` (CSVDiff vs the
`nc_sloped` golds -- bit-for-bit equal to the hand-written input; covers variables,
kernels, geometry, constant `phi_bar`/`K_up`, `ConstantFluidProperties`, sharp
`relperm_uo`, default gravity/eos), `physics/ve_flow_{fe,fv}_capillary` (PetscJacobianTester),
`physics/ve_flow_dissolution_{fe,fv}` + `physics/ve_flow_anisotropic` (CSVDiff), and the
three negative tests below.

The block has more options than that exercises. Suite to close the matrix --
**prefer mirroring an existing hand-written test and CSVDiff/Exodiff the Physics version
against the SAME gold** (equivalence), use PetscJacobianTester for branches that change
the Jacobian, and RunException for the validation guards. All new inputs live in
`test/tests/physics/` with one `[Tests]` block appended to `physics/tests`.
**[DONE]** marks cases implemented and passing (full suite now 63); the rest are still TODO.

### Equivalence (CSVDiff/Exodiff vs the existing gold)

| New input | Option exercised | Mirror | Status |
|---|---|---|---|
| `ve_flow_dissolution_{fe,fv}.i` | `enable_dissolution = true` | `dissolution/constant_rate{,_fv}.i` | **DONE** |
| `ve_flow_anisotropic.i` | `K_up_xy` off-diagonal forwarding | `permeability/anisotropic_permeability.i` | **DONE** |
| `ve_flow_fv_user_geometry.i` | `define_geometry_variables = false` (user-declared FV z_top/z_bottom) | own gold (= nc_sloped_fv) | **DONE** |
| `ve_flow_relperm_table.i` | `relperm_uo` = tabulated UO (not sharp) | `buckley_leverett/bl_flat_table.i` | TODO |
| `ve_flow_exodus_fields.i` | `phi_bar`/`K_up`/`z_top` as coupled Exodus **fields** (exercises the coupled-var path, vs the constant path the nc tests use) | `exodus/exodus_fields.i` | TODO |
| `ve_flow_var_names.i` | `pressure_variable`/`saturation_variable` renamed | `buckley_leverett/bl_flat.i` (rename + matching PPs) | TODO |

Equivalence golds: the DONE cases keep their own gold in `physics/gold/` (verified once to
equal the mirrored hand-written gold), matching the `ve_flow_{fe,fv}` pattern.

### Jacobian (PetscJacobianTester -- branches that alter the Jacobian)

| New input | Option exercised | Mirror | Status |
|---|---|---|---|
| `ve_flow_eos_jacobian.i` | real EOS: `fp_nw = CO2FluidProperties`, `fp_w = SimpleBrineFluidProperties`, `temperature` | `fluidproperties/fv_eos_jacobian.i` | TODO |
| `ve_flow_eos_interface_{fe,fv}.i` | `eos_reference_depth = interface` + `S_wr` | `fluidproperties/eos_interface_{fe,fv}_jacobian.i` | TODO |
| `ve_flow_{fe,fv}_capillary_fringe.i` | `capillary = capillary_fringe` + `pc_uo` (VECapillaryPressureTableUO) | own inputs (FE: cap_equilibrium-style non-uniform; FV: flowing uniform) | **DONE** |

### Negative (RunException -- the validation guards) -- all **DONE**

| Test | Guard | Expected error |
|---|---|---|
| `err_fv_geometry_fe_variable` (`define_geometry_variables=false`, z_top declared FE LAGRANGE; input `ve_flow_fv_bad_geometry.i`) | `VEFlowFV::checkGeometryVariableType` | "is not a finite-volume variable" |
| `err_dissolution_requires_flux` (`enable_dissolution = true` without `dissolution_flux`) | `VEFlowPhysicsBase` ctor guard | "dissolution_flux is required when enable_dissolution = true" |
| `err_interface_requires_swr` (`eos_reference_depth = interface` without `S_wr`) | `VEFlowPhysicsBase` ctor guard | "S_wr is required when eos_reference_depth = interface" |
| `err_missing_field` (`phi_bar` names a variable not present) | `checkRequiredFields` (rebuild-Exodus guard) | "not present in the problem ... ghost_field" |
| `err_geometry_double_declared` (define=true + user `[AuxVariables]` z_top) | `checkGeometryNotUserDeclared` (ActionWarehouse) | "z_top' has already been declared ... define_geometry_variables = true" |

`err_dissolution_requires_flux` / `err_interface_requires_swr` / `err_missing_field` reuse
`ve_flow_fe.i` via `cli_args`. (The earlier `err_fv_geometry_constant` was dropped: z_top is now a name param,
not a coupled-var-or-constant, so a constant is no longer expressible.)

### Caveats / things the writer must get right

- **Dissolution is now fully wired by the toggle** (no longer kernel-only): `enable_dissolution
  = true` adds `VEDissolution` + `c_diss` aux var + `VEDissolvedCO2Aux` + `VE*DissolutionSink`.
  The user passes the model params on the Physics (`dissolution_flux` required, plus optional
  `dissolution_s_ref` / `dissolution_c_cap`); no `[Materials]`/`[AuxKernels]`/`[AuxVariables]`
  dissolution boilerplate -- the `ve_flow_dissolution_{fe,fv}.i` inputs reduce to the Physics
  block + ICs/BCs and still match the `constant_rate` golds.
- **Interface EOS must use `SimpleFluidProperties`, not CO2.** Gotcha D (see "EOS reference
  depth" section): interface mode amplifies the value-vs-analytic-derivative inconsistency
  by `|g|*h`, so `CO2FluidProperties` trips the strict 1e-7 ratio. The real-EOS Jacobian
  test (`ve_flow_eos_jacobian`) can use CO2 in `top_surface` mode; the interface tests
  mirror `eos_interface_*_jacobian.i` which deliberately use `SimpleFluidProperties`.
- **`pc_entry` and non-default `gravity` are still globally untested** (see audit). Fold a
  `pc_entry > 0` variant into one capillary test rather than a standalone case.
- **Geometry thickness-mode -- REMOVED.** `VEGeometry`'s H-supplied "thickness" mode was
  deleted: it was never plumbed through the Physics action, was incompatible with
  `eos_reference_depth = interface` and with FV, and `ex2ve` emits `z_top`/`z_bottom`
  directly. `VEGeometry` now requires `z_top` + `z_bottom` and computes `ve_H = z_top - z_bottom`.
  The dedicated `geometry/thickness_mode.i` test was removed; `exodus/exodus_fields*` and the
  `anticline_a1` benchmark were migrated to supply `z_bottom` (derived `z_top - H` where the
  source carried `H`).
- **`capillary = capillary_fringe` -- DONE (this session).** The fringe flux path is wired in
  both FE and FV; the ctor `paramError` guard is gone. New Physics param `pc_uo` (required in
  fringe) threads one `PorousFlowCapillaryPressure` UO (e.g. `VECapillaryPressureTableUO`) to BOTH
  `VEPlumeReconstruction` (mode=capillary_fringe) and the cap-pressure material. Negative test
  is now `err_capillary_fringe_requires_pc_uo`. See "Capillary-fringe flux path" below for the
  design and gotchas.

## Other open threads
- FV-vs-FE numerical-diffusion comparison for the verification writeup.
- `topography_demo` repo decision: parked on the local `topography-demo` branch
  (commit 224d05c), off `main`. Push if it should live remotely.
- **grad(H) in capillary_fringe mode -- RESOLVED (this session).** The fringe partial is
  `ve_dPcup_dH = delta_rho*g*sat_n/S_n(h)` (and `ve_dPcup_dsatn = delta_rho*g*H/S_n(h)`),
  with `S_n(h) = 1 - Sw(delta_rho*g*h)`; sharp is the `S_n = 1 - S_wr` special case. See
  "Capillary-fringe flux path".

## Conventions / environment
- ASCII-only in `.h`/`.C` (no Unicode). Verify: `grep -Prn '[^\x00-\x7F]' src include`.
- `.gitignore` already excludes run artifacts.
- Gold-file naming: default `file_base` + a `[csv]` block → `<input>_csv.csv`;
  explicit `Outputs/file_base=X` → `X.csv`.
- User commits straight to `main` (solo dev); branch only if asked.
- Build/test not runnable in the assistant's sandbox — the user runs
  `./run_tests` and reports.
