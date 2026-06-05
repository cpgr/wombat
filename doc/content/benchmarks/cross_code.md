# Cross-code benchmark: MRST co2lab <-> wombat VE

Verification-ladder step 3 (see `CLAUDE.md`). The first two ladder rungs are
analytical (Buckley-Leverett with gravity; Nordbotten-Celia sloped aquifer). This
rung is a *cross-code* comparison against an established, independent VE
implementation: the MATLAB **MRST co2lab** module, the canonical reference VE code.

It builds directly on the Exodus real-field input pipeline
(`test/tests/exodus/`): the two codes exchange a shared 2D top-surface mesh, so
that any disagreement is physics, not two different grids.

## Guiding principle: control the confounders

A cross-code comparison is only as sharp as its controllability. Jumping straight
to a full Johansen run (real EOS + dissolution + black-oil PVT on an unstructured
grid) makes any mismatch uninterpretable -- is a 15% footprint error the
sloped-surface flux kernel, the density reference depth, the trapping model, or
just two different Span-Wagner tables?

So the benchmark is **three tiers**, each adding exactly one new source of
disagreement, on the *same* mesh and with the *same* fluid data wherever the two
codes allow it.

## The non-negotiable rule: share the mesh

The two codes must run on the *same* 2D top-surface grid, or physics error can
never be separated from discretization error. Workflow:

1. MRST builds the VE grid and upscales (`z_T`, `H`, `phi_bar`, `K_up`) via
   `topSurfaceGrid` / `averageRock`.
2. MRST exports cell- and node-centred fields to CSV (or `.mat`).
3. A small converter (clone of the `test/tests/exodus/` generator pattern) writes
   a 2D **Exodus** carrying nodal `z_top`, `z_bottom` (LAGRANGE FIRST -- needed for
   non-zero `grad(z_T)`, `grad(H)`; `ve_H = z_top - z_bottom`) and elemental
   `phi_bar`, `K_up_xx/yy/xy` (CONSTANT MONOMIAL).
4. wombat reads it through the existing pipeline (`FileMesh` +
   `initial_from_file_var` -> `VEGeometry` / `VEPorosity` / `VEPermeability`).

For the conceptual tier (below) the grid is analytic, so step 1-3 collapse into a
single MOOSE generator input that *both* codes can reproduce exactly -- no MRST
installation needed to stand the test up.

## Tier A -- conceptual anticline, stripped physics (the workhorse)

Exercises the VE *core* -- the sloped-surface buoyancy drive in genuine 2D plus
structural trapping -- with zero EOS ambiguity. More discriminating than
`nordbotten_celia` (which only checks 1D nose velocity) because the plume migrates
updip and pools in a structural trap.

**Geometry** -- a tilted plane carrying a Gaussian dome (a structural culmination,
i.e. a local crest that forms the trap):

```
domain    10 km x 10 km, 100 x 100 cells (+ a 200 x 200 refinement run)
z_T(x,y)  = -800 + tan(1.5 deg)*x + 40*exp(-((x-7000)^2 + (y-5000)^2)/(2*1200^2))
            (regional dip: +x is updip/shallower; +40 m dome crest near (7km,5km))
H         = 50 m constant  (then an H(x,y)-varying variant to exercise grad(H))
phi_bar   = 0.20
K_up      = 200 mD isotropic  (then an anisotropic Kxx != Kyy variant -- subtlety 6)
```

**Physics (deliberately minimal, identical in both codes):**

| Knob        | Value                                              |
|-------------|----------------------------------------------------|
| Interface   | sharp                                              |
| Densities   | constant: rho_CO2 = 700, rho_brine = 1000 kg/m^3   |
| Viscosities | constant: mu_CO2 = 6e-5, mu_brine = 8e-4 Pa.s      |
| Rel-perm    | sharp/linear, S_wr = 0.2                           |
| Capillary   | off (then a non-zero-fringe variant)               |
| Dissolution | off                                                |
| Trapping    | **off** in A1; **Land** (S_gr = 0.2, matched C) in A2 |
| Far boundary| open hydrostatic brine outflow on the updip edge; no-flow elsewhere |

Note: with constant (incompressible) densities the domain cannot be fully closed --
injected mass needs an outlet or the pressure problem is overconstrained. Tier A
therefore uses an open hydrostatic brine boundary on the updip edge (matched in both
codes). A fully closed, pressure-confined formation is a compressible-EOS case
(Tier C).

**Schedule:** inject 1 Mt/yr for 10 yr at a downdip well (~ (2 km, 5 km)), then
490 yr of pure migration. In wombat the well is a boundary/areal source; match the
total kg/s and location exactly to co2lab's rate well.

**Why two sub-cases.** A1 (no trapping) isolates the advective flux + `grad(z_T)`
drive. A2 turns on hysteretic trapping (the Land/Killough path) against co2lab's
`'residual'` trapping option -- the part of v1 that has so far been checked only
analytically (`column_cycle`), never cross-code.

## Tier B -- Johansen, constant density first

Same physics as Tier A, on the real Johansen sector grid that co2lab ships
(`getAtlasGrid('Johansen')` / `makeJohansenVEgrid`). Still constant rho, mu. This
isolates **grid realism** -- irregular `z_T`, spatially varying `H`, `phi`,
anisotropic `K` -- from EOS. If A passes and B fails, the suspect is the
geometry/permeability plumbing on an unstructured field, not the kernels.

## Tier C -- Johansen, real EOS

Switch wombat's `VEFluidStateBrineCO2` from constant density to the
`PorousFlowBrineCO2` wrap, and co2lab to its `CO2props` PVT.

De-risk the EOS confounder with an intermediate **C0**: feed wombat a tabulated
`rho(p,T)`, `mu(p,T)` exported from co2lab so both codes use byte-identical fluid
tables. Any remaining gap is then purely the depth-integrated flux + reference-depth
convention (subtlety 2), not two Span-Wagner implementations. Add convective
dissolution last, with a matched `q0` / correlation.

## What to compare (benchmark observables)

All already available as wombat post-processors / aux fields:

1. Mobile CO2 mass vs t -- `VEMobileCO2MassPostprocessor`.
2. Trapped / dissolved mass vs t (tiers with those physics on) --
   `VETrappedCO2MassPostprocessor`, `VEDissolvedCO2Aux` integral.
3. Plume footprint area vs t -- `VEPlumeFootprintPostprocessor`.
4. Updip nose position vs t -- max-x of plume > threshold (LineValueSampler).
5. Plume thickness profile `h(x, y0)` at snapshots (t = 10, 50, 100, 500 yr) --
   `VEPlumeHeightAux`.
6. Max pressure rise (compressible tiers only).
7. End-state inventory partition (mobile / trapped / dissolved fractions).

co2lab's matching series are exported once and committed as
`reference/*_reference_*.csv`, read via `data_file` (with `axis = x` for spatial
profiles -- the PiecewiseLinear-defaults-to-time trap flagged in the
McWhorter-Sunada notes).

## Pass criteria

Cross-code agreement is **not** machine precision -- two discretizations of the
same PDE. Tolerances are percent-level and documented as such (this is
verification, not regression):

| Metric                                | Tolerance                       |
|---------------------------------------|---------------------------------|
| Mass conservation *within* each code  | machine zero (hard requirement) |
| Mobile-mass-vs-t curve, rel-L2        | <= 2-3 %                        |
| Footprint area                        | <= 5 %                          |
| Nose position                         | within 1-2 grid cells           |
| h(x) profile, rel-L2                  | <= 5 %, tightening on refinement|
| End inventory partition               | <= 3 % absolute per fraction    |

The refinement check (Tier A at 100^2 and 200^2) is what turns "the curves look
close" into "they converge toward each other" -- the real verification statement.

## Layout

This benchmark lives in the top-level `benchmarks/` tree, NOT under `test/tests/` --
it is a long, field-scale comparison case, not a short regression test, so it is not
run by `./run_tests` by default. Run it explicitly when validating against co2lab.

```
benchmarks/cross_code/
  anticline/    Tier A -- analytic Exodus generator + reader/solver,
                CSVDiff vs reference/anticline_reference_*.csv
  johansen/     Tiers B/C -- Johansen Exodus + reference series (later)
  reference/    committed MRST exports + provenance README
                (co2lab script, MRST version, git SHA)
```

It reuses the two-input generator->reader `prereq` pattern from
`test/tests/exodus/tests`; the `.e` is a gitignored run artifact, the MRST
reference CSVs are committed.

## Pitfalls to nail down before running

Each silently breaks cross-code agreement:

- **Grid identity** -- non-negotiable; share the Exodus.
- **Far-boundary BC** -- co2lab often defaults to a hydrostatic pressure boundary;
  wombat tests default closed. Pick one, set both.
- **Well model** -- co2lab rate well vs a wombat boundary/areal source: match total
  kg/s, cell, and the injection-vs-migration split.
- **Density reference depth** -- subtlety 2 (interface `z_T + h` vs `z_T`); match
  co2lab's convention via the `VEFluidStateBrineCO2` strategy switch.
- **Trapping parameters** -- Land `C`, `S_wr`, `S_gr` identical; co2lab and wombat
  parameterize Land slightly differently -- convert, don't copy.
- **Total simulated time / output cadence** -- agree on snapshot times so profiles
  are compared at equal physical time.
- **Discretization** -- the FE continuous-Galerkin path oscillates on
  advection-dominated migration (HANDOVER, hysteresis Phase 4b note (2)). The FV
  path (+ TVD limiter) is the robust choice for the full migration run; the FE path
  is fine for short/diffusive smoke checks.

## FV input gotchas (getting it to converge)

The FV reader exposed two failure modes that both present as a **dead linear solve**
(LU makes no progress, the linear residual sits perfectly flat) because each leaves
the Jacobian singular. Neither is a hard-nonlinear-problem symptom -- they are
plumbing bugs. Both are fixed in `anticline/anticline_a1.i`; record them here so the
Johansen tiers do not rediscover them.

1. **Every field a flux kernel touches on a face must be a `MooseVariableFVReal`,
   never an FE auxvariable.** `ComputeFVFluxThread` only calls `computeFaceValues` on
   FV variables, so a Material that couples an FE auxvar reads **zero on faces**.
   `VEFVAdvectiveFlux` evaluates `ve_K_up` on faces: with an FE `K_up`, `ve_K_up = 0`
   on every face, the brine flux and its `pp_top` Laplacian Jacobian vanish, `pp_top`
   loses its diagonal, and the matrix is singular. The trap is that the *elemental*
   storage kernel reads the same FE field correctly, so the field looks right in the
   output -- only the face evaluation is zero. Fix: read `phi_bar`, `K_up` (and any
   real geometry, for Tier B/C) from Exodus into **FV** variables. The inverse also
   bites: `VEGeometry` feeds only the elemental mass PP (`ve_H`) but computes
   `ve_grad_z_top` via `coupledGradient`, whose array is **empty for an FV variable**
   (size-0 `MooseArray` -> out-of-bounds crash). So `VEGeometry` must couple **FE**
   (gradient-capable) geometry. The lesson: match the variable family to the
   consumer -- FV for anything evaluated on a face, FE for anything that needs a
   continuous `coupledGradient`.
2. **FV functor geometry must be live before the first Jacobian.** Setting
   `z_top`/`z_bottom` with an `INITIAL` `AuxKernel` is too late -- the functor reads
   stale zeros during the first assembly, `H = z_top - z_bottom` collapses to 0, and
   the `sat_n` storage diagonal `(phi*H*rho/dt)` vanishes. Set them with `FunctionIC`
   (Tier A, analytic geometry) or `initial_from_file_var` (Tier B/C, real geometry) --
   both are part of solution init, so the values are closed/ghosted before assembly.

The proven references for both rules are `test/tests/nordbotten_celia/nc_sloped_fv.i`
and `test/tests/mesh/dip/dip_injection_fv.i`. Pair every Jacobian/convergence change
with a physics check -- a singular matrix is loud, but a *wrong-but-consistent* one is
not (HANDOVER Gotchas A/B).

## Suggested first step

Tier **A1** end-to-end: anticline Exodus generator -> reader -> one mobile-mass
curve vs a hand-built constant-density co2lab run. It reuses everything already
green (the Exodus pipeline, `VEAdvectiveFluxS`/`VEFVAdvectiveFlux`, footprint/mass
PPs) and adds only the MRST export + converter. Get one curve agreeing to a few
percent before touching Johansen or real EOS.

The scaffold under `benchmarks/cross_code/anticline/` stands this up to the point
where the only missing piece is the MRST reference CSV (currently a placeholder --
see `reference/README.md`).

## Step-2 result: Nordbotten-Celia sloped aquifer (nc_sloped)

Before the anticline, the 1D sloped aquifer was run against an MRST
`CO2VEBlackOilTypeModel` mirror. The inputs live in
`benchmarks/cross_code/nordbotten_celia/`:
`nc_sloped_mrst.m` (MRST), `nc_sloped_fv_crosscode.i` (wombat, capillary on, the
comparison input), and `nc_buoyancy_slump.i` (the pure-buoyancy diagnostic). The clean
**capillary-off** regression test stays at `test/tests/nordbotten_celia/nc_sloped_fv.i`.
This is the cleanest rung -- 1D, constant properties, no trap, `S_wr = 0`, continuous
injection -- with three hand-checkable anchors at `t = 2e6 s`: injected mass 5600 kg,
integral(sat_n dx) = 2.0 m, gravity nose `v_N = K*drho*g*sin(theta)/(phi*mu_n)` -> 21.9 m.

**Outcome:** after forcing the MRST CO2 to be incompressible (see below), the two
codes agree to ~0.1% on mass and injection -- mass integral 2.002 vs 2.0, injected
mass 5605 vs 5600 kg. The CO2 *migration*, however, differed substantially, and
chasing it down produced the key lesson of this rung.

**The migration difference is the hydrostatic spreading term `Δρ·g·∂h/∂x`** (the
two-pressure capillary flux `Pc^up = Δρ·g·h`). wombat's `nc_sloped_fv.i` runs with
**capillary off**, so its CO2 feels only the slope drive -- it migrates at the
analytic nose velocity `v_N` and barely spreads. MRST's `CO2VEBlackOilTypeModel`
includes `p_g - p_w = Δρ·g·h` **intrinsically** and cannot disable it, so it spreads
strongly. The mismatch is therefore **saturation-dependent** (`∂h/∂x` grows with
column thickness): ~2x at the thin injection inlet (s≈0.09), ~10x for a thick
gravity-slumping block (s=0.5).

This was proven with a **pure-buoyancy block-slump** (CO2 block in x<20, no
injection, t = 5e5 s):

| run | nose advance | regime |
|-----|--------------|--------|
| wombat, capillary **off** | 4.75 m (= `v_N·t`) | slope-only, barely spreads |
| wombat, capillary **on**  | 29.25 m | strongly spreading |
| MRST | 49.75 m | strongly spreading |

and on the injection inlet saturation: wombat cap-off 0.0913 -> cap-on 0.0595 -> MRST
0.046 (same direction, ~60% of the gap closed by enabling the term). The residual
(0.0595 vs 0.046) is second-order -- exact spreading coefficient, numerical diffusion,
time step -- not a regime difference. The earlier characterisation of this as an
"inherent formulation difference" was wrong: it is a **configuration mismatch** in
which flux terms are active.

**Therefore, for a fair cross-code comparison, enable wombat's two-pressure capillary
flux** -- add `capillary = true` to the CO2 `VEFVAdvectiveFlux` and a `VEFVCapPressure`
functor material (`sat_n`, `z_top`, `z_bottom`, `S_wr`, `gravity = '0 0 -9.81'`); see
`test/tests/capillary_equilibrium/cap_equilibrium_fv.i` for the wiring. Do this in
**both** the nc_sloped and anticline wombat inputs before comparing against MRST.
(Geometry was ruled out as a cause -- `Gt.cells.z` span 3.47 = 0.0349·99.5 and cell
`dz/dx = -0.0349` both correct; the buoyancy *coefficient* matches. It was the
spreading term, not the slope term, that differed.)

### MRST setup gotchas (learned the hard way)

- **`makeVEFluid` CO2 is a real compressible EOS**, not constant density: rho = 700
  only near ~93 bar; at 0.3 bar it collapses to ~1.8 kg/m3. wombat's
  `ConstantFluidProperties` is pressure-independent. Running MRST at wombat's
  `p = 0`-referenced pressures makes the same injected mass occupy ~60-70x the pore
  volume and flood the domain (distortion scales as 1/pressure). To match a
  constant-density verification case, force incompressible:
  `fluid.bG = @(p,varargin) 1+0*p; fluid.bW = @(p,varargin) 1+0*p;`. For real field
  cases keep the EOS and run at true depth pressure.
- **Phase order is black-oil:** `s(:,1) = water, s(:,2) = CO2`; all-brine IC is
  `initResSol(Gt, p, [1,0])` and you must set `state0.sGmax = state0.s(:,2)`.
- **`poreVolume(Gt, rock2D)` for a top-surface grid is `poro*AREA` and omits H.**
  CO2 volume = `sum(s .* poreVolume(Gt,rock2D) .* Gt.cells.H)`.
- **Verify injection from `wellSols`:** `qGs`/`ComponentTotalFlux` are the true
  surface rate/mass; the `flux` field is reservoir volume, so `flux(2)/qGs = 1/bG`
  flags an over-expanded (low-pressure) CO2.
- **Wells can't change count between control periods** -- shut injection with
  `W.val = 0`, not `W = []`.

## The broader lesson: matching a reference is a physics-configuration problem

The hardest part of this benchmark was not the solver, the mesh, or the numerics --
it was discovering *which physics each code had switched on*. The migration mismatch
came down to a single flux term (`Δρ·g·∂h/∂x`) that one code includes by construction
and the other treats as optional. Cross-code comparison is, in large part, the work of
reconciling these modelling choices, and they are rarely written on the tin. Expect to
iterate: match mass first, then isolate each physics term with a stripped-down case
(the pure-buoyancy slump here) until the remaining difference is explained.

This cuts in a way that is worth stating explicitly, because it is a genuine strength of
building on MOOSE rather than a nuisance:

**wombat can dial individual physics terms in and out, so the *same* model can be made
to match very different references.** The Nordbotten-Celia case has two legitimate
"correct" answers depending on what you are comparing against:

- **To match the analytic solution**, capillary must be **off**. The closed-form nose
  velocity `v_N = K·Δρ·g·sin(theta)/(phi·mu_n)` is the sharp-tongue, slope-only result;
  it deliberately omits hydrostatic spreading. wombat with `capillary = false`
  reproduces it (nose advance = `v_N·t` to within first-order-upwind diffusion). This is
  the regression test `test/tests/nordbotten_celia/nc_sloped_fv.i`.
- **To match MRST**, capillary must be **on**, because MRST's VE model carries
  `Pc = Δρ·g·h` intrinsically and cannot remove it. wombat with `capillary = true` enters
  the same spreading-dominated regime.

Same governing code, same mesh, opposite settings of one switch -- each correct for its
target. A monolithic code with the spreading term hard-wired (either always on, like
MRST, or always off) could match only one of these references; it could never be
verified against the analytic `v_N` *and* cross-checked against MRST without editing
source. The fine-grained physics control that makes the matching take a few iterations
is exactly what lets wombat sit on every rung of the verification ladder -- analytic at
the bottom, cross-code in the middle, full field physics at the top -- using one model
with different terms enabled. The cost is that you must know, and document, which terms
are active for each comparison; that is what these benchmark inputs and this page record.
