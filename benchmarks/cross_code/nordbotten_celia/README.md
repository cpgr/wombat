# Nordbotten-Celia sloped aquifer: MRST co2lab <-> wombat cross-code comparison

Verification-ladder **step 2/3**: a 1D sloped aquifer, constant properties, no trap,
`S_wr = 0`. The simplest rung that exercises the sloped-surface buoyancy drive, used
to shake out unit conventions and physics-configuration mismatches before the 2D
anticline (Tier A).

## Files

| file | what it is |
|------|------------|
| `nc_sloped_mrst.m` | MRST `CO2VEBlackOilTypeModel` mirror (run in MATLAB after `mrstModule add co2lab ad-core ad-blackoil`). Writes `nc_mrst_sat_*.csv`. |
| `nc_sloped_fv_crosscode.i` | wombat FV input, **capillary ON**, configured to match MRST. The cross-code comparison input. |
| `nc_buoyancy_slump.i` | pure-buoyancy gravity-slump diagnostic (CO2 block, no injection) that isolates the spreading term. |

The clean regression test (capillary OFF, matches the analytic nose velocity) lives
separately at `test/tests/nordbotten_celia/nc_sloped_fv.i` -- see the note below on why
the two inputs differ.

## Running

```
# wombat side
./wombat-opt -i benchmarks/cross_code/nordbotten_celia/nc_sloped_fv_crosscode.i
./wombat-opt -i benchmarks/cross_code/nordbotten_celia/nc_buoyancy_slump.i Executioner/end_time=5e5

# MRST side (in MATLAB)
run('benchmarks/cross_code/nordbotten_celia/nc_sloped_mrst.m')
```

## The one thing that matters: turn capillary ON to match MRST

MRST's `CO2VEBlackOilTypeModel` includes the two-pressure capillary flux
`Pc^up = delta_rho * g * h` **intrinsically** -- it cannot be switched off. Its
gradient `delta_rho * g * dh/dx` is the hydrostatic **spreading** term, and it
dominates plume migration. wombat exposes this as an option (`capillary = true` on the
CO2 `VEFVAdvectiveFlux` + a `VEFVCapPressure` material), and the regression test leaves
it **off** so it reproduces the analytic sharp-tongue nose velocity
`v_N = K*drho*g*sin(theta)/(phi*mu_n)`. **To compare against MRST you must turn it on**,
or the codes are running different physics:

| | injection inlet sat (t=2e6) | block-slump nose advance (t=5e5) |
|---|---|---|
| wombat, capillary **off** | 0.0913 (= analytic v_N regime) | 4.75 m (= v_N*t) |
| wombat, capillary **on**  | 0.0595 | 29.25 m |
| MRST | 0.046 | 49.75 m |

The mismatch is saturation-dependent (`dh/dx` grows with column thickness): ~2x at the
thin inlet, ~10x for a thick plume. Enabling capillary moves wombat into MRST's
spreading regime; the residual (0.0595 vs 0.046) is second-order (exact spreading
coefficient, numerical diffusion, time step), not a regime difference.

After the incompressible-fluid fix on the MRST side, the codes agree to ~0.1% on mass
(integral sat_n dx = 2.0 m, injected mass 5600 kg). See
`doc/content/benchmarks/cross_code.md` for the full account and the MRST setup gotchas.
