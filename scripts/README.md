# wombat scripts

## ve_upscale_curves.py -- functional (kr/Pc) upscaling over the thickness

The VE solver consumes DEPTH-INTEGRATED relative permeability `kr^up(sat_n)` and an
upscaled column `Sw(Pc)` curve. The geometry/petrophysics upscaling that builds the 2D
mesh handles `H`, `phi_bar`, `K_up`; this script handles the *functional* part -- it
turns fine-scale (vertically resolved) rock `kr(s)` / `Pc(s)` curves into the tables
`VERelPermTableUO` and `VECapillaryPressureTableUO` read. Preprocessing, like
`ve_reconstruct_vertical.py`; not part of the solve.

For each plume height `h` it integrates the capillary-fringe VE equilibrium profile over
the column (depth `d` from the top, contact at `h`, `Pc(d) = delta_rho*g*(h-d)`,
`Sw = Sw(Pc)`), `k`/`phi`-weighted:

    sat_n(h)   = int_0^h phi (1 - Sw) dd / int_0^H phi dd
    kr_n^up(h) = int_0^H k kr_n(Sw) dd / int_0^H k dd     (kr_n = 0 below the plume)
    kr_w^up(h) = int_0^H k kr_w(Sw) dd / int_0^H k dd

Sweeping `h` over `[0, H]` traces `kr_n^up`, `kr_w^up` vs `sat_n`. `--mode sharp` uses the
box profile `S_n = 1 - swr` (reproduces `VERelPermSharpUO` on a uniform column; still
non-trivial with a layered `k(d)`). The emitted `Sw(Pc)` table is the fine-scale rock
curve -- feed the same `pc_uo` to `VEPlumeReconstruction` so the in-solve reconstruction
matches the curve the `kr^up` table was built from.

Fine-scale curves: parametric van Genuchten (`--pc vg`, matching
`PorousFlowCapillaryPressureVG`) or Brooks-Corey (`--pc bc`) Pc; Corey relperm
(`--kr corey`); or tables (`--pc-table Pc,Sw`, `--kr-table Sw,krw,krn`). Vertical
heterogeneity via optional `--k-profile`/`--phi-profile` CSVs (`depth_below_top, value`);
default uniform.

Pure `numpy` (no MOOSE / SEACAS), so it runs in any environment:

```
# van Genuchten Pc + Corey relperm, uniform 30 m column
python scripts/ve_upscale_curves.py --out-prefix sleipner \
    --thickness 30 --rho-n 700 --rho-w 1000 \
    --pc vg --vg-alpha 5e-4 --vg-m 0.5 --swr 0.2 \
    --kr corey --krn-max 1 --krw-max 1 --corey-nn 2 --corey-nw 2

# real-field facies tables + a layered permeability profile
python scripts/ve_upscale_curves.py --out-prefix field \
    --thickness 50 --rho-n 650 --rho-w 1020 \
    --pc table --pc-table rock_pc.csv --kr table --kr-table rock_kr.csv \
    --k-profile k_of_depth.csv --swr 0.15

python scripts/ve_upscale_curves.py -h   # all options
```

Outputs `<prefix>_upscaled.i` (a `[UserObjects]` block with `VERelPermTableUO` and -- in
fringe mode -- `VECapillaryPressureTableUO`, ready to `!include` or paste) plus
`<prefix>_relperm.csv` and `<prefix>_pc.csv` for inspection/plotting. The emitted tables
satisfy the UOs' monotonicity requirements by construction (`sat_n` strictly increasing,
`pc` increasing from 0, `sw` strictly decreasing, `sat_lr = swr <= min(sw)`).

### Notes

- The `kr^up` table is indexed by `sat_n`, so it stays consistent with the solve regardless
  of `phi(d)`; the `phi`-weighting only sets the (pore-volume) `sat_n <-> h` mapping, which
  is exactly the average the `sat_n` variable represents.
- Multi-facies (different rock curves stacked vertically) is a defined next step; v1 is a
  single rock curve with optional layered `k(d)`/`phi(d)`. For a multi-facies column the
  effective `Sw(Pc)` is only approximate (CLAUDE.md key-subtlety 9).
- `--mode sharp` emits no Pc table (the sharp reconstruction needs no `pc_uo`).

## ve_reconstruct_vertical.py -- 3D vertical saturation reconstruction

The VE solver runs on a 2D mesh (the formation top surface) and tracks the plume only
through depth-averaged fields. Under the VE assumption each vertical column is in
instantaneous capillary-gravity equilibrium, so the full 3D CO2 saturation can be
reconstructed *after* a run as a pointwise function of vertical position. This script
reads a 2D VE Exodus time series and writes a transient 3D Exodus on a static extruded
column mesh (no solve -- pure evaluate-and-write), one frame per stored timestep.

- `sharp` mode: box column, `S_n = 1 - S_wr` within the plume thickness `h`.
- `fringe` mode: graded column, `S_n(zeta) = 1 - Sw(delta_rho*g*zeta)` -- the same
  integrand as `VEPlumeReconstruction::satNBar`, evaluated pointwise. `Sw(Pc)` comes
  from a table (`--pc table`, recommended -- export the curve the run used) or analytic
  van Genuchten (`--pc vg`).

### What the 2D run must output

As AuxVariables (see `test/tests/reconstruction/ve_reconstruction_demo.i` for the full
wiring): `z_top` and `z_bottom` (or `H`) for the column geometry; `h_plume`
(`VEPlumeHeightAux`, the mobile-plume thickness `ve_h`); and -- for fringe mode -- the
two phase densities `rho_n`, `rho_w` (two `ADMaterialStdVectorAux` on `ve_density`,
index 0 = non-wetting, 1 = wetting). Output both densities, not a precomputed
`delta_rho`, so the file stays self-documenting; the script forms the difference. For a
trapped tail, also output `sat_n_max`.

### Running it (needs the SEACAS exodus.py)

The Python `exodus`/`exodus3` bindings ship with the conda `moose-seacas` package, which
is NOT part of the `moose` env and is NOT needed to build or test MOOSE (the test
harness uses the `exodiff` binary). Install it into a small dedicated env -- the script
only needs `exodus` + `numpy`:

```
mamba create -n seacas -c https://conda.software.inl.gov/public moose-seacas numpy
conda activate seacas
export PYTHONPATH=$CONDA_PREFIX/seacas/lib   # the env's bindings (NOT a pkgs/ cache path)
python -c "import exodus3; print('ok')"
```

The `moose-seacas` package adds no activation hook for `PYTHONPATH`, so you must point it
at `$CONDA_PREFIX/seacas/lib` yourself (where the env's `exodus3.py` + `libexodus.dylib`
live). If `import exodus3` instead loads a `.../pkgs/...` path and fails on a missing
`libhdf5`, a stale `PYTHONPATH` is shadowing the env copy -- `unset PYTHONPATH` and set it
again, or open a fresh shell.

Then, from the directory holding the 2D output:

```
# sharp interface, footprint clipped to the real plume
python /path/to/scripts/ve_reconstruct_vertical.py run_out.e plume3d_sharp.e \
    --mode sharp --s-wr 0.0 --n-layers 40 --vertical up --sat-min 1e-2 \
    --z-top-var z_top --z-bottom-var z_bottom --h-var h_plume

# capillary fringe (analytic VG curve here; use --pc table on real fields)
python /path/to/scripts/ve_reconstruct_vertical.py run_out.e plume3d_fringe.e \
    --mode fringe --pc vg --vg-alpha 1.0e-3 --vg-m 0.5 --s-wr 0.0 \
    --n-layers 40 --vertical up --sat-min 1e-2 \
    --z-top-var z_top --z-bottom-var z_bottom --h-var h_plume \
    --rho-n-var rho_n --rho-w-var rho_w

# drape onto the ORIGINAL 3D geology (em2ex grid) instead of an extruded column
# mesh -- exact faults/layering/topography, with porosity/permeability carried through
python /path/to/scripts/ve_reconstruct_vertical.py run_out.e plume3d_orig.e \
    --original-mesh field_model.e --mode sharp --s-wr 0.0 --vertical up --sat-min 1e-2 \
    --z-top-var z_top --z-bottom-var z_bottom --h-var h_plume
```

Open `plume3d_*.e` in ParaView and colour by `sat_co2` (`co2_region`: 0 brine / 1 mobile
/ 2 trapped). `--sat-min` drops the thin numerical-diffusion halo so only the real plume
footprint is drawn; `--vertical up`/`down` sets the z convention (`up` -> bottom at
`z_top - H`). `python ve_reconstruct_vertical.py -h` lists all options.

By default the output is a synthetic extruded column mesh (`--n-layers` controls the
vertical resolution). With **`--original-mesh ORIG.e`** the reconstruction is draped onto
the original 3D mesh instead: each original node is mapped to its nearest 2D VE column and
`sat_co2` is evaluated at the node's real depth below `z_top`, so the plume sits in the
exact geology. The original mesh's elemental fields (porosity, permeability, ...) are
carried into the output for context (`--no-carry` to omit them). `--n-layers` is ignored
in this mode (the vertical resolution is the original mesh's).

### Notes

- Works for FE (nodal) and FV (elemental) output -- elemental fields are averaged to
  nodes for the reconstruction.
- The vertical profile is instantaneous (each frame depends only on that step's fields),
  except the trapped tail, which the 2D solve already accumulated into `sat_n_max`.
- For a fully fringe-consistent result, run the model itself in `capillary_fringe` mode
  with the same `pc_uo` and reconstruct with the matching curve; reconstructing a fringe
  on a sharp run is a visualisation overlay (the sharp reconstruction is the
  mass-consistent one for a sharp run).
