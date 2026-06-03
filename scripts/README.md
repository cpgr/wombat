# wombat scripts

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
```

Open `plume3d_*.e` in ParaView and colour by `sat_co2` (`co2_region`: 0 brine / 1 mobile
/ 2 trapped). `--sat-min` drops the thin numerical-diffusion halo so only the real plume
footprint is drawn; `--vertical up`/`down` sets the z convention (`up` -> bottom at
`z_top - H`). `python ve_reconstruct_vertical.py -h` lists all options.

### Notes

- Works for FE (nodal) and FV (elemental) output -- elemental fields are averaged to
  nodes for the reconstruction.
- The vertical profile is instantaneous (each frame depends only on that step's fields),
  except the trapped tail, which the 2D solve already accumulated into `sat_n_max`.
- For a fully fringe-consistent result, run the model itself in `capillary_fringe` mode
  with the same `pc_uo` and reconstruct with the matching curve; reconstructing a fringe
  on a sharp run is a visualisation overlay (the sharp reconstruction is the
  mass-consistent one for a sharp run).
