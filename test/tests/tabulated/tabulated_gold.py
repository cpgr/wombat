#!/usr/bin/env python3
"""Generate the independent gold for tabulated_uo.i.

tabulated_uo.i (gold/tabulated_uo.csv) feeds the ve_upscale_curves.py output tables (gold/curves_*.csv,
emitted as the curves_upscaled.i [UserObjects] block) into the table-reading
UserObjects and samples them at sat_n = 0.1, 0.2, ..., 0.6:

  VERelPermTableUO  -> kr_n(sat_n), kr_w(sat_n)   = linear interpolation of the
                        (sat_n, kr) table  == numpy np.interp here.
  VECapillaryPressureTableUO -> consumed by VEPlumeReconstruction (capillary_fringe)
                        to give ve_h(sat_n). VECapillaryPressureTableUO.saturation(pc)
                        == np.interp(pc, pc_points, sw_points), and the Newton
                        inversion below replicates VEPlumeReconstruction exactly
                        (64-point trapezoidal satNBar, 20 Newton iterations).

So this gold is an INDEPENDENT reference: if MOOSE matches it, the UOs read and
interpolate the script's emitted tables correctly (and the format is valid).

Run from this directory after regenerating gold/curves_*.csv:
    python3 tabulated_gold.py
"""

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Must match tabulated_uo.i.
H = 20.0
RHO_N, RHO_W = 700.0, 1000.0
DRHO = RHO_W - RHO_N
G = 9.81
SWR = 0.2
SAT_N = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])   # = 0.1 * t, t = 1..6

rel = np.loadtxt(os.path.join(HERE, "gold", "curves_relperm.csv"), delimiter=",", skiprows=1)
pcd = np.loadtxt(os.path.join(HERE, "gold", "curves_pc.csv"), delimiter=",", skiprows=1)
satn_tbl, krn_tbl, krw_tbl = rel[:, 0], rel[:, 1], rel[:, 2]
pc_pts, sw_pts = pcd[:, 0], pcd[:, 1]


def sw_of_pc(pc):
    # VECapillaryPressureTableUO.saturation(pc) == LinearInterpolation.sample(pc),
    # which clamps to the end values outside the table -- exactly np.interp.
    return np.interp(pc, pc_pts, sw_pts)


def sat_n_bar(h):
    # Replicates VEPlumeReconstruction::satNBar (N = 64 trapezoid).
    N = 64
    dz = h / N
    total = 0.0
    for i in range(N + 1):
        sn = 1.0 - sw_of_pc(DRHO * G * i * dz)
        total += 0.5 * sn if (i == 0 or i == N) else sn
    return total * dz / H


def ve_h(sat_n):
    # Replicates VEPlumeReconstruction capillary_fringe Newton inversion (value only).
    h = max(0.0, min(sat_n * H / (1.0 - SWR), H))   # sharp warm start
    for _ in range(20):
        f = sat_n_bar(h) - sat_n
        if abs(f) < 1.0e-12:
            break
        sn_top = 1.0 - sw_of_pc(DRHO * G * h)
        if sn_top < 1.0e-14:
            break
        h -= f * H / sn_top
        h = max(0.0, min(h, H))
    return h


rows = []
for t, s in enumerate(SAT_N, start=1):
    rows.append((float(t),
                 float(np.interp(s, satn_tbl, krn_tbl)),
                 float(np.interp(s, satn_tbl, krw_tbl)),
                 ve_h(s)))

out = os.path.join(HERE, "gold", "tabulated_uo.csv")
with open(out, "w") as f:
    f.write("time,krn,krw,ve_h\n")
    for r in rows:
        f.write(f"{r[0]:.6f},{r[1]:.12e},{r[2]:.12e},{r[3]:.12e}\n")
print("wrote", out)
for r in rows:
    print(f"  t={r[0]:.0f} sat_n={0.1*r[0]:.1f}  krn={r[1]:.5f} krw={r[2]:.5f} ve_h={r[3]:.4f}")
