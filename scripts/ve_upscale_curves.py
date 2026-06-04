#!/usr/bin/env python3
"""Upscale fine-scale rock kr/Pc curves into VE depth-integrated curves.

The VE solver consumes DEPTH-AVERAGED relative permeability kr^up(sat_n) and an
upscaled column Sw(Pc) curve. This script generates both by integrating
fine-scale (vertically resolved) rock curves over the formation thickness under
the vertical-equilibrium assumption -- the "functional upscaling" step that sits
in preprocessing alongside ve_reconstruct_vertical.py, NOT in the solve. It is
the curve-generation companion to the geometry/petrophysics upscaling that
produces the 2D mesh.

PHYSICS (capillary-fringe VE equilibrium, Nordbotten & Dahle 2011)
------------------------------------------------------------------
Depth d is measured downward from the formation top; the buoyant CO2 plume
occupies the top h metres. At depth d inside the plume the height above the
CO2-brine contact is zeta = h - d, so the (hydrostatic) capillary pressure is
Pc(d) = delta_rho * g * (h - d) and the fine-scale water saturation is
Sw(Pc(d)) from the rock Pc curve; below the contact the column is fully brine.
For a column with vertical profiles k(d), phi(d) the upscaled quantities are the
transmissivity/pore-volume-weighted vertical averages

    sat_n(h)    = integral_0^h phi (1 - Sw(Pc(d))) dd  /  integral_0^H phi dd
    kr_n^up(h)  = integral_0^H k kr_n(Sw(Pc(d))) dd    /  integral_0^H k dd
    kr_w^up(h)  = integral_0^H k kr_w(Sw(Pc(d))) dd    /  integral_0^H k dd

(kr_n vanishes below the plume because Sw = 1 there.) Sweeping h over [0, H]
traces kr_n^up, kr_w^up against sat_n. 'sharp' mode uses the box profile
S_n = 1 - swr within h instead of the graded fringe -- with a uniform column it
reproduces VERelPermSharpUO, but with a layered k(d) it is still a non-trivial
k-weighted curve.

The upscaled Sw(Pc) table is the fine-scale rock curve itself (a single facies
column); feed it to VECapillaryPressureTableUO / VEPlumeReconstruction so the
in-solve reconstruction uses the same curve the kr^up table was built from.

OUTPUT
------
  <prefix>_upscaled.i  -- a [UserObjects] block (VERelPermTableUO + the optional
                          VECapillaryPressureTableUO) ready to !include or paste.
  <prefix>_relperm.csv -- sat_n, kr_n_up, kr_w_up  (inspection / plotting).
  <prefix>_pc.csv      -- pc, sw                    (inspection / plotting).

USAGE
-----
  # van Genuchten Pc + Corey relperm, uniform column
  python ve_upscale_curves.py --out-prefix sleipner \\
      --thickness 30 --rho-n 700 --rho-w 1000 \\
      --pc vg --vg-alpha 5e-4 --vg-m 0.5 --swr 0.2 \\
      --kr corey --krn-max 1.0 --krw-max 1.0 --corey-nn 2 --corey-nw 2

  # real-field facies curves from tables + a layered permeability profile
  python ve_upscale_curves.py --out-prefix field \\
      --thickness 50 --rho-n 650 --rho-w 1020 \\
      --pc table --pc-table rock_pc.csv --swr 0.15 \\
      --kr table --kr-table rock_kr.csv \\
      --k-profile k_of_depth.csv

  python ve_upscale_curves.py -h   # all options

Input CSVs are header-optional (lines starting with # ignored):
  --pc-table : two columns  Pc[Pa], Sw[-]      (Pc ascending, Sw descending)
  --kr-table : three columns Sw[-], krw[-], krn[-]
  --k-profile / --phi-profile : two columns  depth_below_top[m], value[-]
"""

import argparse
import sys

import numpy as np

# np.trapz was deprecated in NumPy 2.0 and renamed np.trapezoid (which does not
# exist on NumPy < 2.0). Use whichever the installed version provides.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


# ----------------------------------------------------------------------------------
# Fine-scale capillary-pressure curves: each provides sw_of_pc(pc) and pc_of_sw(sw)
# on TRUE water saturation (residual swr folded in), so the emitted table is
# directly consumable by VECapillaryPressureTableUO (sat_lr = swr).
# ----------------------------------------------------------------------------------
class PcVanGenuchten:
    """Sw(Pc) matching PorousFlowCapillaryPressureVG:
    Se = (1 + (alpha*pc)^(1/(1-m)))^(-m),  Sw = swr + (1-swr) Se."""

    def __init__(self, alpha, m, swr):
        if alpha is None or m is None:
            sys.exit("error: --pc vg requires --vg-alpha and --vg-m")
        if not (0.0 < m < 1.0):
            sys.exit("error: --vg-m must be in (0, 1)")
        self.alpha, self.m, self.swr = alpha, m, swr

    def sw_of_pc(self, pc):
        pc = np.asarray(pc, dtype=float)
        se = np.ones_like(pc)
        pos = pc > 0.0
        se[pos] = (1.0 + (self.alpha * pc[pos]) ** (1.0 / (1.0 - self.m))) ** (-self.m)
        return self.swr + (1.0 - self.swr) * se

    def pc_of_sw(self, sw):
        se = np.clip((np.asarray(sw, dtype=float) - self.swr) / (1.0 - self.swr), 1e-12, 1.0)
        return (1.0 / self.alpha) * (se ** (-1.0 / self.m) - 1.0) ** (1.0 - self.m)


class PcBrooksCorey:
    """Sw(Pc) = swr + (1-swr) (Pc/pd)^(-lambda) for Pc >= pd, else 1 (entry pd)."""

    def __init__(self, pd, lam, swr):
        if pd is None or lam is None:
            sys.exit("error: --pc bc requires --bc-pd and --bc-lambda")
        if pd <= 0.0 or lam <= 0.0:
            sys.exit("error: --bc-pd and --bc-lambda must be positive")
        self.pd, self.lam, self.swr = pd, lam, swr

    def sw_of_pc(self, pc):
        pc = np.asarray(pc, dtype=float)
        se = np.ones_like(pc)
        above = pc > self.pd
        se[above] = (pc[above] / self.pd) ** (-self.lam)
        return self.swr + (1.0 - self.swr) * se

    def pc_of_sw(self, sw):
        se = np.clip((np.asarray(sw, dtype=float) - self.swr) / (1.0 - self.swr), 1e-12, 1.0)
        return self.pd * se ** (-1.0 / self.lam)


class PcTable:
    """Sw(Pc) from a CSV (Pc[Pa], Sw[-]); piecewise-linear, clamped at the ends."""

    def __init__(self, path, swr):
        data = _read_csv(path, 2, "--pc-table", "Pc,Sw")
        order = np.argsort(data[:, 0])
        self.pc = data[order, 0]
        self.sw = data[order, 1]
        self.swr = swr
        if np.any(np.diff(self.pc) <= 0):
            sys.exit("error: --pc-table Pc column must be strictly increasing")
        if np.any(np.diff(self.sw) >= 0):
            sys.exit("error: --pc-table Sw column must be strictly decreasing")

    def sw_of_pc(self, pc):
        return np.interp(pc, self.pc, self.sw)

    def pc_of_sw(self, sw):
        # sw descending -> reverse for np.interp (x must be increasing).
        return np.interp(sw, self.sw[::-1], self.pc[::-1])


# ----------------------------------------------------------------------------------
# Fine-scale relative permeability: krw(sw), krn(sw) on true water saturation.
# ----------------------------------------------------------------------------------
class KrCorey:
    """Corey curves: Se = (Sw - swr)/(1 - swr - snr),
    kr_w = krw_max Se^nw,  kr_n = krn_max (1 - Se)^nn."""

    def __init__(self, krn_max, krw_max, nn, nw, swr, snr):
        self.krn_max, self.krw_max = krn_max, krw_max
        self.nn, self.nw = nn, nw
        self.swr, self.snr = swr, snr

    def _se(self, sw):
        return np.clip((np.asarray(sw, dtype=float) - self.swr) / (1.0 - self.swr - self.snr),
                       0.0, 1.0)

    def krw(self, sw):
        return self.krw_max * self._se(sw) ** self.nw

    def krn(self, sw):
        return self.krn_max * (1.0 - self._se(sw)) ** self.nn


class KrTable:
    """krw(sw), krn(sw) from a CSV (Sw[-], krw[-], krn[-]); piecewise-linear."""

    def __init__(self, path):
        data = _read_csv(path, 3, "--kr-table", "Sw,krw,krn")
        order = np.argsort(data[:, 0])
        self.sw = data[order, 0]
        self.krw_pts = data[order, 1]
        self.krn_pts = data[order, 2]
        if np.any(np.diff(self.sw) <= 0):
            sys.exit("error: --kr-table Sw column must be strictly increasing")

    def krw(self, sw):
        return np.interp(sw, self.sw, self.krw_pts)

    def krn(self, sw):
        return np.interp(sw, self.sw, self.krn_pts)


# ----------------------------------------------------------------------------------
# Vertical k(d) / phi(d) profiles (d = depth below the top). Default: uniform.
# ----------------------------------------------------------------------------------
class Profile1D:
    def __init__(self, path, name):
        if path is None:
            self.d = self.v = None
            return
        data = _read_csv(path, 2, name, "depth_below_top,value")
        order = np.argsort(data[:, 0])
        self.d = data[order, 0]
        self.v = data[order, 1]
        if np.any(self.v <= 0):
            sys.exit(f"error: {name} values must be positive")

    def __call__(self, d):
        if self.d is None:
            return np.ones_like(np.asarray(d, dtype=float))
        return np.interp(d, self.d, self.v)


def _read_csv(path, ncol, flag, cols):
    try:
        data = np.loadtxt(path, delimiter=",", comments="#")
    except Exception as exc:  # noqa: BLE001
        sys.exit(f"error: cannot read {flag} '{path}': {exc}")
    data = np.atleast_2d(data)
    if data.shape[1] < ncol:
        sys.exit(f"error: {flag} '{path}' must have {ncol} columns: {cols}")
    return data


# ----------------------------------------------------------------------------------
# Upscaling
# ----------------------------------------------------------------------------------
def upscale(args, pc_curve, kr_curve, k_prof, phi_prof):
    """Return (sat_n, kr_n_up, kr_w_up) sampled over plume heights h in [0, H]."""
    H = args.thickness
    drho = args.rho_w - args.rho_n
    if drho <= 0.0:
        sys.exit("error: rho_w must exceed rho_n (delta_rho > 0)")
    nq = args.n_quad

    # Column-total weights (denominators), evaluated once on a fine full-column grid.
    d_full = np.linspace(0.0, H, nq)
    k_denom = _trapezoid(k_prof(d_full), d_full)
    phi_denom = _trapezoid(phi_prof(d_full), d_full)

    hs = np.linspace(0.0, H, args.n_points)
    sat_n = np.zeros_like(hs)
    kr_n_up = np.zeros_like(hs)
    kr_w_up = np.zeros_like(hs)

    for i, h in enumerate(hs):
        if h <= 0.0:
            # Empty column: no CO2, brine fully mobile over the whole thickness.
            sat_n[i] = 0.0
            kr_n_up[i] = 0.0
            kr_w_up[i] = kr_curve.krw(1.0)
            continue

        d = np.linspace(0.0, h, nq)           # within-plume depths
        k = k_prof(d)
        phi = phi_prof(d)
        if args.mode == "sharp":
            sw = np.full_like(d, args.swr)    # box: S_n = 1 - swr
        else:
            pc = drho * args.gravity * (h - d)  # zeta = h - d above the contact
            sw = pc_curve.sw_of_pc(pc)

        sn = 1.0 - sw
        sat_n[i] = _trapezoid(phi * sn, d) / phi_denom
        kr_n_up[i] = _trapezoid(k * kr_curve.krn(sw), d) / k_denom
        # Brine: mobile in the plume (at sw) AND fully saturated below it (kr_w(1)).
        k_below = k_denom - _trapezoid(k_prof(d), d)
        kr_w_up[i] = (_trapezoid(k * kr_curve.krw(sw), d) + kr_curve.krw(1.0) * k_below) / k_denom

    # VERelPermTableUO needs strictly increasing sat_n. sat_n(h) is monotone; guard
    # against tiny non-monotone steps from quadrature noise by keeping a strict subset.
    keep = [0]
    for i in range(1, len(sat_n)):
        if sat_n[i] > sat_n[keep[-1]] + 1e-12:
            keep.append(i)
    keep = np.array(keep)
    return sat_n[keep], kr_n_up[keep], kr_w_up[keep]


def pc_table(args, pc_curve):
    """Return (pc_points, sw_points) for VECapillaryPressureTableUO: pc strictly
    increasing from 0, sw strictly decreasing, spanning [0, delta_rho*g*H]."""
    drho = args.rho_w - args.rho_n
    pc_max = drho * args.gravity * args.thickness

    if isinstance(pc_curve, PcTable):
        if pc_curve.pc[-1] < pc_max - 1e-9:
            print(f"warning: --pc-table tops out at Pc = {pc_curve.pc[-1]:.0f} Pa but the "
                  f"buoyancy head over H is {pc_max:.0f} Pa; Sw is clamped to "
                  f"{pc_curve.sw[-1]:.3f} above the table. Extend the table to cover the "
                  "full column.", file=sys.stderr)
        # Use the supplied table, restricted to [0, pc_max] (plus a clamp point).
        pc = pc_curve.pc[pc_curve.pc <= pc_max]
        sw = pc_curve.sw[: len(pc)]
        if pc.size == 0 or pc[0] > 0.0:
            pc = np.concatenate(([0.0], pc))
            sw = np.concatenate(([1.0], sw))
        return pc, sw

    # Parametric: sample Sw uniformly from 1 down to Sw(pc_max) and invert to Pc(Sw).
    # Capping at Sw(pc_max) -- NOT swr -- keeps every point inside [0, pc_max] (the
    # only range the column ever sees) and lands the last point exactly on the full
    # buoyancy head; sampling to swr instead would send low-Sw points past pc_max,
    # where clamping/dedup would drop them and leave the high-Pc end unresolved.
    # Uniform-in-Sw naturally clusters Pc near 0, where Sw(Pc) is steepest. Pc=0 is
    # anchored (Sw=1); a Brooks-Corey entry-pressure plateau collapses to that point.
    sw_end = float(pc_curve.sw_of_pc(np.array([pc_max]))[0])
    sw = np.linspace(1.0, sw_end, args.n_points)
    pc = pc_curve.pc_of_sw(sw)
    if pc[0] > 1e-9:  # Brooks-Corey: pc_of_sw(1) = pd > 0 -> add the Pc=0 anchor
        pc = np.concatenate(([0.0], pc))
        sw = np.concatenate(([1.0], sw))
    else:
        pc[0] = 0.0
    # Drop any points that collapse (strict monotonicity for VECapillaryPressureTableUO).
    keep = [0]
    for i in range(1, len(pc)):
        if pc[i] > pc[keep[-1]] + 1e-9 and sw[i] < sw[keep[-1]] - 1e-9:
            keep.append(i)
    keep = np.array(keep)
    return pc[keep], sw[keep]


# ----------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------
def _fmt(vec):
    """Space-separated values wrapped inside a single MOOSE quoted string."""
    parts = [f"{v:.8g}" for v in vec]
    lines, cur = [], "    "
    for p in parts:
        if len(cur) + len(p) + 1 > 96:
            lines.append(cur)
            cur = "    "
        cur += p + " "
    lines.append(cur.rstrip())
    return "\n".join(lines)


def write_outputs(args, sat_n, kr_n_up, kr_w_up, pc, sw):
    rp = f"{args.out_prefix}_relperm.csv"
    np.savetxt(rp, np.column_stack((sat_n, kr_n_up, kr_w_up)), delimiter=",",
               header="sat_n,kr_n_up,kr_w_up", comments="", fmt="%.10e")
    pcf = None
    if args.mode == "fringe":
        pcf = f"{args.out_prefix}_pc.csv"
        np.savetxt(pcf, np.column_stack((pc, sw)), delimiter=",",
                   header="pc,sw", comments="", fmt="%.10e")

    inc = f"{args.out_prefix}_upscaled.i"
    with open(inc, "w") as f:
        f.write("# Upscaled VE curves generated by ve_upscale_curves.py\n")
        f.write(f"# mode={args.mode}  thickness={args.thickness} m  "
                f"delta_rho={args.rho_w - args.rho_n} kg/m3  swr={args.swr}\n")
        f.write("# !include this file, or paste the blocks into [UserObjects].\n")
        f.write("[UserObjects]\n")
        f.write("  [relperm_uo]\n")
        f.write("    type = VERelPermTableUO\n")
        f.write(f"    sat_n_points = '\n{_fmt(sat_n)}'\n")
        f.write(f"    krn_points = '\n{_fmt(kr_n_up)}'\n")
        f.write(f"    krw_points = '\n{_fmt(kr_w_up)}'\n")
        f.write("  []\n")
        if args.mode == "fringe":
            f.write("  [pc_uo] # for capillary = capillary_fringe (VEPlumeReconstruction / cap pressure)\n")
            f.write("    type = VECapillaryPressureTableUO\n")
            f.write(f"    sat_lr = {args.swr:.8g}\n")
            f.write(f"    pc_points = '\n{_fmt(pc)}'\n")
            f.write(f"    sw_points = '\n{_fmt(sw)}'\n")
            f.write("  []\n")
        f.write("[]\n")
    return inc, rp, pcf


# ----------------------------------------------------------------------------------
def build_curves(args):
    # The Pc curve drives the fringe column integral and the Pc table; 'sharp' mode
    # uses the box profile S_n = 1 - swr and emits no Pc table, so it needs no Pc curve.
    if args.mode == "sharp":
        pc_curve = None
    elif args.pc == "vg":
        pc_curve = PcVanGenuchten(args.vg_alpha, args.vg_m, args.swr)
    elif args.pc == "bc":
        pc_curve = PcBrooksCorey(args.bc_pd, args.bc_lambda, args.swr)
    else:
        if not args.pc_table:
            sys.exit("error: --pc table requires --pc-table CSV")
        pc_curve = PcTable(args.pc_table, args.swr)

    if args.kr == "corey":
        kr_curve = KrCorey(args.krn_max, args.krw_max, args.corey_nn, args.corey_nw,
                           args.swr, args.snr)
    else:
        if not args.kr_table:
            sys.exit("error: --kr table requires --kr-table CSV")
        kr_curve = KrTable(args.kr_table)

    return pc_curve, kr_curve


def main():
    p = argparse.ArgumentParser(
        description="Upscale fine-scale rock kr/Pc curves into VE depth-integrated "
                    "VERelPermTableUO / VECapillaryPressureTableUO tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--out-prefix", required=True,
                   help="output file prefix (<prefix>_upscaled.i, _relperm.csv, _pc.csv)")
    p.add_argument("--mode", choices=["fringe", "sharp"], default="fringe",
                   help="vertical saturation model integrated over the column")

    g = p.add_argument_group("column")
    g.add_argument("--thickness", type=float, required=True, help="formation thickness H [m]")
    g.add_argument("--rho-n", type=float, default=700.0, help="CO2 density [kg/m3]")
    g.add_argument("--rho-w", type=float, default=1000.0, help="brine density [kg/m3]")
    g.add_argument("--gravity", type=float, default=9.81, help="|g| [m/s2]")
    g.add_argument("--swr", type=float, default=0.0, help="residual water saturation [-]")
    g.add_argument("--snr", type=float, default=0.0,
                   help="residual CO2 saturation [-] (Corey relperm only)")
    g.add_argument("--k-profile", help="CSV depth_below_top[m],k[-] (default uniform)")
    g.add_argument("--phi-profile", help="CSV depth_below_top[m],phi[-] (default uniform)")

    g = p.add_argument_group("capillary pressure")
    g.add_argument("--pc", choices=["vg", "bc", "table"], default="vg",
                   help="fine-scale Pc model")
    g.add_argument("--vg-alpha", type=float, help="van Genuchten alpha [1/Pa]")
    g.add_argument("--vg-m", type=float, help="van Genuchten m in (0,1)")
    g.add_argument("--bc-pd", type=float, help="Brooks-Corey entry pressure pd [Pa]")
    g.add_argument("--bc-lambda", type=float, help="Brooks-Corey pore-size index lambda")
    g.add_argument("--pc-table", help="CSV Pc[Pa],Sw[-] for --pc table")

    g = p.add_argument_group("relative permeability")
    g.add_argument("--kr", choices=["corey", "table"], default="corey",
                   help="fine-scale relperm model")
    g.add_argument("--krn-max", type=float, default=1.0, help="CO2 endpoint kr [-]")
    g.add_argument("--krw-max", type=float, default=1.0, help="brine endpoint kr [-]")
    g.add_argument("--corey-nn", type=float, default=2.0, help="CO2 Corey exponent")
    g.add_argument("--corey-nw", type=float, default=2.0, help="brine Corey exponent")
    g.add_argument("--kr-table", help="CSV Sw[-],krw[-],krn[-] for --kr table")

    g = p.add_argument_group("sampling")
    g.add_argument("--n-points", type=int, default=50,
                   help="number of table points (plume-height / saturation samples)")
    g.add_argument("--n-quad", type=int, default=400,
                   help="vertical quadrature points per column integral")

    args = p.parse_args()
    if args.swr < 0 or args.swr >= 1:
        sys.exit("error: --swr must be in [0, 1)")
    if args.n_points < 3:
        sys.exit("error: --n-points must be >= 3")

    pc_curve, kr_curve = build_curves(args)
    k_prof = Profile1D(args.k_profile, "--k-profile")
    phi_prof = Profile1D(args.phi_profile, "--phi-profile")

    sat_n, kr_n_up, kr_w_up = upscale(args, pc_curve, kr_curve, k_prof, phi_prof)
    pc, sw = pc_table(args, pc_curve) if args.mode == "fringe" else (None, None)
    inc, rp, pcf = write_outputs(args, sat_n, kr_n_up, kr_w_up, pc, sw)

    print(f"mode               = {args.mode}")
    print(f"sat_n range        = [{sat_n[0]:.4f}, {sat_n[-1]:.4f}]  ({len(sat_n)} points)")
    print(f"kr_n_up range      = [{kr_n_up[0]:.4f}, {kr_n_up[-1]:.4f}]")
    print(f"kr_w_up range      = [{kr_w_up[0]:.4f}, {kr_w_up[-1]:.4f}]")
    if args.mode == "fringe":
        print(f"Pc table           = {len(pc)} points, Pc in [0, {pc[-1]:.1f}] Pa, "
              f"Sw in [{sw[-1]:.3f}, 1]")
    print(f"wrote {inc}")
    print(f"wrote {rp}")
    if pcf:
        print(f"wrote {pcf}")


if __name__ == "__main__":
    main()
