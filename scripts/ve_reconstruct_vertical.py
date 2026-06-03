#!/usr/bin/env python3
"""Reconstruct the vertical CO2 saturation profile from a 2D VE Exodus result.

The wombat VE solver runs on a 2D mesh (the formation top surface) and carries the
plume only through depth-averaged fields: the mobile-plume thickness ve_h, the
density contrast delta_rho, the geometry z_top / H, and -- if hysteresis is on -- the
historical maximum saturation sat_n_max.  Under the VE assumption each vertical column
is in instantaneous capillary-gravity equilibrium, so the full 3D saturation field can
be reconstructed *after the fact* as a pointwise function of vertical position.  This
script does that for every timestep stored in the 2D Exodus file and writes a transient
3D Exodus on a static extruded column mesh (no solve -- pure evaluate-and-write).

Vertical profile (let xi be depth below the top surface, 0 <= xi <= H; the mobile
CO2-brine contact is at xi = h = ve_h; zeta = h - xi is height above the contact):

  sharp_interface : S_n(xi) = 1 - S_wr           for 0 <= xi <= h, else 0
  capillary_fringe: S_n(xi) = 1 - Sw(delta_rho*g*zeta)  for 0 <= xi <= h, else 0

with delta_rho = rho_w - rho_n.  The fringe form is exactly the integrand of
VEPlumeReconstruction::satNBar -- here it is evaluated pointwise instead of integrated,
so no Newton inversion is needed (the model already computed h).  Sw(Pc) is supplied
either as a table (recommended: export the same curve the run used, e.g. the
VECapillaryPressureTable CSV) or as analytic van Genuchten.

The density contrast is taken from the TWO phase-density fields (rho_n, rho_w) rather
than a precomputed delta_rho, so the 2D file stays self-documenting -- a user inspecting
it sees both physical densities, not an opaque difference.  Output them from the run with
two ADMaterialStdVectorAux on ve_density (index 0 = non-wetting, 1 = wetting).

Optional trapped tail (post-injection): where the plume has receded, residual CO2 is
left behind down to the historical maximum extent h_max, at a constant trapped
saturation S_trap.  h_max is taken from --h-max-var if present, else estimated from
--sat-n-max-var with the sharp formula h_max = sat_n_max*H/(1-S_wr).

Requires the SEACAS exodus Python module (exodus.py / exodus3), the same one shipped
with MOOSE's PETSc/SEACAS.  meshio is deliberately NOT used.

Example:
  ./ve_reconstruct_vertical.py run_out.e plume3d.e \\
      --mode fringe --pc table --pc-table swpc.csv \\
      --thickness-var H --h-var ve_h --rho-n-var rho_n --rho-w-var rho_w \\
      --s-wr 0.2 --n-layers 40
"""

import argparse
import os
import sys

import numpy as np


# --------------------------------------------------------------------------------------
# exodus.py import (module name varies across SEACAS versions)
# --------------------------------------------------------------------------------------
def _import_exodus():
    for mod in ("exodus3", "exodus"):
        try:
            m = __import__(mod)
            return getattr(m, "exodus")
        except (ImportError, AttributeError):
            continue
    sys.exit(
        "error: could not import the SEACAS exodus module (tried 'exodus3' and "
        "'exodus').\nMake sure SEACAS is on PYTHONPATH, e.g.\n"
        "  export PYTHONPATH=$SEACAS_DIR/lib:$PYTHONPATH"
    )


Exodus = _import_exodus()


# --------------------------------------------------------------------------------------
# Sw(Pc) providers (capillary_fringe mode only)
# --------------------------------------------------------------------------------------
class SwTable:
    """Piecewise-linear Sw(Pc) from a 2-column CSV (Pc[Pa], Sw[-]), Pc ascending.

    Mirrors VECapillaryPressureTable: export the same curve the run consumed for an
    exact reconstruction.
    """

    def __init__(self, csv_path):
        data = np.loadtxt(csv_path, delimiter=",", comments="#")
        if data.ndim != 2 or data.shape[1] < 2:
            sys.exit("error: --pc-table must have two columns: Pc,Sw")
        order = np.argsort(data[:, 0])
        self.pc = data[order, 0]
        self.sw = data[order, 1]

    def __call__(self, pc):
        # np.interp clamps to the end values outside the table range, which is the
        # desired behaviour: below entry pressure Sw=1, above the last point Sw->S_wr.
        return np.interp(pc, self.pc, self.sw)


class SwVanGenuchten:
    """Analytic van Genuchten effective wetting saturation, matching
    PorousFlowCapillaryPressureVG.saturation(pc): Se = (1 + (alpha*pc)^(1/(1-m)))^(-m)
    for pc > 0, else 1.  For verification cases; prefer a table on real fields.
    """

    def __init__(self, alpha, m):
        if not (0.0 < m < 1.0):
            sys.exit("error: --vg-m must be in (0, 1)")
        self.alpha = alpha
        self.m = m

    def __call__(self, pc):
        pc = np.asarray(pc, dtype=float)
        se = np.ones_like(pc)
        pos = pc > 0.0
        n = 1.0 / (1.0 - self.m)
        se[pos] = (1.0 + (self.alpha * pc[pos]) ** n) ** (-self.m)
        return se


# --------------------------------------------------------------------------------------
# 2D reader helpers
# --------------------------------------------------------------------------------------
class Mesh2D:
    """Reads the 2D mesh and exposes nodal field access (elemental fields are averaged
    to nodes so the reconstruction can run entirely on nodes).
    """

    def __init__(self, path):
        self.exo = Exodus(path, mode="r", array_type="numpy")
        self.x, self.y, _z = self.exo.get_coords()
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        self.nn = self.x.size
        self.times = np.asarray(self.exo.get_times(), dtype=float)
        self.node_vars = list(self.exo.get_node_variable_names() or [])
        self.elem_vars = list(self.exo.get_element_variable_names() or [])

        # Per-block CORNER connectivity (0-based), and an output element type.
        self.blocks = []  # list of dict(id, out_type, corners[Nelem, ncorner])
        for blk_id in self.exo.get_elem_blk_ids():
            elem_type, n_elem, n_npe, _n_attr = self._blk_info(blk_id)
            conn, _ne, _npe = self.exo.get_elem_connectivity(blk_id)
            conn = np.asarray(conn, dtype=np.int64).reshape(n_elem, n_npe) - 1  # 0-based
            et = elem_type.upper()
            if et.startswith("QUAD"):
                ncorner, out_type = 4, "HEX8"
            elif et.startswith("TRI"):
                ncorner, out_type = 3, "WEDGE6"
            else:
                sys.exit(
                    f"error: unsupported 2D element type '{elem_type}'. "
                    "Only QUAD* and TRI* are handled."
                )
            self.blocks.append(
                dict(id=blk_id, out_type=out_type, corners=conn[:, :ncorner])
            )

        # Node incidence (for elemental -> nodal averaging), built once.
        self._inc_count = np.zeros(self.nn)
        for blk in self.blocks:
            c = blk["corners"]
            for col in range(c.shape[1]):
                np.add.at(self._inc_count, c[:, col], 1.0)
        self._inc_count[self._inc_count == 0.0] = 1.0

    def _blk_info(self, blk_id):
        # exodus.py reader is elem_blk_info (no get_ prefix); returns
        # (elem_type, num_blk_elems, num_elem_nodes, num_elem_attrs). elem_type may
        # come back as bytes (e.g. b'QUAD4') in some builds -> normalise to str.
        info = self.exo.elem_blk_info(blk_id)
        elem_type = info[0]
        if isinstance(elem_type, bytes):
            elem_type = elem_type.decode()
        return elem_type, int(info[1]), int(info[2]), int(info[3])

    def has(self, name):
        return name in self.node_vars or name in self.elem_vars

    def nodal(self, name, step):
        """Nodal values of `name` at 1-based time `step` (elemental averaged to nodes)."""
        if name in self.node_vars:
            return np.asarray(
                self.exo.get_node_variable_values(name, step), dtype=float
            )
        if name in self.elem_vars:
            acc = np.zeros(self.nn)
            for blk in self.blocks:
                vals = np.asarray(
                    self.exo.get_element_variable_values(blk["id"], name, step),
                    dtype=float,
                )
                c = blk["corners"]
                for col in range(c.shape[1]):
                    np.add.at(acc, c[:, col], vals)
            return acc / self._inc_count
        sys.exit(f"error: variable '{name}' is not in the 2D file")

    def close(self):
        self.exo.close()


# --------------------------------------------------------------------------------------
# Reconstruction
# --------------------------------------------------------------------------------------
def reconstruct_column_saturation(xi, H, h, drho, args, sw_curve, h_max):
    """CO2 saturation at depth `xi` (array, per node) given per-node column fields.

    Returns (sat, region): region 0=brine, 1=mobile CO2, 2=residually trapped.
    """
    sat = np.zeros_like(xi)
    region = np.zeros_like(xi)

    # Strict inequality: a zero-thickness column (h = 0) has NO mobile nodes, so the
    # top layer (xi = 0) is not painted in empty columns. With xi <= h the test
    # 0 <= 0 would mark every column's top node, sheeting the whole top surface.
    mobile = xi < h
    if args.mode == "sharp":
        sat[mobile] = 1.0 - args.s_wr
    else:
        zeta = np.clip(h - xi, 0.0, None)  # height above the contact
        pc = drho * args.gravity * zeta
        sat[mobile] = (1.0 - sw_curve(pc))[mobile]
    region[mobile] = 1.0

    if h_max is not None and args.s_trap > 0.0:
        tail = (~mobile) & (xi < h_max)
        sat[tail] = args.s_trap
        region[tail] = 2.0

    # Never exceed the column-fillable fraction (numerical safety on the curve tails).
    np.clip(sat, 0.0, 1.0 - args.s_wr if args.mode == "sharp" else 1.0, out=sat)
    return sat, region


def build_curve(args, mesh):
    if args.mode == "sharp":
        return None
    if args.pc == "table":
        if not args.pc_table:
            sys.exit("error: --mode fringe --pc table requires --pc-table CSV")
        return SwTable(args.pc_table)
    if args.pc == "vg":
        if args.vg_alpha is None or args.vg_m is None:
            sys.exit("error: --mode fringe --pc vg requires --vg-alpha and --vg-m")
        return SwVanGenuchten(args.vg_alpha, args.vg_m)
    sys.exit("error: --pc must be 'table' or 'vg'")


def geometry(mesh, args):
    """Per-node z_top and thickness H (read at the first step; geometry is static)."""
    z_top = mesh.nodal(args.z_top_var, 1)
    if mesh.has(args.thickness_var):
        H = mesh.nodal(args.thickness_var, 1)
    elif args.z_bottom_var and mesh.has(args.z_bottom_var):
        H = np.abs(z_top - mesh.nodal(args.z_bottom_var, 1))
    else:
        sys.exit(
            "error: need formation thickness -- provide --thickness-var "
            f"(have '{args.thickness_var}'?) or --z-bottom-var"
        )
    if np.any(H <= 0.0):
        sys.exit("error: non-positive formation thickness H found")
    return z_top, H


def extrude_nodes(mesh, z_top, H, args):
    """Static 3D node coordinates: N+1 layers, layer k at depth xi = (k/N)*H per column.

    z increases upward by default (vertical=up -> bottom at z_top - H); pass
    vertical=down for a depth-positive convention (bottom at z_top + H).
    """
    N = args.n_layers
    nn2 = mesh.nn
    sign = -1.0 if args.vertical == "up" else 1.0
    frac = (np.arange(N + 1) / N)[:, None]  # (N+1, 1)
    xi = frac * H[None, :]  # (N+1, nn2) depth below top per layer/column
    x3 = np.tile(mesh.x, N + 1)
    y3 = np.tile(mesh.y, N + 1)
    z3 = (z_top[None, :] + sign * xi).reshape(-1)
    return x3, y3, z3, xi  # xi shape (N+1, nn2)


def extrude_connectivity(mesh, args):
    """Per-block 3D connectivity (1-based, flat) by stacking N element layers.

    HEX8/WEDGE6 ordering: lower-z face first, then higher-z face, both in the 2D
    CCW corner order -> positive volume.  Which layer is lower-z depends on `vertical`.
    """
    N = args.n_layers
    nn2 = mesh.nn

    def gid(node0, layer):  # 1-based 3D node id
        return layer * nn2 + node0 + 1

    out = []  # list of (block_id, out_type, conn_flat, n_elem3)
    for blk in mesh.blocks:
        corners = blk["corners"]  # (nE, nc) 0-based
        ne, nc = corners.shape
        rows = []
        for k in range(N):
            lower, upper = (k + 1, k) if args.vertical == "up" else (k, k + 1)
            low = np.column_stack([gid(corners[:, c], lower) for c in range(nc)])
            up = np.column_stack([gid(corners[:, c], upper) for c in range(nc)])
            rows.append(np.hstack([low, up]))  # (nE, 2*nc)
        conn3 = np.vstack(rows)  # (N*nE, 2*nc)
        out.append((blk["id"], blk["out_type"], conn3.reshape(-1), N * ne))
    return out


def density_contrast(mesh, args, step):
    """Per-node delta_rho = rho_w - rho_n from the two phase-density fields, or a
    constant fallback.  Only needed in capillary_fringe mode.
    """
    if mesh.has(args.rho_n_var) and mesh.has(args.rho_w_var):
        drho = mesh.nodal(args.rho_w_var, step) - mesh.nodal(args.rho_n_var, step)
        if np.any(drho <= 0.0):
            sys.stderr.write(
                "warning: rho_w - rho_n <= 0 somewhere; the fringe profile assumes "
                "buoyant CO2 (delta_rho > 0)\n"
            )
        return drho
    if args.delta_rho is not None:
        return np.full(mesh.nn, args.delta_rho)
    sys.exit(
        "error: capillary_fringe needs the phase densities --rho-n-var and "
        f"--rho-w-var (have '{args.rho_n_var}', '{args.rho_w_var}'?), "
        "or an explicit constant --delta-rho"
    )


def trapped_extent(mesh, H, args, step):
    """Historical maximum mobile-plume thickness h_max per node, or None."""
    if not args.trapped:
        return None
    if args.h_max_var and mesh.has(args.h_max_var):
        return mesh.nodal(args.h_max_var, step)
    if args.sat_n_max_var and mesh.has(args.sat_n_max_var):
        return mesh.nodal(args.sat_n_max_var, step) * H / (1.0 - args.s_wr)
    sys.exit(
        "error: --trapped needs --h-max-var or --sat-n-max-var present in the file"
    )


def main():
    args = parse_args()
    mesh = Mesh2D(args.input)
    if mesh.times.size == 0:
        sys.exit("error: the 2D file has no time steps to reconstruct")

    sw_curve = build_curve(args, mesh)
    z_top, H = geometry(mesh, args)
    x3, y3, z3, xi = extrude_nodes(mesh, z_top, H, args)
    blocks3 = extrude_connectivity(mesh, args)

    nn3 = x3.size
    ne3 = sum(b[3] for b in blocks3)
    out_vars = [args.sat_name, args.region_name]

    # exodus.py refuses to open an existing file for write ("Cowardly not opening");
    # remove it so re-runs just overwrite.
    if os.path.exists(args.output):
        os.remove(args.output)

    out = Exodus(
        args.output,
        mode="w",
        array_type="numpy",
        title="VE vertical reconstruction",
        numDims=3,
        numNodes=nn3,
        numElems=ne3,
        numBlocks=len(blocks3),
        numNodeSets=0,
        numSideSets=0,
    )
    out.put_coord_names(["x", "y", "z"])
    out.put_coords(x3, y3, z3)
    for blk_id, out_type, conn_flat, ne in blocks3:
        npe = 8 if out_type == "HEX8" else 6
        out.put_elem_blk_info(blk_id, out_type, ne, npe, 0)
        out.put_elem_connectivity(blk_id, conn_flat)
    out.set_node_variable_number(len(out_vars))
    for i, name in enumerate(out_vars, start=1):
        out.put_node_variable_name(name, i)

    N = args.n_layers
    for step in range(1, mesh.times.size + 1):
        h = mesh.nodal(args.h_var, step)
        # Clip the numerical-noise halo: columns whose depth-averaged saturation is
        # below --sat-min are rendered as brine (h -> 0), so only the real plume
        # footprint is drawn. FE advection leaves a small non-zero sat_n ahead of the
        # front, which would otherwise reconstruct as a thin sheet over the surface.
        if args.sat_min > 0.0:
            sat_bar = mesh.nodal(args.sat_var, step)
            h = np.where(sat_bar < args.sat_min, 0.0, h)
        drho = (
            density_contrast(mesh, args, step)
            if args.mode == "fringe"
            else np.zeros(mesh.nn)
        )
        h_max = trapped_extent(mesh, H, args, step)

        sat3 = np.empty(nn3)
        reg3 = np.empty(nn3)
        for k in range(N + 1):
            sl = slice(k * mesh.nn, (k + 1) * mesh.nn)
            s, r = reconstruct_column_saturation(
                xi[k], H, h, drho, args, sw_curve,
                None if h_max is None else h_max,
            )
            sat3[sl] = s
            reg3[sl] = r

        out.put_time(step, float(mesh.times[step - 1]))
        out.put_node_variable_values(args.sat_name, step, sat3)
        out.put_node_variable_values(args.region_name, step, reg3)

    out.close()
    mesh.close()
    print(
        f"wrote {args.output}: {nn3} nodes, {ne3} elems, "
        f"{mesh.times.size} steps, mode={args.mode}"
    )


def parse_args():
    p = argparse.ArgumentParser(
        description="Reconstruct the 3D vertical CO2 saturation from a 2D VE Exodus run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="2D VE Exodus result (with time steps)")
    p.add_argument("output", help="3D Exodus to write")

    p.add_argument("--mode", choices=["sharp", "fringe"], default="sharp",
                   help="vertical saturation model")
    p.add_argument("--n-layers", type=int, default=20,
                   help="number of vertical element layers")
    p.add_argument("--vertical", choices=["up", "down"], default="up",
                   help="z convention: 'up' -> bottom at z_top-H; 'down' -> z_top+H")

    # field names in the 2D file
    p.add_argument("--z-top-var", default="z_top", help="top-surface elevation field")
    p.add_argument("--thickness-var", default="H", help="formation thickness field")
    p.add_argument("--z-bottom-var", default="z_bottom",
                   help="bottom elevation field (used if --thickness-var absent)")
    p.add_argument("--h-var", default="ve_h", help="mobile-plume thickness field")
    p.add_argument("--rho-n-var", default="rho_n",
                   help="non-wetting (CO2) phase density field (fringe mode)")
    p.add_argument("--rho-w-var", default="rho_w",
                   help="wetting (brine) phase density field (fringe mode)")
    p.add_argument("--delta-rho", type=float, default=None,
                   help="constant rho_w - rho_n fallback if the density fields are absent")

    p.add_argument("--s-wr", type=float, default=0.0,
                   help="residual water saturation in the CO2 zone")
    p.add_argument("--gravity", type=float, default=9.81, help="|g| [m/s2]")

    # plume-footprint clip
    p.add_argument("--sat-min", type=float, default=0.0,
                   help="render only columns whose depth-averaged saturation (--sat-var) "
                        "is >= this; drops the thin numerical-noise halo (0 = no clip)")
    p.add_argument("--sat-var", default="sat_n",
                   help="depth-averaged saturation field used for the --sat-min clip")

    # capillary curve (fringe)
    p.add_argument("--pc", choices=["table", "vg"], default="table",
                   help="Sw(Pc) source (fringe mode)")
    p.add_argument("--pc-table", help="CSV Pc,Sw for --pc table")
    p.add_argument("--vg-alpha", type=float, help="van Genuchten alpha [1/Pa]")
    p.add_argument("--vg-m", type=float, help="van Genuchten m in (0,1)")

    # optional trapped tail
    p.add_argument("--trapped", action="store_true",
                   help="add a residually-trapped tail below the mobile plume")
    p.add_argument("--h-max-var", help="historical max plume thickness field")
    p.add_argument("--sat-n-max-var", default="sat_n_max",
                   help="historical max saturation (sharp h_max if --h-max-var absent)")
    p.add_argument("--s-trap", type=float, default=0.0,
                   help="residual trapped CO2 saturation in the tail")

    # output variable names
    p.add_argument("--sat-name", default="sat_co2",
                   help="output saturation variable name")
    p.add_argument("--region-name", default="co2_region",
                   help="output region flag (0 brine, 1 mobile, 2 trapped)")
    return p.parse_args()


if __name__ == "__main__":
    main()
