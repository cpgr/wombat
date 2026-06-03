#!/usr/bin/env python3
"""Generate dual_field.e: a flat VE mesh carrying z_top/z_bottom as BOTH nodal and
elemental fields under the same names.

This stands in for the output of `ex2ve --z-discretization both` (the upscaling
preprocessor). Its whole purpose is to pin a MOOSE framework contract that wombat
relies on: when a field name exists in BOTH the nodal and the elemental Exodus
namespaces, MOOSE resolves which to read from the receiving variable's type --
a LAGRANGE variable reads the nodal field, a MooseVariableFVReal reads the
elemental one -- so a single mesh drives the FE and FV solvers. read_both.i checks
that, and if a future MOOSE change breaks the disambiguation the CSVDiff fails.

MOOSE itself cannot write this fixture (two AuxVariables can't share the name
'z_top'), so it is generated here with netCDF4 directly. The mesh is a flat 10x10
QUAD4 grid over [0,100]^2 (z = 0) with:
    z_top    = -900 + 0.5*x + 0.3*y   (planar dip: 0.5 in x, 0.3 in y)
    z_bottom = z_top - 100            (uniform H = 100)
written nodally (at nodes) and elementally (at cell centroids).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

LEN_STRING = 33
NX = NY = 10
LX = LY = 100.0


def _ztop(x, y):
    return -900.0 + 0.5 * x + 0.3 * y


def _write_names(var, names):
    # S1 char array without stringtochar (broken in netCDF4 >= 1.7 + numpy >= 2.0).
    var[:, :] = np.array(
        [list(s.ljust(LEN_STRING)[:LEN_STRING]) for s in names], dtype="S1"
    )


def main() -> int:
    nxp, nyp = NX + 1, NY + 1
    dx, dy = LX / NX, LY / NY

    # Nodes (i fastest), flat (z = 0).
    ii, jj = np.meshgrid(np.arange(nxp), np.arange(nyp), indexing="xy")
    xs = (ii * dx).ravel()
    ys = (jj * dy).ravel()
    num_nodes = xs.size

    def nid(i, j):
        return j * nxp + i  # 0-indexed

    # QUAD4 connectivity, CCW from above.
    conn = []
    cxs, cys = [], []
    for j in range(NY):
        for i in range(NX):
            conn.append([nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)])
            cxs.append((i + 0.5) * dx)
            cys.append((j + 0.5) * dy)
    conn = np.array(conn, dtype=np.int32) + 1  # 1-indexed
    cxs = np.array(cxs); cys = np.array(cys)
    num_elem = conn.shape[0]

    zt_nod = _ztop(xs, ys)
    zb_nod = zt_nod - 100.0
    zt_elem = _ztop(cxs, cys)
    zb_elem = zt_elem - 100.0

    out = Path(__file__).with_name("dual_field.e")
    ds = Dataset(str(out), "w", format="NETCDF3_64BIT_OFFSET")
    try:
        ds.title = "wombat FE/FV dual-read fixture"
        ds.api_version = np.float32(4.98)
        ds.version = np.float32(4.98)
        ds.floating_point_word_size = np.int32(8)
        ds.file_size = np.int32(1)

        ds.createDimension("len_string", LEN_STRING)
        ds.createDimension("four", 4)
        ds.createDimension("num_dim", 3)
        ds.createDimension("num_nodes", num_nodes)
        ds.createDimension("num_elem", num_elem)
        ds.createDimension("num_el_blk", 1)
        ds.createDimension("num_el_in_blk1", num_elem)
        ds.createDimension("num_nod_per_el1", 4)
        ds.createDimension("num_elem_var", 2)
        ds.createDimension("num_nod_var", 2)
        ds.createDimension("time_step", None)

        ds.createVariable("coordx", "f8", ("num_nodes",))[:] = xs
        ds.createVariable("coordy", "f8", ("num_nodes",))[:] = ys
        ds.createVariable("coordz", "f8", ("num_nodes",))[:] = np.zeros(num_nodes)
        _write_names(ds.createVariable("coor_names", "S1", ("num_dim", "len_string")),
                     ["x", "y", "z"])

        eb_prop1 = ds.createVariable("eb_prop1", "i4", ("num_el_blk",))
        eb_prop1.setncattr("name", "ID")
        eb_prop1[:] = np.array([1], dtype=np.int32)
        ds.createVariable("eb_status", "i4", ("num_el_blk",))[:] = np.array([1], dtype=np.int32)
        _write_names(ds.createVariable("eb_names", "S1", ("num_el_blk", "len_string")),
                     ["ve_layer"])

        connect1 = ds.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
        connect1.elem_type = "QUAD4"
        connect1[:, :] = conn

        ds.createVariable("time_whole", "f8", ("time_step",))[0] = 0.0

        # Elemental z_top / z_bottom.
        _write_names(ds.createVariable("name_elem_var", "S1", ("num_elem_var", "len_string")),
                     ["z_top", "z_bottom"])
        ds.createVariable("elem_var_tab", "i4", ("num_el_blk", "num_elem_var"))[:, :] = 1
        ds.createVariable("vals_elem_var1eb1", "f8", ("time_step", "num_el_in_blk1"))[0, :] = zt_elem
        ds.createVariable("vals_elem_var2eb1", "f8", ("time_step", "num_el_in_blk1"))[0, :] = zb_elem

        # Nodal z_top / z_bottom (same names, different namespace).
        _write_names(ds.createVariable("name_nod_var", "S1", ("num_nod_var", "len_string")),
                     ["z_top", "z_bottom"])
        ds.createVariable("vals_nod_var1", "f8", ("time_step", "num_nodes"))[0, :] = zt_nod
        ds.createVariable("vals_nod_var2", "f8", ("time_step", "num_nodes"))[0, :] = zb_nod
    finally:
        ds.close()

    print(f"wrote {out} ({num_nodes} nodes, {num_elem} QUAD4; z_top/z_bottom nodal+elemental)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
