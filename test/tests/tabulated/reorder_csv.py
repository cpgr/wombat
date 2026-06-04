#!/usr/bin/env python3
"""Write column-permuted copies of the committed gold kr/Pc tables.

This supports the column-order-independence test: VERelPermTableUO and
VECapillaryPressureTableUO match CSV columns by header name (not position), so a
file whose columns are in a different order -- yet carries the same data under
the same headers -- must produce identical results.

Reads the committed gold tables and writes, in this directory:
    curves_relperm_reordered.csv  columns kr_w_up, sat_n, kr_n_up
    curves_pc_reordered.csv       columns sw, pc

Both are scrambled relative to the ve_upscale_curves.py native ordering. The
generated files are run artifacts (gitignored); only gold/ is committed.
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def reorder(src, dst, order):
    with open(os.path.join(HERE, src), newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    index = {name: i for i, name in enumerate(header)}
    missing = [c for c in order if c not in index]
    if missing:
        sys.exit(f"FAIL: {src} is missing expected column(s) {missing} (header={header})")
    with open(os.path.join(HERE, dst), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(order)
        for r in data:
            w.writerow([r[index[c]] for c in order])


def main():
    reorder(os.path.join("gold", "curves_relperm.csv"),
            "curves_relperm_reordered.csv",
            ["kr_w_up", "sat_n", "kr_n_up"])
    reorder(os.path.join("gold", "curves_pc.csv"),
            "curves_pc_reordered.csv",
            ["sw", "pc"])


if __name__ == "__main__":
    main()
