#!/usr/bin/env python3
"""Driver for the ve_upscale_curves.py regression / config RunCommand tests.

Exit 0 on success, non-zero (with a FAIL message) on mismatch. Usage:
    runcheck.py regression    # script output matches gold/curves_*.csv
    runcheck.py equivalence   # --config and the equivalent CLI give identical output
    runcheck.py override      # a CLI flag overrides the config value (CLI wins)

All runs happen in this directory; generated files (curves_*, from_*, ovr_*,
ref_*) are run artifacts (gitignored), except the committed gold/.
"""

import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "scripts", "ve_upscale_curves.py"))

# CLI flags equivalent to upscale.yaml -- kept in sync by the equivalence test.
BASE = dict(out_prefix="curves", mode="fringe", thickness="20", rho_n="700",
            rho_w="1000", swr="0.2", pc="vg", vg_alpha="5e-4", vg_m="0.5",
            kr="corey", krn_max="1", krw_max="1", corey_nn="2", corey_nw="2",
            n_points="20", n_quad="400")


def cli(**overrides):
    opts = dict(BASE)
    opts.update({k: str(v) for k, v in overrides.items()})
    args = []
    for key, value in opts.items():
        args += ["--" + key.replace("_", "-"), value]
    return args


def run(args):
    subprocess.run([sys.executable, SCRIPT] + args, check=True, cwd=HERE)


def load(name):
    return np.loadtxt(os.path.join(HERE, name), delimiter=",", skiprows=1)


def identical(a, b):
    return a.shape == b.shape and np.allclose(a, b, rtol=1e-9, atol=1e-12)


def fail(msg):
    sys.exit("FAIL: " + msg)


def main():
    case = sys.argv[1]

    if case == "regression":
        run(["--config", "upscale.yaml"])
        for name in ("curves_relperm.csv", "curves_pc.csv"):
            if not identical(load(name), load(os.path.join("gold", name))):
                fail(f"{name} differs from gold (ve_upscale_curves.py regressed)")

    elif case == "equivalence":
        run(["--config", "upscale.yaml", "--out-prefix", "from_cfg"])
        run(cli(out_prefix="from_cli"))
        for kind in ("relperm", "pc"):
            if not identical(load(f"from_cfg_{kind}.csv"), load(f"from_cli_{kind}.csv")):
                fail(f"config and CLI {kind} tables differ (config plumbing wrong)")

    elif case == "override":
        # --n-points on the CLI must override the config's 20: the result must equal
        # a pure-CLI run at the overridden value, and differ from the base config run.
        run(["--config", "upscale.yaml", "--n-points", "12", "--out-prefix", "ovr"])
        run(cli(out_prefix="ref", n_points="12"))
        for kind in ("relperm", "pc"):
            ovr, ref = load(f"ovr_{kind}.csv"), load(f"ref_{kind}.csv")
            if not identical(ovr, ref):
                fail(f"CLI override of n_points not applied to {kind} (CLI did not win)")
        if load("ovr_relperm.csv").shape[0] == load(os.path.join("gold",
                                                                 "curves_relperm.csv")).shape[0]:
            fail("override produced the same point count as the config (override was a no-op)")

    else:
        fail(f"unknown case '{case}'")


if __name__ == "__main__":
    main()
