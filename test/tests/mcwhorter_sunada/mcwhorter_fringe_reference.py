#!/usr/bin/env python3
"""
McWhorter-Sunada semi-analytical reference for the CAPILLARY-FRINGE VE flux.

Fringe counterpart of mcwhorter_reference.py. Same 1-D counter-current
capillary-imbibition setup (per-phase-constant densities, flat formation,
closed far boundary => total volume flux q_t = 0), so the depth-averaged CO2
saturation S obeys

    phi dS/dt = d/dx [ D(S) dS/dx ].

The ONLY difference from the sharp case is the slope of the upscaled capillary
pressure that enters D(S). Sharp interface uses the constant
C = (rho_w - rho_n) g H / (1 - S_wr). Capillary fringe uses the S-dependent

    dPc^up/dS = (rho_w - rho_n) g H / S_n(h(S)),

where h(S) is the column thickness that reproduces the depth-averaged
saturation S (Newton/bisection inverse of the trapezoidal integral
S = (1/H) int_0^h [1 - Sw(Pc(zeta))] dzeta, Pc(zeta) = delta_rho g zeta), and
S_n(h) = 1 - Sw(delta_rho g h) is the CO2 saturation at the plume top -- exactly
what VEPlumeReconstruction + VEUpscaledCapPressure / VEFVCapPressure compute.
With linear relperm (kr_n = S, kr_w = 1 - S, S_wr = 0) and M = 1 this gives

    D_fringe(S) = (K/mu) S (1 - S) * (rho_w-rho_n) g H / S_n(h(S))
                = D_sharp(S) / S_n(h(S)).

Because S_n(h(S)) < 1, the fringe column imbibes FASTER than the sharp one -- a
clearly different profile, so this verifies the fringe S_n(h) factor (sign AND
magnitude), which the linear-table no-flow equilibrium (sat_n*H = const, shared
by both models) cannot.

The Sw(Pc) table MUST match the VECapillaryPressureTable in mcwhorter_fringe_fe.i.
Run from this directory:  python3 mcwhorter_fringe_reference.py
"""

import numpy as np
from scipy.integrate import solve_bvp

# --- Physical parameters (MUST match mcwhorter_fringe_fe.i) -----------------
K = 1.0e-10        # upscaled permeability [m^2]
phi = 0.2          # porosity [-]
mu = 1.0e-3        # viscosity of both phases [Pa.s] (M = 1)
H = 20.0           # formation thickness [m] (constant -> no grad(H) term)
g = 9.81           # gravity magnitude [m/s^2]
rho_n = 700.0      # CO2 density [kg/m^3]
rho_w = 1000.0     # brine density [kg/m^3]
delta_rho = rho_w - rho_n

# Linear Sw(Pc) table (must equal the VECapillaryPressureTable in the .i).
PC_POINTS = np.array([0.0, 15696.0, 31392.0, 47088.0, 62784.0])
SW_POINTS = np.array([1.0, 0.8, 0.6, 0.4, 0.2])

# Inlet / far-field CO2 saturation. Both representable by the table on an H = 20 m
# column (max column-average ~0.375 here) and chosen so S_n(h(S)) stays well away
# from zero (S_n ~ 0.27..0.61 over [Si, S0]) -- the fringe D ~ sqrt(S) degenerates
# as S -> 0, so Si is kept clear of it. The fringe coefficient is then 1.6..3.6x
# the sharp one over this band, a strong discriminator.
S0 = 0.30          # inlet CO2 saturation [-]
Si = 0.04          # initial / far-field CO2 saturation [-]

t_compare = 4.0e4  # time at which the profile is compared [s]. Chosen so the
                   # transition spans ~[0, 63] m: comfortably inside the 100 m
                   # domain (far field exactly Si, no boundary pile-up) yet wide
                   # enough that (a) the FE mesh resolves it well and (b) the
                   # fringe vs sharp-coefficient profiles separate by L2 ~ 0.022,
                   # which the gold abs_zero = 0.012 discriminates.
xmax = 100.0       # domain length [m]
nx = 400           # number of cells (sample at nx+1 node positions)

_SN_FLOOR = 1.0e-6


def sw(pc):
    """True water saturation from the (linear) Pc table; clamped at the ends."""
    return np.interp(pc, PC_POINTS, SW_POINTS)


# Precompute the monotone S(h) map and its inverse h(S) on a fine grid.
_h_grid = np.linspace(0.0, H, 4001)
_pc_grid = delta_rho * g * _h_grid
_sn_grid = 1.0 - sw(_pc_grid)                      # S_n(zeta) along the column
# Cumulative trapezoid of S_n over [0, h], divided by H -> depth-averaged S(h).
_cumS = np.concatenate(([0.0], np.cumsum(0.5 * (_sn_grid[1:] + _sn_grid[:-1]) *
                                         np.diff(_h_grid)))) / H
_Smax = _cumS[-1]


def h_of_S(S):
    """Plume thickness h reproducing depth-averaged saturation S (interp inverse)."""
    return np.interp(S, _cumS, _h_grid)


def Sn_top(S):
    """S_n(h(S)) = 1 - Sw(delta_rho g h(S)), floored to keep D finite."""
    sn = 1.0 - sw(delta_rho * g * h_of_S(S))
    return np.maximum(sn, _SN_FLOOR)


def D(S):
    """Fringe capillary-dispersion coefficient D_sharp(S) / S_n(h(S))."""
    S = np.clip(S, 1.0e-9, 1.0 - 1.0e-9)
    dPc_dS = delta_rho * g * H / Sn_top(S)          # S-dependent (vs sharp const)
    return (K / mu) * S * (1.0 - S) * dPc_dS


def odes(eta, y):
    S = np.clip(y[0], 1.0e-9, 1.0 - 1.0e-9)
    w = y[1]
    Dval = D(S)
    return np.vstack((w / Dval, -(phi / 2.0) * eta * w / Dval))


def bc(ya, yb):
    return np.array([ya[0] - S0, yb[0] - Si])


def solve_profile():
    Dmax = D(0.5 * (S0 + Si))
    eta_scale = np.sqrt(Dmax / phi)
    eta_max = 12.0 * eta_scale
    eta = np.linspace(0.0, eta_max, 400)
    S_guess = S0 + (Si - S0) * (eta / eta_max)
    w_guess = -np.full_like(eta, Dmax) * (S0 - Si) / eta_max
    sol = solve_bvp(odes, bc, eta, np.vstack((S_guess, w_guess)),
                    max_nodes=400000, tol=1e-8)
    if not sol.success:
        raise RuntimeError("solve_bvp failed: " + sol.message)
    return sol, eta_max


def mol_crosscheck(sol, eta_max, n=1600):
    """Independent method-of-lines solve of phi dS/dt = d/dx(D dS/dx), for a
    confidence cross-check of the similarity profile (not written out). Explicit
    Euler with a diffusion-stable adaptive dt; S is clamped to the representable
    band so a transient overshoot cannot hit the D ~ sqrt(S) degeneracy."""
    x = np.linspace(0.0, xmax, n + 1)
    dx = x[1] - x[0]
    S = np.full_like(x, Si)
    S[0] = S0
    t = 0.0
    while t < t_compare:
        Sf = 0.5 * (S[1:] + S[:-1])
        Df = D(Sf)
        # Effective diffusivity is D/phi (phi dS/dt = d/dx(D dS/dx)); the explicit
        # stability limit is dt < dx^2 / (2 D/phi).
        dt = min(0.4 * dx * dx * phi / max(Df.max(), 1e-30), t_compare - t)
        flux = Df * (S[1:] - S[:-1]) / dx
        dSdt = np.zeros_like(S)
        dSdt[1:-1] = (flux[1:] - flux[:-1]) / dx / phi
        S = np.clip(S + dt * dSdt, 0.0, _Smax)
        S[0] = S0
        S[-1] = S[-2]
        t += dt
    xe = x / np.sqrt(t_compare)
    S_sim = np.where(xe <= eta_max, sol.sol(np.clip(xe, 0.0, eta_max))[0], Si)
    return np.sqrt(np.mean((S - S_sim) ** 2))


def main():
    sol, eta_max = solve_profile()

    x = np.linspace(0.0, xmax, nx + 1)
    eta = x / np.sqrt(t_compare)
    S_ref = np.where(eta <= eta_max, sol.sol(np.clip(eta, 0.0, eta_max))[0], Si)
    S_ref = np.clip(S_ref, Si, S0)

    np.savetxt("mcwhorter_fringe_reference.csv", np.column_stack((x, S_ref)),
               delimiter=",", header="", comments="", fmt="%.10e")

    print(f"max representable S (h=H) = {_Smax:.4f}  (S0={S0}, Si={Si})")
    print(f"D(S0)   = {float(D(np.array([S0]))[0]):.4e} m^2/s  "
          f"S_n(h(S0)) = {float(Sn_top(np.array([S0]))[0]):.4f}")
    print(f"D(Si)   = {float(D(np.array([Si]))[0]):.4e} m^2/s")
    print(f"eta_max = {eta_max:.4e}")
    print(f"MOL cross-check L2(sim-similarity) = {mol_crosscheck(sol, eta_max):.3e}")
    print(f"S(x=0)  = {S_ref[0]:.4f} (expect {S0})")
    print(f"S(x=L)  = {S_ref[-1]:.4f} (expect ~{Si})")
    print(f"wrote mcwhorter_fringe_reference.csv with {len(x)} rows")


if __name__ == "__main__":
    main()
