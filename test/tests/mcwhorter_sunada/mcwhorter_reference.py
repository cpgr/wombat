#!/usr/bin/env python3
"""
McWhorter-Sunada semi-analytical reference for the VE capillary-diffusion term.

PROBLEM
-------
1-D counter-current capillary imbibition of CO2 into a brine-saturated VE column.
With per-phase-constant densities, a flat formation (grad z_T = 0) and a closed
far boundary, the total VOLUME flux q_t = q_n + q_w is uniform in x and pinned to
zero. Under q_t = 0 the depth-averaged CO2 saturation S = S_n obeys the nonlinear
diffusion equation

    phi dS/dt = d/dx [ D(S) dS/dx ],

with the VE capillary-dispersion coefficient

    D(S) = K * C * lambda_w(S) * lambda_n(S) / (lambda_w(S) + lambda_n(S)),
    C    = (rho_w - rho_n) * g * H / (1 - S_wr),     (= d Pc^up / dS)
    lambda_n = kr_n(S) / mu_n,   lambda_w = kr_w(S) / mu_w.

C is the slope of the upscaled capillary pressure Pc^up = (rho_w-rho_n) g h,
h = S H / (1 - S_wr), exactly as VEUpscaledCapPressure / VEFVCapPressure compute it.

For the verification case we use S_wr = 0, linear relperm (kr_n = S, kr_w = 1 - S)
and equal viscosities mu_n = mu_w = mu (mobility ratio M = 1), so

    D(S) = (K * C / mu) * S * (1 - S).

SELF-SIMILAR SOLUTION
---------------------
Boltzmann variable eta = x / sqrt(t) turns the PDE into the ODE BVP

    S'(eta)  = w / D(S)
    w'(eta)  = -(phi/2) * eta * w / D(S)          where w = D(S) S'

with S(0) = S0 (inlet) and S(eta_max) = Si (far field). Solved with solve_bvp.
The physical profile at time t_compare is S(x) = S_sim(x / sqrt(t_compare)).

OUTPUT
------
Writes a CSV (x [m], S_ref [-]) sampled at the simulation node positions for the
chosen t_compare, for use as a PiecewiseLinear MOOSE Function in an ElementL2Error
postprocessor. Run from this directory:

    python3 mcwhorter_reference.py

Keep the parameters in sync with mcwhorter_fe.i (and mcwhorter_fv.i).
"""

import numpy as np
from scipy.integrate import solve_bvp

# --- Physical parameters (MUST match the .i files) -------------------------
K = 1.0e-10        # upscaled permeability [m^2]
phi = 0.2          # porosity [-]
mu = 1.0e-3        # viscosity of both phases [Pa.s] (M = 1)
H = 10.0           # formation thickness [m]
g = 9.81           # gravity magnitude [m/s^2]
rho_n = 700.0      # CO2 density [kg/m^3]
rho_w = 1000.0     # brine density [kg/m^3]
S_wr = 0.0         # residual water saturation [-]

S0 = 0.7           # inlet CO2 saturation [-]
Si = 0.05          # initial / far-field CO2 saturation [-]

t_compare = 1.0e5  # time at which the profile is compared [s]
xmax = 100.0       # domain length [m]
nx = 200           # number of cells (node positions sampled at cell centres+ends)

dPc_dS = (rho_w - rho_n) * g * H / (1.0 - S_wr)   # C [Pa]


def D(S):
    """Capillary-dispersion coefficient D(S) = (K C / mu) S (1 - S)."""
    return (K * dPc_dS / mu) * S * (1.0 - S)


def odes(eta, y):
    """y[0] = S, y[1] = w = D(S) S'. Returns [S', w']."""
    S = np.clip(y[0], 1.0e-9, 1.0 - 1.0e-9)
    w = y[1]
    Dval = D(S)
    dS = w / Dval
    dw = -(phi / 2.0) * eta * w / Dval
    return np.vstack((dS, dw))


def bc(ya, yb):
    """S(0) = S0, S(eta_max) = Si."""
    return np.array([ya[0] - S0, yb[0] - Si])


def solve_profile():
    # eta_max: a few times the diffusion length sqrt(D/phi). Generous so the
    # far-field clamp S=Si is well inside the flat region.
    Dmax = D(0.5)
    eta_scale = np.sqrt(Dmax / phi)
    eta_max = 12.0 * eta_scale

    eta = np.linspace(0.0, eta_max, 400)
    # Initial guess: linear S from S0 to Si, w small negative (S decreasing).
    S_guess = S0 + (Si - S0) * (eta / eta_max)
    w_guess = -np.full_like(eta, Dmax) * (S0 - Si) / eta_max
    y_guess = np.vstack((S_guess, w_guess))

    sol = solve_bvp(odes, bc, eta, y_guess, max_nodes=200000, tol=1e-8)
    if not sol.success:
        raise RuntimeError("solve_bvp failed: " + sol.message)
    return sol, eta_max


def main():
    sol, eta_max = solve_profile()

    # Sample at the FE node positions (LAGRANGE FIRST on nx cells => nx+1 nodes).
    x = np.linspace(0.0, xmax, nx + 1)
    eta = x / np.sqrt(t_compare)
    eta_clipped = np.clip(eta, 0.0, eta_max)
    S_ref = sol.sol(eta_clipped)[0]
    # Beyond eta_max the solution is the far field Si.
    S_ref = np.where(eta <= eta_max, S_ref, Si)
    S_ref = np.clip(S_ref, Si, S0)

    out = np.column_stack((x, S_ref))
    # Header-less two-column CSV: x, S_ref. Matches MOOSE PiecewiseLinear
    # 'format = columns' (DelimitedFileReader expects numeric rows).
    np.savetxt("mcwhorter_reference.csv", out, delimiter=",",
               header="", comments="", fmt="%.10e")

    # Diagnostics.
    print(f"dPc/dS (C)      = {dPc_dS:.6e} Pa")
    print(f"D(0.5)          = {D(0.5):.6e} m^2/s")
    print(f"eta_max         = {eta_max:.6e}")
    print(f"front x (eta~scale) at t_compare ~ {np.sqrt(D(0.5)/phi)*np.sqrt(t_compare):.3f} m")
    print(f"S at x=0        = {S_ref[0]:.4f} (expect {S0})")
    print(f"S at x=xmax     = {S_ref[-1]:.4f} (expect ~{Si})")
    print(f"wrote mcwhorter_reference.csv with {len(x)} rows")


if __name__ == "__main__":
    main()
