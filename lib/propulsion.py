"""Propeller thrust model — the single source of truth for prop physics.

This lives on the DESIGN side on purpose (separation of concerns): the optimizer
owns the propeller model, and downstream consumers (e.g. the JSBSim FDM built by
the sibling AircraftSim) receive only its RESULT — a thrust-vs-airspeed table baked
into ``soln.json`` — never a re-implementation of the math. If the prop fit changes,
it changes here and flows out through the artifact; nothing downstream needs editing.

Two entry points:
  * ``prop_ct_cp(J, pd)`` — the APC 16in thin-electric C_T/C_P polynomial fit. Pure
    arithmetic (only +,-,*,**), so it evaluates on CasADi symbolics during the solve
    AND on plain floats / numpy arrays when baking the artifact.
  * ``thrust_curve_from_soln(soln)`` — the airspeed thrust-lapse curve derived from
    that map, ready to serialize into ``soln.json`` under ``Propulsion.thrust_curve``.
"""

from __future__ import annotations

import numpy as _np


def prop_ct_cp(J, pd):
    """APC thin-electric 16in-family thrust/power coefficients (degree-3 fit).

    C_T(J,P/D) falls as advance ratio J rises -> at fixed RPM, thrust drops as the
    aircraft speeds up. Fit RMS ~1.6% (C_T) / ~3% (C_P) in the cruise core. Pure
    polynomial so it is C1-smooth for the optimizer and array-safe for baking.
    """
    C_T = (0.022025980141242307
           + 0.10182937396297323 * pd + 0.14703356855390493 * pd**2 - 0.22041866553079092 * pd**3
           - 0.10987881067866853 * J - 0.007686670565597137 * J*pd + 0.383584672564028 * J*pd**2
           - 0.09683894973386901 * J**2 - 0.34512648094474313 * J**2*pd + 0.15150046362483677 * J**3)
    C_P = (-0.0025623311770923388
           + 0.043004718882882834 * pd + 0.0779853604331156 * pd**2 - 0.10716202981227332 * pd**3
           + 0.014350796115607196 * J - 0.07374815736967645 * J*pd + 0.29351716144758694 * J*pd**2
           - 0.04782487322689711 * J**2 - 0.18781875948418347 * J**2*pd + 0.02743631110560021 * J**3)
    return C_T, C_P


def prop_coeff_table(soln: dict, n_points: int = 16) -> dict | None:
    """Raw propeller coefficient map C_T(J), C_P(J) + geometry, for a real JSBSim
    FGPropeller in the sim (higher fidelity than the normalized thrust_curve: the
    sim's electric engine + propeller then solve RPM and advance-ratio thrust
    themselves). Coefficients follow the standard convention T = C_T*rho*n^2*D^4,
    P = C_P*rho*n^3*D^5 (n rev/s) -- the same JSBSim uses -- so they drop straight
    into a <propeller> C_THRUST/C_POWER table. Per-prop values; the consumer scales
    by n_props if it lumps them into one engine. None if the soln lacks prop fields.
    """
    prop = soln.get("Propulsion", {}) or {}
    try:
        pd = float(prop["prop_pitch_over_diam"])
        D = float(prop["propeller_diameter"])
        n_props = int(round(float(prop.get("propeller_n", 1) or 1)))
    except (KeyError, TypeError, ValueError):
        return None
    if not (D > 0):
        return None
    J = _np.linspace(0.0, 1.3, n_points)
    CT, CP = prop_ct_cp(J, pd)
    return {
        "J": [round(float(j), 4) for j in J],
        "CT": [round(float(c), 6) for c in CT],
        "CP": [round(float(c), 6) for c in CP],
        "diameter_m": round(D, 5),
        "numblades": 2,
        "pitch_over_diam": round(pd, 4),
        "n_props": n_props,
    }


def thrust_curve_from_soln(soln: dict, n_points: int = 14) -> dict | None:
    """Airspeed thrust-lapse table for one solved design, normalized to 1.0 at cruise.

    Models the FULL-THROTTLE prop: solve for the RPM n_max that produces the design's
    own max thrust (``thrust_climb``) at cruise speed, then tabulate
    ``frac(V) = T_total(V, n_max) / thrust_climb`` over an airspeed grid. A consumer
    reproduces the actual thrust as ``throttle * max_thrust * frac(V)`` -- at cruise
    frac == 1, above cruise thrust lapses, below it rises (clamped so a near-static
    point can't blow up or go negative from a windmilling C_T<0).

    Returns a JSON-safe dict {"V_mps", "frac", "model"} or None if the soln lacks the
    prop fields (then the consumer keeps constant thrust). Reads plain floats, so call
    it AFTER the raw report has been resolved to numbers.
    """
    prop = soln.get("Propulsion", {}) or {}
    perf = soln.get("Performance", {}) or {}
    env = soln.get("Environment", {}) or {}
    try:
        pd = float(prop["prop_pitch_over_diam"])
        D = float(prop["propeller_diameter"])
        n_prop = float(prop.get("propeller_n", 1) or 1)
        derate = float(prop.get("prop_thrust_derate", 1.0) or 1.0)
        rpm_cruise = float(prop["rpm_cruise"])
        Vc = float(perf["airspeed"])
        T_max = float(perf["thrust_climb"])       # N, design max thrust (all motors)
        rho = float(env.get("density", 1.0))
    except (KeyError, TypeError, ValueError):
        return None
    if not (D > 0 and Vc > 0 and T_max > 0 and rho > 0 and rpm_cruise > 0):
        return None

    def T_total(V, n):   # n in rev/s; V scalar or array
        CT, _ = prop_ct_cp(V / (n * D), pd)
        return n_prop * derate * _np.maximum(CT, 0.0) * rho * n**2 * D**4

    # Solve T_total(Vc, n_max) = T_max by bisection (T rises monotonically with n:
    # higher n -> lower J -> higher C_T, and n^2). Bracket from cruise RPM upward.
    lo, hi = rpm_cruise / 60.0, 6.0 * rpm_cruise / 60.0
    if T_total(Vc, hi) < T_max:      # design max thrust unreachable on this prop -> skip
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if T_total(Vc, mid) < T_max:
            lo = mid
        else:
            hi = mid
    n_max = 0.5 * (lo + hi)

    V = _np.linspace(0.15 * Vc, 2.4 * Vc, n_points)
    frac = _np.clip(T_total(V, n_max) / T_max, 0.02, 1.8)
    return {
        "V_mps": [round(float(v), 5) for v in V],
        "frac": [round(float(f), 6) for f in frac],
        "model": "apc_ctj_fixed_rpm",
    }
