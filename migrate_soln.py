"""Backfill a legacy soln.json to the current (opti_rev7) schema.

Runs made before opti_rev7 (e.g. run_h6ue, written by the old opti.py) lack the
stability & control derivatives and a few geometry fields the flight sim reads, so
they fall back to raw/assumed values and fly wrong. This recomputes everything the
run's OWN geometry determines and writes it back, so a legacy run flies faithfully
in the same pipeline as a fresh opti_rev7 run -- no re-optimization, same aircraft.

What is RECOVERED (from the run's real geometry, via AeroBuildup):
  Lateral_Directional_Stability: Cnb, Clb, Clp, Cnr, Clr, Cnp
  Aerodynamics: Cmq
What is DERIVED (documented analytic model, since the field was never stored):
  Main Wing: chord_tip / taper_ratio  (from S_w, b_w, root chord -- matches area)
What is ASSUMED (the legacy design never sized these -- flagged in the output):
  Main Wing: dihedral_angle (defaults 0 = as-designed flat unless --dihedral given)
  Cl_da (aileron), Cm_de (elevator): conventional-surface analytic estimate

Usage:
  python migrate_soln.py output/run_h6ue            # migrate one run (writes a .bak)
  python migrate_soln.py output/run_h6ue --dihedral 3
  python migrate_soln.py --all                      # migrate every legacy run under output/
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import aerosandbox as asb
from aerosandbox.atmosphere import Atmosphere

# Conventional control-surface sizing assumptions, matching opti_rev7's defaults.
# Used ONLY when the legacy run never designed the surface (flagged in output).
AILERON_Y_START_FRAC = 0.55
AILERON_Y_END_FRAC = 0.95
AILERON_CF_C = 0.25
ELEVATOR_CE_C = 0.30
TAIL_EFFICIENCY = 0.9
SECTION_LIFT_SLOPE = 2 * np.pi


def _airplane(soln: dict, chord_tip: float, dihedral_deg: float) -> asb.Airplane:
    """Build the AeroSandbox airplane from the soln geometry (matches the sim's
    reconstruction: real taper, airfoils, dihedral)."""
    mw, hs, vs, geo = soln["Main Wing"], soln["HStab"], soln["V Stab"], soln["Geometry"]
    span, chord_root, twist = mw["wingspan"], mw["chordlen"], mw["struct_defined_aoa"]
    boom_length, cg = geo["boom_length"], geo["cg_le_dist"]
    tip_z = (span / 2) * np.tan(np.radians(dihedral_deg))
    wf = str(mw.get("wing_airfoil", "s4110"))
    wing = asb.Wing(name="MainWing", symmetric=True, xsecs=[
        asb.WingXSec(xyz_le=[0, 0, 0], chord=chord_root, twist=twist, airfoil=asb.Airfoil(wf)),
        asb.WingXSec(xyz_le=[0, span / 2, tip_z], chord=chord_tip, twist=twist, airfoil=asb.Airfoil(wf)),
    ])
    hf = str(hs.get("hstab_airfoil", "naca0010"))
    hstab = asb.Wing(name="HStab", symmetric=True, xsecs=[
        asb.WingXSec(xyz_le=[0, 0, 0], chord=hs["hstab_chordlen"], twist=hs["hstab_aoa"], airfoil=asb.Airfoil(hf)),
        asb.WingXSec(xyz_le=[0, hs["hstab_span"] / 2, 0], chord=hs["hstab_chordlen"], twist=hs["hstab_aoa"], airfoil=asb.Airfoil(hf)),
    ]).translate([boom_length, 0, 0])
    vf = str(vs.get("vstab_airfoil", "naca0010"))
    vstab = asb.Wing(name="VStab", symmetric=False, xsecs=[
        asb.WingXSec(xyz_le=[0, 0, 0], chord=vs["vstab_root_chord"], twist=0, airfoil=asb.Airfoil(vf)),
        asb.WingXSec(xyz_le=[0, 0, vs["vstab_span"]], chord=hs["hstab_chordlen"], twist=0, airfoil=asb.Airfoil(vf)),
    ]).translate([boom_length, 0, 0])
    return asb.Airplane(name="rev6", wings=[wing, hstab, vstab], xyz_ref=[cg, 0, 0],
                        s_ref=float(mw["S_w"]), c_ref=float(mw["MAC_w"]), b_ref=float(mw["b_w"]))


def migrate(run_dir: Path, dihedral_deg: float | None = None, force: bool = False) -> bool:
    soln_path = run_dir / "soln.json"
    soln = json.loads(soln_path.read_text())
    if "Lateral_Directional_Stability" in soln and not force:
        print(f"  {run_dir.name}: already migrated (has Lateral_Directional_Stability) — skipping. Use --force to redo.")
        return False

    mw = soln["Main Wing"]
    S_w, b_w, chord_root = float(mw["S_w"]), float(mw["b_w"]), float(mw["chordlen"])

    # chord_tip: recover from stored area if not present. For a linear-taper wing
    # S = b*(c_root+c_tip)/2  ->  c_tip = 2S/b - c_root. This reproduces the design's
    # true area/MAC that a constant-chord reconstruction (c_tip=c_root) would inflate.
    if "chord_tip" in mw:
        chord_tip, tip_src = float(mw["chord_tip"]), "stored"
    else:
        chord_tip = max(0.02, 2.0 * S_w / b_w - chord_root)
        tip_src = "derived from S_w/b_w/root"
    taper = chord_tip / chord_root

    # dihedral: legacy runs didn't store one. Default flat (as-designed) unless given.
    if dihedral_deg is not None:
        dih, dih_src = float(dihedral_deg), "override"
    elif "dihedral_angle" in mw:
        dih, dih_src = float(mw["dihedral_angle"]), "stored"
    else:
        dih, dih_src = 0.0, "ASSUMED flat (run never designed dihedral)"

    # Stability derivatives from AeroBuildup on the real geometry, at cruise alpha.
    plane = _airplane(soln, chord_tip, dih)
    V = float(soln["Performance"]["airspeed"])
    atm = Atmosphere(altitude=float(soln["Mission"]["operating_altitude"]), temperature_deviation=34.0)
    op = asb.OperatingPoint(atmosphere=atm, velocity=V, alpha=float(mw["struct_defined_aoa"]))
    d = asb.AeroBuildup(airplane=plane, op_point=op).run_with_stability_derivatives(
        alpha=True, beta=True, p=True, q=True, r=True)

    def g(k):
        return float(np.asarray(d[k]).item())

    Cnb, Clb, Clp, Cnr, Clr, Cnp, Cmq = (g("Cnb"), g("Clb"), g("Clp"), g("Cnr"),
                                         g("Clr"), g("Cnp"), g("Cmq"))

    # Control-surface effectiveness: analytic, conventional assumptions (never sized).
    y1, y2 = AILERON_Y_START_FRAC * b_w / 2, AILERON_Y_END_FRAC * b_w / 2
    aileron_tau = 1.233 * AILERON_CF_C + 0.077
    Cl_da = (SECTION_LIFT_SLOPE * aileron_tau * chord_root / (S_w * b_w)) * (y2**2 - y1**2)
    hstab_AR = float(soln["HStab"]["hstab_AR"])
    a_hstab = 2 * np.pi * hstab_AR / (2 + hstab_AR)
    V_H = float(soln["HStab"].get("V_H_actual", 0.4))
    Cm_de = TAIL_EFFICIENCY * V_H * a_hstab * np.sqrt(ELEVATOR_CE_C)

    # Write back. Backup once.
    bak = soln_path.with_suffix(".json.bak")
    if not bak.exists():
        shutil.copy(soln_path, bak)

    mw["chord_tip"] = chord_tip
    mw["taper_ratio"] = taper
    mw["dihedral_angle"] = dih
    soln["Aerodynamics"]["Cmq"] = Cmq
    soln["Aerodynamics"]["Cm_de"] = Cm_de
    soln["Lateral_Directional_Stability"] = {
        "Cnb": Cnb, "Clb": Clb, "Clp": Clp, "Cnr": Cnr, "Clr": Clr, "Cnp": Cnp,
        "spiral_criterion": Clb * Cnr - Cnb * Clr,
        "dihedral_angle": dih,
        "Cl_da": Cl_da,
    }
    soln.setdefault("Meta", {})["migrated_from_legacy"] = True
    soln_path.write_text(json.dumps(soln, indent=2))

    print(f"  {run_dir.name}: migrated (backup -> {bak.name})")
    print(f"    chord_tip = {chord_tip:.4f} m  ({tip_src}); taper = {taper:.3f}")
    print(f"    dihedral  = {dih:.2f} deg  ({dih_src})")
    print(f"    Cnb={Cnb:+.4f}  Clb={Clb:+.4f}  Clp={Clp:+.4f}  Cnr={Cnr:+.4f}  Clr={Clr:+.4f}  Cnp={Cnp:+.4f}  Cmq={Cmq:+.4f}")
    print(f"    Cl_da={Cl_da:.4f}  Cm_de={Cm_de:.4f}   (ASSUMED conventional surfaces — run never sized them)")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", nargs="?", help="a run directory under output/ (e.g. output/run_h6ue)")
    ap.add_argument("--all", action="store_true", help="migrate every legacy run under output/")
    ap.add_argument("--dihedral", type=float, default=None, help="assume this wing dihedral (deg)")
    ap.add_argument("--force", action="store_true", help="re-migrate even if already migrated")
    args = ap.parse_args()

    if args.all:
        runs = sorted(p.parent for p in Path("output").glob("run_*/soln.json"))
    elif args.run_dir:
        runs = [Path(args.run_dir)]
    else:
        ap.error("give a run directory or --all")

    n = 0
    for r in runs:
        if migrate(r, dihedral_deg=args.dihedral, force=args.force):
            n += 1
    print(f"\nmigrated {n} run(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
