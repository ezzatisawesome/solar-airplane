"""Reconstruct the full design airplane from a solved ``soln.json``.

The optimizer builds the true airframe (tapered/polyhedral wing, in-plane H-stab,
twin STRADDLING fins at +/- boom_y, twin pods + twin booms) at solve time. This
rebuilds that same topology from the exported scalar fields, so anything downstream
-- the GLB/CAD export here, or the sibling sim's aero model -- shares ONE geometry
definition instead of each inventing its own stripped surrogate.

Kept field-for-field faithful to opti_rev7's construction: same surfaces, same
positions, same straddling-fin placement. The wing is reconstructed with a
root -> boom-break -> tip trapezoid (taper + polyhedral outboard of the break),
which reproduces the design planform's area and shape for meshing and aero.
"""

from __future__ import annotations

import aerosandbox as asb
import aerosandbox.numpy as np


def build_airplane_from_soln(soln: dict) -> asb.Airplane:
    mw = soln["Main Wing"]
    hs = soln["HStab"]
    vs = soln["V Stab"]
    geo = soln["Geometry"]

    b = float(mw["wingspan"])
    chord_root = float(mw["chordlen"])
    chord_tip = float(mw.get("chord_tip", chord_root))
    twist = float(mw.get("struct_defined_aoa", 0.0))
    dihedral_deg = float(mw.get("dihedral_angle", 3.0))
    wing_foil = asb.Airfoil(str(mw.get("wing_airfoil", "s4310")))

    boom_y = float(geo["boom_y"])
    boom_length = float(geo["boom_length"])
    boom_radius = float(geo.get("boom_radius", 0.01))
    fuselage_length = float(geo.get("fuselage_length", 0.98))
    fuselage_radius = float(geo.get("fuselage_radius", 0.0281))
    cg_le_dist = float(geo.get("cg_le_dist", 0.0))

    # Wing: root -> boom break (full chord, flat) -> tip (tapered, raised by dihedral
    # outboard of the break). Matches opti_rev7's polyhedral-outboard planform.
    s = b / 2
    z_tip = np.sind(dihedral_deg) * max(s - boom_y, 0.0)
    main_wing = asb.Wing(
        name="Main Wing", symmetric=True,
        xsecs=[
            asb.WingXSec(xyz_le=[0, 0, 0], chord=chord_root, twist=twist, airfoil=wing_foil),
            asb.WingXSec(xyz_le=[0, boom_y, 0], chord=chord_root, twist=twist, airfoil=wing_foil),
            asb.WingXSec(xyz_le=[0, s, z_tip], chord=chord_tip, twist=twist, airfoil=wing_foil),
        ],
    )

    hstab_span = float(hs["hstab_span"])
    hstab_chordlen = float(hs["hstab_chordlen"])
    hstab_aoa = float(hs.get("hstab_aoa", 0.0))
    tail_foil = asb.Airfoil(str(hs.get("hstab_airfoil", "naca0010")))
    hor_stabilizer = asb.Wing(
        name="Horizontal Stabilizer", symmetric=True,
        xsecs=[
            asb.WingXSec(xyz_le=[0, 0, 0], chord=hstab_chordlen, twist=hstab_aoa, airfoil=tail_foil),
            asb.WingXSec(xyz_le=[0, hstab_span / 2, 0], chord=hstab_chordlen, twist=hstab_aoa, airfoil=tail_foil),
        ],
    ).translate([boom_length, 0, 0])

    vstab_span = float(vs["vstab_span"])
    vstab_root_chord = float(vs["vstab_root_chord"])
    vstab_foil = asb.Airfoil(str(vs.get("vstab_airfoil", "naca0010")))

    def _make_fin(name, sign):
        # Straddling fin: -h/2 (below) to +h/2 (above) the H-stab plane, at the boom,
        # LE aft of the H-stab TE. Same as opti_rev7's _make_fin.
        return asb.Wing(
            name=name, symmetric=False,
            xsecs=[
                asb.WingXSec(xyz_le=[0, 0, -0.5 * vstab_span], chord=vstab_root_chord, twist=0, airfoil=vstab_foil),
                asb.WingXSec(xyz_le=[0, 0, 0.5 * vstab_span], chord=vstab_root_chord, twist=0, airfoil=vstab_foil),
            ],
        ).translate([boom_length + hstab_chordlen, sign * boom_y, 0])

    vert_stabilizer_L = _make_fin("Vertical Stabilizer L", +1)
    vert_stabilizer_R = _make_fin("Vertical Stabilizer R", -1)

    def _make_boom(name, sign):
        return asb.Fuselage(name=name, xsecs=[
            asb.FuselageXSec(xyz_c=[boom_length * xi, sign * boom_y, -0.02], radius=boom_radius)
            for xi in np.cosspace(0, 1, 20)
        ])

    def _make_pod(name, sign):
        return asb.Fuselage(name=name, xsecs=[
            asb.FuselageXSec(
                xyz_c=[fuselage_length * xi - fuselage_length, sign * boom_y, -0.02],
                radius=fuselage_radius * min(1.0, 6.0 * xi * (1.0 - 0.4 * xi)),
            )
            for xi in np.cosspace(0, 1, 30)
        ])

    return asb.Airplane(
        name=str(soln.get("Meta", {}).get("run_id", "aircraft")),
        xyz_ref=[cg_le_dist, 0, 0],
        wings=[main_wing, hor_stabilizer, vert_stabilizer_L, vert_stabilizer_R],
        fuselages=[_make_pod("Left Fuse", +1), _make_pod("Right Fuselage", -1),
                   _make_boom("Boom L", +1), _make_boom("Boom R", -1)],
    )
