import json
import sys
from pathlib import Path
from typing import Dict, Tuple

import aerosandbox as asb
import numpy as onp
from aerosandbox.weights import MassProperties
from lib.artifacts import latest_run_dir, load_layout, state_path, update_state


def _mp_point(mass_kg: float, xyz: Tuple[float, float, float]) -> MassProperties:
    return MassProperties(mass=float(mass_kg), x_cg=float(xyz[0]), y_cg=float(xyz[1]), z_cg=float(xyz[2]))


def _mp_box(mass_kg: float, xyz_cg: Tuple[float, float, float], Lx: float, Ly: float, Lz: float) -> MassProperties:
    m = float(mass_kg)
    Lx = float(Lx)
    Ly = float(Ly)
    Lz = float(Lz)
    return MassProperties(
        mass=m,
        x_cg=float(xyz_cg[0]),
        y_cg=float(xyz_cg[1]),
        z_cg=float(xyz_cg[2]),
        Ixx=m / 12 * (Ly**2 + Lz**2),
        Iyy=m / 12 * (Lx**2 + Lz**2),
        Izz=m / 12 * (Lx**2 + Ly**2),
        Ixy=0.0,
        Iyz=0.0,
        Ixz=0.0,
    )


def _mp_cylinder_x(mass_kg: float, xyz_cg: Tuple[float, float, float], radius_m: float, length_m: float) -> MassProperties:
    m = float(mass_kg)
    r = float(radius_m)
    L = float(length_m)
    return MassProperties(
        mass=m,
        x_cg=float(xyz_cg[0]),
        y_cg=float(xyz_cg[1]),
        z_cg=float(xyz_cg[2]),
        Ixx=0.5 * m * r**2,
        Iyy=m / 12 * (3 * r**2 + L**2),
        Izz=m / 12 * (3 * r**2 + L**2),
        Ixy=0.0,
        Iyz=0.0,
        Ixz=0.0,
    )


def _mp_disk_x(mass_kg: float, xyz_cg: Tuple[float, float, float], radius_m: float) -> MassProperties:
    m = float(mass_kg)
    r = float(radius_m)
    return MassProperties(
        mass=m,
        x_cg=float(xyz_cg[0]),
        y_cg=float(xyz_cg[1]),
        z_cg=float(xyz_cg[2]),
        Ixx=0.5 * m * r**2,
        Iyy=0.25 * m * r**2,
        Izz=0.25 * m * r**2,
        Ixy=0.0,
        Iyz=0.0,
        Ixz=0.0,
    )


def _mp_to_dict(mp: MassProperties) -> Dict[str, float]:
    return {
        "mass": float(mp.mass),
        "x_cg": float(mp.x_cg),
        "y_cg": float(mp.y_cg),
        "z_cg": float(mp.z_cg),
        "Ixx": float(mp.Ixx),
        "Iyy": float(mp.Iyy),
        "Izz": float(mp.Izz),
        "Ixy": float(mp.Ixy),
        "Iyz": float(mp.Iyz),
        "Ixz": float(mp.Ixz),
    }


def default_positions(soln: dict) -> Dict[str, dict]:
    """Computes the default center-of-gravity position of every component from the solution.

    Returns a mapping ``{name: {"xyz": [x, y, z], "kind": str, "draggable": bool, "category": str}}``.
    Point-mass components (avionics, motors, ESCs, navlights, interfaces, power boards) and the
    battery pack are marked draggable so they can be relocated live in the UI; structural surfaces
    (wings, booms, fuselages, stabilizers) are geometry-tied and stay fixed.
    """
    main_wing = soln.get("Main Wing", {})
    geometry = soln.get("Geometry", {})
    hstab = soln.get("HStab", {})
    vstab = soln.get("V Stab", {})

    b_w = float(main_wing.get("wingspan", main_wing.get("b_w", 0)))
    c_root = float(main_wing.get("chordlen", 0))
    boom_len = float(geometry.get("boom_length", 0))
    boom_y = float(geometry.get("boom_y", 0))
    x_cg_assumed = float(geometry.get("cg_le_dist", 0))
    vstab_span = float(vstab.get("vstab_span", 0))
    vstab_root_chord = float(vstab.get("vstab_root_chord", 0))
    hstab_chord = float(hstab.get("hstab_chordlen", 0))

    avionics_x = x_cg_assumed + 0.33 * c_root  # Third chord
    hstab_x = boom_len + vstab_root_chord / 4 + 0.25 * hstab_chord
    hstab_z = vstab_span

    def entry(xyz, kind, draggable, category):
        return {"xyz": [float(xyz[0]), float(xyz[1]), float(xyz[2])],
                "kind": kind, "draggable": bool(draggable), "category": category}

    positions: Dict[str, dict] = {
        # Structural surfaces (geometry-tied, not draggable)
        "Main wing": entry((x_cg_assumed, 0.0, 0.0), "box", False, "structure"),
        "Horizontal stabilizer": entry((hstab_x, 0.0, hstab_z), "box", False, "structure"),
        "Vertical stabilizer L": entry((boom_len + 0.25 * vstab_root_chord, boom_y, 0.5 * vstab_span), "box", False, "structure"),
        "Vertical stabilizer R": entry((boom_len + 0.25 * vstab_root_chord, -boom_y, 0.5 * vstab_span), "box", False, "structure"),
        "Boom L": entry((0.5 * boom_len, boom_y, -0.02), "cylinder", False, "structure"),
        "Boom R": entry((0.5 * boom_len, -boom_y, -0.02), "cylinder", False, "structure"),
        "Fuselage L": entry((-0.25, boom_y, -0.02), "cylinder", False, "structure"),
        "Fuselage R": entry((-0.25, -boom_y, -0.02), "cylinder", False, "structure"),
        "Solar cells": entry((x_cg_assumed + 0.25 * c_root, 0.0, 0.001), "box", False, "structure"),
        # Battery pack (box, but relocatable for CG tuning)
        "Batteries": entry((x_cg_assumed + 0.25 * c_root, 0.0, 0.0), "box", True, "power"),
        # Avionics point masses (third chord, centerline)
        "FC": entry((avionics_x, 0.0, 0.0), "point", True, "avionics"),
        "GPS": entry((avionics_x, 0.0, 0.0), "point", True, "avionics"),
        "Telemetry": entry((avionics_x, 0.0, 0.0), "point", True, "avionics"),
        "Receiver": entry((avionics_x, 0.0, 0.0), "point", True, "avionics"),
        # Navlights
        "Navlight Port": entry((avionics_x, -b_w / 4, 0.0), "point", True, "avionics"),
        "Navlight Starboard": entry((avionics_x, b_w / 4, 0.0), "point", True, "avionics"),
        "Navlight Center": entry((avionics_x, 0.0, 0.0), "point", True, "avionics"),
        "Navlight Hstab": entry((hstab_x + 0.33 * hstab_chord, 0.0, hstab_z), "point", True, "avionics"),
        # Propulsion point masses
        "Motor L": entry((-0.5, boom_y, -0.02), "point", True, "propulsion"),
        "Motor R": entry((-0.5, -boom_y, -0.02), "point", True, "propulsion"),
        "ESC L": entry((-0.25, boom_y, -0.01), "point", True, "propulsion"),
        "ESC R": entry((-0.25, -boom_y, -0.01), "point", True, "propulsion"),
        "Prop L": entry((-0.5, boom_y, 0.0), "disk", False, "propulsion"),
        "Prop R": entry((-0.5, -boom_y, 0.0), "disk", False, "propulsion"),
        # Structural interfaces & power boards (point masses)
        "Superstructure L": entry((x_cg_assumed, boom_y, 0.0), "point", True, "structure"),
        "Superstructure R": entry((x_cg_assumed, -boom_y, 0.0), "point", True, "structure"),
        "Power board L": entry((x_cg_assumed, boom_y, 0.0), "point", True, "power"),
        "Power board R": entry((x_cg_assumed, -boom_y, 0.0), "point", True, "power"),
        "Boom_vstab interface L": entry((boom_len, boom_y, -0.02), "point", True, "structure"),
        "Boom_vstab interface R": entry((boom_len, -boom_y, -0.02), "point", True, "structure"),
        "Vstab_hstab interface L": entry((boom_len + vstab_root_chord / 4, boom_y, vstab_span), "point", True, "structure"),
        "Vstab_hstab interface R": entry((boom_len + vstab_root_chord / 4, -boom_y, vstab_span), "point", True, "structure"),
    }
    return positions


def compute_mass_properties(soln: dict, layout: Dict[str, list] = None) -> dict:
    """Computes full mass properties for the aircraft.

    ``layout`` optionally overrides the center position ``[x, y, z]`` of individual (draggable)
    components; anything not present falls back to :func:`default_positions`.
    """
    layout = layout or {}
    positions = default_positions(soln)
    # Apply user layout overrides (only for known, draggable components)
    for name, xyz in layout.items():
        if name in positions and positions[name]["draggable"] and xyz is not None:
            positions[name]["xyz"] = [float(xyz[0]), float(xyz[1]), float(xyz[2])]

    def pos(name):
        return tuple(positions[name]["xyz"])

    # Access data from solution structure
    main_wing = soln.get("Main Wing", {})
    geometry = soln.get("Geometry", {})
    hstab = soln.get("HStab", {})
    vstab = soln.get("V Stab", {})
    masses = soln.get("Masses", {})
    power = soln.get("Power", {})
    propulsion = soln.get("Propulsion", {})
    airfoils = soln.get("airfoils", {})

    b_w = float(main_wing.get("wingspan", main_wing.get("b_w", 0)))
    c_root = float(main_wing.get("chordlen", 0))
    boom_len = float(geometry.get("boom_length", 0))
    boom_y = float(geometry.get("boom_y", 0))
    x_cg_assumed = float(geometry.get("cg_le_dist", 0))
    vstab_span = float(vstab.get("vstab_span", 0))
    vstab_root_chord = float(vstab.get("vstab_root_chord", 0))
    hstab_chord = float(hstab.get("hstab_chordlen", 0))
    hstab_span = float(hstab.get("hstab_span", 0))

    wing_t_over_c = float(airfoils.get("wing_t_over_c", 0.12))
    tail_t_over_c = float(airfoils.get("tail_t_over_c", 0.10))
    t_wing = wing_t_over_c * c_root
    t_tail = tail_t_over_c * max(hstab_chord, 1e-6)

    boom_radius = float(geometry.get("boom_radius", 0.01))
    solar_panel_side_length = float(power.get("solar_panel_side_length", 0.125))
    
    # Fuselage dimensions from opti.py
    # Fuselage length: 0.5m (from x=-0.5 to x=0)
    # Fuselage radius: 0.4 * dae51 airfoil max thickness
    fuselage_length = 0.5
    fuselage_radius = 0.4 * asb.Airfoil("dae51").max_thickness()
    
    # Battery box dimensions for moment of inertia calculation
    batt_box = [0.20, 0.10, 0.03]  # [Lx, Ly, Lz] in meters

    # Split masses (left/right or per-instance)
    vstab_each_mass = 0.5 * float(masses.get("mass_vstab", 0))
    boom_each_mass = 0.5 * float(masses.get("mass_boom", 0))
    fuselage_each_mass = 0.5 * float(masses.get("mass_fuselages", 0))

    # Solar cells - distributed across full span; mean chord covered sets the box Lx
    solar_panels_n = float(power.get("solar_panels_n", 0))
    solar_area = solar_panels_n * solar_panel_side_length**2
    solar_mean_chord_covered = solar_area / max(b_w, 1e-6)

    # Propulsion
    prop_radius = 0.5 * float(propulsion.get("propeller_diameter", 0.4))

    # Individual avionics component masses
    mass_fc = float(masses.get("mass_fc", 0.08))
    mass_gps = float(masses.get("mass_gps", 0.02))
    mass_telemtry = float(masses.get("mass_telemtry", 0.03))
    mass_receiver = float(masses.get("mass_receiver", 0.02))
    mass_navlight_each = float(masses.get("mass_navlights", 0.04)) / 4.0
    power_board_each_mass = 0.5 * float(masses.get("mass_power_board", 0.15))
    superstructure_each_mass = 0.5 * float(masses.get("mass_superstructures", 0))
    boom_vstab_each_mass = 0.5 * float(masses.get("mass_boom_vstab_interfaces", 0))
    vstab_hstab_each_mass = 0.5 * float(masses.get("mass_vstab_stab_interfaces", 0))

    # Build mass properties components. All center positions come from `positions`
    # (defaults + any user layout overrides applied above), keyed by component name.
    mp_components = {
        "Main wing": _mp_box(float(masses.get("mass_main_wing", 0)), pos("Main wing"), Lx=0.75 * c_root, Ly=b_w, Lz=max(t_wing, 0.005)),
        "Horizontal stabilizer": _mp_box(float(masses.get("mass_hstab", 0)), pos("Horizontal stabilizer"), Lx=hstab_chord, Ly=hstab_span, Lz=max(t_tail, 0.003)),
        "Vertical stabilizer L": _mp_box(
            vstab_each_mass,
            pos("Vertical stabilizer L"),
            Lx=max(0.5 * (vstab_root_chord + hstab_chord), 1e-6),
            Ly=max(t_tail, 0.003),
            Lz=max(vstab_span, 1e-6),
        ),
        "Vertical stabilizer R": _mp_box(
            vstab_each_mass,
            pos("Vertical stabilizer R"),
            Lx=max(0.5 * (vstab_root_chord + hstab_chord), 1e-6),
            Ly=max(t_tail, 0.003),
            Lz=max(vstab_span, 1e-6),
        ),
        "Boom L": _mp_cylinder_x(boom_each_mass, pos("Boom L"), radius_m=boom_radius, length_m=boom_len),
        "Boom R": _mp_cylinder_x(boom_each_mass, pos("Boom R"), radius_m=boom_radius, length_m=boom_len),
        "Fuselage L": _mp_cylinder_x(fuselage_each_mass, pos("Fuselage L"), radius_m=fuselage_radius, length_m=fuselage_length),
        "Fuselage R": _mp_cylinder_x(fuselage_each_mass, pos("Fuselage R"), radius_m=fuselage_radius, length_m=fuselage_length),
        # Solar cells: distributed across full span (Ly = wingspan)
        "Solar cells": _mp_box(float(masses.get("mass_solar_cells", 0)), pos("Solar cells"), Lx=max(solar_mean_chord_covered, 0.05), Ly=b_w, Lz=0.002),
        # Batteries: distributed across inboard wing (Ly = 2 * boom_y)
        "Batteries": _mp_box(float(masses.get("mass_batteries", 0)), pos("Batteries"), Lx=batt_box[0], Ly=2 * boom_y, Lz=batt_box[2]),
        # Wires: Disregard (removed)
        # Avionics components at third chord
        "FC": _mp_point(mass_fc, pos("FC")),
        "GPS": _mp_point(mass_gps, pos("GPS")),
        "Telemetry": _mp_point(mass_telemtry, pos("Telemetry")),
        "Receiver": _mp_point(mass_receiver, pos("Receiver")),
        # Navlights
        "Navlight Port": _mp_point(mass_navlight_each, pos("Navlight Port")),
        "Navlight Starboard": _mp_point(mass_navlight_each, pos("Navlight Starboard")),
        "Navlight Center": _mp_point(mass_navlight_each, pos("Navlight Center")),
        "Navlight Hstab": _mp_point(mass_navlight_each, pos("Navlight Hstab")),
        "Motor L": _mp_point(0.5 * float(masses.get("mass_motors_mounted", 0)), pos("Motor L")),
        "Motor R": _mp_point(0.5 * float(masses.get("mass_motors_mounted", 0)), pos("Motor R")),
        "ESC L": _mp_point(0.5 * float(masses.get("mass_escs", masses.get("mass_esc", 0))), pos("ESC L")),
        "ESC R": _mp_point(0.5 * float(masses.get("mass_escs", masses.get("mass_esc", 0))), pos("ESC R")),
        "Prop L": _mp_disk_x(0.5 * float(masses.get("mass_propellers", 0)), pos("Prop L"), radius_m=max(prop_radius, 0.01)),
        "Prop R": _mp_disk_x(0.5 * float(masses.get("mass_propellers", 0)), pos("Prop R"), radius_m=max(prop_radius, 0.01)),
        # Structural interfaces
        "Superstructure L": _mp_point(superstructure_each_mass, pos("Superstructure L")),
        "Superstructure R": _mp_point(superstructure_each_mass, pos("Superstructure R")),
        "Power board L": _mp_point(power_board_each_mass, pos("Power board L")),
        "Power board R": _mp_point(power_board_each_mass, pos("Power board R")),
        "Boom_vstab interface L": _mp_point(boom_vstab_each_mass, pos("Boom_vstab interface L")),
        "Boom_vstab interface R": _mp_point(boom_vstab_each_mass, pos("Boom_vstab interface R")),
        "Vstab_hstab interface L": _mp_point(vstab_hstab_each_mass, pos("Vstab_hstab interface L")),
        "Vstab_hstab interface R": _mp_point(vstab_hstab_each_mass, pos("Vstab_hstab interface R")),
    }

    mp_total = sum(mp_components.values(), MassProperties(mass=0.0))
    
    # Verify total mass matches solution
    total_mass_calculated = mp_total.mass
    total_mass_expected = float(masses.get("total_mass", 0))
    if abs(total_mass_calculated - total_mass_expected) > 0.01:  # 10g tolerance
        print(f"[WARNING] Total mass mismatch: calculated={total_mass_calculated:.6f} kg, expected={total_mass_expected:.6f} kg")
    
    mass_properties = {
        "components": {k: _mp_to_dict(v) for k, v in mp_components.items()},
        "total": _mp_to_dict(mp_total),
        "total_mass_verification": {
            "calculated": total_mass_calculated,
            "expected": total_mass_expected,
            "difference": total_mass_calculated - total_mass_expected,
        },
        # Per-component centers + editability, so the UI knows what can be dragged.
        "positions": positions,
        "applied_overrides": {
            k: v for k, v in layout.items() if k in positions and positions[k]["draggable"] and v is not None
        },
    }
    return mass_properties


# Backwards-compatible alias for the previous public name.
compute_and_embed_mass_properties = compute_mass_properties


def main(soln_path: Path):
    run_dir = soln_path.parent
    with open(soln_path, "r") as f:
        soln = json.load(f)

    # Apply any persisted point-mass layout (dragged positions) from the run's state.json.
    layout = load_layout(run_dir)
    mass_properties = compute_mass_properties(soln, layout)

    # Consolidated per-run artifact: mass properties live inside state.json alongside `layout`.
    update_state(run_dir, mass=mass_properties)

    print(f"[postprocess] Mass properties written to: {state_path(run_dir)}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_dir = Path(sys.argv[1])
        soln_path = run_dir / "soln.json"
        if not soln_path.exists():
            raise FileNotFoundError(f"Solution file not found: {soln_path}")
    else:
        base_dir = Path(__file__).resolve().parent
        output_dir = base_dir / "output"
        run_dir = latest_run_dir(output_dir)
        soln_path = run_dir / "soln.json"
    
    main(soln_path)
