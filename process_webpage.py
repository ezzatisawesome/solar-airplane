#!/usr/bin/env python3
"""
Generate a standalone Three.js-based aircraft mass viewer HTML file.

Key features:
- Embeds soln.json + mass_properties.json directly into a single HTML file
- Renders basic 3D shapes for components (boxes, cylinders, spheres)
- Optional wireframe from an AeroSandbox .aero airplane file
- Sidebar aircraft specs with tabs
- Legend toggles component visibility
- Hover tooltip with mass + inertia info
- Optional Chart.js battery/power plot if battery_states are present
- Energy balance chart FIXED: embeds energy_balance_data.csv into HTML (no runtime fetch)
"""

from __future__ import annotations

import csv
import html
import io
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

import aerosandbox as asb  # type: ignore

try:
    from lib.artifacts import latest_run_dir  # type: ignore
except (ModuleNotFoundError, ImportError):  # pragma: no cover
    def latest_run_dir(output_dir_base: Path) -> Path:
        """Fallback implementation if optional dependencies in lib/ aren't installed."""
        if not output_dir_base.exists():
            raise FileNotFoundError(f"Output directory not found: {output_dir_base}")
        run_dirs = [p for p in output_dir_base.iterdir() if p.is_dir() and p.name.startswith("run_")]
        if not run_dirs:
            raise FileNotFoundError(f"No run directories found under: {output_dir_base}")
        run_dirs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return run_dirs[0]

JSONDict = Dict[str, Any]


# ----------------------------
# Small HTML helpers
# ----------------------------

def _h(text: str) -> str:
    """HTML-escape."""
    return html.escape(text, quote=True)


def build_top_overlay_html(run_dir: Path) -> str:
    """Builds the top-right download overlay HTML for the viewer.

    Links are generated relative to the viewer HTML file (which lives inside `run_dir`).
    Only links to files that exist are shown.
    """
    candidates: List[Tuple[str, str]] = [
        (".step", "airplane.step"),
        (".aero", "airplane.aero"),
        (".xml", "aircraft.xml"),
        (".xfl", "XFLR/aircraft.xfl"),
    ]

    link_htmls: List[str] = []
    for label, rel_path in candidates:
        if (run_dir / rel_path).exists():
            # `download` is best-effort; browsers may ignore it depending on context.
            link_htmls.append(
                f'<a href="{_h(rel_path)}" target="_blank" rel="noopener" download>{_h(label)}</a>'
            )

    links = "\n        ".join(link_htmls)

    return f"""
    <div id="top-overlay" aria-label="AircraftView header">
      <div class="left"></div>
      <div class="right" aria-label="Downloads">
        <div class="app-title">AircraftView</div>
        {links}
      </div>
    </div>
""".rstrip()


def _num(v: Any, default: float = 0.0) -> float:
    """Best-effort numeric conversion."""
    try:
        if v is None:
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def _row(label: str, value: str) -> str:
    return f"<tr><td><b>{_h(label)}:</b></td><td>{value}</td></tr>"


def _row_multi(columns: List[str]) -> str:
    """Create a table row with multiple columns."""
    cells = "".join(f"<td>{col}</td>" for col in columns)
    return f"<tr>{cells}</tr>"


def _fmt(v: Any, fmt: str, default: float = 0.0, suffix: str = "") -> str:
    """Format numeric value with fallback."""
    return f"{_num(v, default):{fmt}}{suffix}"


def _section(title: str, rows: Iterable[str]) -> str:
    rows_html = "".join(rows)
    return (
        f"<h3>{_h(title)}</h3>"
        "<table style='width: 100%; border-collapse: collapse;'>"
        f"{rows_html}"
        "</table>"
    )


# ----------------------------
# Energy balance CSV -> embedded JSON
# ----------------------------

def load_energy_balance_csv(run_dir: Path) -> Optional[JSONDict]:
    """
    Load energy balance CSV and return series data for plotting.

    Supports header-based CSVs like:
      time_hours, solar_flux_W_per_m2, power_generated_W, power_used_W, net_energy_Wh, battery_state_Wh

    If period_type is missing, infer it from solar flux:
      flux == 0 -> night
      0 < flux < 50 -> dawn_dusk
      flux >= 50 -> day
    """
    candidates = [
        run_dir / "energy_balance_data.csv",
        run_dir / "analysis" / "energy_balance_data.csv",
    ]
    csv_path = next((p for p in candidates if p.exists()), None)
    if csv_path is None:
        return None

    text = csv_path.read_text(encoding="utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text))
    rows = list(reader)
    if not rows:
        return None

    def _get_float(row, *keys, default=None):
        for k in keys:
            if k in row and row[k] not in (None, ""):
                try:
                    return float(row[k])
                except ValueError:
                    pass
        return default

    def _get_str(row, *keys, default=None):
        for k in keys:
            if k in row and row[k] not in (None, ""):
                return str(row[k])
        return default

    time_hr: List[float] = []
    solar_flux: List[float] = []
    period_type: List[str] = []
    p_gen: List[float] = []
    p_used: List[float] = []
    batt_Wh: List[float] = []

    for r in rows:
        t = _get_float(r, "time_hr", "time_hours", "time", default=None)
        if t is None:
            continue

        flux = _get_float(r, "solar_flux_W_per_m2", "solar_flux", default=0.0) or 0.0
        g = _get_float(r, "power_generated_W", "p_gen", "power_generated", default=None)
        u = _get_float(r, "power_used_W", "p_used", "power_used", default=None)
        b = _get_float(r, "battery_state_Wh", "battery_Wh", "battery_state", default=None)

        if g is None or u is None or b is None:
            continue

        period = _get_str(r, "period_type", "period", default=None)
        if period is None:
            if flux <= 0.0:
                period = "night"
            elif flux < 50.0:
                period = "dawn_dusk"
            else:
                period = "day"

        time_hr.append(float(t))
        solar_flux.append(float(flux))
        period_type.append(period)
        p_gen.append(float(g))
        p_used.append(float(u))
        batt_Wh.append(float(b))

    if not time_hr:
        return None

    return {
        "source_csv": str(csv_path.name),
        "time_hr": time_hr,
        "solar_flux_W_per_m2": solar_flux,
        "period_type": period_type,
        "power_generated_W": p_gen,
        "power_used_W": p_used,
        "battery_state_Wh": batt_Wh,
    }


def load_sensitivity_analysis_json(run_dir: Path) -> Optional[JSONDict]:
    """
    Load sensitivity analysis JSON and return data for embedding.
    
    Looks for sensitivity_analysis.json first, then falls back to energy_balance_summary.json.
    """
    candidates = [
        run_dir / "sensitivity_analysis.json",
        run_dir / "analysis" / "sensitivity_analysis.json",
    ]
    json_path = next((p for p in candidates if p.exists()), None)
    
    if json_path is None:
        # Try loading from energy_balance_summary.json
        summary_candidates = [
            run_dir / "energy_balance_summary.json",
            run_dir / "analysis" / "energy_balance_summary.json",
        ]
        summary_path = next((p for p in summary_candidates if p.exists()), None)
        if summary_path:
            try:
                with summary_path.open("r", encoding="utf-8") as f:
                    summary = json.load(f)
                sensitivity_data = summary.get("sensitivity_analysis")
                if sensitivity_data:
                    return sensitivity_data
            except Exception:
                pass
        return None
    
    try:
        with json_path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


# ----------------------------
# XFLR5 loads (distribution TXT) -> embedded JSON
# ----------------------------

_FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _first_float(s: str) -> Optional[float]:
    m = _FLOAT_RE.search(s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def _parse_header_value(lines: List[str], key: str) -> Optional[float]:
    """
    Best-effort parse of scalar values from XFLR text headers.
    Examples:
      QInf  =   12.120000 m/s
      Alpha =    0.000000
      Freestream speed :  12.120 m/s
    """
    key_lower = key.lower()
    for ln in lines[:60]:
        l = ln.strip()
        if not l:
            continue
        if key_lower == "qinf":
            if l.lower().startswith("qinf"):
                return _first_float(l)
            if "freestream speed" in l.lower():
                return _first_float(l)
        elif l.lower().startswith(key_lower):
            return _first_float(l)
    return None


def parse_xflr_distribution_file(txt_path: Path) -> Optional[JSONDict]:
    """
    Parse an XFLR5 distribution export file (per-surface station tables).
    Returns:
      {
        "source_file": "...",
        "alpha_deg": float|None,
        "qinf_mps": float|None,
        "surfaces": {
          "<surface name>": [
             {"y_span_m":..., "chord_m":..., "Cl":..., "CmGeom":..., "CmAirf_c4":..., "BM":...},
             ...
          ],
        },
      }
    """
    if not txt_path.exists():
        return None
    text = txt_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    alpha = _parse_header_value(lines, "Alpha")
    qinf = _parse_header_value(lines, "QInf")

    surfaces: Dict[str, List[JSONDict]] = {}

    i = 0
    n = len(lines)
    while i < n - 2:
        name_line = lines[i].rstrip("\n")
        name = name_line.strip()
        # Surface name lines appear unindented (e.g., "Main Wing"),
        # followed immediately by a header line starting with "y-span".
        if (
            name
            and (name_line[:1] != " ")  # unindented
            and ("cp coefficients" not in name.lower())
        ):
            header = lines[i + 1].strip()
            if ("y-span" in header) and ("Chord" in header) and ("Cl" in header) and ("BM" in header):
                # Parse table rows until blank line
                rows: List[JSONDict] = []
                j = i + 2
                while j < n:
                    row_line = lines[j].strip()
                    if not row_line:
                        break
                    parts = row_line.split()
                    # Expect: y-span, Chord, Ai, Cl, PCd, ICd, CmGeom, CmAirf@c/4, XTrtop, XTrBot, XCP, BM
                    if len(parts) < 12:
                        break
                    try:
                        y_span = float(parts[0])
                        chord_m = float(parts[1])
                        cl = float(parts[3])
                        cm_geom = float(parts[6])
                        cm_airf = float(parts[7])
                        bm = float(parts[11])
                    except ValueError:
                        break
                    rows.append(
                        {
                            "y_span_m": y_span,
                            "chord_m": chord_m,
                            "Cl": cl,
                            "CmGeom": cm_geom,
                            "CmAirf_c4": cm_airf,
                            "BM": bm,
                        }
                    )
                    j += 1
                if rows:
                    surfaces[name] = rows
                i = j + 1
                continue
        i += 1

    if not surfaces:
        return None

    return {
        "source_file": txt_path.name,
        "alpha_deg": alpha,
        "qinf_mps": qinf,
        "surfaces": surfaces,
    }


def _to_half_span(rows: List[JSONDict], y_round_decimals: int = 4) -> List[JSONDict]:
    """
    Convert +/- span data to a single-sided representation by grouping on abs(y),
    and averaging duplicate entries.
    """
    buckets: Dict[float, List[JSONDict]] = {}
    for r in rows:
        y = float(r.get("y_span_m", 0.0))
        key = round(abs(y), y_round_decimals)
        buckets.setdefault(key, []).append(r)

    out: List[JSONDict] = []
    for y_key, rs in buckets.items():
        if not rs:
            continue
        def avg(field: str) -> float:
            vals = [float(x.get(field, 0.0)) for x in rs]
            return sum(vals) / len(vals) if vals else 0.0

        out.append(
            {
                "station_y_m": float(y_key),
                "chord_m": avg("chord_m"),
                "Cl": avg("Cl"),
                "CmGeom": avg("CmGeom"),
                "CmAirf_c4": avg("CmAirf_c4"),
                "BM": avg("BM"),
            }
        )

    out.sort(key=lambda d: float(d.get("station_y_m", 0.0)))
    return out


def _compute_station_loads(
    stations: List[JSONDict],
    rho_kg_m3: Optional[float],
    v_mps: Optional[float],
) -> List[JSONDict]:
    rho = float(rho_kg_m3) if rho_kg_m3 is not None else float("nan")
    v = float(v_mps) if v_mps is not None else float("nan")
    q = 0.5 * rho * v * v if (math.isfinite(rho) and math.isfinite(v)) else float("nan")

    out: List[JSONDict] = []
    for s in stations:
        c = float(s.get("chord_m", 0.0) or 0.0)
        cl = float(s.get("Cl", 0.0) or 0.0)
        cmg = float(s.get("CmGeom", 0.0) or 0.0)
        cma = float(s.get("CmAirf_c4", 0.0) or 0.0)

        lprime = q * c * cl if math.isfinite(q) else float("nan")  # N/m
        mprime_g = q * c * c * cmg if math.isfinite(q) else float("nan")  # N*m per m
        mprime_a = q * c * c * cma if math.isfinite(q) else float("nan")  # N*m per m

        out.append(
            {
                **s,
                "Lprime_N_per_m": lprime,
                "Mprime_CmGeom_Nm_per_m": mprime_g,
                "Mprime_CmAirf_Nm_per_m": mprime_a,
            }
        )
    return out


def load_xflr_loads(run_dir: Path, soln: Mapping[str, Any]) -> Optional[JSONDict]:
    """
    Loads per-station distributions from XFLR export TXT files and converts to physical loads.
    """
    xflr_dir = run_dir / "XFLR"
    if not xflr_dir.exists():
        return None

    # Use soln environment density if present.
    env = soln.get("Environment", {}) or {}
    rho = env.get("density", None)

    # Prefer distribution files with the naming convention: *_a=<alpha>_v=<V>ms.txt
    rx = re.compile(r"_a=[-+]?\d+(?:\.\d+)?_v=\d+(?:\.\d+)?ms\.txt$")
    candidates = [p for p in xflr_dir.glob("*.txt") if rx.search(p.name)]
    candidates.sort(key=lambda p: p.name)
    if not candidates:
        return None

    cases: List[JSONDict] = []
    for p in candidates:
        parsed = parse_xflr_distribution_file(p)
        if not parsed:
            continue

        alpha_deg = parsed.get("alpha_deg", None)
        qinf_mps = parsed.get("qinf_mps", None)
        if qinf_mps is None:
            perf = soln.get("Performance", {}) or {}
            qinf_mps = perf.get("airspeed", None)

        surfaces_out: Dict[str, Any] = {}
        for surf_name, rows in (parsed.get("surfaces", {}) or {}).items():
            half = _to_half_span(rows)
            with_loads = _compute_station_loads(half, rho_kg_m3=rho, v_mps=qinf_mps)
            surfaces_out[surf_name] = {"stations": with_loads}

        # Case label for UI
        label_parts = []
        if alpha_deg is not None:
            label_parts.append(f"α={float(alpha_deg):.2f}°")
        if qinf_mps is not None:
            label_parts.append(f"V={float(qinf_mps):.2f} m/s")
        label = " | ".join(label_parts) if label_parts else p.stem

        cases.append(
            {
                "id": p.stem,
                "label": label,
                "source_file": p.name,
                "alpha_deg": alpha_deg,
                "qinf_mps": qinf_mps,
                "rho_kg_m3": rho,
                "surfaces": surfaces_out,
            }
        )

    if not cases:
        return None

    return {
        "cases": cases,
        "units": {
            "station_y_m": "m",
            "chord_m": "m",
            "Lprime_N_per_m": "N/m",
            "Mprime_Nm_per_m": "N·m/m",
            "BM": "N·m",
        },
    }


# ----------------------------
# Mass Override / Build Helpers
# ----------------------------

def load_mass_updates(run_dir: Path) -> Dict[str, float]:
    """Legacy loader for mass_updates.json (kept for backward compatibility)."""
    update_path = run_dir / "mass_updates.json"
    if not update_path.exists():
        return {}
    
    try:
        with open(update_path, "r", encoding="utf-8") as f:
            updates = json.load(f)
        # Validate and convert to float, filtering out invalid values
        result: Dict[str, float] = {}
        for key, value in updates.items():
            try:
                val = float(value)
                if val >= 0 and math.isfinite(val):
                    result[key] = val
            except (TypeError, ValueError):
                continue
        return result
    except (json.JSONDecodeError, IOError) as e:
        print(f"[warning] Could not load mass updates from {update_path}: {e}")
        return {}


def load_builds(run_dir: Path) -> Tuple[Optional[str], Dict[str, float]]:
    """
    Load mass builds from builds.json (preferred) or fall back to mass_updates.json.
    
    Expected new format (either of these):
      1) [
           { "name": "Mojave 1", "masses": { "mass_solar_cells": 0.62, ... } }
         ]
      2) { "builds": [ { "name": "Mojave 1", "masses": { ... } } ] }
    """
    builds_path = run_dir / "builds.json"
    build_name: Optional[str] = None
    updates: Dict[str, float] = {}

    if builds_path.exists():
        try:
            with open(builds_path, "r", encoding="utf-8") as f:
                data = json.load(f)

            # Support both top-level list and {"builds": [...]} container
            if isinstance(data, dict):
                builds = data.get("builds") or []
            else:
                builds = data or []

            if isinstance(builds, dict):
                builds = [builds]

            if isinstance(builds, list) and builds:
                first = builds[0] or {}
                if isinstance(first, dict):
                    build_name = str(first.get("name") or "Mojave 1")
                    raw_updates = first.get("masses") or first.get("mass_updates") or {}
                    # If masses block missing, allow mass_* keys at top level
                    if not raw_updates:
                        raw_updates = {k: v for k, v in first.items() if k.startswith("mass_")}

                    for key, value in raw_updates.items():
                        try:
                            val = float(value)
                            if val >= 0 and math.isfinite(val):
                                updates[key] = val
                        except (TypeError, ValueError):
                            continue
        except (json.JSONDecodeError, IOError) as e:
            print(f"[warning] Could not load builds from {builds_path}: {e}")

        return build_name, updates

    # Fallback: legacy mass_updates.json behaviour
    legacy_updates = load_mass_updates(run_dir)
    if legacy_updates:
        build_name = "Mojave 1"
        updates = legacy_updates
    return build_name, updates


def calculate_mass_changes(original: Dict[str, float], updates: Dict[str, float]) -> Dict[str, Dict[str, float]]:
    """Calculate deltas and percentages for mass components."""
    changes = {}
    for key, orig_val in original.items():
        update_val = updates.get(key, orig_val)
        delta = update_val - orig_val
        if orig_val != 0:
            pct_change = (delta / orig_val) * 100.0
        else:
            pct_change = 0.0 if delta == 0 else float("inf") if delta > 0 else float("-inf")
        
        changes[key] = {
            "original": orig_val,
            "update": update_val,
            "delta": delta,
            "pct_change": pct_change,
        }
    return changes


# ----------------------------
# Specs HTML
# ----------------------------

def format_specs_html(soln: Mapping[str, Any], mass_properties: Mapping[str, Any], run_dir: Optional[Path] = None) -> str:
    """Format aircraft specifications as HTML."""
    meta = soln.get("Meta", {}) or {}
    perf = soln.get("Performance", {}) or {}
    main_wing = soln.get("Main Wing", {}) or {}
    geom = soln.get("Geometry", {}) or {}
    total = (mass_properties.get("total", {}) or {})
    aero = soln.get("Aerodynamics", {}) or {}
    mission = soln.get("Mission", {}) or {}
    env = soln.get("Environment", {}) or {}
    hstab = soln.get("HStab", {}) or {}
    vstab = soln.get("V Stab", {}) or soln.get("V_Stab", {}) or {}
    power = soln.get("Power", {}) or {}
    prop = soln.get("Propulsion", {}) or {}
    masses = soln.get("Masses", {}) or {}

    tab_contents: Dict[str, List[str]] = {
        "overview": [],
        "performance": [],
        "aerodynamics": [],
        "wings": [],
        "geometry": [],
        "mass": [],
        "power": [],
        "propulsion": [],
        "structure": [],
        "loads": [],
    }

    # Overview tab
    if meta:
        meta_rows = [_row("Run ID", _h(str(meta.get("run_id", "N/A"))))]
        if meta.get("timestamp_utc"):
            meta_rows.append(_row("Timestamp (UTC)", _h(str(meta.get("timestamp_utc")))))
        tab_contents["overview"].append(_section("Meta", meta_rows))

    if mission:
        tab_contents["overview"].append(_section("Mission", [
            _row("Mission Date", _h(str(mission.get("mission_date", "N/A")))),
            _row("Operating Latitude", _fmt(mission.get("operating_lat"), ".4f", suffix="°")),
            _row("Operating Altitude", _fmt(mission.get("operating_altitude"), ".0f", suffix=" m")),
        ]))

    if env:
        tab_contents["overview"].append(_section("Environment", [
            _row("Temperature (K)", _fmt(env.get("temperature_K"), ".2f", suffix=" K")),
            _row("Temperature (F)", _fmt(env.get("temperature_F"), ".1f", suffix=" °F")),
            _row("Pressure", _fmt(env.get("pressure"), ".0f", suffix=" Pa")),
            _row("Density", _fmt(env.get("density"), ".4f", suffix=" kg/m³")),
            _row("Speed of Sound", _fmt(env.get("speed_of_sound"), ".2f", suffix=" m/s")),
        ]))

    # Performance tab
    tab_contents["performance"].append(_section("Performance", [
        _row("Total Mass", _fmt(perf.get("total_mass"), ".3f", suffix=" kg")),
        _row("Airspeed", _fmt(perf.get("airspeed"), ".2f", suffix=" m/s")),
        _row("Thrust (Cruise)", _fmt(perf.get("thrust_cruise"), ".2f", suffix=" N")),
        _row("Thrust (Climb)", _fmt(perf.get("thrust_climb"), ".2f", suffix=" N")),
        _row("Power (Cruise, All Motors)", _fmt(perf.get("power_cruise (all motors)"), ".2f", suffix=" W")),
        _row("Power (Cruise, One Motor)", _fmt(perf.get("power_cruise (one motor)"), ".2f", suffix=" W")),
        _row("Power (Shaft Cruise, One Motor)", _fmt(perf.get("power_shaft_cruise (one motor)"), ".2f", suffix=" W")),
        _row("Power Out Max", _fmt(perf.get("power_out_max"), ".2f", suffix=" W")),
        _row("Min Climb Angle", _fmt(perf.get("min_climb_angle"), ".1f", suffix="°")),
        _row("L/D Ratio", _fmt(perf.get("L_over_D"), ".2f")),
        _row("Stall Speed", _fmt(perf.get("stall_speed"), ".2f", suffix=" m/s")),
        _row("Takeoff Gross Weight", _fmt(perf.get("togw"), ".2f", suffix=" N")),
    ]))

    # Aerodynamics tab
    tab_contents["aerodynamics"].append(_section("Aerodynamics", [
        _row("CL", _fmt(aero.get("CL"), ".3f")),
        _row("CD", _fmt(aero.get("CD"), ".4f")),
        _row("Cm", _fmt(aero.get("Cm"), ".4f")),
        _row("Lift", _fmt(aero.get("L"), ".2f", suffix=" N")),
        _row("Drag", _fmt(aero.get("D"), ".2f", suffix=" N")),
        _row("Static Margin", _fmt(aero.get("static_margin"), ".3f")),
        _row("Neutral Point X", _fmt(aero.get("x_np"), ".3f", suffix=" m")),
    ]))

    # Wings tab
    main_wing_rows = [
        _row("Wingspan", _fmt(main_wing.get("wingspan"), ".3f", suffix=" m")),
        _row("Chord", _fmt(main_wing.get("chordlen"), ".3f", suffix=" m")),
        _row("Area", _fmt(main_wing.get("S_w"), ".3f", suffix=" m²")),
        _row("Aspect Ratio", _fmt(main_wing.get("main_wing_AR"), ".2f")),
    ]
    if main_wing.get("struct_defined_aoa") is not None:
        main_wing_rows.append(_row("Structural AOA", _fmt(main_wing.get("struct_defined_aoa"), ".2f", suffix="°")))
    if main_wing.get("MAC_w") is not None:
        main_wing_rows.append(_row("Mean Aerodynamic Chord", _fmt(main_wing.get("MAC_w"), ".3f", suffix=" m")))
    if main_wing.get("wing_airfoil"):
        main_wing_rows.append(_row("Airfoil", _h(str(main_wing.get("wing_airfoil", "N/A")))))
    tab_contents["wings"].append(_section("Main Wing", main_wing_rows))

    if hstab:
        tab_contents["wings"].append(_section("Horizontal Stabilizer", [
            _row("Span", _fmt(hstab.get("hstab_span"), ".3f", suffix=" m")),
            _row("Chord", _fmt(hstab.get("hstab_chordlen"), ".3f", suffix=" m")),
            _row("Area", _fmt(hstab.get("hstab_area"), ".3f", suffix=" m²")),
            _row("Aspect Ratio", _fmt(hstab.get("hstab_AR"), ".2f")),
            _row("Angle of Attack", _fmt(hstab.get("hstab_aoa"), ".2f", suffix="°")),
            _row("Volume Coefficient", _fmt(hstab.get("V_H_actual"), ".3f")),
            _row("Airfoil", _h(str(hstab.get("hstab_airfoil", "N/A")))),
        ]))

    if vstab:
        tab_contents["wings"].append(_section("Vertical Stabilizer", [
            _row("Span (Height)", _fmt(vstab.get("vstab_span"), ".3f", suffix=" m")),
            _row("Root Chord", _fmt(vstab.get("vstab_root_chord"), ".3f", suffix=" m")),
            _row("Area (each)", _fmt(vstab.get("vstab_area"), ".3f", suffix=" m²")),
            _row("Total Area", _fmt(vstab.get("vstab_area_total"), ".3f", suffix=" m²")),
            _row("Aspect Ratio", _fmt(vstab.get("vstab_AR"), ".2f")),
            _row("Volume Coefficient", _fmt(vstab.get("V_V_actual"), ".4f")),
            _row("Airfoil", _h(str(vstab.get("vstab_airfoil", "N/A")))),
        ]))

    # Geometry tab
    tab_contents["geometry"].append(_section("Geometry", [
        _row("Boom Length", _fmt(geom.get("boom_length"), ".3f", suffix=" m")),
        _row("Boom Y Position", _fmt(geom.get("boom_y"), ".3f", suffix=" m")),
        _row("CG Location (from LE)", _fmt(geom.get("cg_le_dist"), ".3f", suffix=" m")),
    ]))

    # Mass tab
    tab_contents["mass"].append(_section("Mass Properties", [
        _row("Total Mass", _fmt(total.get("mass"), ".3f", suffix=" kg")),
        _row("CG X", _fmt(total.get("x_cg"), ".3f", suffix=" m")),
        _row("CG Y", _fmt(total.get("y_cg"), ".3f", suffix=" m")),
        _row("CG Z", _fmt(total.get("z_cg"), ".3f", suffix=" m")),
        _row("Ixx", _fmt(total.get("Ixx"), ".3f", suffix=" kg·m²")),
        _row("Iyy", _fmt(total.get("Iyy"), ".3f", suffix=" kg·m²")),
        _row("Izz", _fmt(total.get("Izz"), ".3f", suffix=" kg·m²")),
    ]))

    build_name: Optional[str] = None
    if masses:
        # Load mass builds/updates if run_dir is provided
        updates: Dict[str, float] = {}
        changes: Dict[str, Dict[str, float]] = {}
        if run_dir:
            build_name, updates = load_builds(run_dir)
            # Convert masses dict to float dict for calculation
            original_masses = {k: _num(v) for k, v in masses.items() if k.startswith("mass_")}
            changes = calculate_mass_changes(original_masses, updates)
        
        # Mass component mapping: display name -> key in soln.json
        mass_components = [
            ("Solar Cells", "mass_solar_cells"),
            ("Power Board", "mass_power_board"),
            ("Batteries", "mass_batteries"),
            ("Wires", "mass_wires"),
            ("Avionics", "mass_avionics"),
            ("Servos", "mass_servos"),
            ("Motors (Mounted)", "mass_motors_mounted"),
            ("ESCs", "mass_escs"),
            ("Propellers", "mass_propellers"),
            ("Main Wing", "mass_main_wing"),
            ("Horizontal Stabilizer", "mass_hstab"),
            ("Vertical Stabilizer", "mass_vstab"),
            ("Boom", "mass_boom"),
            ("Fuselages", "mass_fuselages"),
            ("Superstructures", "mass_superstructures"),
            ("Boom-VStab Interfaces", "mass_boom_vstab_interfaces"),
            ("VStab-HStab Interfaces", "mass_vstab_stab_interfaces"),
        ]
        
        # Build table with header
        run_name = _h(str(meta.get("run_id", "Original"))) if meta else "Original"
        build_label = _h(build_name or "Mojave 1")
        table_html = "<h3>Masses Breakdown</h3>"
        table_html += "<table style='width: 100%; border-collapse: collapse;'>"
        table_html += "<thead style='background: rgba(255,255,255,0.06);'>"
        table_html += _row_multi([
            "<b>Component</b>",
            f"<b>{run_name} (kg)</b>",
            f"<b>{build_label} (kg)</b>",
            "<b>Delta (kg)</b>",
            "<b>% Change</b>"
        ])
        table_html += "</thead>"
        table_html += "<tbody>"
        
        for display_name, key in mass_components:
            orig_val = _num(masses.get(key, 0))
            if key in changes:
                change = changes[key]
                update_val = change["update"]
                delta = change["delta"]
                pct = change["pct_change"]
            else:
                update_val = orig_val
                delta = 0.0
                pct = 0.0
            
            # Format delta with color coding
            if abs(delta) < 0.001:
                delta_str = _fmt(delta, ".3f", suffix=" kg")
            else:
                color = "rgb(244, 67, 54)" if delta > 0 else "rgb(76, 175, 80)"
                delta_str = f"<span style='color: {color};'>{_fmt(delta, '.3f', suffix=' kg')}</span>"
            
            # Format percentage with color coding
            if abs(pct) < 0.01:
                pct_str = _fmt(pct, ".2f", suffix="%")
            else:
                color = "rgb(244, 67, 54)" if pct > 0 else "rgb(76, 175, 80)"
                sign = "+" if pct > 0 else ""
                pct_str = f"<span style='color: {color};'>{sign}{_fmt(pct, '.2f', suffix='%')}</span>"
            
            table_html += _row_multi([
                f"<b>{_h(display_name)}</b>",
                _fmt(orig_val, ".3f", suffix=" kg"),
                _fmt(update_val, ".3f", suffix=" kg"),
                delta_str,
                pct_str
            ])
        
        # Total mass row
        total_orig = _num(masses.get("total_mass", 0))
        if changes:
            # Recalculate total from component updates
            component_keys = [key for _, key in mass_components]
            total_update = sum(changes.get(k, {}).get("update", _num(masses.get(k, 0))) 
                               for k in component_keys)
            total_delta = total_update - total_orig
            total_pct = (total_delta / total_orig * 100.0) if total_orig != 0 else 0.0
        else:
            total_update = total_orig
            total_delta = 0.0
            total_pct = 0.0
        
        if abs(total_delta) < 0.001:
            total_delta_str = _fmt(total_delta, ".3f", suffix=" kg")
        else:
            color = "rgb(244, 67, 54)" if total_delta > 0 else "rgb(76, 175, 80)"
            total_delta_str = f"<span style='color: {color};'>{_fmt(total_delta, '.3f', suffix=' kg')}</span>"
        
        if abs(total_pct) < 0.01:
            total_pct_str = _fmt(total_pct, ".2f", suffix="%")
        else:
            color = "rgb(244, 67, 54)" if total_pct > 0 else "rgb(76, 175, 80)"
            sign = "+" if total_pct > 0 else ""
            total_pct_str = f"<span style='color: {color};'>{sign}{_fmt(total_pct, '.2f', suffix='%')}</span>"
        
        table_html += _row_multi([
            f"<b>{_h('Total Mass')}</b>",
            _fmt(total_orig, ".3f", suffix=" kg"),
            _fmt(total_update, ".3f", suffix=" kg"),
            total_delta_str,
            total_pct_str
        ])
        
        table_html += "</tbody>"
        table_html += "</table>"
        tab_contents["mass"].append(table_html)
        
        # Add dual pie chart section: baseline run vs build masses
        tab_contents["mass"].append(
            "<div style='margin-top: 20px; padding-top: 15px; border-top: 2px solid #444;'>"
            "<h4 style='margin-bottom: 12px; color: #4CAF50;'>Mass Distribution</h4>"
            "<div style='display:flex; flex-wrap:wrap; gap:16px;'>"
            "  <div style='flex:1 1 260px; min-width:240px;'>"
            "    <div id='massPieLabelOriginal' style='color:#e0e0e0; font-size:13px; margin-bottom:6px;'></div>"
            "    <div style='position: relative; height: 260px; width: 100%; max-height: 260px; overflow: hidden;'>"
            "      <canvas id='massPieChartOriginal' style='max-height: 260px;'></canvas>"
            "    </div>"
            "  </div>"
            "  <div style='flex:1 1 260px; min-width:240px;'>"
            "    <div id='massPieLabelBuild' style='color:#e0e0e0; font-size:13px; margin-bottom:6px;'></div>"
            "    <div style='position: relative; height: 260px; width: 100%; max-height: 260px; overflow: hidden;'>"
            "      <canvas id='massPieChartBuild' style='max-height: 260px;'></canvas>"
            "    </div>"
            "  </div>"
            "</div>"
            "</div>"
        )

    # Power tab
    if power:
        power_rows = [
            _row("Solar Panels", _fmt(power.get("solar_panels_n"), ".1f")),
            _row("Solar Panels (Rows)", _fmt(power.get("solar_panels_n_rows"), ".0f")),
            _row("Battery Capacity", _fmt(power.get("battery_capacity"), ".1f", suffix=" Wh")),
            _row("Battery Voltage", _fmt(power.get("battery_voltage"), ".1f", suffix=" V")),
            _row("Number of Packs", _fmt(power.get("num_packs"), ".1f")),
        ]
        if power.get("solar_panel_side_length") is not None:
            power_rows.append(_row("Panel Size", _fmt(_num(power.get("solar_panel_side_length")) * 1000, ".1f", suffix=" mm")))
        if power.get("solar_encapsulation_eff_hit") is not None:
            power_rows.append(_row("Solar Encapsulation Efficiency (HIT)", _fmt(power.get("solar_encapsulation_eff_hit"), ".3f")))
        if power.get("solar_cell_efficiency") is not None:
            power_rows.append(_row("Solar Cell Efficiency", _fmt(power.get("solar_cell_efficiency"), ".3f")))
        if power.get("energy_generation_margin") is not None:
            power_rows.append(_row("Energy Generation Margin", _fmt(power.get("energy_generation_margin"), ".3f")))
        tab_contents["power"].append(_section("Power System", power_rows))

        battery_states = power.get("battery_states") or []
        if isinstance(battery_states, list) and len(battery_states) > 0:
            tab_contents["power"].append(
                "<div style='margin-top: 15px;'>"
                "<h4 style='margin-bottom: 10px; color: #4CAF50;'>Power Generation & Battery State</h4>"
                "<div style='position: relative; height: 200px; width: 100%; max-height: 200px; overflow: hidden;'>"
                "<canvas id='powerChart' style='max-height: 200px;'></canvas>"
                "</div>"
                "</div>"
            )

    # Sensitivity analysis: single chart with variable + delta dropdowns
    tab_contents["power"].append(
        "<div style='margin-top: 20px; padding-top: 15px; border-top: 2px solid #444;'>"
        "<h4 style='margin-bottom: 10px; color: #4CAF50;'>Sensitivity Analysis (48h)</h4>"
        "<div style='color:#aaa; font-size:12px; margin-bottom:10px;'>"
        "Select which variable to perturb and the magnitude of the change, then the 48-hour energy balance "
        "will be re-plotted using the pre-computed sensitivity data."
        "</div>"
        "<div style='display:flex; flex-wrap:wrap; gap:10px; align-items:center; margin-bottom:10px;'>"
        "  <label style='display:flex; flex-direction:column; gap:4px; font-size:12px;'>"
        "    <span style='color:#e0e0e0; font-weight:600;'>Variable</span>"
        "    <select id='sensVarSelect' style='padding:6px 8px; background:#222; color:#e0e0e0; border:1px solid #444; border-radius:6px; min-width:180px;'></select>"
        "  </label>"
        "  <label style='display:flex; flex-direction:column; gap:4px; font-size:12px;'>"
        "    <span id='sensDeltaLabel' style='color:#e0e0e0; font-weight:600;'>Change</span>"
        "    <select id='sensDeltaSelect' style='padding:6px 8px; background:#222; color:#e0e0e0; border:1px solid #444; border-radius:6px; min-width:140px;'></select>"
        "  </label>"
        "</div>"
        "<div style='position:relative; height:340px; width:100%; margin-top:5px; overflow:hidden;'>"
        "  <canvas id='sensChart_main'></canvas>"
        "</div>"
        "<div style='color:#888; font-size:11px; margin-top:6px;'>"
        "Note: CG shift is modeled as an electrical power penalty proportional to |shift|; this can be tuned in the JS."
        "</div>"
        "</div>"
    )

    # Propulsion tab
    if prop:
        prop_rows: List[str] = []
        if prop.get("propeller_diameter") is not None:
            prop_rows.append(_row("Propeller Diameter", _fmt(prop.get("propeller_diameter"), ".3f", suffix=" m")))
        if prop.get("propeller_n") is not None:
            prop_rows.append(_row("Number of Motors", _fmt(prop.get("propeller_n"), ".0f")))
        if prop.get("advanced_ratio") is not None:
            prop_rows.append(_row("Advanced Ratio", _fmt(prop.get("advanced_ratio"), ".3f")))
        if prop.get("motor_kv") is not None:
            prop_rows.append(_row("Motor KV", _fmt(prop.get("motor_kv"), ".2f", suffix=" rpm/V")))
        if prop.get("rpm_cruise") is not None:
            prop_rows.append(_row("RPM (Cruise)", _fmt(prop.get("rpm_cruise"), ".1f", suffix=" rpm")))
        if prop.get("i_cruise") is not None:
            prop_rows.append(_row("Current (Cruise)", _fmt(prop.get("i_cruise"), ".3f", suffix=" A")))
        if prop_rows:
            tab_contents["propulsion"].append(_section("Propulsion", prop_rows))

    # Structure tab
    if geom.get("boom_radius") is not None:
        tab_contents["structure"].append(_section("Structural Details", [
            _row("Boom Radius", _fmt(_num(geom.get("boom_radius")) * 1000, ".1f", suffix=" mm")),
            _row("Boom Spacing Fraction", _fmt(geom.get("boom_spacing_frac"), ".2f")),
        ]))

    # Loads tab (XFLR distributions + exports)
    tab_contents["loads"].append(
        "<div style='margin-top: 10px;'>"
        "<h3 style='margin-bottom: 8px;'>XFLR Loads</h3>"
        "<div style='display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 12px;'>"
        "  <label style='display:flex; flex-direction:column; gap:4px;'>"
        "    <span style='color:#aaa; font-size:12px;'>Case</span>"
        "    <select id='xflrCaseSelect' style='padding:6px 8px; background:#222; color:#e0e0e0; border:1px solid #444; border-radius:6px; min-width: 220px;'></select>"
        "  </label>"
        "  <label style='display:flex; flex-direction:column; gap:4px;'>"
        "    <span style='color:#aaa; font-size:12px;'>Surface</span>"
        "    <select id='xflrSurfaceSelect' style='padding:6px 8px; background:#222; color:#e0e0e0; border:1px solid #444; border-radius:6px; min-width: 220px;'></select>"
        "  </label>"
        "  <div style='flex:1; min-width: 200px;'></div>"
        "  <button id='xflrCopyCsvBtn' style='padding:8px 10px; background:#2e7d32; color:white; border:none; border-radius:6px; cursor:pointer;'>Copy CSV</button>"
        "  <button id='xflrDownloadCsvBtn' style='padding:8px 10px; background:#1565c0; color:white; border:none; border-radius:6px; cursor:pointer;'>Download CSV</button>"
        "</div>"
        "<div id='xflrLoadsMessage' style='margin: 6px 0 12px 0; color:#aaa; font-size: 13px;'></div>"
        "<div style='display: grid; grid-template-columns: 1fr; gap: 12px;'>"
        "  <div style='position: relative; height: 220px; width: 100%; max-height: 220px; overflow: hidden; background: rgba(0,0,0,0.12); border:1px solid #333; border-radius: 8px; padding: 8px;'>"
        "    <canvas id='xflrLiftChart' style='max-height: 200px;'></canvas>"
        "  </div>"
        "  <div style='position: relative; height: 220px; width: 100%; max-height: 220px; overflow: hidden; background: rgba(0,0,0,0.12); border:1px solid #333; border-radius: 8px; padding: 8px;'>"
        "    <canvas id='xflrMomentChart' style='max-height: 200px;'></canvas>"
        "  </div>"
        "  <div style='position: relative; height: 220px; width: 100%; max-height: 220px; overflow: hidden; background: rgba(0,0,0,0.12); border:1px solid #333; border-radius: 8px; padding: 8px;'>"
        "    <canvas id='xflrBendingChart' style='max-height: 200px;'></canvas>"
        "  </div>"
        "</div>"
        "<div style='margin-top: 14px;'>"
        "  <h4 style='margin: 0 0 8px 0; color:#4CAF50;'>Per-station values</h4>"
        "  <div style='overflow-x:auto; border:1px solid #333; border-radius: 8px;'>"
        "    <table id='xflrLoadsTable' style='width:100%; border-collapse: collapse; font-size: 12px;'>"
        "      <thead style='background: rgba(255,255,255,0.06);'>"
        "        <tr>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>y (m)</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>chord (m)</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>Cl</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>CmGeom</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>CmAirf@c/4</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>L' (N/m)</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>M' CmGeom (N·m/m)</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>M' CmAirf (N·m/m)</th>"
        "          <th style='text-align:left; padding:8px; border-bottom:1px solid #333;'>BM</th>"
        "        </tr>"
        "      </thead>"
        "      <tbody></tbody>"
        "    </table>"
        "  </div>"
        "  <div style='margin-top: 10px; color:#888; font-size: 12px;'>Tip: After clicking <b>Copy CSV</b>, paste directly into Google Sheets (cell A1).</div>"
        "  <textarea id='xflrCsvFallback' readonly style='width:100%; margin-top:10px; height:120px; background:#111; color:#ddd; border:1px solid #333; border-radius:8px; padding:8px; display:none;'></textarea>"
        "</div>"
        "</div>"
    )

    # Build tab HTML structure
    tab_names = {
        "overview": "Overview",
        "performance": "Performance",
        "aerodynamics": "Aerodynamics",
        "wings": "Wings",
        "geometry": "Geometry",
        "mass": "Mass",
        "power": "Power",
        "propulsion": "Propulsion",
        "structure": "Structure",
        "loads": "Loads",
    }

    parts: List[str] = []
    parts.append("<div style='font-family: Arial, sans-serif; padding: 0; margin-bottom: 40px;'>")
    parts.append("<h2 style='margin-top: 0;'>Aircraft Specifications</h2>")

    parts.append("<div class='tabs-container'>")
    parts.append("<div class='tabs'>")
    for tab_id, tab_name in tab_names.items():
        active_class = " active" if tab_id == "overview" else ""
        parts.append(f"<button class='tab-button{active_class}' data-tab='{tab_id}'>{_h(tab_name)}</button>")
    parts.append("</div>")

    for tab_id in tab_names:
        active_class = " active" if tab_id == "overview" else ""
        parts.append(f"<div class='tab-content{active_class}' id='tab-{tab_id}'>")
        parts.extend(tab_contents[tab_id])
        parts.append("</div>")

    parts.append("</div>")  # tabs-container
    parts.append("</div>")  # wrapper
    return "".join(parts)


# ----------------------------
# AeroSandbox geometry extraction
# ----------------------------

def extract_airplane_geometry(airplane_path: Path) -> JSONDict:
    """Extract geometry from an AeroSandbox .aero airplane file for wireframe visualization."""
    if not airplane_path.exists():
        return {}

    if asb is None:
        print("[warning] AeroSandbox not installed; skipping airplane wireframe extraction.")
        return {}

    try:
        airplane = asb.load(str(airplane_path))
    except Exception as e:
        print(f"[warning] Could not load airplane file {airplane_path}: {e}")
        return {}

    geometry: JSONDict = {"wings": [], "fuselages": []}

    for wing in getattr(airplane, "wings", []) or []:
        wing_data: JSONDict = {
            "name": getattr(wing, "name", "Unknown"),
            "symmetric": bool(getattr(wing, "symmetric", False)),
            "xsecs": [],
        }
        for xsec in getattr(wing, "xsecs", []) or []:
            xyz_le = getattr(xsec, "xyz_le", [0, 0, 0])
            if hasattr(xyz_le, "tolist"):
                xyz_le = xyz_le.tolist()
            wing_data["xsecs"].append({
                "xyz_le": list(xyz_le),
                "chord": float(getattr(xsec, "chord", 0.0) or 0.0),
            })
        geometry["wings"].append(wing_data)

    for fuselage in getattr(airplane, "fuselages", []) or []:
        fuselage_data: JSONDict = {"name": getattr(fuselage, "name", "Unknown"), "xsecs": []}
        for xsec in getattr(fuselage, "xsecs", []) or []:
            xyz_c = getattr(xsec, "xyz_c", [0, 0, 0])
            if hasattr(xyz_c, "tolist"):
                xyz_c = xyz_c.tolist()
            fuselage_data["xsecs"].append({
                "xyz_c": list(xyz_c),
                "radius": float(getattr(xsec, "radius", 0.0) or 0.0),
            })
        geometry["fuselages"].append(fuselage_data)

    xyz_ref = getattr(airplane, "xyz_ref", [0, 0, 0])
    if hasattr(xyz_ref, "tolist"):
        xyz_ref = xyz_ref.tolist()
    geometry["xyz_ref"] = list(xyz_ref)
    return geometry


# ----------------------------
# HTML template
# ----------------------------

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>AircraftView</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
  <style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: Arial, sans-serif; overflow-x: hidden; overflow-y: auto; background: #1a1a1a; color: #e0e0e0; }}

    #canvas-container {{
      width: 100vw;
      height: 60vh;
      position: relative;
      background: #1a1a1a;
    }}

    #top-overlay {{
      position: absolute;
      top: 10px;
      left: 10px;
      right: 10px;
      z-index: 1000;
      user-select: none;
      display: flex;
      align-items: center;
      justify-content: space-between;
      color: #cfcfcf;
      font-size: 12px;
      font-weight: 400;
      letter-spacing: 0.1px;
      pointer-events: none;
    }}
    #top-overlay .app-title {{
      font-size: 14px;
      padding: 4px 8px;
      border: 1px solid #ffffff;
      border-radius: 0px;
      color: #ffffff;
    }}
    #top-overlay .left,
    #top-overlay .right {{
      display: flex;
      align-items: center;
      gap: 10px;
      pointer-events: auto;
    }}
    #top-overlay a {{
      color: #cfcfcf;
      text-decoration: underline;
      text-underline-offset: 2px;
      opacity: 0.9;
    }}
    #top-overlay a:hover {{
      opacity: 1.0;
    }}
    #top-overlay .sep {{
      opacity: 0.5;
    }}

    #content-area {{
      padding: 0;
      width: 100%;
      margin: 0;
    }}

    .tabs-container {{
      margin-top: 20px;
    }}

    .tabs {{
      display: flex;
      flex-wrap: wrap;
      gap: 5px;
      margin-bottom: 20px;
      border-bottom: 2px solid #444;
    }}

    .tab-button {{
      background: rgba(40,40,40,0.8);
      color: #b0b0b0;
      border: none;
      padding: 10px 20px;
      cursor: pointer;
      font-size: 14px;
      border-radius: 5px 5px 0 0;
      transition: all 0.3s ease;
      border-bottom: 2px solid transparent;
      margin-bottom: -2px;
    }}

    .tab-button:hover {{
      background: rgba(50,50,50,0.9);
      color: #e0e0e0;
    }}

    .tab-button.active {{
      background: rgba(30,30,30,0.95);
      color: #4CAF50;
      border-bottom: 2px solid #4CAF50;
    }}

    .tab-content {{
      display: none;
    }}

    .tab-content.active {{
      display: block;
    }}

    #legend {{
      position: absolute;
      top: 10px;
      left: 10px;
      background: rgba(30,30,30,0.95);
      padding: 15px;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.5);
      min-width: 200px;
      max-width: 250px;
      max-height: calc(60vh - 20px);
      overflow-y: auto;
      z-index: 1000;
      border: 1px solid #444;
    }}
    #legend h3 {{ margin-top: 0; margin-bottom: 10px; color: #4CAF50; font-size: 16px; }}
    .legend-item {{ display: flex; align-items: center; padding: 5px 0; cursor: pointer; user-select: none; color: #e0e0e0; }}
    .legend-item:hover {{ background: rgba(255,255,255,0.1); border-radius: 4px; }}
    .legend-color {{ width: 20px; height: 20px; border: 1px solid #555; margin-right: 8px; flex-shrink: 0; border-radius: 3px; }}
    .legend-item.hidden {{ opacity: 0.3; }}

    #sidebar {{
      flex: 1;
      background: rgba(30,30,30,0.95);
      padding: 20px;
      border-radius: 0px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.5);
      min-width: 400px;
    }}
    #sidebar h1 {{ margin-top: 0; color: #4CAF50; font-size: 24px; }}
    #sidebar h2 {{ margin-top: 0; color: #4CAF50; font-size: 20px; }}
    #sidebar h3 {{ margin-top: 15px; margin-bottom: 10px; color: #4CAF50; font-size: 16px; }}
    #sidebar h4 {{ margin-top: 10px; margin-bottom: 8px; color: #4CAF50; font-size: 14px; }}
    #sidebar p {{ color: #aaa; }}
    #sidebar table {{ width: 100%; border-collapse: collapse; margin: 10px 0; font-size: 13px; }}
    #sidebar table td {{ padding: 8px; border-bottom: 1px solid #444; color: #e0e0e0; }}
    #sidebar table td:first-child {{ font-weight: bold; width: 60%; color: #b0b0b0; }}

    #info {{
      position: fixed;
      bottom: 10px;
      left: 10px;
      background: rgba(0,0,0,0.8);
      color: #e0e0e0;
      padding: 10px 15px;
      border-radius: 5px;
      font-size: 12px;
      z-index: 1000;
      border: 1px solid #444;
    }}

    #tooltip {{
      position: fixed;
      background: rgba(20,20,20,0.95);
      color: #e0e0e0;
      padding: 10px 15px;
      border-radius: 5px;
      font-size: 12px;
      pointer-events: none;
      z-index: 2000;
      display: none;
      max-width: 320px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.5);
      border: 1px solid #444;
    }}
    #tooltip.visible {{ display: block; }}
    #tooltip h4 {{ margin: 0 0 5px 0; color: #4CAF50; font-size: 14px; }}
    #tooltip p {{ margin: 3px 0; color: #e0e0e0; }}
  </style>
</head>
<body>
  <div id="canvas-container">
{top_overlay_html}
    <div id="legend">
      <h3>Components</h3>
      <div id="legend-content"></div>
    </div>
  </div>

  <div id="content-area">
    <div id="sidebar">
      {specs_html}
    </div>
  </div>

  <div id="info">
    Left click: Rotate | Right click: Pan | Scroll: Zoom
  </div>

  <div id="tooltip"></div>

  <script type="importmap">
    {{
      "imports": {{
        "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
        "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"
      }}
    }}
  </script>

  <script type="module">
    import * as THREE from 'three';
    import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';

    const soln = {soln_json};
    const massProperties = {mass_properties_json};
    const airplaneGeometry = {airplane_geometry_json};
    const energyBalanceData = {energy_balance_data_json};
    const sensitivityAnalysisData = {sensitivity_analysis_json};
    const xflrLoadsData = {xflr_loads_data_json};
    const massChartData = {mass_chart_data_json};

    const batteryStates = soln.Power?.battery_states || null;
    const batteryCapacity = soln.Power?.battery_capacity || 0;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x1a1a1a);

    // Align to expected aircraft view: rotate scene instead of remapping every component
    scene.rotation.z = Math.PI;
    scene.rotation.x = -Math.PI / 2;

    const canvasContainer = document.getElementById('canvas-container');
    const canvasHeight = Math.floor(window.innerHeight * 0.6);
    const canvasWidth = window.innerWidth;

    const camera = new THREE.PerspectiveCamera(75, canvasWidth / canvasHeight, 0.1, 1000);
    camera.position.set(5, 5, 5);
    camera.lookAt(0, 0, 0);

    const renderer = new THREE.WebGLRenderer({{ antialias: true }});
    renderer.setSize(canvasWidth, canvasHeight);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.shadowMap.enabled = true;
    canvasContainer.appendChild(renderer.domElement);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;

    scene.add(new THREE.AmbientLight(0xffffff, 0.6));
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight.position.set(10, 10, 5);
    directionalLight.castShadow = true;
    scene.add(directionalLight);

    function getColorForCategory(name) {{
      const n = (name || '').toLowerCase();
      if (n.includes('wing') || n.includes('stabilizer')) return 0x4A90E2;
      if (n.includes('boom') || n.includes('fuselage') || n.includes('superstructure') || n.includes('interface')) return 0x8B7355;
      if (n.includes('motor') || n.includes('esc') || n.includes('prop')) return 0xE67E22;
      if (n.includes('batter')) return 0xE74C3C;
      if (n.includes('solar')) return 0xF1C40F;
      if (n.includes('fc') || n.includes('gps') || n.includes('telemetry') || n.includes('receiver') || n.includes('navlight') || n.includes('power board')) return 0x27AE60;
      return 0x9B59B6;
    }}

    // Lightweight materials + geometry (reduced segments)
    function createMaterial(color, opacity=0.7) {{
      return new THREE.MeshPhongMaterial({{
        color,
        opacity,
        transparent: opacity < 1.0,
        side: THREE.DoubleSide,
        emissive: 0x000000,
        emissiveIntensity: 0.0
      }});
    }}

    function createBox(center, Lx, Ly, Lz, color, userData) {{
      const geometry = new THREE.BoxGeometry(Lx, Ly, Lz); // no subdivisions
      const material = createMaterial(color, 0.7);
      const mesh = new THREE.Mesh(geometry, material);
      mesh.position.set(center[0], center[1], center[2]);
      mesh.castShadow = true;
      mesh.receiveShadow = true;
      mesh.userData = userData;
      return mesh;
    }}

    function createCylinder(center, radius, length, axis, color, userData) {{
      const geometry = new THREE.CylinderGeometry(radius, radius, length, 24);
      const material = createMaterial(color, 0.7);
      const mesh = new THREE.Mesh(geometry, material);

      if (axis === 'x') mesh.rotation.z = Math.PI / 2;
      else if (axis === 'z') mesh.rotation.x = Math.PI / 2;

      mesh.position.set(center[0], center[1], center[2]);
      mesh.castShadow = true;
      mesh.receiveShadow = true;
      mesh.userData = userData;
      return mesh;
    }}

    function createSphere(center, radius, color, userData) {{
      const geometry = new THREE.SphereGeometry(radius, 16, 12);
      const material = createMaterial(color, 0.8);
      const mesh = new THREE.Mesh(geometry, material);
      mesh.position.set(center[0], center[1], center[2]);
      mesh.castShadow = true;
      mesh.receiveShadow = true;
      mesh.userData = userData;
      return mesh;
    }}

    function reconstructGeometry(soln, massProperties) {{
      const mainWing = soln['Main Wing'] || soln.Main_Wing || {{}};
      const geometry = soln.Geometry || {{}};
      const hstab = soln.HStab || {{}};
      const vstab = soln['V Stab'] || soln.V_Stab || soln.VStab || {{}};
      const power = soln.Power || {{}};
      const propulsion = soln.Propulsion || {{}};
      const airfoils = soln.airfoils || {{}};

      const b_w = mainWing.wingspan || mainWing.b_w || 0;
      const c_root = mainWing.chordlen || 0;
      const boom_len = geometry.boom_length || 0;
      const boom_y = geometry.boom_y || 0;
      const vstab_span = vstab.vstab_span || 0;
      const vstab_root_chord = vstab.vstab_root_chord || 0;
      const hstab_chord = hstab.hstab_chordlen || 0;
      const hstab_span = hstab.hstab_span || 0;

      const wing_t_over_c = airfoils.wing_t_over_c || 0.12;
      const tail_t_over_c = airfoils.tail_t_over_c || 0.10;
      const t_wing = wing_t_over_c * c_root;
      const t_tail = tail_t_over_c * Math.max(hstab_chord, 1e-6);

      const boom_radius = geometry.boom_radius || 0.01;
      const solar_panel_side_length = power.solar_panel_side_length || 0.125;

      const fuselage_length = 0.5;
      const fuselage_radius = 0.015;
      const batt_box = [0.20, 0.10, 0.03];

      const solar_panels_n = power.solar_panels_n || 0;
      const solar_area = solar_panels_n * solar_panel_side_length * solar_panel_side_length;
      const solar_mean_chord_covered = solar_area / Math.max(b_w, 1e-6);

      const prop_radius = 0.5 * (propulsion.propeller_diameter || 0.4);

      const components = massProperties.components || {{}};
      const geometryMap = {{}};

      const masses = Object.values(components).map(c => c.mass || 0);
      const maxMass = Math.max(1e-9, ...(masses.length ? masses : [1e-9]));

      for (const [name, comp] of Object.entries(components)) {{
        const x_cg = comp.x_cg || 0;
        const y_cg = comp.y_cg || 0;
        const z_cg = comp.z_cg || 0;
        const mass = comp.mass || 0;

        if (name === 'Main wing') {{
          geometryMap[name] = {{ type: 'box', center: [x_cg, y_cg, z_cg], Lx: 0.75 * c_root, Ly: b_w, Lz: Math.max(t_wing, 0.005), mass }};
        }} else if (name === 'Horizontal stabilizer') {{
          geometryMap[name] = {{ type: 'box', center: [x_cg, y_cg, z_cg], Lx: hstab_chord, Ly: hstab_span, Lz: Math.max(t_tail, 0.003), mass }};
        }} else if (name.includes('Vertical stabilizer')) {{
          geometryMap[name] = {{ type: 'box', center: [x_cg, y_cg, z_cg], Lx: Math.max(0.5 * (vstab_root_chord + hstab_chord), 1e-6), Ly: Math.max(t_tail, 0.003), Lz: Math.max(vstab_span, 1e-6), mass }};
        }} else if (name === 'Boom L' || name === 'Boom R') {{
          geometryMap[name] = {{ type: 'cylinder', center: [x_cg, y_cg, z_cg], radius: boom_radius, length: boom_len, axis: 'x', mass }};
        }} else if (name === 'Fuselage L' || name === 'Fuselage R') {{
          geometryMap[name] = {{ type: 'cylinder', center: [x_cg, y_cg, z_cg], radius: fuselage_radius, length: fuselage_length, axis: 'x', mass }};
        }} else if (name === 'Solar cells') {{
          geometryMap[name] = {{ type: 'box', center: [x_cg, y_cg, z_cg], Lx: Math.max(solar_mean_chord_covered, 0.05), Ly: b_w, Lz: 0.002, mass }};
        }} else if (name === 'Batteries') {{
          geometryMap[name] = {{ type: 'box', center: [x_cg, y_cg, z_cg], Lx: batt_box[0], Ly: 2 * boom_y, Lz: batt_box[2], mass }};
        }} else if (name === 'Prop L' || name === 'Prop R') {{
          geometryMap[name] = {{ type: 'cylinder', center: [x_cg, y_cg, z_cg], radius: Math.max(prop_radius, 0.01), length: 0.001, axis: 'x', mass }};
        }} else {{
          const sphereRadius = 0.05 * (mass / maxMass) + 0.02;
          geometryMap[name] = {{ type: 'point', center: [x_cg, y_cg, z_cg], radius: sphereRadius, mass }};
        }}
      }}

      return geometryMap;
    }}

    const geometryMap = reconstructGeometry(soln, massProperties);
    const components = massProperties.components || {{}};

    const componentMeshes = {{}};
    const legendContent = document.getElementById('legend-content');

    for (const [name, geom] of Object.entries(geometryMap)) {{
      const color = getColorForCategory(name);
      const componentData = components[name] || {{}};

      const inertia = (componentData.Ixx !== undefined) ? {{
        Ixx: componentData.Ixx, Iyy: componentData.Iyy, Izz: componentData.Izz,
        Ixy: componentData.Ixy || 0, Iyz: componentData.Iyz || 0, Ixz: componentData.Ixz || 0
      }} : null;

      const userData = {{
        name, type: geom.type, mass: geom.mass, center: geom.center, inertia
      }};

      let mesh = null;
      if (geom.type === 'box') mesh = createBox(geom.center, geom.Lx, geom.Ly, geom.Lz, color, userData);
      else if (geom.type === 'cylinder') mesh = createCylinder(geom.center, geom.radius, geom.length, geom.axis, color, userData);
      else if (geom.type === 'point') mesh = createSphere(geom.center, geom.radius, color, userData);

      if (!mesh) continue;

      scene.add(mesh);
      componentMeshes[name] = mesh;

      const legendItem = document.createElement('div');
      legendItem.className = 'legend-item';
      legendItem.innerHTML = `
        <div class="legend-color" style="background-color: #${{color.toString(16).padStart(6, '0')}}"></div>
        <span>${{name}}</span>
      `;
      legendItem.onclick = () => {{
        mesh.visible = !mesh.visible;
        legendItem.classList.toggle('hidden', !mesh.visible);
      }};
      legendContent.appendChild(legendItem);
    }}

    // Total CG marker
    const total = massProperties.total || {{}};
    const cgGeometry = new THREE.SphereGeometry(0.05, 16, 12);
    const cgMaterial = new THREE.MeshPhongMaterial({{ color: 0xff0000, emissive: 0x000000, emissiveIntensity: 0.0 }});
    const cgMesh = new THREE.Mesh(cgGeometry, cgMaterial);
    cgMesh.position.set(total.x_cg || 0, total.y_cg || 0, total.z_cg || 0);
    cgMesh.userData = {{
      name: 'Total CG', type: 'cg', mass: total.mass || 0,
      center: [total.x_cg || 0, total.y_cg || 0, total.z_cg || 0],
      inertia: null
    }};
    scene.add(cgMesh);
    componentMeshes['Total CG'] = cgMesh;

    // Wireframe tube helper
    function createWireframeTube(start, end, radius=0.003, color=0x808080) {{
      const direction = new THREE.Vector3().subVectors(end, start);
      const length = direction.length();
      const geometry = new THREE.CylinderGeometry(radius, radius, length, 8);
      const material = new THREE.MeshBasicMaterial({{ color, depthWrite: false }});
      const mesh = new THREE.Mesh(geometry, material);

      const midpoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
      mesh.position.copy(midpoint);

      const axis = new THREE.Vector3(0, 1, 0);
      const quaternion = new THREE.Quaternion();
      quaternion.setFromUnitVectors(axis, direction.normalize());
      mesh.quaternion.copy(quaternion);

      mesh.renderOrder = 1000;
      return mesh;
    }}

    function addWireframe(airplaneGeometry) {{
      if (!airplaneGeometry?.wings?.length) return;

      const wireframeColor = 0x808080;
      const wireframeRadius = 0.003;
      const xyz_ref = airplaneGeometry.xyz_ref || [0,0,0];
      const cgX = (massProperties.total || {{}}).x_cg || 0;

      for (const wing of airplaneGeometry.wings) {{
        const xsecs = wing.xsecs || [];
        if (xsecs.length < 2) continue;

        function le(xsec) {{
          return new THREE.Vector3(xsec.xyz_le[0] + xyz_ref[0] - cgX, xsec.xyz_le[1] + xyz_ref[1], xsec.xyz_le[2] + xyz_ref[2]);
        }}
        function te(xsec) {{
          return new THREE.Vector3(xsec.xyz_le[0] + xsec.chord + xyz_ref[0] - cgX, xsec.xyz_le[1] + xyz_ref[1], xsec.xyz_le[2] + xyz_ref[2]);
        }}

        for (let i = 0; i < xsecs.length - 1; i++) {{
          const a = xsecs[i], b = xsecs[i+1];
          scene.add(createWireframeTube(le(a), le(b), wireframeRadius, wireframeColor));
          scene.add(createWireframeTube(te(a), te(b), wireframeRadius, wireframeColor));
          scene.add(createWireframeTube(le(a), te(a), wireframeRadius, wireframeColor));
        }}
        scene.add(createWireframeTube(le(xsecs[xsecs.length-1]), te(xsecs[xsecs.length-1]), wireframeRadius, wireframeColor));

        if (wing.symmetric) {{
          for (let i = 0; i < xsecs.length - 1; i++) {{
            const a = xsecs[i], b = xsecs[i+1];
            const leA = le(a); leA.y = -leA.y + 2*xyz_ref[1];
            const leB = le(b); leB.y = -leB.y + 2*xyz_ref[1];
            const teA = te(a); teA.y = -teA.y + 2*xyz_ref[1];
            const teB = te(b); teB.y = -teB.y + 2*xyz_ref[1];

            scene.add(createWireframeTube(leA, leB, wireframeRadius, wireframeColor));
            scene.add(createWireframeTube(teA, teB, wireframeRadius, wireframeColor));
            scene.add(createWireframeTube(leA, teA, wireframeRadius, wireframeColor));
          }}
          const last = xsecs[xsecs.length-1];
          const leL = le(last); leL.y = -leL.y + 2*xyz_ref[1];
          const teL = te(last); teL.y = -teL.y + 2*xyz_ref[1];
          scene.add(createWireframeTube(leL, teL, wireframeRadius, wireframeColor));
        }}
      }}

      for (const fuselage of airplaneGeometry.fuselages || []) {{
        const xsecs = fuselage.xsecs || [];
        if (xsecs.length < 2) continue;
        for (let i = 0; i < xsecs.length - 1; i++) {{
          const a = xsecs[i], b = xsecs[i+1];
          const c1 = new THREE.Vector3(a.xyz_c[0] + xyz_ref[0] + cgX, a.xyz_c[1] + xyz_ref[1], a.xyz_c[2] + xyz_ref[2]);
          const c2 = new THREE.Vector3(b.xyz_c[0] + xyz_ref[0] - cgX, b.xyz_c[1] + xyz_ref[1], b.xyz_c[2] + xyz_ref[2]);
          scene.add(createWireframeTube(c1, c2, wireframeRadius, wireframeColor));
        }}
      }}
    }}

    addWireframe(airplaneGeometry);

    scene.add(new THREE.AxesHelper(2));
    const gridHelper = new THREE.GridHelper(10, 10);
    gridHelper.rotation.x = -Math.PI / 2;
    scene.add(gridHelper);

    function createAxisLabel(text, position, color) {{
      const canvas = document.createElement('canvas');
      const ctx = canvas.getContext('2d');
      canvas.width = 64; canvas.height = 64;
      ctx.fillStyle = color;
      ctx.font = 'Bold 48px Arial';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(text, 32, 32);

      const texture = new THREE.CanvasTexture(canvas);
      const sprite = new THREE.Sprite(new THREE.SpriteMaterial({{ map: texture }}));
      sprite.position.set(position[0], position[1], position[2]);
      sprite.scale.set(0.5, 0.5, 1);
      return sprite;
    }}

    scene.add(createAxisLabel('X', [2.5, 0, 0], '#ff0000'));
    scene.add(createAxisLabel('Y', [0, 2.5, 0], '#00ff00'));
    scene.add(createAxisLabel('Z', [0, 0, 2.5], '#0000ff'));

    // Hover tooltip + highlight (no material cloning; adjust emissive)
    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2();
    const tooltip = document.getElementById('tooltip');
    const allMeshes = Object.values(componentMeshes);

    let hovered = null;

    function setHover(mesh, enabled) {{
      if (!mesh || !mesh.material) return;
      if (!mesh.material.emissive) return;
      mesh.material.emissive.setHex(enabled ? 0xffffff : 0x000000);
      mesh.material.emissiveIntensity = enabled ? 0.35 : 0.0;
      mesh.material.needsUpdate = true;
    }}

    function onMouseMove(event) {{
      const canvasRect = canvasContainer.getBoundingClientRect();
      mouse.x = ((event.clientX - canvasRect.left) / canvasRect.width) * 2 - 1;
      mouse.y = -((event.clientY - canvasRect.top) / canvasRect.height) * 2 + 1;

      raycaster.setFromCamera(mouse, camera);
      const hits = raycaster.intersectObjects(allMeshes, true);

      if (hits.length) {{
        const obj = hits[0].object;
        if (hovered !== obj) {{
          if (hovered) setHover(hovered, false);
          hovered = obj;
          setHover(obj, true);

          const ud = obj.userData || {{}};
          tooltip.style.left = (event.clientX + 10) + 'px';
          tooltip.style.top = (event.clientY + 10) + 'px';

          let content = `<h4>${{ud.name || 'Component'}}</h4>`;
          content += `<p><b>Type:</b> ${{ud.type}}</p>`;
          content += `<p><b>Mass:</b> ${{(ud.mass ?? 0).toFixed(4)}} kg</p>`;
          if (ud.center) {{
            content += `<p><b>Position:</b> (${{ud.center[0].toFixed(3)}}, ${{ud.center[1].toFixed(3)}}, ${{ud.center[2].toFixed(3)}}) m</p>`;
          }}
          if (ud.inertia) {{
            content += `<p><b>Inertia (kg·m²):</b></p>`;
            content += `<p style="margin-left: 10px;">Ixx: ${{ud.inertia.Ixx.toFixed(6)}}</p>`;
            content += `<p style="margin-left: 10px;">Iyy: ${{ud.inertia.Iyy.toFixed(6)}}</p>`;
            content += `<p style="margin-left: 10px;">Izz: ${{ud.inertia.Izz.toFixed(6)}}</p>`;
            if (ud.inertia.Ixy || ud.inertia.Iyz || ud.inertia.Ixz) {{
              content += `<p style="margin-left: 10px;">Ixy: ${{(ud.inertia.Ixy || 0).toFixed(6)}}</p>`;
              content += `<p style="margin-left: 10px;">Iyz: ${{(ud.inertia.Iyz || 0).toFixed(6)}}</p>`;
              content += `<p style="margin-left: 10px;">Ixz: ${{(ud.inertia.Ixz || 0).toFixed(6)}}</p>`;
            }}
          }}

          tooltip.innerHTML = content;
          tooltip.classList.add('visible');
        }} else {{
          tooltip.style.left = (event.clientX + 10) + 'px';
          tooltip.style.top = (event.clientY + 10) + 'px';
        }}
      }} else {{
        if (hovered) {{
          setHover(hovered, false);
          hovered = null;
        }}
        tooltip.classList.remove('visible');
      }}
    }}

    renderer.domElement.addEventListener('mousemove', onMouseMove);
    renderer.domElement.addEventListener('mouseleave', () => {{
      if (hovered) setHover(hovered, false);
      hovered = null;
      tooltip.classList.remove('visible');
    }});

    function animate() {{
      requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    }}

    window.addEventListener('resize', () => {{
      const h = Math.floor(window.innerHeight * 0.6);
      const w = window.innerWidth;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
      renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    }});

    animate();

    // Battery/power chart
    if (batteryStates && batteryStates.length > 0) {{
      const canvas = document.getElementById('powerChart');
      if (canvas) {{
        canvas.style.height = '200px';
        canvas.style.maxHeight = '200px';

        const timeHours = Array.from({{length: batteryStates.length}}, (_, i) => (i / batteryStates.length) * 24);
        const dt = 24 / batteryStates.length;

        const powerGeneration = [];
        for (let i = 0; i < batteryStates.length - 1; i++) {{
          const dE = batteryStates[i + 1] - batteryStates[i];
          powerGeneration.push(dE / dt);
        }}
        if (powerGeneration.length) powerGeneration.push(powerGeneration[powerGeneration.length - 1]);

        new Chart(canvas, {{
          type: 'line',
          data: {{
            labels: timeHours.map(t => t.toFixed(1)),
            datasets: [
              {{
                label: 'Battery State (Wh)',
                data: batteryStates,
                borderColor: 'rgb(76, 175, 80)',
                backgroundColor: 'rgba(76, 175, 80, 0.2)',
                yAxisID: 'y',
                tension: 0.1,
                pointRadius: 0
              }},
              {{
                label: 'Power Generation (W)',
                data: powerGeneration,
                borderColor: 'rgb(255, 152, 0)',
                backgroundColor: 'rgba(255, 152, 0, 0.2)',
                yAxisID: 'y1',
                tension: 0.1,
                pointRadius: 0
              }}
            ]
          }},
          options: {{
            responsive: true,
            maintainAspectRatio: false,
            interaction: {{ mode: 'index', intersect: false }},
            plugins: {{
              title: {{ display: true, text: '24-Hour Power Generation & Battery State', color: '#4CAF50' }},
              legend: {{ display: true, position: 'top', labels: {{ color: '#e0e0e0' }} }}
            }},
            scales: {{
              x: {{
                display: true,
                title: {{ display: true, text: 'Time (hours)', color: '#e0e0e0' }},
                ticks: {{ color: '#aaa' }},
                grid: {{ color: '#444' }}
              }},
              y: {{
                type: 'linear',
                display: true,
                position: 'left',
                title: {{ display: true, text: 'Battery State (Wh)', color: '#e0e0e0' }},
                ticks: {{ color: '#aaa' }},
                grid: {{ color: '#444' }},
                min: 0,
                max: batteryCapacity ? batteryCapacity * 1.1 : undefined
              }},
              y1: {{
                type: 'linear',
                display: true,
                position: 'right',
                title: {{ display: true, text: 'Power Generation (W)', color: '#e0e0e0' }},
                ticks: {{ color: '#aaa' }},
                grid: {{ drawOnChartArea: false, color: '#444' }}
              }}
            }}
          }}
        }});
      }}
    }}

    // Energy balance chart (embedded data; no fetch)
    let energyBalanceChart = null;

    function initEnergyBalanceChart() {{
      if (energyBalanceChart) return;

      const canvas = document.getElementById('energyBalanceChart');
      if (!canvas) return;

      const tabContent = canvas.closest('.tab-content');
      if (!tabContent || !tabContent.classList.contains('active')) return;

      if (!energyBalanceData && !sensitivityAnalysisData) {{
        const parent = canvas.parentElement;
        parent.innerHTML = "<div style='color:#aaa; font-size:13px; padding:10px; border:1px dashed #555; border-radius:6px;'>No energy balance data found.</div>";
        return;
      }}

      // Use the same helper as the sensitivity charts so the baseline matches the sensitivity 0% case.
      // Prefer an explicit 0-perturbation weight variant if available; otherwise fall back to baseline.
      let series48 = null;
      if (sensitivityAnalysisData && sensitivityAnalysisData.weight && sensitivityAnalysisData.weight["0"]) {{
        series48 = expandTo48h(sensitivityAnalysisData.weight["0"]);
      }} else if (energyBalanceData) {{
        series48 = expandTo48h(energyBalanceData);
      }} else {{
        series48 = {{
          time_hr: [],
          period_type: [],
          solar_flux_W_per_m2: [],
          power_generated_W: [],
          power_used_W: [],
          battery_state_Wh: [],
        }};
      }}

      const timeHours = series48.time_hr || [];
      const periodType = series48.period_type || [];
      const powerGenerated = series48.power_generated_W || [];
      const powerUsed = series48.power_used_W || [];
      const batteryState = series48.battery_state_Wh || [];
      if (!timeHours.length) return;

      canvas.style.height = '300px';
      canvas.style.maxHeight = '300px';

      // Calculate battery state range for y axis with padding
      const validBatteryState = batteryState.filter(v => !isNaN(v) && isFinite(v));
      const battMin = validBatteryState.length > 0 ? Math.min(...validBatteryState) : 0;
      const battMax = validBatteryState.length > 0 ? Math.max(...validBatteryState) : 100;
      const battRange = battMax - battMin;
      const battPadding = Math.max(battRange * 0.1, 5); // 10% padding or at least 5Wh
      const battYMin = Math.max(0, battMin - battPadding);
      const battYMax = battMax + battPadding;

      // Calculate power range for y1 axis with padding
      const allPowerValues = [...powerGenerated, ...powerUsed].filter(v => !isNaN(v) && isFinite(v));
      const powerMin = allPowerValues.length > 0 ? Math.min(...allPowerValues) : 0;
      const powerMax = allPowerValues.length > 0 ? Math.max(...allPowerValues) : 100;
      const powerRange = powerMax - powerMin;
      const powerPadding = Math.max(powerRange * 0.1, 10); // 10% padding or at least 10W
      const powerYMin = Math.min(0, powerMin - powerPadding);
      const powerYMax = powerMax + powerPadding;

      // Create background datasets for period highlighting
      // Use the maximum range that covers both y and y1 axes for full coverage
      const fullYMin = Math.min(battYMin, powerYMin);
      const fullYMax = Math.max(battYMax, powerYMax);
      
      // Create data arrays for period backgrounds (fill from min to max)
      const createPeriodBackground = (periodName, color) => {{
        const data = timeHours.map((t, i) => {{
          const period = periodType[i] || 'day';
          // Return fullYMax where period matches, fullYMin otherwise to create fill effect
          return (period === periodName || (periodName === 'dawn_dusk' && (period === 'dawn_dusk' || period === 'dusk'))) ? fullYMax : fullYMin;
        }});
        return {{
          label: periodName === 'dawn_dusk' ? 'Dawn/Dusk' : periodName.charAt(0).toUpperCase() + periodName.slice(1),
          data: data,
          backgroundColor: color,
          borderColor: 'transparent',
          borderWidth: 0,
          fill: {{ value: fullYMin }},
          pointRadius: 0,
          pointHoverRadius: 0,
          order: -1, // Render behind main datasets
          yAxisID: 'y',
          tension: 0
        }};
      }};

      const backgroundDatasets = [];
      if (periodType.length > 0) {{
        // Night background (dark blue)
        backgroundDatasets.push(createPeriodBackground('night', 'rgba(0, 0, 128, 0.20)'));
        // Dawn/Dusk background (orange)
        backgroundDatasets.push(createPeriodBackground('dawn_dusk', 'rgba(255, 165, 0, 0.18)'));
        // Day background (light yellow) - subtle indication
        backgroundDatasets.push(createPeriodBackground('day', 'rgba(255, 250, 205, 0.15)'));
      }}

      energyBalanceChart = new Chart(canvas, {{
        type: 'line',
        data: {{
          labels: timeHours.map(t => t.toFixed(2)),
          datasets: [
            ...backgroundDatasets,
            {{
              label: 'Power Generated (W)',
              data: powerGenerated,
              borderColor: 'rgb(76, 175, 80)',
              backgroundColor: 'rgba(76, 175, 80, 0.2)',
              yAxisID: 'y1',
              tension: 0.1,
              pointRadius: 0
            }},
            {{
              label: 'Power Used (W)',
              data: powerUsed,
              borderColor: 'rgb(244, 67, 54)',
              backgroundColor: 'rgba(244, 67, 54, 0.2)',
              yAxisID: 'y1',
              tension: 0.1,
              pointRadius: 0
            }},
            {{
              label: 'Battery State (Wh)',
              data: batteryState,
              borderColor: 'rgb(33, 150, 243)',
              backgroundColor: 'rgba(33, 150, 243, 0.2)',
              yAxisID: 'y',
              tension: 0.1,
              pointRadius: 0
            }}
          ]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          interaction: {{ mode: 'index', intersect: false }},
          plugins: {{
            title: {{
              display: true,
              text: `Energy Balance Analysis${{energyBalanceData.source_csv ? " (" + energyBalanceData.source_csv + ")" : ""}}`,
              color: '#4CAF50'
            }},
            legend: {{ display: true, position: 'top', labels: {{ color: '#e0e0e0' }} }}
          }},
          scales: {{
            x: {{
              display: true,
              title: {{ display: true, text: 'Time (hours)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }},
            y: {{
              type: 'linear',
              display: true,
              position: 'left',
              title: {{ display: true, text: 'Battery State (Wh)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }},
              min: battYMin,
              max: battYMax
            }},
            y1: {{
              type: 'linear',
              display: true,
              position: 'right',
              title: {{ display: true, text: 'Power (W)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ drawOnChartArea: false, color: '#444' }},
              min: powerYMin,
              max: powerYMax
            }}
          }}
        }}
      }});
    }}

    // ----------------------------
    // Sensitivity Analysis (48h) - in-browser recompute based on embedded baseline
    // ----------------------------
    const sensCharts = {{}};
    function expandTo48h(series) {{
      const t = series.time_hr || [];
      if (!t.length) return series;

      const tMax = Math.max(...t);
      // If it's already ~48h, do nothing
      if (tMax > 30) return series;

      // Duplicate day 1 -> day 2 (avoid duplicating the first point)
      const start = 1;
      const t2 = t.slice(start).map(v => Number(v) + 24);

      const dup = (arr) => (arr || []).concat((arr || []).slice(start));
      const dupShift = (arr) => (arr || []).concat((arr || []).slice(start));

      return {{
        time_hr: t.concat(t2),
        period_type: dupShift(series.period_type),
        solar_flux_W_per_m2: dup(series.solar_flux_W_per_m2),
        power_generated_W: dup(series.power_generated_W),
        power_used_W: dup(series.power_used_W),
        battery_state_Wh: dup(series.battery_state_Wh),
      }};
    }}

    
    function destroyIfExists(id) {{
      if (sensCharts[id]) {{
        try {{ sensCharts[id].destroy(); }} catch (e) {{}}
        delete sensCharts[id];
      }}
    }}


    // Builds a chart visually consistent with your existing energy balance plot,
    // plus an optional Solar Flux dataset on a third axis (to match your "bottom of page" 48h style).
    // If baselineSeries is provided, it will be drawn alongside the main series for comparison.
    function buildEnergyBalanceChart(canvasId, series, titleText, baselineSeries=null) {{
      const ZERO_PAD_WH = 20;   // ensures space past 0 Wh
      const ZERO_PAD_W  = 30;   // ensures space past 0 W
      const canvas = document.getElementById(canvasId);
      if (!canvas) return null;

      const timeHours = series.time_hr || [];
      const periodType = series.period_type || [];
      const powerGenerated = series.power_generated_W || [];
      const powerUsed = series.power_used_W || [];
      const batteryState = series.battery_state_Wh || [];
      const solarFlux = series.solar_flux_W_per_m2 || [];

      const basePowerGenerated = baselineSeries ? (baselineSeries.power_generated_W || []) : [];
      const basePowerUsed = baselineSeries ? (baselineSeries.power_used_W || []) : [];
      const baseBatteryState = baselineSeries ? (baselineSeries.battery_state_Wh || []) : [];
      const baseSolarFlux = baselineSeries ? (baselineSeries.solar_flux_W_per_m2 || []) : [];

      if (!timeHours.length) return null;

      // Axis ranges (same approach as initEnergyBalanceChart), including baseline if present
      const allBatteryValues = [...batteryState, ...baseBatteryState].filter(v => !isNaN(v) && isFinite(v));
      const validBatteryState = allBatteryValues;
      const battMin = validBatteryState.length > 0 ? Math.min(...validBatteryState) : 0;
      const battMax = validBatteryState.length > 0 ? Math.max(...validBatteryState) : 100;
      const battRange = battMax - battMin;
      const battPadding = Math.max(battRange * 0.15, ZERO_PAD_WH);
      const battYMin = Math.min(0 - ZERO_PAD_WH, battMin - battPadding);
      const battYMax = battMax + battPadding;
      const allPowerValues = [...powerGenerated, ...powerUsed, ...basePowerGenerated, ...basePowerUsed].filter(v => !isNaN(v) && isFinite(v));
      const powerMin = allPowerValues.length > 0 ? Math.min(...allPowerValues) : 0;
      const powerMax = allPowerValues.length > 0 ? Math.max(...allPowerValues) : 100;
      const powerRange = powerMax - powerMin;
      const powerPadding = Math.max(powerRange * 0.15, ZERO_PAD_W);
      const powerYMin = Math.min(0 - ZERO_PAD_W, powerMin - powerPadding);
      const powerYMax = powerMax + powerPadding;


      const fullYMin = Math.min(battYMin, powerYMin);
      const fullYMax = Math.max(battYMax, powerYMax);

      const createPeriodBackground = (periodName, color) => {{
        const data = timeHours.map((t, i) => {{
          const period = periodType[i] || 'day';
          return (period === periodName || (periodName === 'dawn_dusk' && (period === 'dawn_dusk' || period === 'dusk')))
            ? fullYMax
            : fullYMin;
        }});
        return {{
          label: periodName === 'dawn_dusk' ? 'Dawn/Dusk' : periodName.charAt(0).toUpperCase() + periodName.slice(1),
          data,
          backgroundColor: color,
          borderColor: 'transparent',
          borderWidth: 0,
          fill: {{ value: fullYMin }},
          pointRadius: 0,
          pointHoverRadius: 0,
          order: -1,
          yAxisID: 'y',
          tension: 0
        }};
      }};

      const backgroundDatasets = [];
      if (periodType.length > 0) {{
        backgroundDatasets.push(createPeriodBackground('night', 'rgba(0, 0, 128, 0.20)'));
        backgroundDatasets.push(createPeriodBackground('dawn_dusk', 'rgba(255, 165, 0, 0.18)'));
        backgroundDatasets.push(createPeriodBackground('day', 'rgba(255, 250, 205, 0.15)'));
      }}

      const datasets = [
        ...backgroundDatasets,
      ];

      // Optional baseline series (drawn underneath with dashed lines).
      // Only battery baseline is drawn to avoid cluttering the legend with extra power/flux entries.
      if (baselineSeries) {{
        datasets.push(
          {{
            label: 'Battery Energy (Base, Wh)',
            data: baseBatteryState,
            borderColor: 'rgba(33, 150, 243, 0.7)',
            backgroundColor: 'rgba(33, 150, 243, 0.1)',
            yAxisID: 'y',
            tension: 0.1,
            pointRadius: 0,
            borderDash: [4, 3]
          }}
        );
      }}

      // Main (selected) series
      datasets.push(
        {{
          label: 'Power Generated (W)',
          data: powerGenerated,
          borderColor: 'rgb(76, 175, 80)',
          backgroundColor: 'rgba(76, 175, 80, 0.2)',
          yAxisID: 'y1',
          tension: 0.1,
          pointRadius: 0
        }},
        {{
          label: 'Power Used (W)',
          data: powerUsed,
          borderColor: 'rgb(244, 67, 54)',
          backgroundColor: 'rgba(244, 67, 54, 0.2)',
          yAxisID: 'y1',
          tension: 0.1,
          pointRadius: 0
        }},
        {{
          label: 'Battery Energy (Wh)',
          data: batteryState,
          borderColor: 'rgb(33, 150, 243)',
          backgroundColor: 'rgba(33, 150, 243, 0.2)',
          yAxisID: 'y',
          tension: 0.1,
          pointRadius: 0
        }},
        // Solar flux (optional) to match your 48h "bottom of page" style
        ...(solarFlux && solarFlux.length ? [{{
          label: 'Solar Flux (W/m²)',
          data: solarFlux,
          borderColor: 'rgb(255, 152, 0)',
          backgroundColor: 'rgba(255, 152, 0, 0.0)',
          yAxisID: 'y2',
          tension: 0.1,
          pointRadius: 0,
          borderDash: [6, 4]
        }}] : [])
      );

      return new Chart(canvas, {{
        type: 'line',
        data: {{
          labels: timeHours.map(t => Number(t).toFixed(2)),
          datasets
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          interaction: {{ mode: 'index', intersect: false }},
          plugins: {{
            title: {{ display: true, text: titleText, color: '#4CAF50' }},
            legend: {{ display: true, position: 'top', labels: {{ color: '#e0e0e0' }} }}
          }},
          scales: {{
            x: {{
              display: true,
              title: {{ display: true, text: 'Time (hours)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }},
            y: {{
              type: 'linear',
              display: true,
              position: 'left',
              title: {{ display: true, text: 'Battery Energy (Wh)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }},
              min: battYMin,
              max: battYMax
            }},
            y1: {{
              type: 'linear',
              display: true,
              position: 'right',
              title: {{ display: true, text: 'Power (W)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ drawOnChartArea: false, color: '#444' }},
              min: powerYMin,
              max: powerYMax
            }},
            y2: {{
              type: 'linear',
              display: true,
              position: 'right',
              offset: true,
              title: {{ display: true, text: 'Solar Flux (W/m²)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ drawOnChartArea: false, color: '#444' }}
            }}
          }}
        }}
      }});
    }}

    // Get pre-computed sensitivity variant series
    function getVariantSeries(kind, selectionValue) {{
      // Prefer explicit sensitivity data if present
      if (sensitivityAnalysisData && sensitivityAnalysisData[kind]) {{
        const kindData = sensitivityAnalysisData[kind] || {{}};
        const variant = kindData[String(selectionValue)];
        if (variant) {{
          // Return pre-computed variant (may need to expand to 48h if only 24h)
          return expandTo48h(variant);
        }}
      }}

      // Fallback: use baseline energyBalanceData if available
      if (energyBalanceData) {{
        return expandTo48h(energyBalanceData);
      }}

      // Last-resort safe empty series (avoids JS crashes if data is missing)
      return {{
        time_hr: [],
        period_type: [],
        solar_flux_W_per_m2: [],
        power_generated_W: [],
        power_used_W: [],
        battery_state_Wh: [],
      }};
    }}

    function fillPercentSelect(selId) {{
      const sel = document.getElementById(selId);
      if (!sel) return;
      sel.innerHTML = '';
      const opts = [
        {{ label: '-20%', v: -0.2 }},
        {{ label: '-10%', v: -0.1 }},
        {{ label: '0%', v: 0.0 }},
        {{ label: '+10%', v: 0.1 }},
        {{ label: '+20%', v: 0.2 }}
      ];
      for (const o of opts) {{
        const el = document.createElement('option');
        el.value = String(o.v);
        el.textContent = o.label;
        if (o.v === 0) el.selected = true;
        sel.appendChild(el);
      }}
    }}

    function fillCgSelect(selId) {{
      const sel = document.getElementById(selId);
      if (!sel) return;
      sel.innerHTML = '';
      for (let cm = -40; cm <= 20; cm += 5) {{
        const el = document.createElement('option');
        el.value = String(cm);
        el.textContent = `${{cm}} cm`;
        if (cm === 0) el.selected = true;
        sel.appendChild(el);
      }}
    }}

    function initSensitivityCharts() {{
      // Only run if power tab is active (same guard style as initEnergyBalanceChart)
      const powerTab = document.getElementById('tab-power');
      if (!powerTab || !powerTab.classList.contains('active')) return;

      // Require sensitivity data; we can operate even if energyBalanceData is missing
      if (!sensitivityAnalysisData) return;

      const varSelect = document.getElementById('sensVarSelect');
      const deltaSelect = document.getElementById('sensDeltaSelect');
      const deltaLabel = document.getElementById('sensDeltaLabel');
      if (!varSelect || !deltaSelect || !deltaLabel) return;

      // Define variables available for sensitivity (must match keys in sensitivityAnalysisData)
      const varDefs = [
        {{ kind: 'weight', label: 'Weight' }},
        {{ kind: 'drag', label: 'Drag' }},
        {{ kind: 'solar', label: 'Solar (Irradiance)' }},
        {{ kind: 'cellEff', label: 'Solar Cell Efficiency' }},
        {{ kind: 'solarPower', label: 'Solar Power (Panel/Area Proxy)' }},
        {{ kind: 'motorEff', label: 'Motor Efficiency' }},
        {{ kind: 'cg', label: 'Center of Gravity Shift' }},
      ];

      // Populate variable select
      varSelect.innerHTML = '';
      for (const v of varDefs) {{
        const el = document.createElement('option');
        el.value = v.kind;
        el.textContent = v.label;
        varSelect.appendChild(el);
      }}
      // Default to weight
      varSelect.value = 'weight';

      function refreshDeltaOptions() {{
        const kind = varSelect.value || 'weight';
        if (kind === 'cg') {{
          // Use CG-specific options (in cm)
          fillCgSelect('sensDeltaSelect');
          deltaLabel.textContent = 'CG Shift (cm)';
        }} else {{
          // Use generic percent deltas
          fillPercentSelect('sensDeltaSelect');
          deltaLabel.textContent = 'Change (%)';
        }}
      }}

      function getCurrentSeries() {{
        const kind = varSelect.value || 'weight';
        const selVal = deltaSelect.value || '0';
        const series = getVariantSeries(kind, selVal);
        return {{ series, kind }};
      }}

      let mainSensChart = null;

      function renderSensitivity() {{
        const {{ series, kind }} = getCurrentSeries();
        const varDef = varDefs.find(v => v.kind === kind) || varDefs[0];
        const title = `Sensitivity: ${{varDef.label}}`;

        // Baseline series: 0-perturbation for the same variable (or fallback inside getVariantSeries)
        const baseSeries = getVariantSeries(kind, 0);

        // Destroy any existing chart so we can rebuild with the new baseline + variant
        if (mainSensChart) {{
          try {{ mainSensChart.destroy(); }} catch (e) {{}}
          mainSensChart = null;
        }}

        mainSensChart = buildEnergyBalanceChart('sensChart_main', series, title, baseSeries);
      }}

      // Initialize controls and chart
      refreshDeltaOptions();
      renderSensitivity();

      // Wire up change handlers
      varSelect.addEventListener('change', () => {{
        refreshDeltaOptions();
        renderSensitivity();
      }});
      deltaSelect.addEventListener('change', () => {{
        renderSensitivity();
      }});
    }}


    // XFLR loads tab (embedded data; no fetch)
    let xflrLoadsInitialized = false;
    let xflrLiftChart = null;
    let xflrMomentChart = null;
    let xflrBendingChart = null;

    function destroyChart(ch) {{
      if (!ch) return null;
      try {{ ch.destroy(); }} catch (e) {{}}
      return null;
    }}

    function fmtNum(v, digits=6) {{
      const x = Number(v);
      if (!isFinite(x)) return '';
      // Avoid scientific notation for common ranges while keeping precision
      return x.toFixed(digits).replace(/\\.?0+$/, '');
    }}

    function setXflrMessage(text) {{
      const el = document.getElementById('xflrLoadsMessage');
      if (el) el.textContent = text || '';
    }}

    function getXflrCases() {{
      const cases = xflrLoadsData && xflrLoadsData.cases ? xflrLoadsData.cases : [];
      return Array.isArray(cases) ? cases : [];
    }}

    function getSelectedXflrCase() {{
      const caseSel = document.getElementById('xflrCaseSelect');
      const cases = getXflrCases();
      if (!caseSel || !cases.length) return null;
      const id = caseSel.value;
      return cases.find(c => c.id === id) || cases[0];
    }}

    function getSelectedXflrSurface(caseObj) {{
      const surfSel = document.getElementById('xflrSurfaceSelect');
      if (!surfSel || !caseObj) return null;
      const surfaces = caseObj.surfaces || {{}};
      const keys = Object.keys(surfaces);
      if (!keys.length) return null;
      const name = surfSel.value;
      return surfaces[name] ? name : keys[0];
    }}

    function populateXflrSelectors() {{
      const caseSel = document.getElementById('xflrCaseSelect');
      const surfSel = document.getElementById('xflrSurfaceSelect');
      const cases = getXflrCases();
      if (!caseSel || !surfSel) return;

      caseSel.innerHTML = '';
      for (const c of cases) {{
        const opt = document.createElement('option');
        opt.value = c.id;
        opt.textContent = c.label || c.id;
        caseSel.appendChild(opt);
      }}

      const c0 = getSelectedXflrCase();
      surfSel.innerHTML = '';
      if (c0 && c0.surfaces) {{
        const names = Object.keys(c0.surfaces);
        names.sort();
        for (const s of names) {{
          const opt = document.createElement('option');
          opt.value = s;
          opt.textContent = s;
          surfSel.appendChild(opt);
        }}
      }}
    }}

    function updateXflrSurfaceOptions() {{
      const caseObj = getSelectedXflrCase();
      const surfSel = document.getElementById('xflrSurfaceSelect');
      if (!surfSel) return;
      const prev = surfSel.value;
      surfSel.innerHTML = '';
      if (!caseObj || !caseObj.surfaces) return;
      const names = Object.keys(caseObj.surfaces);
      names.sort();
      for (const s of names) {{
        const opt = document.createElement('option');
        opt.value = s;
        opt.textContent = s;
        surfSel.appendChild(opt);
      }}
      if (names.includes(prev)) surfSel.value = prev;
    }}

    function buildXflrDelimited(stations, meta, includeMeta=false, delimiter=',') {{
      const lines = [];
      if (includeMeta) {{
        lines.push(`# source_file: ${{meta.source_file || ''}}`);
        if (meta.alpha_deg !== null && meta.alpha_deg !== undefined) lines.push(`# alpha_deg: ${{meta.alpha_deg}}`);
        if (meta.qinf_mps !== null && meta.qinf_mps !== undefined) lines.push(`# qinf_mps: ${{meta.qinf_mps}}`);
        if (meta.rho_kg_m3 !== null && meta.rho_kg_m3 !== undefined) lines.push(`# rho_kg_m3: ${{meta.rho_kg_m3}}`);
        if (meta.surface) lines.push(`# surface: ${{meta.surface}}`);
        lines.push(`# units: station_y_m(m), chord_m(m), Lprime(N/m), Mprime(Nm/m), BM(Nm)`);
      }}

      const header = [
        'station_y_m',
        'chord_m',
        'Cl',
        'CmGeom',
        'CmAirf_c4',
        'Lprime_N_per_m',
        'Mprime_CmGeom_Nm_per_m',
        'Mprime_CmAirf_Nm_per_m',
        'BM'
      ].join(delimiter);
      lines.push(header);

      for (const s of stations) {{
        const row = [
          fmtNum(s.station_y_m, 6),
          fmtNum(s.chord_m, 6),
          fmtNum(s.Cl, 6),
          fmtNum(s.CmGeom, 6),
          fmtNum(s.CmAirf_c4, 6),
          fmtNum(s.Lprime_N_per_m, 6),
          fmtNum(s.Mprime_CmGeom_Nm_per_m, 6),
          fmtNum(s.Mprime_CmAirf_Nm_per_m, 6),
          fmtNum(s.BM, 6),
        ].join(delimiter);
        lines.push(row);
      }}
      return lines.join('\\n');
    }}

    function buildXflrCsv(stations, meta, includeMeta=false) {{
      return buildXflrDelimited(stations, meta, includeMeta, ',');
    }}

    function buildXflrTsv(stations, meta) {{
      // TSV is the most reliable for copy/paste into Google Sheets / Excel
      // (commas sometimes paste into a single column depending on locale/settings)
      return buildXflrDelimited(stations, meta, false, '\\t');
    }}

    async function copyTextToClipboard(text) {{
      try {{
        await navigator.clipboard.writeText(text);
        return true;
      }} catch (e) {{
        return false;
      }}
    }}

    function fallbackShowCsv(text) {{
      const ta = document.getElementById('xflrCsvFallback');
      if (!ta) return;
      ta.value = text;
      ta.style.display = 'block';
      ta.focus();
      ta.select();
      try {{ document.execCommand('copy'); }} catch (e) {{}}
    }}

    function safeFilename(s) {{
      return String(s || 'xflr_loads').replace(/[^a-zA-Z0-9._-]+/g, '_');
    }}

    function downloadCsv(text, filename) {{
      const blob = new Blob([text], {{ type: 'text/csv;charset=utf-8' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      setTimeout(() => URL.revokeObjectURL(url), 2000);
    }}

    function renderXflrLoads() {{
      const caseObj = getSelectedXflrCase();
      if (!caseObj) {{
        setXflrMessage('No XFLR loads found in this run (missing output/<run>/XFLR/*_a=*_v=*ms.txt).');
        return;
      }}
      const surfaceName = getSelectedXflrSurface(caseObj);
      if (!surfaceName) {{
        setXflrMessage('No surfaces found in selected case.');
        return;
      }}
      const surface = (caseObj.surfaces || {{}})[surfaceName] || {{}};
      const stations = surface.stations || [];
      if (!stations.length) {{
        setXflrMessage('No station data for selected surface.');
        return;
      }}

      setXflrMessage(`Source: ${{caseObj.source_file}} | ${{caseObj.label}} | ρ=${{fmtNum(caseObj.rho_kg_m3, 4)}} kg/m³`);

      // Update table
      const tbody = document.querySelector('#xflrLoadsTable tbody');
      if (tbody) {{
        tbody.innerHTML = '';
        for (const s of stations) {{
          const tr = document.createElement('tr');
          const cells = [
            fmtNum(s.station_y_m, 4),
            fmtNum(s.chord_m, 4),
            fmtNum(s.Cl, 6),
            fmtNum(s.CmGeom, 6),
            fmtNum(s.CmAirf_c4, 6),
            fmtNum(s.Lprime_N_per_m, 3),
            fmtNum(s.Mprime_CmGeom_Nm_per_m, 3),
            fmtNum(s.Mprime_CmAirf_Nm_per_m, 3),
            fmtNum(s.BM, 3),
          ];
          for (const c of cells) {{
            const td = document.createElement('td');
            td.textContent = c;
            td.style.padding = '6px 8px';
            td.style.borderBottom = '1px solid #2a2a2a';
            tr.appendChild(td);
          }}
          tbody.appendChild(tr);
        }}
      }}

      // Charts
      xflrLiftChart = destroyChart(xflrLiftChart);
      xflrMomentChart = destroyChart(xflrMomentChart);
      xflrBendingChart = destroyChart(xflrBendingChart);

      const liftCanvas = document.getElementById('xflrLiftChart');
      const momCanvas = document.getElementById('xflrMomentChart');
      const bendCanvas = document.getElementById('xflrBendingChart');
      if (!liftCanvas || !momCanvas || !bendCanvas) return;

      const liftPoints = stations.map(s => ({{ x: s.station_y_m, y: s.Lprime_N_per_m }}));
      const mGeomPoints = stations.map(s => ({{ x: s.station_y_m, y: s.Mprime_CmGeom_Nm_per_m }}));
      const mAirfPoints = stations.map(s => ({{ x: s.station_y_m, y: s.Mprime_CmAirf_Nm_per_m }}));
      const bmPoints = stations.map(s => ({{ x: s.station_y_m, y: s.BM }}));

      xflrLiftChart = new Chart(liftCanvas, {{
        type: 'line',
        data: {{
          datasets: [{{
            label: "Lift per span L' (N/m)",
            data: liftPoints,
            borderColor: 'rgb(76, 175, 80)',
            backgroundColor: 'rgba(76, 175, 80, 0.15)',
            tension: 0.05,
            pointRadius: 0
          }}]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          parsing: false,
          interaction: {{ mode: 'index', intersect: false }},
          plugins: {{
            title: {{ display: true, text: "Lift per span", color: '#4CAF50' }},
            legend: {{ display: true, position: 'top', labels: {{ color: '#e0e0e0' }} }}
          }},
          scales: {{
            x: {{
              type: 'linear',
              title: {{ display: true, text: 'Station y (m)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }},
            y: {{
              title: {{ display: true, text: \"L' (N/m)\", color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }}
          }}
        }}
      }});

      xflrMomentChart = new Chart(momCanvas, {{
        type: 'line',
        data: {{
          datasets: [
            {{
              label: \"Pitching moment per span M' from CmGeom (N·m/m)\",
              data: mGeomPoints,
              borderColor: 'rgb(33, 150, 243)',
              backgroundColor: 'rgba(33, 150, 243, 0.10)',
              tension: 0.05,
              pointRadius: 0
            }},
            {{
              label: \"Pitching moment per span M' from CmAirf@c/4 (N·m/m)\",
              data: mAirfPoints,
              borderColor: 'rgb(255, 152, 0)',
              backgroundColor: 'rgba(255, 152, 0, 0.10)',
              tension: 0.05,
              pointRadius: 0
            }}
          ]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          parsing: false,
          interaction: {{ mode: 'index', intersect: false }},
          plugins: {{
            title: {{ display: true, text: 'Pitching moment per span', color: '#4CAF50' }},
            legend: {{ display: true, position: 'top', labels: {{ color: '#e0e0e0' }} }}
          }},
          scales: {{
            x: {{
              type: 'linear',
              title: {{ display: true, text: 'Station y (m)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }},
            y: {{
              title: {{ display: true, text: \"M' (N·m/m)\", color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }}
          }}
        }}
      }});

      xflrBendingChart = new Chart(bendCanvas, {{
        type: 'line',
        data: {{
          datasets: [{{
            label: 'Bending moment (BM)',
            data: bmPoints,
            borderColor: 'rgb(244, 67, 54)',
            backgroundColor: 'rgba(244, 67, 54, 0.10)',
            tension: 0.05,
            pointRadius: 0
          }}]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          parsing: false,
          interaction: {{ mode: 'index', intersect: false }},
          plugins: {{
            title: {{ display: true, text: 'Bending moment (XFLR BM column)', color: '#4CAF50' }},
            legend: {{ display: true, position: 'top', labels: {{ color: '#e0e0e0' }} }}
          }},
          scales: {{
            x: {{
              type: 'linear',
              title: {{ display: true, text: 'Station y (m)', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }},
            y: {{
              title: {{ display: true, text: 'BM', color: '#e0e0e0' }},
              ticks: {{ color: '#aaa' }},
              grid: {{ color: '#444' }}
            }}
          }}
        }}
      }});
    }}

    function initXflrLoadsTab() {{
      if (xflrLoadsInitialized) return;
      xflrLoadsInitialized = true;

      const caseSel = document.getElementById('xflrCaseSelect');
      const surfSel = document.getElementById('xflrSurfaceSelect');
      const copyBtn = document.getElementById('xflrCopyCsvBtn');
      const dlBtn = document.getElementById('xflrDownloadCsvBtn');

      const cases = getXflrCases();
      if (!cases.length) {{
        setXflrMessage('No XFLR loads found in this run (expected output/<run>/XFLR/*_a=*_v=*ms.txt).');
        if (caseSel) caseSel.disabled = true;
        if (surfSel) surfSel.disabled = true;
        if (copyBtn) copyBtn.disabled = true;
        if (dlBtn) dlBtn.disabled = true;
        return;
      }}

      populateXflrSelectors();

      function currentStations() {{
        const c = getSelectedXflrCase();
        if (!c) return {{ stations: [], meta: {{}} }};
        const sName = getSelectedXflrSurface(c);
        const surf = sName ? (c.surfaces || {{}})[sName] : null;
        return {{
          stations: (surf && surf.stations) ? surf.stations : [],
          meta: {{ ...c, surface: sName }}
        }};
      }}

      if (caseSel) {{
        caseSel.addEventListener('change', () => {{
          updateXflrSurfaceOptions();
          renderXflrLoads();
        }});
      }}
      if (surfSel) {{
        surfSel.addEventListener('change', () => {{
          renderXflrLoads();
        }});
      }}

      if (copyBtn) {{
        copyBtn.addEventListener('click', async () => {{
          const {{ stations, meta }} = currentStations();
          const tsv = buildXflrTsv(stations, meta);
          const ok = await copyTextToClipboard(tsv);
          if (ok) {{
            setXflrMessage('Copied (TSV). Paste into Google Sheets (cell A1).');
          }} else {{
            setXflrMessage('Clipboard blocked by browser. Data shown below—select all and copy.');
            fallbackShowCsv(tsv);
          }}
        }});
      }}

      if (dlBtn) {{
        dlBtn.addEventListener('click', () => {{
          const {{ stations, meta }} = currentStations();
          const csv = buildXflrCsv(stations, meta, false);
          const fname = safeFilename(`xflr_loads_${{meta.id}}_${{meta.surface}}.csv`);
          downloadCsv(csv, fname);
        }});
      }}

      renderXflrLoads();
    }}

    // ----------------------------
    // Mass Pie Charts (baseline vs build)
    // ----------------------------
    let massPieCharts = [];
    function destroyMassPieChart() {{
      if (Array.isArray(massPieCharts) && massPieCharts.length) {{
        massPieCharts.forEach(chart => {{
          try {{ chart.destroy(); }} catch (e) {{}}
        }});
      }}
      massPieCharts = [];
    }}

    function createMassPieChart(canvasId, labelId, series, titleText) {{
      const canvas = document.getElementById(canvasId);
      if (!canvas || !Array.isArray(series) || series.length === 0) return;

      const totalMass = series.reduce((sum, item) => sum + (item.mass || 0), 0);
      if (totalMass === 0) return;

      const labelElem = labelId ? document.getElementById(labelId) : null;
      if (labelElem) {{
        labelElem.textContent = titleText;
      }}

      // Generate colors for pie chart segments
      const colors = [
        'rgb(76, 175, 80)',   // Green
        'rgb(33, 150, 243)',  // Blue
        'rgb(255, 152, 0)',   // Orange
        'rgb(244, 67, 54)',   // Red
        'rgb(156, 39, 176)',  // Purple
        'rgb(0, 188, 212)',   // Cyan
        'rgb(255, 235, 59)',  // Yellow
        'rgb(121, 85, 72)',   // Brown
        'rgb(96, 125, 139)',  // Blue Grey
        'rgb(255, 87, 34)',   // Deep Orange
        'rgb(63, 81, 181)',   // Indigo
        'rgb(233, 30, 99)',   // Pink
        'rgb(205, 220, 57)',  // Lime
        'rgb(0, 150, 136)',   // Teal
        'rgb(139, 195, 74)',  // Light Green
        'rgb(3, 169, 244)',   // Light Blue
        'rgb(255, 193, 7)',   // Amber
      ];

      const chartData = series.map((item, idx) => {{
        const mass = item.mass || 0;
        const pct = (mass / totalMass * 100).toFixed(1);
        return {{
          label: `${{item.label}}: ${{mass.toFixed(3)}} kg (${{pct}}%)`,
          value: mass,
          color: colors[idx % colors.length]
        }};
      }});

      const chart = new Chart(canvas, {{
        type: 'pie',
        data: {{
          labels: chartData.map(d => d.label.split(':')[0]),
          datasets: [{{
            data: chartData.map(d => d.value),
            backgroundColor: chartData.map(d => d.color),
            // Remove slice borders for a cleaner look
            borderColor: 'transparent',
            borderWidth: 0
          }}]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          plugins: {{
            // Use the HTML header / label as the only visible title to avoid duplicates
            title: {{
              display: false,
              text: titleText,
              color: '#4CAF50',
              font: {{ size: 16 }}
            }},
            legend: {{
              // Hide legend for the mass pie charts
              display: false
            }},
            tooltip: {{
              callbacks: {{
                label: function(context) {{
                  const label = context.label || '';
                  const value = context.parsed || 0;
                  const pct = ((value / totalMass) * 100).toFixed(1);
                  return `${{label}}: ${{value.toFixed(3)}} kg (${{pct}}%)`;
                }}
              }},
              backgroundColor: 'rgba(0, 0, 0, 0.8)',
              titleColor: '#e0e0e0',
              bodyColor: '#e0e0e0',
              borderColor: '#444',
              borderWidth: 1
            }}
          }}
        }}
      }});

      massPieCharts.push(chart);
    }}

    function initMassPieChart() {{
      destroyMassPieChart();
      if (!massChartData || typeof massChartData !== 'object') {{
        return;
      }}

      const originalSeries = Array.isArray(massChartData.original) ? massChartData.original : [];
      const buildSeries = Array.isArray(massChartData.build) ? massChartData.build : [];
      const runName = massChartData.runName || 'Original';
      const buildName = massChartData.buildName || 'Build';

      if (!originalSeries.length && !buildSeries.length) {{
        return;
      }}

      if (originalSeries.length) {{
        createMassPieChart(
          'massPieChartOriginal',
          'massPieLabelOriginal',
          originalSeries,
          `${{runName}}`
        );
      }}

      if (buildSeries.length) {{
        createMassPieChart(
          'massPieChartBuild',
          'massPieLabelBuild',
          buildSeries,
          `${{buildName}}`
        );
      }}
    }}

    // Tab switching
    document.querySelectorAll('.tab-button').forEach(button => {{
      button.addEventListener('click', () => {{
        const tabId = button.getAttribute('data-tab');

        document.querySelectorAll('.tab-button').forEach(btn => btn.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(content => content.classList.remove('active'));

        button.classList.add('active');
        const tabContent = document.getElementById(`tab-${{tabId}}`);
        if (tabContent) {{
          tabContent.classList.add('active');
          if (tabId === 'power') {{
            // allow layout to settle before Chart.js init
            setTimeout(() => {{
              initEnergyBalanceChart();
              initSensitivityCharts();
            }}, 50);
          }}
          if (tabId === 'loads') {{
            // allow layout to settle before Chart.js init
            setTimeout(() => initXflrLoadsTab(), 50);
          }}
          if (tabId === 'mass') {{
            // allow layout to settle before Chart.js init
            setTimeout(() => initMassPieChart(), 50);
          }}
        }}
      }});
    }});
  </script>
</body>
</html>
"""


# ----------------------------
# Viewer generator
# ----------------------------

def create_threejs_viewer(
    soln_path: Path,
    mass_properties_path: Path,
    airplane_path: Optional[Path] = None,
    output_path: Optional[Path] = None,
) -> Path:
    """Create a Three.js-based viewer HTML file. Returns the output HTML path."""
    if output_path is None:
        output_path = soln_path.parent / "aircraft_mass_viewer.html"

    soln = json.loads(soln_path.read_text(encoding="utf-8"))
    mass_properties = json.loads(mass_properties_path.read_text(encoding="utf-8"))

    run_dir = soln_path.parent
    specs_html = format_specs_html(soln, mass_properties, run_dir)
    top_overlay_html = build_top_overlay_html(run_dir)
    
    # Calculate mass data for pie charts (baseline vs build)
    masses = soln.get("Masses", {}) or {}
    meta = soln.get("Meta", {}) or {}
    run_name = str(meta.get("run_id", "Original"))
    build_name, updates = load_builds(run_dir)
    mass_components = [
        ("Solar Cells", "mass_solar_cells"),
        ("Power Board", "mass_power_board"),
        ("Batteries", "mass_batteries"),
        ("Wires", "mass_wires"),
        ("Avionics", "mass_avionics"),
        ("Servos", "mass_servos"),
        ("Motors (Mounted)", "mass_motors_mounted"),
        ("ESCs", "mass_escs"),
        ("Propellers", "mass_propellers"),
        ("Main Wing", "mass_main_wing"),
        ("Horizontal Stabilizer", "mass_hstab"),
        ("Vertical Stabilizer", "mass_vstab"),
        ("Boom", "mass_boom"),
        ("Fuselages", "mass_fuselages"),
        ("Superstructures", "mass_superstructures"),
        ("Boom-VStab Interfaces", "mass_boom_vstab_interfaces"),
        ("VStab-HStab Interfaces", "mass_vstab_stab_interfaces"),
    ]
    mass_chart_original: List[Dict[str, float]] = []
    mass_chart_build: List[Dict[str, float]] = []
    for display_name, key in mass_components:
        orig_val = float(masses.get(key, 0) or 0)
        if orig_val > 0:
            mass_chart_original.append({
                "label": display_name,
                "mass": orig_val,
            })

        build_val = float(updates.get(key, orig_val))
        if build_val > 0:
            mass_chart_build.append({
                "label": display_name,
                "mass": build_val,
            })

    mass_chart_data = {
        "runName": run_name,
        "buildName": build_name or "Mojave 1",
        "original": mass_chart_original,
        "build": mass_chart_build,
    }

    airplane_geometry: JSONDict = {}
    if airplane_path and airplane_path.suffix.lower() == ".aero" and airplane_path.exists():
        airplane_geometry = extract_airplane_geometry(airplane_path)

    energy_balance_data = load_energy_balance_csv(run_dir)
    xflr_loads_data = load_xflr_loads(run_dir, soln)
    sensitivity_analysis_data = load_sensitivity_analysis_json(run_dir)

    # Compact JSON to keep HTML size reasonable
    soln_json = json.dumps(soln, separators=(",", ":"))
    mass_properties_json = json.dumps(mass_properties, separators=(",", ":"))
    airplane_geometry_json = json.dumps(airplane_geometry, separators=(",", ":"))
    energy_balance_data_json = json.dumps(energy_balance_data, separators=(",", ":")) if energy_balance_data else "null"
    xflr_loads_data_json = json.dumps(xflr_loads_data, separators=(",", ":")) if xflr_loads_data else "null"
    sensitivity_analysis_json = json.dumps(sensitivity_analysis_data, separators=(",", ":")) if sensitivity_analysis_data else "null"
    mass_chart_data_json = json.dumps(mass_chart_data, separators=(",", ":"))

    html_content = HTML_TEMPLATE.format(
        specs_html=specs_html,
        top_overlay_html=top_overlay_html,
        soln_json=soln_json,
        mass_properties_json=mass_properties_json,
        airplane_geometry_json=airplane_geometry_json,
        energy_balance_data_json=energy_balance_data_json,
        xflr_loads_data_json=xflr_loads_data_json,
        sensitivity_analysis_json=sensitivity_analysis_json,
        mass_chart_data_json=mass_chart_data_json,
    )

    output_path.write_text(html_content, encoding="utf-8")
    print(f"[viewer] Saved Three.js visualization: {output_path}")
    return output_path


# ----------------------------
# CLI
# ----------------------------

def resolve_inputs(arg: Optional[str]) -> Tuple[Path, Path, Optional[Path], Path]:
    """
    Resolve file paths from CLI argument.
    - If arg is a directory: use that directory as run_dir.
    - If arg is a name: interpret as output/<name> under this script folder.
    - If arg is omitted: use latest run dir under output/.
    """
    base_dir = Path(__file__).resolve().parent
    output_dir_base = base_dir / "output"

    if arg is None:
        run_dir = latest_run_dir(output_dir_base)
    else:
        p = Path(arg)
        if p.exists() and p.is_dir():
            run_dir = p.resolve()
        else:
            run_dir = (output_dir_base / arg).resolve()

    soln_path = run_dir / "soln.json"
    mass_properties_path = run_dir / "mass_properties.json"

    airplane_path = run_dir / "airplane.aero"
    if not airplane_path.exists():
        airplane_path = None

    return soln_path, mass_properties_path, airplane_path, run_dir


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    soln_path, mass_properties_path, airplane_path, run_dir = resolve_inputs(arg)

    if not soln_path.exists():
        raise FileNotFoundError(f"Solution file not found: {soln_path}")
    if not mass_properties_path.exists():
        raise FileNotFoundError(f"Mass properties file not found: {mass_properties_path}")

    print(f"[viewer] Loading solution from: {soln_path}")
    print(f"[viewer] Loading mass properties from: {mass_properties_path}")
    if airplane_path:
        print(f"[viewer] Loading airplane geometry from: {airplane_path}")

    out_html = create_threejs_viewer(
        soln_path,
        mass_properties_path,
        airplane_path,
        run_dir / "aircraft_mass_viewer.html",
    )
    print(f"[viewer] Complete! Open {out_html} in a browser.")


if __name__ == "__main__":
    main()