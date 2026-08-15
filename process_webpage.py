#!/usr/bin/env python3
"""
Assemble a run's consolidated state as JSON-able data for the static site.

Provides the data helpers consumed by ``export_site.py`` (and the mass CLI):
- ``gather_run_state`` — soln + live mass properties + energy/xflr/sensitivity + mass chart + specs
- ``build_specs_dict`` — structured, tabbed aircraft specifications (data, not HTML)
- ``build_mass_chart_data`` / ``run_metrics`` / ``load_all_builds`` — chart, index, and as-built data
- loaders for energy-balance CSV, XFLR5 distribution files, and sensitivity JSON
"""

from __future__ import annotations

import csv
import io
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

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


def load_all_builds(run_dir: Path) -> List[Dict[str, Any]]:
    """Load every as-built airframe/flight from builds.json (not just the first).

    Returns a list of ``{"name": str, "masses": {mass_*: float}, "meta": {...}}``. Falls back to
    the legacy single build via :func:`load_builds` when builds.json is absent.
    """
    builds_path = run_dir / "builds.json"
    result: List[Dict[str, Any]] = []
    if builds_path.exists():
        try:
            with open(builds_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            raw = data.get("builds") if isinstance(data, dict) else data
            if isinstance(raw, dict):
                raw = [raw]
            for entry in (raw or []):
                if not isinstance(entry, dict):
                    continue
                masses_raw = entry.get("masses") or entry.get("mass_updates") or {
                    k: v for k, v in entry.items() if k.startswith("mass_")
                }
                masses: Dict[str, float] = {}
                for k, v in masses_raw.items():
                    try:
                        val = float(v)
                        if val >= 0 and math.isfinite(val):
                            masses[k] = val
                    except (TypeError, ValueError):
                        continue
                meta = {k: v for k, v in entry.items() if k not in ("name", "masses", "mass_updates")
                        and not k.startswith("mass_")}
                result.append({"name": str(entry.get("name") or "Build"), "masses": masses, "meta": meta})
        except (json.JSONDecodeError, IOError) as e:
            print(f"[warning] Could not load builds from {builds_path}: {e}")
    if not result:
        name, updates = load_builds(run_dir)
        if name and updates:
            result.append({"name": name, "masses": updates, "meta": {}})
    return result


def build_specs_dict(soln: Mapping[str, Any], mass_properties: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Structured specs (tabs -> sections -> [label, value] rows) for the React SpecsPanel.

    This is the data-first replacement for the HTML-string ``format_specs_html`` (no canvases,
    no embedded chart markup), so the UI can render clean tabs and compare values.
    """
    g = lambda d, k, default=None: (soln.get(d, {}) or {}).get(k, default)
    total = mass_properties.get("total", {}) or {}

    def fmt(v, spec=".3f", suffix=""):
        if v is None:
            return "—"
        try:
            return f"{float(v):{spec}}{suffix}"
        except (TypeError, ValueError):
            return f"{v}{suffix}"

    def section(title, rows):
        return {"title": title, "rows": [[lbl, val] for lbl, val in rows if val is not None]}

    def mission_date(v):
        """soln stores mission date as a day-of-year integer; show a readable calendar date."""
        try:
            from datetime import datetime, timedelta
            doy = int(float(v))
            d = datetime(2025, 1, 1) + timedelta(days=doy - 1)
            return f"{d.strftime('%B %-d')} (day {doy})"
        except (TypeError, ValueError):
            return None if v is None else str(v)

    tabs: List[Dict[str, Any]] = [
        {"id": "overview", "label": "Overview", "sections": [
            section("Mission", [
                ["Run ID", str(g("Meta", "run_id", "N/A"))],
                ["Mission Date", mission_date(g("Mission", "mission_date"))],
                ["Operating Latitude", fmt(g("Mission", "operating_lat"), ".4f", "°")],
                ["Operating Altitude (MSL)", fmt(g("Mission", "operating_altitude"), ".0f", " m")],
            ]),
            section("Environment", [
                ["Temperature", fmt(g("Environment", "temperature_K"), ".2f", " K")],
                ["Pressure", fmt(g("Environment", "pressure"), ".0f", " Pa")],
                ["Density", fmt(g("Environment", "density"), ".4f", " kg/m³")],
            ]),
            section("Performance", [
                ["Total Mass", fmt(g("Performance", "total_mass"), ".3f", " kg")],
                ["Airspeed", fmt(g("Performance", "airspeed"), ".2f", " m/s")],
                ["Thrust (Cruise)", fmt(g("Performance", "thrust_cruise"), ".2f", " N")],
                ["Power (Cruise, all motors)", fmt(g("Performance", "power_cruise (all motors)"), ".2f", " W")],
                ["L/D Ratio", fmt(g("Performance", "L_over_D"), ".2f")],
                ["Stall Speed", fmt(g("Performance", "stall_speed"), ".2f", " m/s")],
            ]),
            section("Aerodynamics", [
                ["CL", fmt(g("Aerodynamics", "CL"), ".3f")],
                ["CD", fmt(g("Aerodynamics", "CD"), ".4f")],
                ["Cm", fmt(g("Aerodynamics", "Cm"), ".4f")],
                ["Static Margin", fmt(g("Aerodynamics", "static_margin"), ".3f")],
                ["Neutral Point X", fmt(g("Aerodynamics", "x_np"), ".3f", " m")],
            ]),
        ]},
        {"id": "geometry", "label": "Geometry & Mass", "sections": [
            section("Main Wing", [
                ["Wingspan", fmt(g("Main Wing", "wingspan"), ".3f", " m")],
                ["Chord", fmt(g("Main Wing", "chordlen"), ".3f", " m")],
                ["Area", fmt(g("Main Wing", "S_w"), ".3f", " m²")],
                ["Aspect Ratio", fmt(g("Main Wing", "main_wing_AR"), ".2f")],
                ["Airfoil", g("Main Wing", "wing_airfoil")],
            ]),
            section("Stabilizers", [
                ["HStab Span", fmt(g("HStab", "hstab_span"), ".3f", " m")],
                ["HStab Chord", fmt(g("HStab", "hstab_chordlen"), ".3f", " m")],
                ["HStab Volume Coeff", fmt(g("HStab", "V_H_actual"), ".3f")],
                ["VStab Span", fmt(g("V Stab", "vstab_span"), ".3f", " m")],
                ["VStab Root Chord", fmt(g("V Stab", "vstab_root_chord"), ".3f", " m")],
                ["VStab Volume Coeff", fmt(g("V Stab", "V_V_actual"), ".4f")],
            ]),
            section("Geometry", [
                ["Boom Length", fmt(g("Geometry", "boom_length"), ".3f", " m")],
                ["Boom Y Position", fmt(g("Geometry", "boom_y"), ".3f", " m")],
                ["CG Location (from LE)", fmt(g("Geometry", "cg_le_dist"), ".3f", " m")],
            ]),
            section("Mass Properties", [
                ["Total Mass", fmt(total.get("mass"), ".3f", " kg")],
                ["CG X", fmt(total.get("x_cg"), ".3f", " m")],
                ["Ixx", fmt(total.get("Ixx"), ".3f", " kg·m²")],
                ["Iyy", fmt(total.get("Iyy"), ".3f", " kg·m²")],
                ["Izz", fmt(total.get("Izz"), ".3f", " kg·m²")],
            ]),
        ]},
        {"id": "power", "label": "Power", "sections": [
            section("Power Budget", [
                ["Cruise Power (all motors)", fmt(g("Performance", "power_cruise (all motors)"), ".1f", " W")],
                ["Cruise Power (one motor)", fmt(g("Performance", "power_cruise (one motor)"), ".1f", " W")],
                ["Shaft Power (one motor)", fmt(g("Performance", "power_shaft_cruise (one motor)"), ".1f", " W")],
                ["Max Power Out", fmt(g("Performance", "power_out_max"), ".1f", " W")],
                ["Thrust (Cruise)", fmt(g("Performance", "thrust_cruise"), ".2f", " N")],
                ["Thrust (Climb)", fmt(g("Performance", "thrust_climb"), ".2f", " N")],
            ]),
            section("Power", [
                ["Solar Panels", fmt(g("Power", "solar_panels_n"), ".0f")],
                ["Solar Cell Efficiency", fmt(g("Power", "solar_cell_efficiency"), ".3f")],
                ["Battery Capacity", fmt(g("Power", "battery_capacity"), ".1f", " Wh")],
                ["Battery Voltage", fmt(g("Power", "battery_voltage"), ".1f", " V")],
                ["Battery Packs", fmt(g("Power", "num_packs"), ".0f")],
                ["Generation Margin", fmt(g("Power", "energy_generation_margin"), ".3f")],
            ]),
            section("Propulsion", [
                ["Propellers", fmt(g("Propulsion", "propeller_n"), ".0f")],
                ["Propeller Diameter", fmt(g("Propulsion", "propeller_diameter"), ".3f", " m")],
                ["Motor Kv", fmt(g("Propulsion", "motor_kv"), ".0f", " rpm/V")],
                ["Cruise RPM", fmt(g("Propulsion", "rpm_cruise"), ".0f")],
                ["Cruise Current", fmt(g("Propulsion", "i_cruise"), ".2f", " A")],
                ["Advance Ratio", fmt(g("Propulsion", "advanced_ratio"), ".3f")],
            ]),
        ]},
    ]
    # Drop empty sections/tabs so the UI never shows blank panels.
    for tab in tabs:
        tab["sections"] = [s for s in tab["sections"] if s["rows"]]
    return [t for t in tabs if t["sections"]]


def run_metrics(soln: Mapping[str, Any], mass_properties: Mapping[str, Any]) -> Dict[str, Any]:
    """Compact headline metrics used by the gallery/compare index."""
    total = mass_properties.get("total", {}) or {}
    perf = soln.get("Performance", {}) or {}
    aero = soln.get("Aerodynamics", {}) or {}
    main_wing = soln.get("Main Wing", {}) or {}
    return {
        "total_mass": total.get("mass"),
        "x_cg": total.get("x_cg"),
        "static_margin": aero.get("static_margin"),
        "L_over_D": perf.get("L_over_D"),
        "wingspan": main_wing.get("wingspan", main_wing.get("b_w")),
        "airspeed": perf.get("airspeed"),
    }


# ----------------------------
# Viewer generator
# ----------------------------

def build_mass_chart_data(soln: Mapping[str, Any], run_dir: Path) -> JSONDict:
    """Builds the baseline-vs-build mass pie-chart data for the viewer."""
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
            mass_chart_original.append({"label": display_name, "mass": orig_val})
        build_val = float(updates.get(key, orig_val))
        if build_val > 0:
            mass_chart_build.append({"label": display_name, "mass": build_val})
    return {
        "runName": run_name,
        "buildName": build_name or "Mojave 1",
        "original": mass_chart_original,
        "build": mass_chart_build,
    }


def gather_run_state(run_dir: Path, *, persist: bool = True) -> JSONDict:
    """Assembles the full consolidated state for a run and returns it as a JSON-able dict.

    Computes live mass properties from ``soln.json`` + the persisted ``layout`` in state.json and
    bundles every dataset the UI needs (energy balance, xflr loads, sensitivity, mass chart, specs
    HTML). When ``persist`` is True the fresh mass block is written back into state.json; the static
    exporter passes ``persist=False`` since it serializes the full state elsewhere and should not
    dirty the run directory.
    """
    from process_mass import compute_mass_properties  # local import to avoid heavy import at load
    from lib.artifacts import load_layout, update_state

    run_dir = Path(run_dir)
    soln = json.loads((run_dir / "soln.json").read_text(encoding="utf-8"))
    layout = load_layout(run_dir)
    mass_properties = compute_mass_properties(soln, layout)
    if persist:
        update_state(run_dir, mass=mass_properties)

    return {
        "run": run_dir.name,
        "soln": soln,
        "layout": layout,
        "massProperties": mass_properties,
        "massChart": build_mass_chart_data(soln, run_dir),
        "energyBalance": load_energy_balance_csv(run_dir),
        "xflrLoads": load_xflr_loads(run_dir, soln),
        "sensitivity": load_sensitivity_analysis_json(run_dir),
        "specs": build_specs_dict(soln, mass_properties),
    }