"""Build-time exporter for the static design site.

Walks ``output/run_*``, computes each run's full state (geometry, mass properties, energy,
loads, specs) plus every as-built airframe/flight, and writes plain JSON into ``web/public/`` so
the Vite/React site can be published as a fully static bundle (no Python at runtime).

Usage:
    python export_site.py            # export every run under output/
    python export_site.py run_31sn   # export a single run
"""

import json
import sys
from pathlib import Path

import aerosandbox as asb

from lib.artifacts import to_jsonable
from process_webpage import gather_run_state, load_all_builds, run_metrics

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "output"
PUBLIC_DIR = BASE_DIR / "web" / "public"

# Precompute the one airfoil-derived constant the mass model needs, so the browser port of the
# mass math (lib/mass.ts) requires no aerosandbox. Matches process_mass.py:210.
FUSELAGE_RADIUS = 0.4 * asb.Airfoil("dae51").max_thickness()


def _airfoil_coords(name: str) -> list | None:
    """Return an airfoil's [x, y] outline (Selig order) for the 2D viewer, or None if unknown."""
    try:
        coords = asb.Airfoil(str(name)).coordinates
        return [[round(float(x), 5), round(float(y), 5)] for x, y in coords]
    except Exception:
        return None


def collect_airfoils(soln: dict) -> list:
    """Collect the distinct airfoils used by the aircraft with their outline coordinates."""
    wanted = [
        ("Main Wing", (soln.get("Main Wing", {}) or {}).get("wing_airfoil")),
        ("H Stabilizer", (soln.get("HStab", {}) or {}).get("hstab_airfoil")),
        ("V Stabilizer", (soln.get("V Stab", {}) or {}).get("vstab_airfoil")),
    ]
    out, seen = [], set()
    for surface, name in wanted:
        if not name or name in seen:
            continue
        seen.add(name)
        coords = _airfoil_coords(name)
        if coords:
            out.append({"surface": surface, "name": str(name), "coords": coords})
    return out


def _write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(to_jsonable(data), f, separators=(",", ":"))


def export_run(run_dir: Path) -> dict:
    """Exports one run to web/public/runs/<id>.json and returns its manifest entry."""
    state = gather_run_state(run_dir, persist=False)
    state["builds"] = load_all_builds(run_dir)
    state["constants"] = {"fuselage_radius": FUSELAGE_RADIUS}
    state["airfoils"] = collect_airfoils(state["soln"])
    metrics = run_metrics(state["soln"], state["massProperties"])

    # Sensitivity analysis can be many MB; write it as a sidecar the UI fetches only when needed.
    sensitivity = state.pop("sensitivity", None)
    state["hasSensitivity"] = bool(sensitivity)
    if sensitivity:
        _write_json(PUBLIC_DIR / "runs" / f"{run_dir.name}.sensitivity.json", sensitivity)

    _write_json(PUBLIC_DIR / "runs" / f"{run_dir.name}.json", state)

    # Emit a standalone mass_properties.json beside soln.json. The sibling
    # AircraftSim (AircraftView's Simulate button, ?aircraft=<run_id>) derives the
    # flown airframe from soln.json + this file; without it a run isn't flyable.
    if (state.get("massProperties") or {}).get("total"):
        _write_json(run_dir / "mass_properties.json", state["massProperties"])

    meta = state["soln"].get("Meta", {}) or {}
    return {
        "id": run_dir.name,
        "label": str(meta.get("run_id", run_dir.name)),
        "date": meta.get("timestamp_utc"),
        "metrics": metrics,
        "builds": [b["name"] for b in state["builds"]],
    }


def main(arg: str | None) -> None:
    if arg:
        p = Path(arg)
        run_dir = p if p.is_absolute() else OUTPUT_DIR / arg
        if not (run_dir / "soln.json").exists():
            raise SystemExit(f"Run not found (no soln.json): {run_dir}")
        run_dirs = [run_dir]
    else:
        run_dirs = sorted(p for p in OUTPUT_DIR.glob("run_*") if (p / "soln.json").exists())

    if not run_dirs:
        raise SystemExit("No runs found to export.")

    # On a FULL export, prune orphaned run JSONs in web/public/runs/ whose output/<run> dir no
    # longer exists, so the site's run list always matches what's actually under output/.
    if not arg:
        valid_ids = {rd.name for rd in run_dirs}
        runs_public = PUBLIC_DIR / "runs"
        if runs_public.exists():
            for f in runs_public.glob("run_*.json"):
                rid = f.name[:-len(".sensitivity.json")] if f.name.endswith(".sensitivity.json") else f.stem
                if rid not in valid_ids:
                    f.unlink()
                    print(f"[export] pruned orphan {f.name} (no output/{rid})")

    runs = []
    for rd in run_dirs:
        entry = export_run(rd)
        runs.append(entry)
        print(f"[export] {rd.name}: {len(entry['builds'])} build(s)")

    # Merge into any existing manifest so a single-run export doesn't drop the others.
    manifest_path = PUBLIC_DIR / "manifest.json"
    existing = {}
    if manifest_path.exists() and arg:
        try:
            for r in json.loads(manifest_path.read_text()).get("runs", []):
                existing[r["id"]] = r
        except (json.JSONDecodeError, OSError):
            pass
    for r in runs:
        existing[r["id"]] = r
    ordered = sorted(existing.values(), key=lambda r: (r.get("date") or "", r["id"]))
    _write_json(manifest_path, {"runs": ordered})
    print(f"[export] wrote manifest with {len(ordered)} run(s) -> {manifest_path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
