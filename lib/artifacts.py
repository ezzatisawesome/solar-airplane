import json
import random
import string
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as onp


def run_id_utc(prefix: str = "run") -> str:
    """Returns a filesystem-safe UTC timestamped run id."""
    return datetime.now(timezone.utc).strftime(f"{prefix}_%Y%m%dT%H%M%SZ")


def run_id_random(length: int = 4, prefix: str = "run") -> str:
    """Returns a filesystem-safe random alphanumeric run id."""
    chars = string.ascii_lowercase + string.digits
    random_suffix = ''.join(random.choice(chars) for _ in range(length))
    return f"{prefix}_{random_suffix}"


def to_jsonable(x: Any):
    """Converts common numpy-ish objects into JSON-serializable types."""
    if isinstance(x, (float, int, str, bool)) or x is None:
        return x

    if isinstance(x, (onp.floating, onp.integer)):
        return x.item()

    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]

    if isinstance(x, dict):
        return {str(k): to_jsonable(v) for k, v in x.items()}

    try:
        arr = onp.array(x)
        if arr.shape == ():
            return float(arr)
        return arr.tolist()
    except Exception:
        return str(x)


def write_json(path: Path, data: Dict[str, Any], *, indent: int = 2) -> None:
    """Writes JSON atomically (temp file + os.replace) so a crash mid-write can't corrupt the file."""
    import os

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        json.dump(to_jsonable(data), f, indent=indent)
    os.replace(tmp, path)


STATE_FILENAME = "state.json"


def state_path(run_dir: Path) -> Path:
    """Path to the consolidated per-run state file."""
    return Path(run_dir) / STATE_FILENAME


def load_state(run_dir: Path) -> Dict[str, Any]:
    """Loads the consolidated state.json for a run, or an empty dict if absent/corrupt."""
    p = state_path(run_dir)
    if not p.exists():
        return {}
    try:
        with open(p, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}


def load_layout(run_dir: Path) -> Dict[str, Any]:
    """Loads just the persisted point-mass ``layout`` (dragged positions) for a run."""
    layout = load_state(run_dir).get("layout", {})
    return layout if isinstance(layout, dict) else {}


def update_state(run_dir: Path, **fields: Any) -> Dict[str, Any]:
    """Merges ``fields`` into state.json (preserving existing keys) and writes it back."""
    state = load_state(run_dir)
    state.update(fields)
    write_json(state_path(run_dir), state)
    return state


def latest_run_dir(output_dir: Path, *, prefix: str = "run_") -> Path:
    """Returns the most-recent run directory inside `output_dir` (by name sort)."""
    output_dir = Path(output_dir)
    run_dirs = sorted([p for p in output_dir.glob(f"{prefix}*") if p.is_dir()])
    if not run_dirs:
        raise FileNotFoundError(f"No run dirs found in {output_dir} matching {prefix}*")
    return run_dirs[-1]


def _try_sol_value(sol: Any, x: Any) -> Tuple[bool, Any]:
    """Attempts to evaluate an Opti variable/expression via a solution object."""
    if sol is None:
        return False, x

    # AeroSandbox solutions are typically callable: sol(x)
    try:
        if callable(sol):
            return True, sol(x)
    except Exception:
        pass

    # Some solutions expose a `.value()` method.
    try:
        value_fn = getattr(sol, "value", None)
        if callable(value_fn):
            return True, value_fn(x)
    except Exception:
        pass

    return False, x


def process_raw_values(raw: Any, sol: Optional[Any] = None) -> Any:
    """Processes 'raw' values into JSON-safe python types.

    - Recursively walks dicts/lists/tuples
    - If `sol` is provided, tries to evaluate Opti symbols/expressions to numeric values
    - Converts numpy-ish scalars/arrays to JSON-safe types
    """
    if isinstance(raw, dict):
        return {str(k): process_raw_values(v, sol) for k, v in raw.items()}

    if isinstance(raw, (list, tuple)):
        return [process_raw_values(v, sol) for v in raw]

    _, v = _try_sol_value(sol, raw)
    return to_jsonable(v)
