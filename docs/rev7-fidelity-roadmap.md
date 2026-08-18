# opti_rev7 — Fidelity Roadmap

Master inventory of modeling fidelity: what's implemented, what it's simplified from,
and what remains. Origin codes:

- **Ours** — our own addition (driven by the rev6 stability/control notes), not in Adam's repo.
- **Ported** — adopted closely from Adam's repo (`/tmp/solar-uav-design`).
- **Surrogate** — Adam's *concept*, re-implemented here as a smooth/analytic optimizer-friendly simplification.
- **Adam-only** — exists in Adam's repo, not yet ported.
- **Neither** — identified gap, in neither codebase.

Status: ✅ implemented · 🟡 partial/surrogate · ❌ missing

---

## Implemented

| # | Feature | Domain | Status | Origin | Notes |
|---|---------|--------|--------|--------|-------|
| 1 | Aileron roll authority (strip-theory roll rate) | S&C | ✅ | Ours | roll_rate ≥ 30°/s |
| 2 | Weathercock / vertical-tail volume floor `V_V` | S&C | ✅ | Ours | yaw stiffness proxy |
| 3 | Dihedral effect floor | S&C | ✅ | Ours | spiral/roll coupling |
| 4 | Elevator pitch **trim** authority | S&C | ✅ | Ours | Cm_de·δe ≥ SM·dCL_trim (static only) |
| 5 | Static-margin bounds (SM_min/SM_max) | S&C | ✅ | Ours | |
| 6 | Wing taper (flat LE) + outboard polyhedral | Geometry | ✅ | Ours | |
| 7 | Faster loiter speed target | Mission | ✅ | Ours | |
| 8 | Amprius SA03 battery data (378 Wh/kg) | Power | ✅ | Ported | verified vs datasheet |
| 9 | Solar loss chain (cell η × encap × wiring × MPPT) | Power | ✅ | Ported | |
| 10 | NOCT cell-temp derate | Power | ✅ | Ported | assumes 75 °C noon (FLAGGED) |
| 11 | ASHRAE IAM incidence loss | Power | ✅ | Ported | |
| 12 | Row-aware panel packing (taper sheds rows) | Power | ✅ | Ported | from `_cell_quads` |
| 13 | MPPT string sizing (cells/string by topology) | Power | ✅ | Ported | single string length only |
| 14 | Solar-panel visualizer (web) | Viz | ✅ | Ported | |
| 15 | Real mass-weighted CG (fixed-point) | Mass | ✅ | Ported | |
| 16 | #3 Battery OCV + CV taper + I²R | Power | 🟡 | Surrogate | **linear** OCV vs Adam's 11-pt table |
| 17 | #4 Propeller η(J) map | Propulsion | 🟡 | Surrogate | **parabola** vs Adam's APC/Mejzlik tables |
| 18 | #2 Bus-voltage motor feasibility | Propulsion | 🟡 | Surrogate | V=iR+rpm/Kv ≤ bus; direct-drive only |
| 19 | #1 Climb energy debit (dawn) | Mission | 🟡 | Surrogate | refill half dropped intentionally |
| 20 | Prop specs surfaced in soln.json + web | Viz | ✅ | Ours | pitch, P/D, η, DxP label |

---

## Missing — Adam has it, not yet ported

| # | Feature | Domain | Status | Origin | Adam value / location | Priority |
|---|---------|--------|--------|--------|-----------------------|----------|
| 21 | Crosswind airspeed floor | Mission | ❌ | Adam-only | `MIN_AIRSPEED_WIND_MS = 7.7` | **T1** |
| 22 | Banked-loiter load factor | Aero/Struct | ❌ | Adam-only | `LOITER_STALL_FACTOR = 1.3` | **T1** |
| 23 | Climb-**rate** requirement | Mission | ❌ | Adam-only | `CLIMB_RATE_REQ_MS = 1.5` | **T1** |
| 24 | Pitch-acceleration authority | S&C (dyn) | ❌ | Adam-only | `PITCH_ACCEL_MIN_RAD_S2 = 0.30` | **T1** |
| 25 | Yaw/rudder unwind authority | S&C (dyn) | ❌ | Adam-only | "5°/s unwind" spec | T2 |
| 26 | Mass moment of inertia (pitch/yaw) | S&C (dyn) | ❌ | Adam-only | `pitch/yaw_inertia_kgm2` + mesh | **T1** (unlocks 24/25) |
| 27 | Geared motor + gearbox efficiency | Propulsion | ❌ | Adam-only | `motor.py` A40-12L kv410 | T2 |
| 28 | Discrete prop × free-Kv matcher | Propulsion | ❌ | Adam-only | `optimize_kv.py` | T2 |
| 29 | pvlib clear-sky irradiance (site/date/alt) | Power | ❌ | Adam-only | `environment.py` | T2 |
| 30 | Lifting-line washout Oswald factor | Aero | ❌ | Adam-only | `lifting_line.py` (we use VLM) | T3 |
| 31 | Mixed-string / per-bay MPPT | Power | ❌ | Adam-only | `run_mixed_strings.py` | T2 |
| 32 | NeuralFoil polars + airfoil selection | Aero | ❌ | Adam-only | `components/airfoils.py` | T2 |
| 33 | Prop mass from geometry | Mass | ❌ | Adam-only | `estimate_prop_mass_kg` | T3 |
| 34 | Monte-Carlo UQ + tornado | Process | ❌ | Adam-only | `monte_carlo.py` | **T2** (rigor) |
| 35 | Validation pipeline (bench ingest) | Process | ❌ | Adam-only | `validation.py` | T3 |
| 36 | ASB conservative-check discipline | Process | ❌ | Adam-only | `asb_physics.py` | T3 |
| 37 | Alternate empennage (π-tail / split) | Config | ❌ | Adam-only | `pi_tail_*/`, `split_empennage/` | T3 |

---

## Missing — neither codebase (identified gaps)

| # | Feature | Domain | Status | Origin | Priority |
|---|---------|--------|--------|--------|----------|
| 38 | Prop **diameter** as a variable (currently frozen 16") | Propulsion | ❌ | Neither | **T1** (quick) |
| 39 | Raise P/D bound (pinned at 0.85) | Propulsion | ❌ | Neither | **T1** (quick) |
| 40 | Structural spar mass from bending/gust loads | Structures | ❌ | Neither | **T1** (biggest optimism) |
| 41 | Reynolds-dependent airfoil polars across mission | Aero | ❌ | Neither | T2 |
| 42 | Trim-drag accounting (hstab download) | Aero | ❌ | Neither | T2 |
| 43 | Parasite/interference drag build-up | Aero | ❌ | Neither | T2 |
| 44 | Dynamic modes (phugoid / dutch-roll / spiral eigenvalues) | S&C (dyn) | ❌ | Neither | T2 |
| 45 | Servo/actuator hinge-moment sizing | Systems | ❌ | Neither | T3 |
| 46 | Loiter altitude as a variable | Mission | ❌ | Neither | T2 |
| 47 | Atmospheric density vs altitude | Environment | ❌ | Neither | T2 |
| 48 | Takeoff / landing / launch phase | Mission | ❌ | Neither | T2 |
| 49 | Payload power duty cycle | Mission | ❌ | Neither | T3 |
| 50 | Partial-shading solar derate | Power | ❌ | Neither | T2 (see 31) |
| 51 | MPPT part-load efficiency curve | Power | ❌ | Neither | T3 |
| 52 | Harness voltage drop & mass by current | Power | ❌ | Neither | T3 |
| 53 | Torsional stiffness → flutter/divergence | Structures | ❌ | Neither | T3 |
| 54 | Twin-boom bending / tail loads | Structures | ❌ | Neither | T3 |
| 55 | Joint / fitting / fastener mass | Mass | ❌ | Neither | T3 |
| 56 | Integrality (num_packs, cell count discrete) | Numerics | ❌ | Neither | T2 |
| 57 | Multi-start (global-minimum confidence) | Numerics | ❌ | Neither | T3 |
| 58 | Close the loop: fly soln.json in AircraftSim/JSBSim | Validation | ❌ | Neither | **T2** |

---

## Recommended execution order (physics-adding)

1. **Banked-loiter load factor (22)** + **crosswind floor (21)** — near-free, both under-counted today.
2. **Prop diameter variable + P/D bound (38, 39)** — 20-min win, we're pinned at a bound.
3. **Structural spar mass from loads (40)** — the single biggest change to the answer.
4. **Inertia (26) → dynamic pitch/yaw authority (24, 25) + climb rate (23)**.
5. **Mixed-string / partial-shading solar (31, 50)**.
6. **Monte-Carlo UQ (34)** — tells us which FLAGGED knob to attack next.
