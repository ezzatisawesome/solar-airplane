#!/usr/bin/env python3
"""
48-hour energy balance analysis for a solar aircraft in circular flight.

Reads aircraft/mission parameters from an artifact folder containing soln.json,
computes solar generation vs. power consumption (including turn load factor),
integrates battery state, classifies day/dawn-dusk/night, and exports plots + data.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from aerosandbox.library import power_solar

try:
    import aerosandbox as asb
    from aerosandbox.atmosphere import Atmosphere
except ImportError:
    asb = None
    Atmosphere = None

from lib.artifacts import latest_run_dir


# Workaround for some AeroSandbox installs where `aerosandbox.library.power_solar` forgets to import Atmosphere.
# (Seen as: NameError: name 'Atmosphere' is not defined inside power_solar.airmass()).
try:
    from aerosandbox.atmosphere import Atmosphere as _Atmosphere
    power_solar.Atmosphere = _Atmosphere
except Exception:
    # If AeroSandbox internals change, fail gracefully; downstream call will surface a clearer error.
    pass


# ----------------------------
# Data model / configuration
# ----------------------------

@dataclass(frozen=True)
class ArtifactData:
    # Mission
    mission_date: int = 100
    operating_lat: float = 37.398928
    operating_altitude: float = 1200.0

    # Performance
    airspeed: float = 10.0  # m/s
    power_cruise_all_motors: float = 50.0  # W

    # Power system
    solar_panels_n: int = 50
    battery_capacity_Wh: float = 400.0
    solar_panel_side_length_m: float = 0.125
    solar_cell_efficiency: float = 0.2187
    energy_generation_margin: float = 1.05
    allowable_battery_depth_of_discharge: float = 0.85

    # Avionics
    avionics_power_W: float = 9.0  # 8W avionics + 1W navlights


def _get(d: Dict[str, Any], key: str, default: Any) -> Any:
    """Helper: get key from dict with fallback."""
    v = d.get(key, default)
    return default if v is None else v


def load_artifact_data(artifact_path: Path) -> ArtifactData:
    soln_path = artifact_path / "soln.json"
    if not soln_path.exists():
        raise FileNotFoundError(f"Solution file not found: {soln_path}")

    with soln_path.open("r") as f:
        soln = json.load(f)

    mission = soln.get("Mission", {})
    performance = soln.get("Performance", {})
    power = soln.get("Power", {})

    return ArtifactData(
        mission_date=int(_get(mission, "mission_date", 100)),
        operating_lat=float(_get(mission, "operating_lat", 37.398928)),
        operating_altitude=float(_get(mission, "operating_altitude", 1200.0)),
        airspeed=float(_get(performance, "airspeed", 10.0)),
        power_cruise_all_motors=float(_get(performance, "power_cruise (all motors)", 50.0)),
        solar_panels_n=int(_get(power, "solar_panels_n", 50)),
        battery_capacity_Wh=float(_get(power, "battery_capacity", 400.0)),
        solar_panel_side_length_m=float(_get(power, "solar_panel_side_length", 0.125)),
        solar_cell_efficiency=float(_get(power, "solar_cell_efficiency", 0.2187)),
        energy_generation_margin=float(_get(power, "energy_generation_margin", 1.05)),
        allowable_battery_depth_of_discharge=float(_get(power, "allowable_battery_depth_of_discharge", 0.85)),
        avionics_power_W=9.0,
    )


# ----------------------------
# Physics / computation
# ----------------------------

def calculate_turn_power_W(power_straight_W: float, airspeed_mps: float, turn_radius_m: float) -> float:
    """
    Coordinated turn load factor n = sqrt(1 + (V^2/(gR))^2)
    Power ~ n^1.5 (approx drag scaling).
    """
    g = 9.81
    if turn_radius_m <= 0:
        raise ValueError("turn_radius_m must be > 0")

    n = np.sqrt(1.0 + (airspeed_mps**2 / (g * turn_radius_m)) ** 2)
    return power_straight_W * (n ** 1.5)


def classify_periods(solar_flux_Wm2: np.ndarray) -> np.ndarray:
    """
    Classify time points into 'night', 'dawn_dusk', 'day' based on solar flux.
    """
    max_flux = float(np.max(solar_flux_Wm2)) if solar_flux_Wm2.size else 0.0
    night_threshold = max(1.0, 0.01 * max_flux)
    dawn_dusk_threshold = 0.20 * max_flux

    # np.select is clearer than nested np.where for 3+ cases
    return np.select(
        [
            solar_flux_Wm2 < night_threshold,
            solar_flux_Wm2 < dawn_dusk_threshold,
        ],
        [
            "night",
            "dawn_dusk",
        ],
        default="day",
    )


def compute_energy_balance(
    a: ArtifactData,
    circular_diameter_m: float,
    n_points: int = 360,
    duration_hours: float = 48.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """
    Returns:
        time_hours, solar_flux, power_generated, power_used, net_energy_Wh,
        battery_state_Wh, period_type, battery_min_threshold_Wh, battery_capacity_Wh
    """
    if n_points < 2:
        raise ValueError("n_points must be >= 2")
    if circular_diameter_m <= 0:
        raise ValueError("circular_diameter_m must be > 0")

    # Time grid
    time_s = np.linspace(0.0, duration_hours * 3600.0, n_points)
    time_hours = time_s / 3600.0
    dt_hours = (time_s[1] - time_s[0]) / 3600.0

    # Turn power
    turn_radius_m = 0.5 * circular_diameter_m
    turn_power_W = calculate_turn_power_W(a.power_cruise_all_motors, a.airspeed, turn_radius_m)
    power_used_W = float(turn_power_W + a.avionics_power_W)  # constant in this scenario

    # Solar area
    solar_area_m2 = a.solar_panels_n * (a.solar_panel_side_length_m ** 2)

    # Solar flux vectorization:
    # Compute time-of-day in seconds for each point (repeat each 24h)
    time_of_day_s = np.mod(time_s, 24.0 * 3600.0).astype(float)

    # Call Aerosandbox solar_flux for each time point.
    # If solar_flux supports numpy arrays for time, this will be fast;
    # if it only supports scalars, np.vectorize will still simplify code.
    try:
        flux = power_solar.solar_flux(
            latitude=a.operating_lat,
            day_of_year=a.mission_date,
            time=time_of_day_s,
            altitude=a.operating_altitude,
            panel_azimuth_angle=0,
            panel_tilt_angle=0,
        )
        solar_flux_Wm2 = np.asarray(flux, dtype=float)
    except Exception:
        # Fallback: safe vectorization for scalar-only APIs
        solar_flux_Wm2 = np.vectorize(
            lambda t: power_solar.solar_flux(
                latitude=a.operating_lat,
                day_of_year=a.mission_date,
                time=float(t),
                altitude=a.operating_altitude,
                panel_azimuth_angle=0,
                panel_tilt_angle=0,
            ),
            otypes=[float],
        )(time_of_day_s)

    solar_flux_Wm2 = np.nan_to_num(solar_flux_Wm2, nan=0.0, posinf=0.0, neginf=0.0)
    solar_flux_Wm2 = np.maximum(solar_flux_Wm2, 0.0)

    # Power generated
    power_generated_W = solar_flux_Wm2 * solar_area_m2 * a.solar_cell_efficiency / a.energy_generation_margin

    # Net energy per step (Wh)
    power_used_series_W = np.full_like(power_generated_W, power_used_W)
    net_energy_Wh = (power_generated_W - power_used_series_W) * dt_hours

    # Battery: +7% oversizing, clamp at max only
    battery_capacity_Wh = a.battery_capacity_Wh * 1.07
    # Start at full capacity, then accumulate net energy changes
    battery_state_Wh = np.zeros_like(net_energy_Wh)
    battery_state_Wh[0] = battery_capacity_Wh
    for i in range(1, len(battery_state_Wh)):
        battery_state_Wh[i] = battery_state_Wh[i-1] + net_energy_Wh[i-1]
        # Only clamp to maximum, allow going below minimum
        battery_state_Wh[i] = min(battery_state_Wh[i], battery_capacity_Wh)

    battery_min_threshold_Wh = battery_capacity_Wh * (1.0 - a.allowable_battery_depth_of_discharge)

    period_type = classify_periods(solar_flux_Wm2)

    return (
        time_hours,
        solar_flux_Wm2,
        power_generated_W,
        power_used_series_W,
        net_energy_Wh,
        battery_state_Wh,
        period_type,
        float(battery_min_threshold_Wh),
        float(battery_capacity_Wh),
    )


# ----------------------------
# Battery integration with limits (for sensitivity analysis)
# ----------------------------

def estimate_battery_power_limits(
    time_hours: np.ndarray,
    battery_state_Wh: np.ndarray,
) -> Dict[str, float]:
    """
    Estimate effective charge/discharge power limits from a baseline energy trace.
    Uses finite differences dE/dt (Wh/hr == W) and returns conservative (95th percentile) limits.
    """
    if time_hours.size < 3 or battery_state_Wh.size != time_hours.size:
        return {"maxChargeW": float("inf"), "maxDischargeW": float("inf")}

    rates = []
    for i in range(time_hours.size - 1):
        dt = time_hours[i + 1] - time_hours[i]
        if dt <= 1e-9:
            continue
        dE = battery_state_Wh[i + 1] - battery_state_Wh[i]
        if np.isfinite(dE):
            rates.append(dE / dt)  # W

    if not rates:
        return {"maxChargeW": float("inf"), "maxDischargeW": float("inf")}

    rates = np.array(rates)
    pos = rates[rates > 0]
    neg = -rates[rates < 0]  # magnitudes

    def percentile(arr: np.ndarray, p: float) -> float:
        if arr.size == 0:
            return float("inf")
        idx = int(np.clip(p * (arr.size - 1), 0, arr.size - 1))
        return float(np.sort(arr)[idx])

    max_charge_W = percentile(pos, 0.95)
    max_discharge_W = percentile(neg, 0.95)

    return {
        "maxChargeW": max_charge_W if np.isfinite(max_charge_W) and max_charge_W > 0 else float("inf"),
        "maxDischargeW": max_discharge_W if np.isfinite(max_discharge_W) and max_discharge_W > 0 else float("inf"),
    }


def integrate_battery_with_limits(
    time_hours: np.ndarray,
    power_generated_W: np.ndarray,
    power_used_W: np.ndarray,
    E0_Wh: float,
    capacity_Wh: float,
    limits: Dict[str, float],
) -> np.ndarray:
    """
    Robust integration of battery energy (Wh) from power balance (W).
    - Trapezoidal rule (reduces stair-step artifacts)
    - Optional charge/discharge power limits (inferred from baseline curve)
    - Efficiency + saturation behavior without "energy injection" spikes
    """
    n = time_hours.size
    if n == 0:
        return np.array([])

    eta_charge = 0.98  # simple lumped charge efficiency
    eta_discharge = 0.98  # lumped discharge efficiency
    max_charge_W = limits.get("maxChargeW", float("inf"))
    max_discharge_W = limits.get("maxDischargeW", float("inf"))

    E = np.zeros(n)
    E[0] = np.clip(E0_Wh, 0, capacity_Wh)

    def net_power(i: int) -> float:
        gen = power_generated_W[i] if i < power_generated_W.size else 0.0
        use = power_used_W[i] if i < power_used_W.size else 0.0
        net = gen - use  # +charges battery, -discharges battery
        # Apply power limits (pre-efficiency)
        return np.clip(net, -max_discharge_W, max_charge_W)

    for i in range(n - 1):
        t0 = time_hours[i]
        t1 = time_hours[i + 1]
        dt = t1 - t0
        if dt <= 0:
            E[i + 1] = E[i]
            continue

        net0 = net_power(i)
        net1 = net_power(i + 1)
        # Trapezoid average power
        net_avg = 0.5 * (net0 + net1)

        # Efficiency: charging stores less than electrical surplus; discharging removes more than load deficit.
        if net_avg >= 0:
            net_avg = net_avg * eta_charge
        else:
            net_avg = net_avg / eta_discharge

        # Enforce saturation *continuously* (don't allow a single step to "teleport" past cap)
        next_E = E[i] + net_avg * dt

        if next_E > capacity_Wh:
            # If we're at/near full, surplus should be curtailed (no additional stored energy)
            next_E = capacity_Wh
        elif next_E < 0:
            next_E = 0.0

        E[i + 1] = next_E

    return E


# ----------------------------
# Weight sensitivity: compute power from weight using AeroSandbox
# ----------------------------

def compute_power_for_weight(
    artifact_path: Path,
    total_mass_kg: float,
    airspeed_mps: float,
    altitude_m: float,
    circular_diameter_m: float,
    avionics_power_W: float,
) -> float:
    """
    Compute power required for a given weight using AeroSandbox aerodynamics.
    
    For given weight: computes required lift L = mass × g
    Uses AeroSandbox AeroBuildup to find trimmed angle of attack (alpha) that produces required lift
    Computes drag at trimmed state
    Applies turn load factor: n = sqrt(1 + (V²/(gR))²), power scales as n^1.5
    Returns power required (including avionics)
    Falls back to linear scaling if airplane file not available
    """
    airplane_path = artifact_path / "airplane.aero"
    
    if not airplane_path.exists() or asb is None:
        # Fallback: linear scaling based on weight ratio
        # Load baseline mass from soln.json
        soln_path = artifact_path / "soln.json"
        if soln_path.exists():
            with soln_path.open("r") as f:
                soln = json.load(f)
            baseline_mass = float(soln.get("Performance", {}).get("total_mass", 9.0))
            baseline_power = float(soln.get("Performance", {}).get("power_cruise (all motors)", 50.0))
            # Scale power by (mass/baseline_mass)^1.5 (approximate drag scaling)
            mass_ratio = total_mass_kg / baseline_mass if baseline_mass > 0 else 1.0
            power_straight = baseline_power * (mass_ratio ** 1.5)
        else:
            # Last resort: assume 50W baseline
            mass_ratio = total_mass_kg / 9.0
            power_straight = 50.0 * (mass_ratio ** 1.5)
        
        turn_radius_m = 0.5 * circular_diameter_m
        turn_power = calculate_turn_power_W(power_straight, airspeed_mps, turn_radius_m)
        return float(turn_power + avionics_power_W)
    
    try:
        # Load airplane
        airplane = asb.load(str(airplane_path))
        
        # Required lift
        g = 9.81
        required_lift_N = total_mass_kg * g
        
        # Get wing area from airplane
        main_wing = None
        for wing in getattr(airplane, "wings", []):
            if hasattr(wing, "name") and ("main" in str(wing.name).lower() or "wing" in str(wing.name).lower()):
                main_wing = wing
                break
        if main_wing is None and len(getattr(airplane, "wings", [])) > 0:
            main_wing = airplane.wings[0]
        
        if main_wing is None:
            raise ValueError("Could not find main wing")
        
        wing_area_m2 = float(main_wing.area())
        
        # Create atmosphere
        try:
            from aerosandbox.atmosphere import Atmosphere as ASB_Atmosphere
            atm = ASB_Atmosphere(altitude_m)
        except:
            # Fallback
            rho = 1.225 * np.exp(-altitude_m / 8400.0)  # Simple exponential model
            atm = None
        
        # Find trimmed alpha (angle of attack that produces required lift)
        # Use binary search or iterative approach
        alpha_deg = 0.0
        tolerance = 0.01  # 0.01 N tolerance
        max_iter = 20
        
        for _ in range(max_iter):
            op_point = asb.OperatingPoint(
                atmosphere=atm if atm else None,
                velocity=airspeed_mps,
                alpha=alpha_deg,
            )
            
            aero = asb.AeroBuildup(
                airplane=airplane,
                op_point=op_point,
            ).run()
            
            lift_N = float(aero.get("L", 0.0))
            
            if abs(lift_N - required_lift_N) < tolerance:
                break
            
            # Adjust alpha (simple proportional control)
            if lift_N < required_lift_N:
                alpha_deg += 1.0
            else:
                alpha_deg -= 1.0
            
            # Clamp alpha to reasonable range
            alpha_deg = np.clip(alpha_deg, -10.0, 15.0)
        
        # Get drag at trimmed state
        drag_N = float(aero.get("D", 0.0))
        
        # Power from drag (straight flight)
        power_straight_W = drag_N * airspeed_mps
        
        # Apply turn load factor
        turn_radius_m = 0.5 * circular_diameter_m
        turn_power_W = calculate_turn_power_W(power_straight_W, airspeed_mps, turn_radius_m)
        
        return float(turn_power_W + avionics_power_W)
        
    except Exception as e:
        # Fallback on any error
        print(f"[warning] AeroSandbox weight sensitivity failed: {e}, using linear scaling")
        soln_path = artifact_path / "soln.json"
        if soln_path.exists():
            with soln_path.open("r") as f:
                soln = json.load(f)
            baseline_mass = float(soln.get("Performance", {}).get("total_mass", 9.0))
            baseline_power = float(soln.get("Performance", {}).get("power_cruise (all motors)", 50.0))
            mass_ratio = total_mass_kg / baseline_mass if baseline_mass > 0 else 1.0
            power_straight = baseline_power * (mass_ratio ** 1.5)
        else:
            mass_ratio = total_mass_kg / 9.0
            power_straight = 50.0 * (mass_ratio ** 1.5)
        
        turn_radius_m = 0.5 * circular_diameter_m
        turn_power = calculate_turn_power_W(power_straight, airspeed_mps, turn_radius_m)
        return float(turn_power + avionics_power_W)


# ----------------------------
# Sensitivity analysis computation
# ----------------------------

def compute_sensitivity_variants(
    baseline_results: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float],
    artifact_data: ArtifactData,
    artifact_path: Path,
    soln: Dict[str, Any],
    circular_diameter_m: float,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """
    Compute sensitivity analysis variants for all parameters.
    
    Returns dictionary with structure:
    {
      "weight": {"-0.2": {...}, "-0.1": {...}, ...},
      "drag": {...},
      ...
    }
    """
    (
        time_hours,
        solar_flux_Wm2,
        power_generated_W,
        power_used_W,
        net_energy_Wh,
        battery_state_Wh,
        period_type,
        battery_min_threshold_Wh,
        battery_capacity_Wh,
    ) = baseline_results
    
    # Get baseline mass
    baseline_mass_kg = float(soln.get("Performance", {}).get("total_mass", 9.0))
    
    # Estimate power limits from baseline
    limits = estimate_battery_power_limits(time_hours, battery_state_Wh)
    
    # Solar area
    solar_area_m2 = artifact_data.solar_panels_n * (artifact_data.solar_panel_side_length_m ** 2)
    
    sensitivity_data: Dict[str, Dict[str, Dict[str, Any]]] = {}
    
    # Weight sensitivity
    print("[sensitivity] Computing weight variants...")
    weight_variants: Dict[str, Dict[str, Any]] = {}
    for delta in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        delta_str = f"{delta:.1f}"
        mass_new = baseline_mass_kg * (1 + delta)
        power_new = compute_power_for_weight(
            artifact_path,
            mass_new,
            artifact_data.airspeed,
            artifact_data.operating_altitude,
            circular_diameter_m,
            artifact_data.avionics_power_W,
        )
        
        # Create power series (constant)
        power_used_variant = np.full_like(power_generated_W, power_new)
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_W - power_used_variant) * dt_hours
        
        # Re-integrate battery
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_W,
            power_used_variant,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        weight_variants[delta_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_W.tolist(),
            "power_used_W": power_used_variant.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["weight"] = weight_variants
    
    # Drag sensitivity (simple scaling)
    print("[sensitivity] Computing drag variants...")
    drag_variants: Dict[str, Dict[str, Any]] = {}
    for delta in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        delta_str = f"{delta:.1f}"
        power_used_variant = power_used_W * (1 + delta)
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_W - power_used_variant) * dt_hours
        
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_W,
            power_used_variant,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        drag_variants[delta_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_W.tolist(),
            "power_used_W": power_used_variant.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["drag"] = drag_variants
    
    # Solar (irradiance) sensitivity
    print("[sensitivity] Computing solar variants...")
    solar_variants: Dict[str, Dict[str, Any]] = {}
    for delta in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        delta_str = f"{delta:.1f}"
        power_generated_variant = power_generated_W * (1 + delta)
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_variant - power_used_W) * dt_hours
        
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_variant,
            power_used_W,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        solar_variants[delta_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_variant.tolist(),
            "power_used_W": power_used_W.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["solar"] = solar_variants
    
    # Cell efficiency sensitivity
    print("[sensitivity] Computing cell efficiency variants...")
    cell_eff_variants: Dict[str, Dict[str, Any]] = {}
    for delta in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        delta_str = f"{delta:.1f}"
        power_generated_variant = power_generated_W * (1 + delta)
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_variant - power_used_W) * dt_hours
        
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_variant,
            power_used_W,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        cell_eff_variants[delta_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_variant.tolist(),
            "power_used_W": power_used_W.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["cellEff"] = cell_eff_variants
    
    # Solar power (panel/area proxy) sensitivity
    print("[sensitivity] Computing solar power variants...")
    solar_power_variants: Dict[str, Dict[str, Any]] = {}
    for delta in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        delta_str = f"{delta:.1f}"
        power_generated_variant = power_generated_W * (1 + delta)
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_variant - power_used_W) * dt_hours
        
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_variant,
            power_used_W,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        solar_power_variants[delta_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_variant.tolist(),
            "power_used_W": power_used_W.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["solarPower"] = solar_power_variants
    
    # Motor efficiency sensitivity
    print("[sensitivity] Computing motor efficiency variants...")
    motor_eff_variants: Dict[str, Dict[str, Any]] = {}
    for delta in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        delta_str = f"{delta:.1f}"
        # Efficiency up => electrical power down. Use reciprocal scaling.
        scale = 1.0 / max(0.2, (1 + delta))  # prevent divide-by-near-zero
        power_used_variant = power_used_W * scale
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_W - power_used_variant) * dt_hours
        
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_W,
            power_used_variant,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        motor_eff_variants[delta_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_W.tolist(),
            "power_used_W": power_used_variant.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["motorEff"] = motor_eff_variants
    
    # CG shift sensitivity
    print("[sensitivity] Computing CG shift variants...")
    cg_variants: Dict[str, Dict[str, Any]] = {}
    for shift_cm in range(-40, 21, 5):  # -40 to 20 in steps of 5
        shift_str = str(shift_cm)
        # Penalty: 0.2% electrical power per cm of |shift|
        penalty_per_cm = 0.002
        scale = 1.0 + penalty_per_cm * abs(shift_cm)
        power_used_variant = power_used_W * scale
        dt_hours = (time_hours[1] - time_hours[0]) if time_hours.size > 1 else 0.0
        net_energy_variant = (power_generated_W - power_used_variant) * dt_hours
        
        battery_variant = integrate_battery_with_limits(
            time_hours,
            power_generated_W,
            power_used_variant,
            battery_capacity_Wh,
            battery_capacity_Wh,
            limits,
        )
        
        cg_variants[shift_str] = {
            "time_hr": time_hours.tolist(),
            "period_type": period_type.tolist(),
            "solar_flux_W_per_m2": solar_flux_Wm2.tolist(),
            "power_generated_W": power_generated_W.tolist(),
            "power_used_W": power_used_variant.tolist(),
            "battery_state_Wh": battery_variant.tolist(),
        }
    sensitivity_data["cg"] = cg_variants
    
    print("[sensitivity] Sensitivity analysis complete")
    return sensitivity_data


# ----------------------------
# Plotting / export
# ----------------------------

def _add_period_shading(ax, time_hours, period_type, y_min, y_max):
    # Masks
    night = period_type == "night"
    dawn = period_type == "dawn_dusk"
    day = period_type == "day"

    # Night
    if np.any(night):
        ax.fill_between(time_hours, y_min, y_max, where=night, alpha=0.20, color="#000080",
                        label="Night", interpolate=True, zorder=0)

    # Dawn/Dusk
    if np.any(dawn):
        ax.fill_between(time_hours, y_min, y_max, where=dawn, alpha=0.18, color="#FFA500",
                        label="Dawn/Dusk", interpolate=True, zorder=0)

    # Day (no label to avoid legend clutter)
    if np.any(day):
        ax.fill_between(time_hours, y_min, y_max, where=day, alpha=0.15, color="#FFFACD",
                        label="", interpolate=True, zorder=0)

    # 24h marker
    ax.axvline(x=24, color="gray", linestyle="--", alpha=0.5, linewidth=1, zorder=10)
    ax.text(
        24, y_max * 0.98, "Day 2", rotation=90, va="top", ha="right", fontsize=9, alpha=0.7, zorder=20,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
    )


def plot_results(
    time_hours: np.ndarray,
    solar_flux_Wm2: np.ndarray,
    power_generated_W: np.ndarray,
    power_used_W: np.ndarray,
    battery_state_Wh: np.ndarray,
    period_type: np.ndarray,
    battery_min_threshold_Wh: float,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(16, 8))
    fig.suptitle("48-Hour (2-Day) Energy Balance Analysis", fontsize=16, fontweight="bold")

    # y-ranges (for shading)
    battery_y_min, battery_y_max = 0.0, float(np.max(battery_state_Wh)) * 1.05
    power_y_min = float(min(np.min(power_generated_W), np.min(power_used_W))) * 0.95
    power_y_max = float(max(np.max(power_generated_W), np.max(power_used_W))) * 1.05
    flux_y_min, flux_y_max = 0.0, float(np.max(solar_flux_Wm2)) * 1.05

    overall_y_min = min(battery_y_min, power_y_min, flux_y_min)
    overall_y_max = max(battery_y_max, power_y_max, flux_y_max)

    _add_period_shading(ax, time_hours, period_type, overall_y_min, overall_y_max)

    # Battery-bad region
    below_min = battery_state_Wh < battery_min_threshold_Wh
    if np.any(below_min):
        ax.fill_between(
            time_hours, overall_y_min, overall_y_max,
            where=below_min, alpha=0.25, color="red",
            label="BELOW MINIMUM (unsafe)", interpolate=True, zorder=1
        )

    # Battery on left axis
    ax.plot(time_hours, battery_state_Wh, "b-", linewidth=2.5, label="Battery Energy (Wh)", zorder=6)
    ax.axhline(battery_min_threshold_Wh, color="r", linestyle="--", alpha=0.5, linewidth=1.5,
               label="Minimum Safe Level", zorder=6)
    ax.set_ylabel("Battery Energy (Wh)", color="blue", fontsize=12, fontweight="bold")
    ax.tick_params(axis="y", labelcolor="blue")
    ax.set_ylim(battery_y_min, battery_y_max)

    # Power + flux on right axis
    ax2 = ax.twinx()
    ax2.plot(time_hours, power_generated_W, "g-", linewidth=2, label="Power Generated (W)", zorder=5)
    ax2.plot(time_hours, power_used_W, "r-", linewidth=2, label="Power Used (W)", zorder=5)

    surplus = power_generated_W >= power_used_W
    deficit = ~surplus
    if np.any(surplus):
        ax2.fill_between(time_hours, power_generated_W, power_used_W, where=surplus, alpha=0.2, color="green",
                         label="Surplus", zorder=3)
    if np.any(deficit):
        ax2.fill_between(time_hours, power_generated_W, power_used_W, where=deficit, alpha=0.2, color="red",
                         label="Deficit", zorder=3)

    ax2.plot(time_hours, solar_flux_Wm2, "orange", linewidth=2, linestyle="--",
             label="Solar Flux (W/m²)", zorder=5)

    ax2.set_ylabel("Power (W) / Solar Flux (W/m²)", color="green", fontsize=12, fontweight="bold")
    ax2.tick_params(axis="y", labelcolor="green")
    ax2.set_ylim(0.0, max(power_y_max, flux_y_max))

    ax.set_xlabel("Time (hours)", fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3, zorder=1)

    # Legend: merge handles and drop empty labels
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    handles, labels = [], []
    for h, l in list(zip(h1, l1)) + list(zip(h2, l2)):
        if l and l not in labels:
            handles.append(h)
            labels.append(l)
    ax.legend(handles, labels, loc="upper left", fontsize=9)

    plt.tight_layout()
    plot_path = output_dir / "energy_balance.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"[plot] Saved plot to: {plot_path}")
    plt.show()


def export_data(
    time_hours: np.ndarray,
    period_type: np.ndarray,
    solar_flux_Wm2: np.ndarray,
    power_generated_W: np.ndarray,
    power_used_W: np.ndarray,
    net_energy_Wh: np.ndarray,
    battery_state_Wh: np.ndarray,
    artifact_data: ArtifactData,
    circular_diameter_m: float,
    battery_capacity_Wh: float,
    output_dir: Path,
    sensitivity_data: Optional[Dict[str, Any]] = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # Summary metrics (trapezoidal integration over hours gives Wh)
    # Use trapezoid (numpy 1.20+) or trapz (older versions)
    if hasattr(np, 'trapezoid'):
        trapz_func = np.trapezoid
    elif hasattr(np, 'trapz'):
        trapz_func = np.trapz
    else:
        # Fallback: manual trapezoidal integration
        def trapz_func(y, x):
            return float(np.sum((y[1:] + y[:-1]) / 2.0 * np.diff(x)))
    
    total_gen_Wh = float(trapz_func(power_generated_W, time_hours))
    total_use_Wh = float(trapz_func(power_used_W, time_hours))
    surplus_Wh = total_gen_Wh - total_use_Wh

    day1 = time_hours <= 24.0
    day2 = time_hours > 24.0

    day1_gen = float(trapz_func(power_generated_W[day1], time_hours[day1]))
    day1_use = float(trapz_func(power_used_W[day1], time_hours[day1]))
    day2_gen = float(trapz_func(power_generated_W[day2], time_hours[day2]))
    day2_use = float(trapz_func(power_used_W[day2], time_hours[day2]))

    # Period stats
    dt_hours = float(time_hours[1] - time_hours[0])
    is_day = period_type == "day"
    is_dawn = period_type == "dawn_dusk"
    is_night = period_type == "night"

    def hours(mask: np.ndarray) -> float:
        return float(np.sum(mask) * dt_hours)

    battery_at_24h = float(battery_state_Wh[np.argmin(np.abs(time_hours - 24.0))])

    summary = {
        "analysis_parameters": {
            "circular_diameter_m": float(circular_diameter_m),
            "turn_radius_m": float(circular_diameter_m / 2.0),
            "n_points": int(time_hours.size),
            "simulation_duration_hours": float(time_hours[-1]),
        },
        "aircraft_parameters": artifact_data.__dict__,
        "period_statistics": {
            "total_daytime_hours": hours(is_day),
            "total_dawn_dusk_hours": hours(is_dawn),
            "total_nighttime_hours": hours(is_night),
            "day1_daytime_hours": hours(is_day & day1),
            "day1_dawn_dusk_hours": hours(is_dawn & day1),
            "day1_nighttime_hours": hours(is_night & day1),
            "day2_daytime_hours": hours(is_day & day2),
            "day2_dawn_dusk_hours": hours(is_dawn & day2),
            "day2_nighttime_hours": hours(is_night & day2),
        },
        "energy_metrics": {
            "total_energy_generated_Wh": total_gen_Wh,
            "total_energy_used_Wh": total_use_Wh,
            "energy_surplus_Wh": float(surplus_Wh),
            "energy_surplus_percent": float(100.0 * surplus_Wh / total_use_Wh) if total_use_Wh > 0 else 0.0,
        },
        "day1_metrics": {
            "energy_generated_Wh": day1_gen,
            "energy_used_Wh": day1_use,
            "energy_surplus_Wh": float(day1_gen - day1_use),
            "energy_surplus_percent": float(100.0 * (day1_gen - day1_use) / day1_use) if day1_use > 0 else 0.0,
        },
        "day2_metrics": {
            "energy_generated_Wh": day2_gen,
            "energy_used_Wh": day2_use,
            "energy_surplus_Wh": float(day2_gen - day2_use),
            "energy_surplus_percent": float(100.0 * (day2_gen - day2_use) / day2_use) if day2_use > 0 else 0.0,
        },
        "battery_metrics": {
            "battery_capacity_Wh": float(battery_capacity_Wh),
            "initial_state_Wh": float(battery_state_Wh[0]),
            "final_state_Wh": float(battery_state_Wh[-1]),
            "min_state_Wh": float(np.min(battery_state_Wh)),
            "max_state_Wh": float(np.max(battery_state_Wh)),
            "state_change_Wh": float(battery_state_Wh[-1] - battery_state_Wh[0]),
            "battery_at_24h_Wh": battery_at_24h,
        },
        "power_metrics": {
            "avg_power_generated_W": float(np.mean(power_generated_W)),
            "max_power_generated_W": float(np.max(power_generated_W)),
            "avg_power_used_W": float(np.mean(power_used_W)),
            "max_power_used_W": float(np.max(power_used_W)),
        },
    }

    json_path = output_dir / "energy_balance_summary.json"
    if sensitivity_data is not None:
        summary["sensitivity_analysis"] = sensitivity_data
    with json_path.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"[export] Saved JSON summary to: {json_path}")

    # Console summary
    print("\n" + "=" * 70)
    print("48-HOUR (2-DAY) ENERGY BALANCE SUMMARY")
    print("=" * 70)
    print(f"Total Energy Generated: {total_gen_Wh:.2f} Wh")
    print(f"Total Energy Used:      {total_use_Wh:.2f} Wh")
    print(f"Energy Surplus:         {surplus_Wh:.2f} Wh ({summary['energy_metrics']['energy_surplus_percent']:.1f}%)")
    print(f"Battery: initial={battery_state_Wh[0]:.2f} Wh, at24h={battery_at_24h:.2f} Wh, "
          f"final={battery_state_Wh[-1]:.2f} Wh, min={np.min(battery_state_Wh):.2f} Wh")
    print("=" * 70)


# ----------------------------
# CLI
# ----------------------------

def main() -> None:
    p = argparse.ArgumentParser(description="48-hour energy balance for solar aircraft in circular flight")
    p.add_argument("--artifact-path", type=str, default=None,
                   help="Path to artifact folder containing soln.json (default: latest run in output/)")
    p.add_argument("--circular-diameter", type=float, default=1000.0,
                   help="Diameter of circular flight pattern in meters (default: 1000)")
    p.add_argument("--n-points", type=int, default=3000,
                   help="Number of time discretization points (default: 3000 for 2 days)")
    args = p.parse_args()

    if args.artifact_path:
        artifact_path = Path(args.artifact_path)
    else:
        artifact_root = Path(__file__).parent / "output"
        artifact_path = latest_run_dir(artifact_root)
        print(f"[config] Using latest run: {artifact_path}")

    if not artifact_path.exists():
        raise FileNotFoundError(f"Artifact path does not exist: {artifact_path}")

    # Output directory is the same as artifact directory
    output_dir = artifact_path

    print(f"[load] Loading artifact data from: {artifact_path}")
    a = load_artifact_data(artifact_path)
    
    # Load soln.json for sensitivity analysis
    soln_path = artifact_path / "soln.json"
    with soln_path.open("r") as f:
        soln = json.load(f)

    print(f"[compute] Computing 48h energy balance (diameter={args.circular_diameter} m, n_points={args.n_points})")
    (
        time_hours,
        solar_flux,
        power_generated,
        power_used,
        net_energy,
        battery_state,
        period_type,
        battery_min_threshold,
        battery_capacity,
    ) = compute_energy_balance(a, args.circular_diameter, n_points=args.n_points)
    
    # Compute sensitivity analysis
    print("[sensitivity] Computing sensitivity variants...")
    baseline_results = (
        time_hours,
        solar_flux,
        power_generated,
        power_used,
        net_energy,
        battery_state,
        period_type,
        battery_min_threshold,
        battery_capacity,
    )
    sensitivity_data = compute_sensitivity_variants(
        baseline_results,
        a,
        artifact_path,
        soln,
        args.circular_diameter,
    )

    print("[plot] Generating plot...")
    plot_results(
        time_hours=time_hours,
        solar_flux_Wm2=solar_flux,
        power_generated_W=power_generated,
        power_used_W=power_used,
        battery_state_Wh=battery_state,
        period_type=period_type,
        battery_min_threshold_Wh=battery_min_threshold,
        output_dir=output_dir,
    )

    print("[export] Exporting CSV and JSON...")
    export_data(
        time_hours=time_hours,
        period_type=period_type,
        solar_flux_Wm2=solar_flux,
        power_generated_W=power_generated,
        power_used_W=power_used,
        net_energy_Wh=net_energy,
        battery_state_Wh=battery_state,
        artifact_data=a,
        circular_diameter_m=args.circular_diameter,
        battery_capacity_Wh=battery_capacity,
        output_dir=output_dir,
        sensitivity_data=sensitivity_data,
    )

    print(f"\n[complete] Analysis complete. Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
