import os
import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.library import power_solar, propulsion_electric, propulsion_propeller
from aerosandbox.atmosphere.atmosphere import Atmosphere
import numpy as onp
from pathlib import Path
from datetime import datetime, timezone

from lib.artifacts import process_raw_values, run_id_random, write_json
from lib.exports import export_xflr5_xml_from_soln
from lib.models import hacker_motor_resistance, hacker_motor_i0



opti = asb.Opti()
opti.solver("ipopt")



### CONSTANTS
## Physics
g = 9.81

## Mission
mission_date = 100
operating_lat = 37.398928
togw_max = 12 # kg
temperature_high = 34 # 34°K hotter than ISA --> leads to ~93 °F operating temperature
operating_altitude = 1200 # in meters
operating_atm = Atmosphere(operating_altitude, temperature_deviation=temperature_high)

## Performance
min_climb_angle = 30 # degrees

## Aerodynamics
# Airfoils
# SPLIT WING (rev7): different airfoil inboard (root) vs outboard (tip), selectable via env.
# Inboard carries the fuselages/batteries/spar loads (thicker, higher-camber section);
# outboard is the tapered efficient panel. AeroSandbox blends between the sections.
root_airfoil = asb.Airfoil(os.environ.get("WING_AIRFOIL_ROOT", "s4310"))
tip_airfoil  = asb.Airfoil(os.environ.get("WING_AIRFOIL_TIP", "s4110"))
wing_airfoil = root_airfoil  # inboard section; used for battery-thickness references
tail_airfoil = asb.Airfoil("naca0010")
# Stall CLmax = the more limiting (lower-CLmax) of the two sections, derated 0.9.
def _clmax(af):
    p = af.get_aero_from_neuralfoil(alpha=onp.linspace(0, 16, 33), Re=200000, mach=0.03)
    return float(onp.max(p["CL"]))
CLmax_cruise = 0.9 * min(_clmax(root_airfoil), _clmax(tip_airfoil))
print(f"[airfoil] root={root_airfoil.name} tip={tip_airfoil.name}: CLmax_cruise (derated) = {CLmax_cruise:.3f}")

# Tail sizing via tail-volume coefficients. V_H fixed; V_V is a design variable
# (grown to meet the weathercock floor below).
V_H = 0.40  # horizontal tail volume coefficient

# Static margin, SM = (x_np - x_cg)/MAC. SM_max is generous because the dihedral +
# larger fins push the neutral point aft and we still need room to trim Cm=0.
SM_min = 0.05
SM_max = 0.70

## Lateral-directional stability. Constrained via the smooth geometry levers below;
## the AeroBuildup per-rad derivatives are also computed and reported as a cross-check.
V_V_floor = 0.015    # min vertical tail volume coeff (weathercock floor)
dihedral_min = 2.0   # deg, min wing dihedral (keeps the spiral mode benign)
Cnb_min = 0.03       # /rad, reference weathercock target (report only)
Clb_max = -0.015     # /rad, reference dihedral-effect target (report only)

## Aileron / roll authority (sized so the aircraft actually rolls at cruise).
roll_rate_target = np.radians(30)         # rad/s steady roll rate at full aileron
aileron_deflection_max = np.radians(20)   # rad, max aileron deflection for sizing
aileron_y_start_frac = 0.55               # aileron inboard edge, frac of semi-span
aileron_y_end_frac = 0.95                 # aileron outboard edge, frac of semi-span
section_lift_slope = 2 * np.pi            # a0 [1/rad], thin-airfoil lift-curve slope

# Stall (AeroBuildup doesn't model stall; CLmax_cruise derived from the airfoil above).
stall_speed_margin = 1.2  # V_cruise >= margin * V_stall

## Structural
structural_mass_markup = 1.2
boom_radius = 0.01 # radius of the boom in meters.
boom_spacing_frac = 0.25 # Fraction that the inboard wing spans.
fuselage_length = 0.6 # m, length of each pod fuselage (rev7: was 0.5 m).

## Propulsion
propeller_n = 2 # Number of motors / propellers.
# (advance ratio is now a design variable in the PROPULSION section, with an eta_prop(J) map)

## Power
battery_voltage = 22.2
N = 180 # Number of discretization points.
time = onp.linspace(0, 24 * 60 * 60, N) # s (numeric; does not need to be symbolic)
dt = onp.diff(time)[0]  # s
solar_panel_side_length = 0.125 # m (12.5cm)
solar_panel_n_rows = 2
# Solar cell model: SunPower C60 (bin I) + full loss stack.
# System eff = cell(0.223) x encapsulation(0.97) x wiring/mismatch/soiling(0.965) x MPPT(0.95).
# The temperature-derate block below scales this by a per-timestep factor.
cell_efficiency_stc = 0.223         # SunPower C60 bin I (datasheet, STC)
encapsulation_transmission = 0.97   # thin fiberglass
wiring_mismatch_soiling_eff = 0.965 # combined array losses
mppt_efficiency = 0.95              # Genasun boost MPPT
# STC system efficiency (before temperature). This is the base the derate block scales.
solar_cell_efficiency = cell_efficiency_stc * encapsulation_transmission * wiring_mismatch_soiling_eff * mppt_efficiency
energy_generation_margin = 1.0      # losses now explicit in the chain above (no extra fudge)
# Battery: Upgrade Energy GOLD V1 6S2P, Amprius SA03 silicon-anode (verified datasheet).
# 483 Wh, 1.278 kg -> 378 Wh/kg NAMEPLATE. Cycle the 20-100% SOC window (80% usable).
allowable_battery_depth_of_discharge = 0.80  # 20-100% SOC window
pack_energy_Wh = 483.0                        # Wh per Amprius SA03 pack (nameplate)
pack_mass_kg = 1.278                          # kg per pack -> 378 Wh/kg
battery_charge_c_rate = 0.25                  # 6 A / pack ~ 0.25C hard charge limit (datasheet)
# Usable-energy derate (Amprius Si-anode). The 450->360 Wh/kg rate-derate is a 1C figure
# and doesn't apply here (~0.06C discharge); 8% covers hot-cell loss + datasheet optimism.
# Cycle life (~125-150 cycles to 80% SOH) is a campaign risk -> budget fresh packs. FLAGGED.
battery_usable_derate = 0.92
# Optional: pin the pack count to an integer for a discrete-reality solve (env FIXED_PACKS).
fixed_packs = os.environ.get("FIXED_PACKS")

# Precompute solar flux profile numerically (robust vs NaNs from symbolic trig/acos edge cases)
solar_flux_profile = onp.array(
    [
        power_solar.solar_flux(
            latitude=operating_lat,
            day_of_year=mission_date,
            time=float(t),
            altitude=operating_altitude,
            panel_azimuth_angle=0,
            panel_tilt_angle=0,
        )
        for t in time
    ],
    dtype=float,
)
solar_flux_profile = onp.nan_to_num(solar_flux_profile, nan=0.0, posinf=0.0, neginf=0.0)
solar_flux_profile = onp.maximum(solar_flux_profile, 0.0)

# FIX: Incidence-Angle Modifier (IAM) reflection loss. The geometric sun-angle projection
# (cos of zenith) is already in solar_flux; IAM is the *additional* surface-reflection loss
# that grows at low sun angles (morning/evening). ASHRAE model on a horizontal panel where
# incidence angle = zenith = 90 - solar_elevation. This is what stretches the battery-night.
_elev = onp.array([power_solar.solar_elevation_angle(operating_lat, mission_date, float(t)) for t in time], dtype=float)
_aoi = onp.clip(90.0 - _elev, 0.0, 89.0)                      # incidence angle on horizontal panel, deg
_iam_b0 = 0.05                                                # ASHRAE reflection coefficient (glass/encapsulant)
iam_profile = onp.clip(1.0 - _iam_b0 * (1.0 / onp.cos(onp.radians(_aoi)) - 1.0), 0.0, 1.0)
iam_profile = onp.where(_elev > 0, iam_profile, 0.0)          # no beam when sun is down
solar_flux_profile = solar_flux_profile * iam_profile         # apply reflection loss
print(
    f"[solar_flux_profile] min={solar_flux_profile.min():.3f} W/m^2, "
    f"max={solar_flux_profile.max():.3f} W/m^2, "
    f"nan_count={onp.isnan(solar_flux_profile).sum()}"
)
print(f"[IAM] noon iam={iam_profile.max():.3f}, mean(daylight)={iam_profile[_elev>0].mean():.3f}")

# Cell-temperature derate. A solar cell in the desert sun is far hotter
# than the air: only ~22% of absorbed sunlight leaves as electricity, the rest is
# heat, so the cell self-heats well above ambient. NOCT-style model:
#     T_cell = T_amb + k * flux,   k calibrated so peak-flux (noon) cell = 75 C.
# Output derates by the datasheet power coefficient gamma (per C above 25 C ref).
# This is deliberately conservative (static, no flight-airflow cooling credited).
cell_gamma_pmp = -0.0032   # 1/C, C60/Maxeon power temp coefficient (datasheet). FLAGGED
cell_T_ref_C = 25.0        # rating reference temperature
cell_T_noon_C = 75.0       # user: peak static cell temp at solstice noon. FLAGGED
# Design-day ambient sinusoid (min at 05:00, peak at 15:30). Cleanup: the afternoon PEAK is
# tied to operating_atm (the SAME hot design-day temperature used for aero density), so there
# is one temperature source, not two. Amplitude is the desert day/night swing. FLAGGED.
_amb_amp = 8.0
_amb_peak_c = operating_atm.temperature() - 273.15    # hot design-day peak (shared with density)
_amb_mean = _amb_peak_c - _amb_amp
_hours = time / 3600.0
T_amb_profile = _amb_mean + _amb_amp * onp.cos(2 * onp.pi * (_hours - 15.5) / 24.0)
_i_noon = int(onp.argmax(solar_flux_profile))
_k_cell = (cell_T_noon_C - T_amb_profile[_i_noon]) / max(solar_flux_profile[_i_noon], 1e-6)
T_cell_profile = T_amb_profile + _k_cell * solar_flux_profile
# Per-timestep effective cell efficiency (never below 0).
solar_cell_efficiency_profile = onp.maximum(
    solar_cell_efficiency * (1 + cell_gamma_pmp * (T_cell_profile - cell_T_ref_C)),
    0.0,
)
print(
    f"[cell_temp] noon T_cell={T_cell_profile[_i_noon]:.1f} C, "
    f"noon derate={1 + cell_gamma_pmp*(T_cell_profile[_i_noon]-cell_T_ref_C):.3f}, "
    f"eff {solar_cell_efficiency:.3f} -> {solar_cell_efficiency_profile[_i_noon]:.3f} at noon"
)



### VARIABLES
## Performance
airspeed = opti.variable(init_guess=10.7, lower_bound=5, upper_bound=15, scale=5, category="airspeed")
togw_design = opti.variable(init_guess=10.5, lower_bound=1, upper_bound=togw_max, scale=1, category="togw_max")
power_out_max = opti.variable(init_guess=500, lower_bound=25*16, scale=100, category="power_out_max")

## Propulsion
propeller_diameter = opti.parameter(value=0.4)  # m (fixed, ~16 in)
motor_kv = opti.variable(init_guess=500, lower_bound=50, upper_bound=2000, scale=100, category="motor_kv")

## Avionics
solar_panel_n = opti.variable(init_guess=88, lower_bound=10, scale=10, category="solar_panel_n")
battery_capacity = opti.variable(init_guess=1450, lower_bound=100, scale=100, category="battery_capacity") # NAMEPLATE installed energy in Wh (drives mass + pack count)
# Deliverable energy after the usable-energy derate. SOC dynamics/bounds use this; mass and
# num_packs stay on nameplate (you carry the full pack mass but only get derate x energy).
battery_capacity_effective = battery_capacity * battery_usable_derate
battery_states = opti.variable(n_vars=N, init_guess=850, scale=100, category="battery_states")
# Seed battery_states with a realistic day/night SOC sawtooth (numeric march at nominal loads).
# A flat 850 Wh seed violates every one of the 179 recursion equalities at once, which leaves
# IPOPT stuck at inf_pr ~ 12 Wh; a physically-shaped seed lets it close feasibility.
_cap_nom = 1450.0 * battery_usable_derate
_area_nom = 88.0 * (0.125 ** 2)
_gen_nom = solar_flux_profile * _area_nom * solar_cell_efficiency_profile  # W
_soc_seed = onp.empty(N)
_soc_seed[0] = 0.85 * _cap_nom
for _i in range(N - 1):
    _net = (_gen_nom[_i] - 45.0) * dt / 3600.0  # 45 W nominal total load
    _soc_seed[_i + 1] = onp.clip(_soc_seed[_i] + _net, 0.2 * _cap_nom, _cap_nom)
opti.set_initial(battery_states, _soc_seed)

## Aerodynamics
# Main wing
wingspan = opti.variable(init_guess=5.9, lower_bound=2, upper_bound=9.0, scale=1, category="wingspan")  # rev7: cap 9 m (heavier realistic tails need a bit more span at 8 packs)
chordlen = opti.variable(init_guess=0.45, lower_bound=0.05, scale=0.1, category="chordlen")
struct_defined_aoa = opti.variable(init_guess=4.7, lower_bound=0, upper_bound=10, scale=1, category="struct_aoa")
# Real CG (rev7): cg_le_dist is now a VARIABLE tied by a fixed-point constraint to the actual
# mass-weighted CG of all components (computed after the mass rollup). The battery position is
# free so the optimizer can place it to balance the aircraft, exactly like a real build.
cg_le_dist = opti.variable(init_guess=0.08, lower_bound=-0.20, upper_bound=0.50, scale=0.1, category="cg_le_dist")
battery_x = opti.variable(init_guess=-0.20, lower_bound=-0.55, upper_bound=0.30, scale=0.1, category="battery_x")

# Wing dihedral (from the root) and taper (flat LE; TE tapers forward, chord shrinks
# from chordlen at the root to taper_ratio*chordlen at the tip).
dihedral_angle = opti.variable(init_guess=3.0, lower_bound=0, upper_bound=8, scale=1, category="dihedral")
taper_ratio = opti.variable(init_guess=0.8, lower_bound=0.4, upper_bound=1.0, scale=0.1, category="taper_ratio")

# Empennage. hstab_AR is NOT here: the H-stab planform is fully set by hstab_span (= boom
# spacing) and hstab_chordlen (= area/span), so its AR is derived after the geometry below.
# The vertical tail is grown by the optimizer to meet the weathercock floor.
vstab_AR = opti.variable(init_guess=2.5, lower_bound=2.0, upper_bound=4.0, scale=1, category="vstab_AR")  # tall+slim fin
V_V = opti.variable(init_guess=0.016, lower_bound=0.01, upper_bound=0.08, scale=0.01, category="V_V")
hstab_aoa = opti.variable(init_guess=-4, lower_bound=-10, upper_bound=5, scale=1, category="hstab_aoa")

# Control-surface chord fractions.
aileron_chord_ratio = opti.variable(init_guess=0.35, lower_bound=0.15, upper_bound=0.40, scale=0.1, category="aileron_cf_c")
elevator_chord_ratio = opti.variable(init_guess=0.30, lower_bound=0.15, upper_bound=0.50, scale=0.1, category="elevator_ce_c")

# Structural
boom_length = opti.variable(init_guess=1.0, lower_bound=0.1, upper_bound=1.5, scale=1, category="boom_length")  # m, capped at 1.5
boom_y = 0.5 * boom_spacing_frac * wingspan



### GEOMETRIES
# Wing planform (rev7): inboard section (root -> boom_y, between the fuselages) is full chord
# and flat. Taper AND polyhedral both start at the inboard break (boom_y); outboard the chord
# tapers to the tip and rises with dihedral. How many solar panels actually fit on the tapered
# wing is handled honestly by the ROW-AWARE panel model in the constraints (panels need chord
# >= 1 pitch for a row, >= 2 pitch for two rows, etc.) -- so taper trades against panel rows.
semi_span = wingspan / 2
_panel_pitch = 0.13                              # m, 12.5 cm cell + gap (matches panel model)
_y_break = boom_y                                # taper + polyhedral break (inboard section end)
_y_mid = boom_y + wingspan * boom_spacing_frac / 2
_y_tip = wingspan / 2
_outboard_span = semi_span - _y_break            # outboard half-span (break -> tip)

def chord_at(y):
    # Full chord (chordlen) inboard of the break; linear taper outboard. np.softmax clamps the
    # inboard side to 0 so this is valid across the whole span (used by the row-aware model).
    return chordlen * (1 - (1 - taper_ratio) * np.softmax((y - _y_break) / _outboard_span, 0, hardness=30))

def z_at(y):
    # Flat inboard of the break; polyhedral rise measured from the break outboard.
    return np.sind(dihedral_angle) * np.softmax(y - _y_break, 0, hardness=50)

chord_root = chordlen
chord_boom = chordlen                             # inboard is full chord
chord_mid = chord_at(_y_mid)
chord_tip = chordlen * taper_ratio

main_wing = asb.Wing(
    name="Main Wing",
    symmetric=True,  # Should this wing be mirrored across the XZ plane?
    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            xyz_le=[0, 0, 0], # leading edge, relative to the wing's leading edge.
            chord=chord_root,
            twist=struct_defined_aoa,
            airfoil=root_airfoil,
        ),
        asb.WingXSec(  # Inboard break (full chord, flat -> polyhedral & taper start here)
            xyz_le=[0.00, _y_break, 0.0],
            chord=chord_boom,
            twist=struct_defined_aoa,
            airfoil=root_airfoil,
        ),
        asb.WingXSec(  # Outboard mid (tapered + polyhedral) -> blends toward tip airfoil
            xyz_le=[0.00, _y_mid, z_at(_y_mid)],
            chord=chord_mid,
            twist=struct_defined_aoa,
            airfoil=tip_airfoil,
        ),
        asb.WingXSec(  # Tip
            xyz_le=[0.00, _y_tip, z_at(_y_tip)],
            chord=chord_tip,
            twist=struct_defined_aoa,
            airfoil=tip_airfoil,
        ),
    ],
)

# Tail sizing from volume coefficients + aspect ratio (rectangular planforms)
S_w = main_wing.area()
MAC_w = main_wing.mean_aerodynamic_chord()
b_w = main_wing.span()

# Horizontal tail
L_H = boom_length - cg_le_dist
hstab_area = V_H * S_w * MAC_w / L_H
hstab_span = boom_y * 2
hstab_chordlen = hstab_area / hstab_span
hstab_AR = hstab_span / hstab_chordlen   # DERIVED real aspect ratio (rectangular H-stab)

# Vertical tail
L_V = boom_length - cg_le_dist
vstab_area_total = V_V * S_w * b_w / L_V
vstab_area = 0.5 * vstab_area_total
vstab_span = np.sqrt(vstab_AR * vstab_area)  # "span" here is height
vstab_root_chord = (2 * vstab_area / vstab_span) - hstab_chordlen  # trapezoid area: S = b*(cr+ct)/2

# Model horizontal and vertical stabilizers.
hor_stabilizer = asb.Wing(
    name="Horizontal Stabilizer",
    symmetric=True,
    xsecs=[
        asb.WingXSec(  # root
            xyz_le=[0, 0, 0],
            chord=hstab_chordlen,
            twist=hstab_aoa,
            airfoil=tail_airfoil,
        ),
        asb.WingXSec(  # tip
            xyz_le=[0.0, hstab_span / 2, 0],
            chord=hstab_chordlen,
            twist=hstab_aoa,
            airfoil=tail_airfoil,
        ),
    ],
).translate([boom_length + vstab_root_chord/4, 0, vstab_span])

# Both fins are identical trapezoids mirrored across XZ — build once, place by sign.
def _make_fin(name, sign):
    return asb.Wing(
        name=name,
        symmetric=False,
        xsecs=[
            asb.WingXSec(xyz_le=[0, 0, 0], chord=vstab_root_chord, twist=0, airfoil=tail_airfoil),
            asb.WingXSec(xyz_le=[vstab_root_chord/4, 0, vstab_span], chord=hstab_chordlen, twist=0, airfoil=tail_airfoil),
        ],
    ).translate([boom_length, sign * boom_y, 0])

vert_stabilizer_L = _make_fin("Vertical Stabilizer L", +1)
vert_stabilizer_R = _make_fin("Vertical Stabilizer R", -1)

# Booms and pods are L/R mirror pairs — build each once, place by sign.
_dae51 = asb.Airfoil("dae51")   # hoisted: one instance, not one per section


def _make_boom(name, sign):
    return asb.Fuselage(
        name=name,
        xsecs=[
            asb.FuselageXSec(xyz_c=[boom_length * xi, sign * boom_y, -0.02], radius=boom_radius)
            for xi in np.cosspace(0, 1, 20)
        ],
    )


def _make_pod(name, sign):
    # pod fuselage of length `fuselage_length`, starting at the LE going forward
    return asb.Fuselage(
        name=name,
        xsecs=[
            asb.FuselageXSec(
                xyz_c=[fuselage_length * xi - fuselage_length, sign * boom_y, -0.02],
                radius=0.4 * _dae51.local_thickness(x_over_c=xi),
            )
            for xi in np.cosspace(0, 1, 30)
        ],
    )


boom_L = _make_boom("Boom L", +1)
boom_R = _make_boom("Boom R", -1)
left_pod = _make_pod("Left Fuse", +1)
right_pod = _make_pod("Right Fuselage", -1)

# Model the airplane.
airplane = asb.Airplane(
    name="Rev 7",
    xyz_ref=[cg_le_dist, 0, 0],  # CG location (measured from main wing LE)
    wings=[main_wing, hor_stabilizer, vert_stabilizer_L, vert_stabilizer_R],
    fuselages=[left_pod, right_pod, boom_L, boom_R],
)



### AERODYNAMICS
vlm = asb.AeroBuildup(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        atmosphere=operating_atm,
        velocity=airspeed,  # m/s
    ),
)
aero = vlm.run_with_stability_derivatives(
    alpha=True,
    beta=True,
    p=True,   # rev7: needed for roll damping Clp (roll-rate authority, note #1)
    q=True,   # pitch damping Cmq -- exported for the JSBSim flight model (AircraftSim)
    r=True,   # rev7: needed for spiral-mode diagnostic (Cnr, Clr)
)

### DRAG BREAKDOWN
# AeroBuildup returns a per-component drag buildup: every wing/fuselage carries a
# PROFILE (skin-friction + form) drag, and the whole aircraft carries a single lift-
# dependent INDUCED drag. We tally these into a table for the report + web display.
# Component order matches the airplane definition:
#   wings     = [main_wing, hor_stabilizer, vert_stabilizer_L, vert_stabilizer_R]
#   fuselages = [left_pod, right_pod, boom_L, boom_R]
_drag_comp_names = [
    "main_wing", "hstab", "vstab_L", "vstab_R",   # wings
    "pod_L", "pod_R", "boom_L", "boom_R",         # fuselages
]
_drag_comps = [*aero["wing_aero_components"], *aero["fuselage_aero_components"]]
assert len(_drag_comps) == len(_drag_comp_names), "drag component names out of sync with airplane"
# qS from the aggregate so CD contributions use AeroBuildup's own reference area.
_qS = aero["D"] / aero["CD"]
drag_breakdown = {
    "components": [
        {"name": name, "D_profile": comp.D, "CD_profile": comp.D / _qS}
        for name, comp in zip(_drag_comp_names, _drag_comps)
    ],
    "D_profile_total": aero["D_profile"],
    "D_induced": aero["D_induced"],
    "D_total": aero["D"],
    "CD_profile_total": aero["D_profile"] / _qS,
    "CD_induced": aero["D_induced"] / _qS,
    "CD_total": aero["CD"],
}



### STABILITY
# --- Longitudinal (unchanged from rev6) ---
static_margin = (aero["x_np"] - cg_le_dist) / main_wing.mean_aerodynamic_chord()

# --- Lateral-directional (rev7) ---
# AeroBuildup returns per-radian derivatives, body axes.
Cnb = aero["Cnb"]  # yaw stiffness / weathercock (>0 stable)
Clb = aero["Clb"]  # dihedral effect (<0 stable)
Clp = aero["Clp"]  # roll damping (<0)
Cnr = aero["Cnr"]  # yaw damping (<0)
Clr = aero["Clr"]  # roll due to yaw rate
Cnp = aero["Cnp"]  # yaw due to roll rate (adverse/proverse yaw)
Cmq = aero["Cmq"]  # pitch damping (<0) -- from q=True above, for the JSBSim model

# Spiral-mode diagnostic (reported, not constrained): spiral is stable when
# Clb * Cnr - Cnb * Clr > 0. Weathercock + dihedral together drive this positive.
spiral_criterion = Clb * Cnr - Cnb * Clr

### ROLL AUTHORITY (rev7 — analytic aileron strip-theory model, note #1)
# AeroBuildup cannot model deflected control surfaces. We use a fully closed-form
# strip-theory steady-roll balance (aileron moment == roll-damping moment) so the
# constraint is smooth for IPOPT (no finite-differenced derivatives in the graph).
#
# Flap effectiveness tau(cf/c): linear fit to the classic Perkins & Hage control-
# surface chart over cf/c in [0.1, 0.4]  (tau ~ 0.38 at 0.25c, ~0.57 at 0.40c).
aileron_tau = 1.233 * aileron_chord_ratio + 0.077
# Aileron spanwise extent as fractions of semi-span.
f1_ail, f2_ail = aileron_y_start_frac, aileron_y_end_frac
y1_ail = f1_ail * (wingspan / 2)
y2_ail = f2_ail * (wingspan / 2)
# Strip-theory steady roll rate (constant-chord wing). Balancing aileron rolling
# moment  ~ a0*tau*c*(y2^2-y1^2)*da  against roll damping  ~ a0*(p/V)*c*b^3/12
# the section lift slope a0 and chord c cancel, leaving:
#   p = 3 * V * tau * (f2^2 - f1^2) * da / b
roll_rate = 3 * airspeed * aileron_tau * (f2_ail**2 - f1_ail**2) * aileron_deflection_max / b_w  # rad/s
# Reported aileron control power (for the report / cross-check only, not constrained):
Cl_da = (section_lift_slope * aileron_tau * chordlen / (S_w * b_w)) * (y2_ail**2 - y1_ail**2)

### ELEVATOR AUTHORITY (rev7 — pitch trim/maneuver). Elevator control power per radian:
#   Cm_de = eta_h * V_H * a_h * tau_e   (tail efficiency x tail volume x tail lift slope x flap eff)
# To trim across a CL range dCL the elevator deflection changes by  dδe = SM * dCL / Cm_de,
# and that must stay within the max deflection: Cm_de * δe_max >= SM * dCL_trim.
elevator_defl_max = np.radians(15)          # max elevator deflection for trim/maneuver
elevator_dCL_trim = 0.4                      # CL range to trim over (cruise + gust/maneuver)
a_hstab = 2 * np.pi * hstab_AR / (2 + hstab_AR)   # hstab 3D lift slope /rad
tail_efficiency = 0.9                        # dynamic-pressure ratio at the tail
elevator_tau = np.sqrt(elevator_chord_ratio) # empirical flap effectiveness (their model)
Cm_de = tail_efficiency * V_H * a_hstab * elevator_tau   # elevator pitch control power /rad



### PROPULSION
# APC-class 2-blade electric prop, low-Re. eta_prop(J) is a parabola peaking at J_peak, and
# J_peak scales with the prop's PITCH (P/D) -- so we let the optimizer SELECT the pitch (i.e.
# pick the right APC prop, e.g. P/D 0.5 = 16x8E, 0.75 = 16x12E) to match where the motor runs,
# keeping the prop near its efficiency peak instead of stuck off-design.
advance_ratio = opti.variable(init_guess=0.6, lower_bound=0.25, upper_bound=0.90, scale=0.1, category="advance_ratio")
prop_pitch_diam = opti.variable(init_guess=0.6, lower_bound=0.40, upper_bound=0.85, scale=0.1, category="prop_P_over_D")
# eta_prop(J) surrogate calibrated to APC-E 2-blade low-Re behavior (FLAGGED: fit, not table).
# Anchored to the prop's PITCH: thrust falls to zero near the windmill advance ratio
# J_windmill ~ 1.1*(P/D) (a bit past geometric pitch), and efficiency peaks around 75% of that.
# So eta -> 0 AT the windmill point (physical), instead of the old symmetric parabola that ran
# past it. The optimizer picks P/D so the operating J sits near the peak.
eta_prop_peak = 0.80              # peak propulsive efficiency for this prop class at low Re
J_windmill = 1.10 * prop_pitch_diam       # zero-thrust advance ratio (scales with pitch)
J_peak = 0.75 * J_windmill                # peak-efficiency J (~75% of windmill)
J_width = J_windmill - J_peak             # high side hits eta=0 exactly at the windmill point
# Floor at 0.15: the parabola goes negative far off-peak, and shaft power = thrust*V/eta_prop
# would blow up (NaN). A real prop below ~0.15 efficiency is unusable anyway; the optimizer
# stays near the peak regardless.
eta_prop = np.maximum(eta_prop_peak * (1 - ((advance_ratio - J_peak) / J_width) ** 2), 0.15)
prop_pitch_in = prop_pitch_diam * (propeller_diameter / 0.0254)  # selected pitch in inches (report)
rpm_cruise = 60 * airspeed / (advance_ratio * propeller_diameter)
# Guard drag positive: AeroBuildup can return D<0 on intermediate iterates, and the motor
# resistance model (power^-0.68) would then NaN. Floor it (restores the original guard).
_drag_per_motor = np.maximum(aero["D"], 0.1) / propeller_n
power_shaft_cruise = _drag_per_motor * airspeed / eta_prop   # shaft power per motor
Q = power_shaft_cruise / (2 * np.pi * rpm_cruise / 60)
Kt = 60 / (2 * np.pi * motor_kv)
# Total motor current = torque current (Q/Kt) + no-load current I0. At night cruise Q is
# tiny, so I0 is a large fraction of the draw -- omitting it understates night bus watts.
motor_i0 = hacker_motor_i0(motor_kv)  # A, FLAGGED (kv410 -> ~1.04 A datasheet anchor)
i_cruise = Q / Kt + motor_i0
motor_resistance = hacker_motor_resistance(power_shaft_cruise * 2.5, motor_kv)
# ESC efficiency: 0.95 is optimistic at ~27% night-cruise duty, so derate. FLAGGED.
esc_efficiency = 0.92
power_cruise = (power_shaft_cruise + (i_cruise**2) * motor_resistance) / esc_efficiency

# Motor terminal-voltage feasibility: the voltage the ESC must supply = back-EMF (rpm/Kv) plus
# the IR drop (i*R). If it exceeds the available bus, that operating point can't exist. Cruise
# must fit within the nominal 6S bus; enforced in CONSTRAINTS below.
bus_v_climb_max = 25.2  # V, fully-charged 6S bus (used for the climb feasibility check)
v_motor_cruise = i_cruise * motor_resistance + rpm_cruise / motor_kv

# Propeller tip mach limit.
propeller_tip_mach = 0.7  # From Dongjoon, 4/30/20
propeller_rads_per_sec = propeller_tip_mach * operating_atm.speed_of_sound() / (propeller_diameter / 2)
propeller_rpm_max = propeller_rads_per_sec * 30 / np.pi

# Climb.
thrust_climb = togw_design * g * np.sin(min_climb_angle * np.pi / 180) + aero["D"]
# Climb motor operating point (same prop rpm, higher torque) -> required terminal voltage must
# fit within the fully-charged 25.2 V bus, else the climb thrust is unachievable.
power_shaft_climb = (thrust_climb / propeller_n) * airspeed / eta_prop   # shaft W per motor
Q_climb = power_shaft_climb / (2 * np.pi * rpm_cruise / 60)
i_climb = Q_climb / Kt + motor_i0
v_motor_climb = i_climb * motor_resistance + rpm_cruise / motor_kv



### POWER
def _soft_min_cap(a, cap, width=20.0):
    # Overflow-safe smooth min(a, cap): a - softplus(a - cap). softplus is evaluated in the
    # numerically-stable max(z,0)+log1p(exp(-|z|)) form so exp() never overflows (unlike
    # np.softmin with hardness on ~1000 Wh args, which under/overflowed to NaN).
    z = (a - cap) / width
    softplus = width * (np.maximum(z, 0.0) + np.log(1.0 + np.exp(-np.abs(z))))
    return a - softplus

charge_power_max = battery_charge_c_rate * battery_capacity  # W, 0.25C constant-current limit
pack_r_internal_ohm = 0.045
r_bank = pack_r_internal_ohm * pack_energy_Wh / battery_capacity  # ohm (packs in parallel)
# FIX #3: Battery OCV vs SOC. The 6S bus sags toward 19.2 V when empty (20% SOC) and rises to
# 25.2 V full (nominal 22.2 V mid-SOC). Linear approx. Used for BOTH the I2R current (a fuller
# battery -> higher V -> lower current -> lower loss) and the CV charge taper below.
bus_v_empty, bus_v_full = 19.2, 25.2
def bus_voltage_at(soc_frac):
    return bus_v_empty + (bus_v_full - bus_v_empty) * (soc_frac - 0.20) / 0.80
# FIX #1: daily climb energy debit -- potential energy to reach loiter altitude divided by the
# drivetrain efficiency (eta_prop * ~0.85 motor+ESC), deducted once at dawn.
climb_altitude_gain = 300.0  # m daily climb to loiter (FLAGGED)
e_climb_wh = togw_design * g * climb_altitude_gain / (eta_prop * 0.85) / 3600.0  # Wh
# Dawn = SOC trough (start of day) -- where the one-time climb debit lands. Taken from the nominal
# SOC march so it tracks the true diurnal phase regardless of the flux array's phasing.
i_dawn = int(onp.argmin(_soc_seed))
for i in range(N-1):
    solar_flux = solar_flux_profile[i]  # W / m^2 (numeric constant)
    solar_area = solar_panel_n * solar_panel_side_length**2 # m^2
    # Per-timestep temperature-derated cell efficiency (hot midday cells generate less).
    power_generated = solar_flux * solar_area * solar_cell_efficiency_profile[i] / energy_generation_margin
    power_used = (power_cruise*propeller_n + 8 + 1)  # 8W avionics, 1W for NavLights
    charge_power = power_generated - power_used                          # W (>0 charging, <0 discharging)
    # Clamp SOC to [0.10, 1.0] so v_bus stays positive and <= bus_v_full at infeasible iterates
    # (hard min/max -> no NaN/overflow). At the solution the SOC constraints keep it in [0.2,1.0].
    soc_frac = np.minimum(np.maximum(battery_states[i] / battery_capacity_effective, 0.10), 1.0)
    v_bus = bus_voltage_at(soc_frac)  # OCV: bus sags when empty -> higher current -> more I2R loss
    # Charge cap = min(0.25C CC limit, CV taper). CV acceptable current = (Vmax - Vbus)/R, so as
    # the pack fills (Vbus -> Vmax) the charge it accepts -> 0. Hard mins (fmin) -> overflow-safe.
    cv_power_limit = v_bus * np.maximum(bus_v_full - v_bus, 0.0) / r_bank
    charge_cap = np.minimum(charge_power_max, cv_power_limit)
    charge_power_capped = np.minimum(charge_power, charge_cap)            # cap charging (unlimited discharge)
    i2r_loss = (charge_power_capped / v_bus) ** 2 * r_bank               # W, always a loss (uses OCV bus)
    climb_debit = e_climb_wh if i == i_dawn else 0.0                     # Wh, one-time dawn climb debit
    net_energy = (charge_power_capped - i2r_loss) * (dt / 3600) - climb_debit  # Wh

    # Smooth, overflow-safe cap at full capacity (surplus solar is dumped once the pack fills).
    battery_update = _soft_min_cap(battery_states[i] + net_energy, battery_capacity_effective)
    opti.subject_to(battery_states[i+1] == battery_update)



### Mass
# Power
mass_solar_cells = 0.008 * solar_panel_n  # C60 ~6.5 g cell + ~1.5 g interconnect (solar-uav-design)
mass_power_board = 0.050 * 2 # 75g estimate for each power board
# MPPT string sizing. The max cells per series string depends on the CONVERTER TOPOLOGY --
# this is a hardware decision, so it is a single explicit switch here:
#   "boost"      : boost-only (e.g. Genasun GVB). String Voc must stay below the DEPLETED bus
#                  (19.2 V) so the converter can always boost -> ~26 cells. Conservative.
#   "buck-boost" : steps voltage up OR down, so the limit is the MPPT input rating -> ~44 cells.
# DECISION NEEDED: set to match the MPPT you will actually fly. Default is the optimistic case.
mppt_topology = "buck-boost"       # FLAGGED: "boost" (26) or "buck-boost" (44)
cell_voc_cold = 0.713              # C60 open-circuit V at cold dawn (10 C) -> worst-case high V
_mppt_string_v_limit = 32.0 if mppt_topology == "buck-boost" else 19.2   # V
max_cells_per_string = onp.floor(_mppt_string_v_limit / cell_voc_cold)   # 44 (buck-boost) / 26 (boost)
mppt_unit_mass = 0.108            # kg per MPPT (Genasun PCB). Series strings: Imp ~5.9A < 8A ok.
# Number of MPPTs scales with cells (>= 1); mass counted in the budget.
n_mppt = np.softmax(solar_panel_n / max_cells_per_string, 1.0, hardness=2.0)
mass_mppt = n_mppt * mppt_unit_mass
# Amprius SA03 pack mass, scaled by nameplate energy (battery_capacity is NAMEPLATE Wh).
mass_batteries = battery_capacity * (pack_mass_kg / pack_energy_Wh)  # 378 Wh/kg
mass_wires = propulsion_electric.mass_wires(
    wire_length=wingspan / 2,
    max_current=power_out_max / battery_voltage,
    allowable_voltage_drop=battery_voltage * 0.01,
    material="aluminum"
)
num_packs = battery_capacity / pack_energy_Wh

# Avionics
mass_fc = .08 # Flight computer --> orange cube w/ carrier board.
mass_gps = 0.02 # GPS module.
mass_telemtry = 0.03 # 915Mhz telemetry module.
mass_receiver = 0.02 # Radio receiver.
mass_navlights = 0.01 * 4 # Four navlights.
mass_pitot = 0.01 # Pitot tube.
mass_avionics = mass_fc + mass_gps + mass_telemtry + mass_receiver + mass_navlights + mass_pitot  # power board counted separately (was double-counted)

# Actuators
mass_servos = .02 * 4 # 20g a servo

# Propulsion (measured/catalog masses per unit, x propeller_n)
mass_motor_raw = 0.275
mass_motors_mounted = mass_motor_raw * 1.1 * propeller_n  # +10% for mounts
mass_escs = 0.030 * propeller_n
mass_propellers = 0.055 * propeller_n

# Structures
mass_main_wing =  main_wing.area("wetted")*1.2 * 0.70 # 700g/m^2 planform (realistic; wing is ~2.5 kg, not buildable lighter)
# rev7: realistic tail masses (from build experience). The vertical fins are light --
# ~100 g each regardless of planform area -- so they're a fixed 0.1 kg each. The hstab
# stays area-based (it carries panels + elevator, so more skin = more mass).
mass_hstab = hor_stabilizer.area("wetted")*1.10 * 0.7
mass_vstabs = 0.1 * 2  # 100 g per fin, two fins (build experience)
mass_empennage = mass_hstab + mass_vstabs
mass_booms = 2 * 0.09 * boom_length  # kg (90 g/m per boom)
mass_fuselages = 0.1 * 2 # 0.1 kg per pod, two pods (rev7)
mass_superstructures = 0.2 # Two superstructures.
mass_boom_vstab_interfaces = 0.1 # Two boom-vstab interfaces.
mass_vstab_hstab_interfaces = 0.06 # Two vstab-stab interfaces.

## Total
total_mass = (
    mass_solar_cells +
    mass_mppt +        # rev7: buck-boost MPPTs (one per <=44-cell string)
    mass_power_board +
    mass_batteries +
    mass_wires + 
    mass_avionics + 
    mass_servos + 
    mass_motors_mounted + 
    mass_escs + 
    mass_propellers + 
    mass_main_wing +
    mass_empennage +   # rev7: floored hstab+vstab (was mass_hstab + mass_vstabs)
    mass_booms +
    mass_fuselages + 
    mass_superstructures +
    mass_boom_vstab_interfaces +
    mass_vstab_hstab_interfaces
)

### REAL CG (rev7). x measured aft (+) from the wing LE. Component x-positions are geometric;
# the battery position (battery_x) is a free variable so the aircraft can be balanced.
# The fixed-point constraint cg_le_dist == (sum m_i x_i)/total_mass makes the CG (and hence
# static margin & tail arms) reflect the ACTUAL mass layout instead of an assumed c/4.
_x_wing = 0.40 * chordlen                    # wing structure CG ~40% root chord
_x_solar = 0.35 * chordlen                   # cells sit forward/mid chord
_x_tail = boom_length + 0.25 * hstab_chordlen # empennage at the tail
_x_boom = 0.5 * boom_length                   # boom CG mid-length
_x_motor = -0.45                              # motors at the pod/boom fronts (forward of LE)
_x_prop = -0.55                               # props ahead of the motors
_x_pod = -0.5 * fuselage_length               # pod structure CG (pods extend forward)
_x_avionics = -0.25                           # avionics/power boards in the pods, forward
cg_moment = (
    mass_main_wing * _x_wing
    + mass_solar_cells * _x_solar
    + mass_mppt * _x_avionics
    + mass_empennage * _x_tail
    + mass_booms * _x_boom
    + mass_motors_mounted * _x_motor
    + mass_propellers * _x_prop
    + mass_escs * _x_motor
    + mass_fuselages * _x_pod
    + mass_avionics * _x_avionics
    + mass_power_board * _x_avionics
    + mass_wires * _x_avionics
    + mass_servos * _x_wing
    + mass_batteries * battery_x            # <- free position (balance ballast)
    + mass_superstructures * (0.5 * chordlen)
    + mass_boom_vstab_interfaces * _x_tail
    + mass_vstab_hstab_interfaces * _x_tail
)
opti.subject_to(cg_le_dist * total_mass == cg_moment)   # fixed-point: real mass-weighted CG



### CONSTRAINTS
# Mission
opti.subject_to(total_mass < togw_design)

# Aerodynamic
opti.subject_to(aero["L"] >= togw_design * g)
opti.subject_to(aero["Cm"] == 0)
opti.subject_to(aero["CL"] <= CLmax_cruise)

# Propulsion
opti.subject_to(v_motor_cruise <= battery_voltage)   # cruise: ESC voltage within nominal 6S bus
opti.subject_to(v_motor_climb <= bus_v_climb_max)    # climb: within fully-charged 6S bus (thrust achievable)
opti.subject_to(propeller_rpm_max >= rpm_cruise)
opti.subject_to(power_out_max >= power_cruise * propeller_n)
opti.subject_to(power_out_max >= thrust_climb * airspeed / (eta_prop * esc_efficiency))  # electrical watts, same domain as cruise

# Stall-speed margin (equivalent to constraining required CL below CLmax at cruise)
V_stall = np.sqrt(2 * (togw_design * g) / (operating_atm.density() * S_w * CLmax_cruise))
opti.subject_to(airspeed >= stall_speed_margin * V_stall)
opti.subject_to(wing_airfoil.max_thickness() * chord_root >= 0.025)  # root must accomodate batteries (20mm)
# ROW-AWARE solar panel model (honest about taper): at each spanwise station the number of
# chordwise panel ROWS that fit is floor(local_chord / panel_pitch), so a full-chord inboard
# station carries 2 rows while a narrow tapered tip carries 1 or 0. We integrate the panel-able
# area across the half-span with a smooth row-counter (sum of soft steps), so taper is paid for
# in lost rows rather than by a rigid rectangular strip.
panel_pitch = 0.13            # m, 12.5 cm cell + gap
panel_packing_eff = 0.85      # spanwise gaps, TE curvature, spar rail
_panel_rows_max = 4           # most rows the widest plausible chord could hold

def _rows_at(c):
    # Smooth floor(c / panel_pitch): +1 for each row-threshold the chord clears.
    return sum(0.5 * (1 + np.tanh((c - r * panel_pitch) * 60.0)) for r in range(1, _panel_rows_max + 1))

# Integrate usable panel area over the half-span (both wings), minus the outboard aileron TE.
_K = 30
_panel_area_available = 0
for _i in range(_K):
    _yf = (_i + 0.5) / _K                                   # station center, fraction of semi-span
    _c = chord_at(_yf * semi_span)
    _c_avail = _c - aileron_chord_ratio * _c * (0.5 * (1 + np.tanh((_yf * semi_span - aileron_y_start_frac * semi_span) * 60.0)))  # aileron eats TE outboard
    _panel_area_available = _panel_area_available + _rows_at(_c_avail) * panel_pitch * (semi_span / _K)
# H-stab carries cells on its FIXED portion only: elevator is a ~30% chord keep-out and a
# ~10% LE margin, so ~0.6 of the hstab planform is panel-able (matches solar-uav-design).
hstab_panel_frac = 0.6
# H-stab panel area is ROW-AWARE too (not raw area): the panel-able band is hstab_panel_frac of
# the hstab chord, and _rows_at() returns ~0 when that band is narrower than one panel pitch --
# so a short-chord hstab honestly carries ZERO cells (its area alone can't hold a 12.5cm panel).
hstab_panel_area = _rows_at(hstab_panel_frac * hstab_chordlen) * panel_pitch * hstab_span
wing_panel_capacity = 2 * _panel_area_available * panel_packing_eff / panel_pitch**2   # cells the wing can hold
hstab_panel_capacity = hstab_panel_area * panel_packing_eff / panel_pitch**2           # cells the hstab can hold
usable_panel_area = (2 * _panel_area_available + hstab_panel_area) * panel_packing_eff
opti.subject_to(solar_panel_n * panel_pitch**2 <= usable_panel_area)  # panels must physically fit (row-aware)
# Split the solved cell count: fill the wing first, overflow onto the hstab (reported numbers).
# Fill the wing first; overflow above wing capacity goes to the hstab, capped at what the hstab
# can physically hold (hstab_panel_capacity, which is ~0 when the hstab chord is too short for a
# panel row). Report-only, so hard min/max is fine. Viewer reads these exact numbers.
_overflow = np.maximum(solar_panel_n - wing_panel_capacity, 0.0)
panels_on_hstab = np.minimum(_overflow, hstab_panel_capacity)
panels_on_wing = solar_panel_n - panels_on_hstab
# Root chord must fit at least the 2 inboard rows (plus a little for the aileron behind them).
opti.subject_to(chordlen >= solar_panel_n_rows * panel_pitch + 0.03)
# Taper must leave the tip at least one panel-row wide (keeps the taper physical).
opti.subject_to(chord_tip >= panel_pitch)
opti.subject_to(hstab_chordlen >= 0.22)  # >=22cm so the 60% panel band (>=0.13m) fits a full row of cells
opti.subject_to(vstab_root_chord >= hstab_chordlen)  # enforce an actual taper (root >= tip) and avoid negative root chord

# Stability — longitudinal
opti.subject_to(L_H >= 1e-3)
opti.subject_to(L_V >= 1e-3)
opti.subject_to(boom_length >= chordlen)
opti.subject_to(static_margin >= SM_min)
opti.subject_to(static_margin <= SM_max)

# Stability — lateral-directional (rev7, the flight-critical rev6 fixes).
# NOTE: constrained via smooth geometry levers, NOT AeroBuildup's finite-differenced
# Cnb/Clb/Clp (those are noisy and make IPOPT's dual infeasibility blow up). The
# AeroBuildup derivatives are still evaluated and reported as an independent check.
opti.subject_to(V_V >= V_V_floor)              # weathercock: fin volume floor (note #2)
opti.subject_to(dihedral_angle >= dihedral_min)  # dihedral effect / spiral (note #3)
opti.subject_to(roll_rate >= roll_rate_target)   # aileron roll authority (note #1)
# Elevator must trim the CL range within +/-15 deg (this pushes back on very high static margin).
opti.subject_to(Cm_de * elevator_defl_max >= static_margin * elevator_dCL_trim)

# Power (SOC window on DELIVERABLE energy)
opti.subject_to(battery_states >= battery_capacity_effective * (1-allowable_battery_depth_of_discharge))
opti.subject_to(battery_states <= battery_capacity_effective)
opti.subject_to(battery_states[0] <= battery_states[N-1])
# NOTE: no explicit daily-refill constraint. Minimizing span (hence solar area + battery mass)
# already drives the smallest pack that is fully cycled -- min SOC floors at 0.20 and the peak
# naturally reaches ~0.95, so the pack does not carry dead capacity. A fixed-index refill target
# only mis-pinned a morning point and inflated span with no physical benefit.
opti.subject_to(num_packs <= 8)
# Optional integer-pack pin for a discrete-reality solve (env FIXED_PACKS=2/3/4...).
if fixed_packs is not None:
    opti.subject_to(battery_capacity == float(fixed_packs) * pack_energy_Wh)
    print(f"[packs] pinned to {fixed_packs} packs = {float(fixed_packs)*pack_energy_Wh:.0f} Wh nameplate")



### SOLVE
opti.minimize(wingspan)

try:
    sol = opti.solve(max_iter=5000)
    _bs = sol.value(battery_states); _cap = sol.value(battery_capacity_effective)
    print(f"[SOC] min={_bs.min()/_cap:.3f}  peak={_bs.max()/_cap:.3f} (i={int(_bs.argmax())})  dawn(i={i_dawn})={_bs[i_dawn]/_cap:.3f}  start={_bs[0]/_cap:.3f} end={_bs[-1]/_cap:.3f}")
except Exception as e:
    print(f"[solve] failed: {e}")
    print("\n[DEBUG] Variable values at failure:")
    # (label, expr, format) — format is applied to opti.debug.value(expr).
    _debug_vars = [
        ("airspeed", airspeed, "{:.3f} m/s"),
        ("togw_design", togw_design, "{:.3f} kg"),
        ("power_out_max", power_out_max, "{:.3f} W"),
        ("solar_panel_n", solar_panel_n, "{:.3f}"),
        ("battery_capacity", battery_capacity, "{:.3f} Wh"),
    ]
    battery_states_vals = opti.debug.value(battery_states)
    _debug_vars += [
        ("battery_states min", float(onp.min(battery_states_vals)), "{:.3f} Wh"),
        ("battery_states max", float(onp.max(battery_states_vals)), "{:.3f} Wh"),
        ("battery_states[0]", float(battery_states_vals[0]), "{:.3f} Wh"),
        ("battery_states[-1]", float(battery_states_vals[-1]), "{:.3f} Wh"),
        ("wingspan", wingspan, "{:.3f} m"),
        ("chordlen", chordlen, "{:.3f} m"),
        ("struct_defined_aoa", struct_defined_aoa, "{:.3f} deg"),
        ("hstab_AR", hstab_AR, "{:.3f}"),
        ("hstab_aoa", hstab_aoa, "{:.3f} deg"),
        ("boom_length", boom_length, "{:.3f} m"),
        ("num_packs", num_packs, "{:.3f}"),
        ("cg_le_dist", cg_le_dist, "{:.3f} m"),
        ("boom_y", boom_y, "{:.3f} m"),
        ("total_mass", total_mass, "{:.3f} kg"),
        ("static_margin", static_margin, "{:.3f}"),
        ("dihedral_angle", dihedral_angle, "{:.3f} deg"),
        ("V_V", V_V, "{:.4f}"),
        ("vstab_AR", vstab_AR, "{:.3f}"),
        ("Cnb", Cnb, f"{{:.4f}} /rad (min {Cnb_min})"),
        ("Clb", Clb, f"{{:.4f}} /rad (max {Clb_max})"),
        ("Clp", Clp, "{:.4f} /rad"),
        ("spiral_criterion", spiral_criterion, "{:.5f}"),
        ("aileron_chord_ratio", aileron_chord_ratio, "{:.3f}"),
        ("Cl_da", Cl_da, "{:.4f} /rad"),
        ("L_H", L_H, "{:.3f} m"),
        ("L_V", L_V, "{:.3f} m"),
        ("mass_wires", mass_wires, "{:.3f} kg"),
        ("power_cruise", power_cruise, "{:.3f} W"),
        ("power_shaft_cruise", power_shaft_cruise, "{:.3f} W"),
        ("aero Cm", aero["Cm"], "{:.6f}"),
    ]
    # `float(...)` entries are already resolved; only CasADi exprs need debug.value.
    for label, expr, fmt in _debug_vars:
        val = expr if isinstance(expr, float) else opti.debug.value(expr)
        print(f"  {label}: " + fmt.format(val))
    # roll_rate needs a degrees conversion around the resolved value, so it stays inline.
    print(f"  roll_rate: {onp.degrees(opti.debug.value(roll_rate)):.2f} deg/s (target {onp.degrees(roll_rate_target):.1f})")
    sol = None



### REPORT
report_raw = {
    "Meta": {
        # run_id is filled in once run_dir is created (near bottom of file).
        "run_id": run_id_random(4, "run"),
        "timestamp_utc": datetime.now(timezone.utc).isoformat()
    },
    "Mission": {
        "mission_date": mission_date,
        "operating_lat": operating_lat,
        "operating_altitude": operating_altitude
    },
    "Environment": {
        "temperature_K": operating_atm.temperature(),
        "temperature_F": operating_atm.temperature() * 1.8 - 459.67,
        "pressure": operating_atm.pressure(),
        "density": operating_atm.density(),
        "speed_of_sound": operating_atm.speed_of_sound(),
    },
    "Performance": {
        "airspeed": airspeed,
        "thrust_cruise": aero["D"],
        "thrust_climb": thrust_climb,
        "power_cruise (all motors)": power_cruise * propeller_n,
        "power_cruise (one motor)": power_cruise,
        "power_shaft_cruise (one motor)": power_shaft_cruise,
        "min_climb_angle": min_climb_angle,
        "power_out_max": power_out_max,
        "total_mass": total_mass,
        "togw": togw_design * g,
        "stall_speed": V_stall,
        "L_over_D": aero["CL"] / aero["CD"],
    },
    "Aerodynamics": {
        "CL": aero["CL"],
        "CD": aero["CD"],
        "Cm": aero["Cm"],
        "L": aero["L"],
        "D": aero["D"],
        "x_np": aero["x_np"],
        "static_margin": static_margin,
        "Cmq": Cmq,                       # pitch damping (<0), for the JSBSim flight model
        "Cm_de": Cm_de,                   # elevator pitch control power /rad (analytic, rev7)
    },
    "Drag": drag_breakdown,
    "Lateral_Directional_Stability": {
        "Cnb": Cnb,                       # weathercock (>0 stable), target >= Cnb_min
        "Cnb_min": Cnb_min,
        "Clb": Clb,                       # dihedral effect (<0 stable), target <= Clb_max
        "Clb_max": Clb_max,
        "Clp": Clp,                       # roll damping
        "Cnr": Cnr,                       # yaw damping
        "Clr": Clr,
        "Cnp": Cnp,                       # yaw due to roll rate (adverse/proverse yaw)
        "spiral_criterion": spiral_criterion,  # >0 spiral-stable
        "dihedral_angle": dihedral_angle,
        "roll_rate_rad_s": roll_rate,
        "roll_rate_deg_s": np.degrees(roll_rate),
        "roll_rate_target_deg_s": np.degrees(roll_rate_target),
        "aileron_chord_ratio": aileron_chord_ratio,
        "aileron_tau": aileron_tau,
        "aileron_y_start": y1_ail,
        "aileron_y_end": y2_ail,
        "Cl_da": Cl_da,
    },
    "Main Wing": {
        "wingspan": wingspan,
        "chordlen": chordlen,
        "taper_ratio": taper_ratio,
        "chord_tip": chord_tip,
        "dihedral_angle": dihedral_angle,
        "struct_defined_aoa": struct_defined_aoa,
        "S_w": S_w,
        "MAC_w": MAC_w,
        "b_w": b_w,
        "main_wing_area": main_wing.area(),
        "main_wing_AR": main_wing.aspect_ratio(),
        "wing_airfoil": str(getattr(wing_airfoil, "name", "s4110")),
    },
    "HStab": {
        "hstab_AR": hstab_AR,
        "hstab_aoa": hstab_aoa,
        "hstab_area": hstab_area,
        "hstab_span": hstab_span,
        "hstab_chordlen": hstab_chordlen,
        "L_H": L_H,
        "V_H_actual": (hstab_area * L_H) / (S_w * MAC_w),
        "hstab_airfoil": str(getattr(tail_airfoil, "name", "naca0010")),
    },
    "V Stab": {
        "vstab_AR": vstab_AR,
        "V_V": V_V,
        "vstab_area": vstab_area,
        "vstab_area_total": vstab_area_total,
        "vstab_span": vstab_span,
        "vstab_root_chord": vstab_root_chord,
        "L_V": L_V,
        "V_V_actual": (vstab_area_total * L_V) / (S_w * b_w),
        "vstab_airfoil": str(getattr(tail_airfoil, "name", "naca0010")),
    },
    "Geometry": {
        "cg_le_dist": cg_le_dist,
        "boom_y": boom_y,
        "boom_length": boom_length,
        "boom_radius": boom_radius,
        "boom_spacing_frac": boom_spacing_frac,
        "fuselage_length": fuselage_length,
    },
    "Power": {
        "solar_panel_n": solar_panel_n,
        "panels_on_wing": panels_on_wing,
        "panels_on_hstab": panels_on_hstab,
        "wing_panel_capacity": wing_panel_capacity,
        "hstab_panel_capacity": hstab_panel_capacity,
        "n_mppt": n_mppt,
        "max_cells_per_string": max_cells_per_string,
        "battery_capacity": battery_capacity,
        "battery_states": battery_states,
        "battery_voltage": battery_voltage,
        "solar_panel_side_length": solar_panel_side_length,
        "solar_panel_n_rows": solar_panel_n_rows,
        "cell_efficiency_stc": cell_efficiency_stc,
        "solar_system_eff_stc": solar_cell_efficiency,
        "energy_generation_margin": energy_generation_margin,
        "num_packs": num_packs,
    },
    "Propulsion": {
        "propeller_n": propeller_n,
        "propeller_diameter": propeller_diameter,
        "advance_ratio": advance_ratio,
        "eta_prop": eta_prop,
        "prop_pitch_over_diam": prop_pitch_diam,
        "prop_pitch_in": prop_pitch_in,
        "motor_kv": motor_kv,
        "rpm_cruise": rpm_cruise,
        "i_cruise": i_cruise
    },
    "Masses": {
        "mass_solar_cells": mass_solar_cells,
        "mass_mppt": mass_mppt,
        "mass_power_board": mass_power_board,
        "mass_batteries": mass_batteries,
        "mass_wires": mass_wires,
        "mass_avionics": mass_avionics,
        "mass_servos": mass_servos,
        "mass_motors_mounted": mass_motors_mounted,
        "mass_escs": mass_escs,
        "mass_propellers": mass_propellers,
        "mass_main_wing": mass_main_wing,
        "mass_hstab": mass_hstab,
        "mass_vstab": mass_vstabs,
        "mass_boom": mass_booms,
        "mass_fuselages": mass_fuselages,
        "mass_superstructures": mass_superstructures,
        "mass_boom_vstab_interfaces": mass_boom_vstab_interfaces,
        "mass_vstab_stab_interfaces": mass_vstab_hstab_interfaces,
        "total_mass": total_mass,
    },
}



### OUTPUT ARTIFACTS
# Only create artifacts if simulation succeeded.
if sol is None:
    print("[artifacts] simulation failed, skipping artifact generation")
else:
    # Convert symbolic/raw Opti values into plain JSON-safe numbers/lists.
    run_dir = Path(__file__).resolve().parent / "output" / report_raw["Meta"]["run_id"]
    run_dir.mkdir(parents=True, exist_ok=False)
    
    soln = process_raw_values(report_raw, sol)
    write_json(run_dir / "soln.json", soln, indent=2)
    
    # Compute solved airplane once for all exports.
    airplane_sol = sol(airplane) if callable(sol) else airplane
    
    # Export an XFLR5 XML (synthesizes a single fin for the twin-fin geometry).
    try:
        export_xflr5_xml_from_soln(airplane_sol=airplane_sol, soln=soln, out_path=run_dir / "aircraft.xml")
    except Exception as e:
        print(f"[export_xflr5_xml_from_soln] skipped due to error: {e}")

    # Draw airplane (needs the optional `pyvista` viewer; skip cleanly if absent).
    try:
        airplane_sol.draw()
    except Exception as e:
        print(f"[draw] skipped ({e})")