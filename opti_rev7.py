import os
import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.library import power_solar, propulsion_electric, propulsion_propeller
from aerosandbox.atmosphere.atmosphere import Atmosphere
import numpy as onp
import pvlib
from pathlib import Path
from datetime import datetime, timezone

from lib.artifacts import process_raw_values, run_id_random, write_json
from lib.exports import export_xflr5_xml_from_soln, export_airplane_glb
from lib.models import hacker_motor_resistance, hacker_motor_i0
from lib.propulsion import prop_ct_cp, thrust_curve_from_soln, prop_coeff_table



opti = asb.Opti()
opti.solver("ipopt")



### CONSTANTS
## Physics
g = 9.81

## Mission
mission_date = 172
operating_lat = 37.398928
togw_max = 12 # kg
temperature_high = 34 # 34°K hotter than ISA --> leads to ~93 °F operating temperature
operating_altitude = 1200 # in meters
operating_atm = Atmosphere(operating_altitude, temperature_deviation=temperature_high)

## Performance
# Climb requirement (rev7, ported from Adam Yang's solar-uav-design): the design must sustain
# a real climb RATE at loiter speed, not an unrealistically steep fixed angle. A solar plane
# does not climb at 30 deg; 1.5 m/s is the operational spec. thrust_climb is derived from this
# rate at the flown airspeed (see CLIMB block). The old 30 deg angle is kept only for the report.
climb_rate_req = 1.5   # m/s, minimum sustained climb rate at loiter speed
min_climb_angle = 30   # degrees (legacy, report-only; superseded by climb_rate_req)

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
stall_speed_margin = 1.2   # V_cruise >= margin * V_stall (straight-line reference)
# Banked-loiter margin (ported from Adam): a station-keeping orbit is flown banked, so the
# operative floor is 1.3*Vs (load factor n~1.3, ~40 deg bank) -- higher than the 1.2 cruise
# margin. This is the binding stall floor in the CONSTRAINTS block.
loiter_stall_factor = 1.3  # V_cruise >= loiter_stall_factor * V_stall
# Crosswind floor (ported from Adam): never fly slower than a 15 kt crosswind, or the aircraft
# cannot make headway / hold the orbit on a windy day.
crosswind_floor_ms = 7.7   # m/s (15 kt)

## Structural
structural_mass_markup = 1.2
boom_radius = 0.01 # radius of the boom in meters.
boom_spacing_frac = 0.25 # Fraction that the inboard wing spans.
fuselage_length = 0.98   # m, per-pod length (Adam: CAD pod is 0.980 m, wing TE -> prop disk)
fuselage_radius = 0.0281 # m, pod radius from Adam's CAD wetted area (173128 mm^2 / (2*pi*0.98 m))

## Propulsion
propeller_n = 2 # Number of motors / propellers.
# (advance ratio is now a design variable in the PROPULSION section, with an eta_prop(J) map)

## Power
battery_voltage = 22.2
N = 180 # Number of discretization points.
time = onp.linspace(0, 24 * 60 * 60, N) # s (numeric; does not need to be symbolic)
dt = onp.diff(time)[0]  # s
solar_panel_side_length = 0.125 # m (12.5cm)
min_panel_rows = 2   # design floor: the wing chord must hold >= this many cell rows (min-chord constraint).
# NOTE: this is NOT the actual row count -- the real chordwise rows = floor(chord/panel_pitch),
# computed and exported as solar_panel_n_rows near the report so the data panel is honest.
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

# SOLAR IRRADIANCE -- pvlib Ineichen clear-sky model (replaces aerosandbox.solar_flux).
# Physically-based: Linke turbidity (aerosol/water-vapor), altitude correction, and a proper
# BEAM/DIFFUSE split so the incidence-angle reflection loss (IAM) is applied per-component --
# beam gets the ASHRAE cos-law loss (large at dawn/dusk), diffuse gets a near-constant sky loss.
# Computed in SOLAR time (noon = 43200 s) to match the model's clock; longitude only shifts noon,
# which the model self-locates from the flux peak, so it isn't needed here.
LINKE_TURBIDITY = 3.0   # clear dry-continental TL (FLAGGED; Adam uses a lat/lon climatology lookup)
IAM_B0 = 0.05           # ASHRAE beam reflection coefficient (glass/encapsulant)
IAM_DIFFUSE = 0.95      # near-constant sky-diffuse incidence-angle modifier
_lat_rad = onp.radians(operating_lat)
_decl = onp.radians(23.45) * onp.sin(2 * onp.pi * (284 + mission_date) / 365.0)   # solar declination
_H = onp.radians(15.0 * (time / 3600.0 - 12.0))                                   # hour angle (noon=12h)
_cos_zen = onp.clip(onp.sin(_lat_rad) * onp.sin(_decl)
                    + onp.cos(_lat_rad) * onp.cos(_decl) * onp.cos(_H), -1.0, 1.0)
_zenith = onp.degrees(onp.arccos(_cos_zen))                    # solar zenith angle, deg
_app_zen = onp.minimum(_zenith, 90.0)
_press = pvlib.atmosphere.alt2pres(operating_altitude)
_airmass = pvlib.atmosphere.get_absolute_airmass(
    pvlib.atmosphere.get_relative_airmass(_app_zen), _press)
_dni_extra = pvlib.irradiance.get_extra_radiation(mission_date)
_cs = pvlib.clearsky.ineichen(_app_zen, _airmass, LINKE_TURBIDITY,
                              altitude=operating_altitude, dni_extra=_dni_extra)
_dni = onp.nan_to_num(onp.asarray(_cs["dni"], dtype=float), nan=0.0)   # beam (direct normal)
_dhi = onp.nan_to_num(onp.asarray(_cs["dhi"], dtype=float), nan=0.0)   # diffuse horizontal
_sun_up = _cos_zen > 0.0
# Per-component IAM on a HORIZONTAL panel: beam AOI = zenith, so cos(AOI) = cos_zen.
_iam_beam = onp.clip(1.0 - IAM_B0 * (1.0 / onp.maximum(_cos_zen, 1e-3) - 1.0), 0.0, 1.0)
solar_flux_profile = onp.where(
    _sun_up, _dni * _cos_zen * _iam_beam + _dhi * IAM_DIFFUSE, 0.0)
solar_flux_profile = onp.maximum(onp.nan_to_num(solar_flux_profile, nan=0.0, posinf=0.0, neginf=0.0), 0.0)
iam_profile = onp.where(_sun_up, _iam_beam, 0.0)   # kept for the [IAM] report line below
_daily_wh_m2 = onp.sum((solar_flux_profile[:-1] + solar_flux_profile[1:]) / 2 * onp.diff(time)) / 3600.0
print(
    f"[solar_flux_profile] pvlib Ineichen (TL={LINKE_TURBIDITY}): "
    f"peak={solar_flux_profile.max():.1f} W/m^2, daily={_daily_wh_m2:.0f} Wh/m^2, "
    f"nan_count={onp.isnan(solar_flux_profile).sum()}"
)
print(f"[IAM] noon iam={iam_profile.max():.3f}, mean(daylight)={iam_profile[_sun_up].mean():.3f}")

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

# PHASE 5 -- DAY/NIGHT DENSITY (ported from Adam). The aircraft loiters at its min-power CL, where
# drag D = W*CD/CL is density-INDEPENDENT but power P = D*V scales as V ~ rho^-0.5. So per timestep
# the cruise power scales with the local air density, which (at fixed altitude) tracks ambient
# temperature: rho_i/rho_peak = T_peak/T_i, hence P_i/P_peak = (T_i/T_peak)^0.5. Cool night air is
# denser -> the plane flies slower -> less power; hot afternoon air is thinner -> more power. The
# stall/climb sizing still uses operating_atm (the hottest, THINNEST air = worst case), so this
# only RELAXES the night bus, which is the binding quantity. Numeric array -> IPOPT-friendly.
_T_peak_K = float(operating_atm.temperature())              # hottest (thinnest) design point
power_density_scale = onp.sqrt((T_amb_profile + 273.15) / _T_peak_K)   # <=1, min at pre-dawn
print(f"[day/night rho] cruise-power scale: night(min)={power_density_scale.min():.3f}  "
      f"day(peak)={power_density_scale.max():.3f}  (1.0 = hot-thin design density)")



### VARIABLES
## Performance
airspeed = opti.variable(init_guess=10.8, lower_bound=5, upper_bound=15, scale=5, category="airspeed")
togw_design = opti.variable(init_guess=11.9, lower_bound=1, upper_bound=togw_max, scale=1, category="togw_max")
power_out_max = opti.variable(init_guess=500, lower_bound=25*16, scale=100, category="power_out_max")

## Propulsion
propeller_diameter = opti.parameter(value=0.4)  # m (fixed, ~16 in)
motor_kv = opti.variable(init_guess=500, lower_bound=50, upper_bound=2000, scale=100, category="motor_kv")

## Avionics
solar_panel_n = opti.variable(init_guess=128, lower_bound=10, scale=10, category="solar_panel_n")
battery_capacity = opti.variable(init_guess=1449, lower_bound=100, scale=100, category="battery_capacity") # NAMEPLATE installed energy in Wh (drives mass + pack count)
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
# Span cap depends on the objective regime: in FREE (margin-max) mode, size is Adam's hard
# constraint (6.2 m -- the tightest span at which this airframe still closes; 6.0 is infeasible
# with honest props). In PINNED-pack (min-span) mode, span is free to find the natural value that
# closes with a fixed integer pack count, so the cap is generous.
# Span cap = the real physical/transport limit.
#   KNOWN-INFEASIBLE at 6.0 m (deliberate): with honest physics + the 2-row chord cap (max_panel_rows),
#   no design closes at a 6.0 m span -- the thin fully-packed wing needs ~7.6 m of span for its lift
#   area, and even a fat-chord wing only reaches ~6.15 m. This matches Adam's optimizer verdict
#   ("nothing closes under 12 kg / 6 m"). The 6.0 m value is kept as the stated hard requirement, so
#   a run will report local infeasibility -- that IS the finding (the 6 m goal needs better tech or a
#   relaxed mission, not a modeling fudge). To get a converging design: raise to ~8.0 (thin packed
#   wing, +101 Wh) or relax max_panel_rows for a fatter 6.15 m wing.
_span_cap_m = 6.0
wingspan = opti.variable(init_guess=min(5.95, _span_cap_m), lower_bound=2, upper_bound=_span_cap_m, scale=1, category="wingspan")
chordlen = opti.variable(init_guess=0.56, lower_bound=0.05, scale=0.1, category="chordlen")
struct_defined_aoa = opti.variable(init_guess=4.7, lower_bound=0, upper_bound=10, scale=1, category="struct_aoa")
# Real CG (rev7): cg_le_dist is now a VARIABLE tied by a fixed-point constraint to the actual
# mass-weighted CG of all components (computed after the mass rollup). The battery position is
# free so the optimizer can place it to balance the aircraft, exactly like a real build.
cg_le_dist = opti.variable(init_guess=0.08, lower_bound=-0.20, upper_bound=0.50, scale=0.1, category="cg_le_dist")
battery_x = opti.variable(init_guess=-0.20, lower_bound=-0.55, upper_bound=0.30, scale=0.1, category="battery_x")

# Wing dihedral (from the root) and taper (flat LE; TE tapers forward, chord shrinks
# from chordlen at the root to taper_ratio*chordlen at the tip).
dihedral_angle = opti.variable(init_guess=2.5, lower_bound=0, upper_bound=4, scale=1, category="dihedral")
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
# Adam-style straddling fin: a rectangular fin split half above / half below the H-stab plane, so
# its area = height * chord (no T-tail trapezoid tip). chord = area / height.
vstab_root_chord = vstab_area / vstab_span

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
).translate([boom_length, 0, 0])   # Adam-style: H-stab in the wing/boom plane (z=0), between the booms

# Both fins are identical trapezoids mirrored across XZ — build once, place by sign.
def _make_fin(name, sign):
    # Adam-style straddling fin: extends from -h/2 (below) to +h/2 (above) the H-stab plane so its
    # side-force acts near z=0 (little roll from yaw). Placed at the boom, LE aft of the H-stab TE.
    return asb.Wing(
        name=name,
        symmetric=False,
        xsecs=[
            asb.WingXSec(xyz_le=[0, 0, -0.5 * vstab_span], chord=vstab_root_chord, twist=0, airfoil=tail_airfoil),
            asb.WingXSec(xyz_le=[0, 0, 0.5 * vstab_span], chord=vstab_root_chord, twist=0, airfoil=tail_airfoil),
        ],
    ).translate([boom_length + hstab_chordlen, sign * boom_y, 0])

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
                radius=fuselage_radius * min(1.0, 6.0 * xi * (1.0 - 0.4 * xi)),   # Adam-style near-constant pod, faired nose/tail
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
# APC-class 2-blade electric prop, low-Re. The optimizer SELECTS the pitch P/D (i.e. which real
# APC prop: 0.5 = 16x8E, 0.625 = 16x10E, 0.75 = 16x12E) and the operating advance ratio J, and the
# real C_T/C_P map below (fit to APC data) sets thrust/power -- no hand-tuned parabola.
advance_ratio = opti.variable(init_guess=0.6, lower_bound=0.25, upper_bound=0.90, scale=0.1, category="advance_ratio")
prop_pitch_diam = opti.variable(init_guess=0.6, lower_bound=0.40, upper_bound=0.75, scale=0.1, category="prop_P_over_D")  # 0.75 = 16x12E, top of the real APC-E fit range
# REAL APC PROP MAP (rev7): C_T(J,P/D) and C_P(J,P/D) are smooth degree-3 polynomial fits to the
# ACTUAL APC thin-electric 16in family (16x6E/8E/10E/12E). The fit itself lives in
# lib/propulsion.py -- the SINGLE source of truth for prop physics -- so the optimizer and the
# baked thrust curve (Propulsion.thrust_curve, consumed by the sibling sim) can never diverge.
# The real propulsive efficiency at an operating point is eta = C_T*J/C_P; the optimizer picks
# P/D (which real APC prop) and J (operating point) to sit near the efficiency peak.
prop_CT, prop_CP = prop_ct_cp(advance_ratio, prop_pitch_diam)
eta_prop = prop_CT * advance_ratio / np.maximum(prop_CP, 1e-3)   # real propulsive efficiency from the map
# Propulsion conservatism (ported from Adam): even the real table is ~8% optimistic on thrust and
# ~8% low on power vs a measured stand -> fold into eta_eff = eta * 0.92/1.08. Downstream shaft-
# power/climb-energy conversions use eta_prop_eff.
prop_thrust_derate = 0.92   # real thrust / table thrust (FLAGGED: calibrate vs thrust stand)
prop_power_inflate = 1.08   # real shaft power / table power (FLAGGED)
eta_prop_eff = eta_prop * prop_thrust_derate / prop_power_inflate
prop_pitch_in = prop_pitch_diam * (propeller_diameter / 0.0254)  # selected pitch in inches (report)
rpm_cruise = 60 * airspeed / (advance_ratio * propeller_diameter)
# Guard drag positive: AeroBuildup can return D<0 on intermediate iterates, and the motor
# resistance model (power^-0.68) would then NaN. Floor it (restores the original guard).
_drag_per_motor = np.maximum(aero["D"], 0.1) / propeller_n
power_shaft_cruise = _drag_per_motor * airspeed / eta_prop_eff   # shaft power per motor (derated)
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

# Climb. Ported from Adam: size for a sustained climb RATE, not a fixed 30 deg angle. At loiter
# speed the required climb angle is gamma = asin(climb_rate/V), so the thrust to hold that climb
# is W*sin(gamma) + D = W*(climb_rate/airspeed) + D. This is the realistic (and much gentler)
# spec for an endurance solar plane than the old ~30 deg angle.
climb_angle_rad = np.arcsin(climb_rate_req / airspeed)   # reported climb angle at loiter speed
thrust_climb = togw_design * g * (climb_rate_req / airspeed) + aero["D"]
# Climb motor operating point (same prop rpm, higher torque) -> required terminal voltage must
# fit within the fully-charged 25.2 V bus, else the climb thrust is unachievable.
power_shaft_climb = (thrust_climb / propeller_n) * airspeed / eta_prop_eff   # shaft W per motor (derated)
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
e_climb_wh = togw_design * g * climb_altitude_gain / (eta_prop_eff * 0.85) / 3600.0  # Wh (derated prop)
# Avionics POWER budget, itemized per component (mirrors the avionics MASS breakdown below), day
# + night standby. Replaces the old unjustified 8 W lump, which was ~2x a bottoms-up estimate of
# these same parts and 1.6x Adam's team-measured 5 W. If a companion computer / payload is added,
# add its draw here explicitly rather than padding the total.
power_fc = 2.5           # W, flight computer (Orange Cube + carrier board)
power_gps = 0.25         # W, GPS module
power_telemetry = 1.0    # W, 915 MHz telemetry (avg incl. TX bursts)
power_receiver = 0.15    # W, radio receiver
power_pitot = 0.10       # W, pitot / airspeed sensor
power_avionics_w = power_fc + power_gps + power_telemetry + power_receiver + power_pitot  # ~4.0 W
power_navlights_w = 0.25 * 4   # W, four nav lights (~1.0 W)  -> total standby ~5.0 W (matches Adam)
# Dawn = SOC trough (start of day) -- where the one-time climb debit lands. Taken from the nominal
# SOC march so it tracks the true diurnal phase regardless of the flux array's phasing.
i_dawn = int(onp.argmin(_soc_seed))
for i in range(N-1):
    solar_flux = solar_flux_profile[i]  # W / m^2 (numeric constant)
    solar_area = solar_panel_n * solar_panel_side_length**2 # m^2
    # Per-timestep temperature-derated cell efficiency (hot midday cells generate less).
    power_generated = solar_flux * solar_area * solar_cell_efficiency_profile[i] / energy_generation_margin
    # Day/night density (Phase 5): scale the propulsion power by the local air density; the
    # itemized avionics + nav-light standby (~5 W total) is density-independent.
    power_used = (power_cruise * propeller_n * power_density_scale[i] + power_avionics_w + power_navlights_w)
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
# PER-BAY MPPT PACKING (ported from Adam). Series strings cannot cross structural breaks, so each
# bay -- inboard wing (between the booms), the two outboard wing panels (L/R), and the H-stab --
# carries its OWN string(s), each capped at max_cells_per_string by the cold-Voc limit. Every
# populated bay rounds its string count UP independently, so the aircraft needs MORE MPPTs (mass)
# than one aggregate string. Bay cell counts are apportioned by row-aware panel area.
def _mppt_rows(c):   # smooth floor(c/pitch): usable panel rows at chord c
    return sum(0.5 * (1 + np.tanh((c - r * _panel_pitch) * 60.0)) for r in range(1, 5))
def _bay_area(y_lo, y_hi, n=16):   # row-aware panel area over one spanwise band
    a = 0.0
    for _j in range(n):
        _y = y_lo + (y_hi - y_lo) * (_j + 0.5) / n
        a = a + _mppt_rows(chord_at(_y)) * _panel_pitch * ((y_hi - y_lo) / n)
    return a
_area_inboard = 2.0 * _bay_area(0.0, boom_y)          # center section spans -boom_y..+boom_y
_area_outboard = _bay_area(boom_y, semi_span)         # one outboard panel (L or R)
_area_hstab = _mppt_rows(0.6 * hstab_chordlen) * _panel_pitch * hstab_span
_area_total = _area_inboard + 2.0 * _area_outboard + _area_hstab + 1e-9
def _bay_strings(cells):   # smooth ceil(cells / max_cells_per_string); 0 when the bay is empty
    return sum(0.5 * (1 + np.tanh((cells - ((k - 1) * max_cells_per_string + 0.5)) * 0.5)) for k in range(1, 6))
_cells_inboard = solar_panel_n * _area_inboard / _area_total
_cells_outboard = solar_panel_n * _area_outboard / _area_total   # per side
_cells_hstab = solar_panel_n * _area_hstab / _area_total
n_mppt = _bay_strings(_cells_inboard) + 2.0 * _bay_strings(_cells_outboard) + _bay_strings(_cells_hstab)
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
mass_main_wing = main_wing.area("wetted") * 0.672  # 0.672 kg/m^2 wetted (Adam Yang's updated
# build estimate, 2026-08-17: 0.8x the old 0.84 kg/m^2 wetted; incl. encapsulation, excl. cells).
# This is the difference-maker that lets the design close at 12 kg with honest (derated) props.
# rev7: realistic tail masses (from build experience). The vertical fins are light --
# ~100 g each regardless of planform area -- so they're a fixed 0.1 kg each. The hstab
# stays area-based (it carries panels + elevator, so more skin = more mass).
mass_hstab = hor_stabilizer.area("wetted") * 0.44  # 0.44 kg/m^2 wetted (Adam's tail areal mass)
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

### PITCH INERTIA & AUTHORITY (Phase 4, ported from Adam). Lumped pitch inertia Iyy about the CG
# from the SAME component masses/positions used for the CG (parallel-axis, point masses), plus the
# wing's own chordwise spread as a conservative pad. Feeds the pitch-acceleration authority check:
# a full elevator deflection at loiter speed must produce >= 0.30 rad/s^2 (~20 deg in 1.5 s from
# rest) -- i.e. the aircraft can actually ROTATE in pitch, not merely trim.
Iyy_pitch = (
    mass_main_wing * (_x_wing - cg_le_dist) ** 2
    + mass_solar_cells * (_x_solar - cg_le_dist) ** 2
    + mass_mppt * (_x_avionics - cg_le_dist) ** 2
    + mass_empennage * (_x_tail - cg_le_dist) ** 2
    + mass_booms * (_x_boom - cg_le_dist) ** 2
    + mass_motors_mounted * (_x_motor - cg_le_dist) ** 2
    + mass_propellers * (_x_prop - cg_le_dist) ** 2
    + mass_escs * (_x_motor - cg_le_dist) ** 2
    + mass_fuselages * (_x_pod - cg_le_dist) ** 2
    + mass_avionics * (_x_avionics - cg_le_dist) ** 2
    + mass_power_board * (_x_avionics - cg_le_dist) ** 2
    + mass_wires * (_x_avionics - cg_le_dist) ** 2
    + mass_servos * (_x_wing - cg_le_dist) ** 2
    + mass_batteries * (battery_x - cg_le_dist) ** 2
    + mass_superstructures * (0.5 * chordlen - cg_le_dist) ** 2
    + mass_boom_vstab_interfaces * (_x_tail - cg_le_dist) ** 2
    + mass_vstab_hstab_interfaces * (_x_tail - cg_le_dist) ** 2
    + mass_main_wing * (0.29 * chordlen) ** 2   # wing own chordwise spread (conservative pad)
)
_q_pitch = 0.5 * operating_atm.density() * airspeed ** 2
pitch_moment_max = _q_pitch * S_w * MAC_w * Cm_de * elevator_defl_max   # N*m at full elevator
pitch_accel = pitch_moment_max / Iyy_pitch                              # rad/s^2
pitch_accel_min = 0.30                                                  # ~20 deg in 1.5 s from rest



### CONSTRAINTS
# Mission
opti.subject_to(total_mass < togw_design)

# Aerodynamic
opti.subject_to(aero["L"] >= togw_design * g)
opti.subject_to(aero["Cm"] == 0)
opti.subject_to(aero["CL"] <= CLmax_cruise)

# Propulsion
# Windmill floor (real map): keep C_T comfortably positive so the operating point stays on the
# thrust-PRODUCING side of the map, below the windmill J where the polynomial fit goes non-physical.
opti.subject_to(prop_CT >= 0.02)
opti.subject_to(v_motor_cruise <= battery_voltage)   # cruise: ESC voltage within nominal 6S bus
opti.subject_to(v_motor_climb <= bus_v_climb_max)    # climb: within fully-charged 6S bus (thrust achievable)
opti.subject_to(propeller_rpm_max >= rpm_cruise)
opti.subject_to(power_out_max >= power_cruise * propeller_n)
opti.subject_to(power_out_max >= thrust_climb * airspeed / (eta_prop_eff * esc_efficiency))  # electrical watts, same domain as cruise
# Phase 3 (ported from Adam): battery DISCHARGE power cap. The datasheet sustained discharge is
# 313 W/pack; peak draw (climb electrical) must stay within the bank's rating. Guard constraint --
# slack at cruise (~78 W) and climb (~400 W) vs 3*313 = 939 W, but binds if a future design draws
# hard. (MPPT 8 A / 210 W per-string caps are already satisfied by the 44-cell cold-Voc string
# limit above: 44 cells x ~5.9 A Imp < 8 A, and 44 x 3.4 W = 150 W < 210 W.)
pack_discharge_max_w = 313.0   # W/pack, datasheet sustained discharge
opti.subject_to(power_out_max <= pack_discharge_max_w * num_packs)

# Stall-speed margin (equivalent to constraining required CL below CLmax at cruise).
# Binding floor is the banked-loiter factor (1.3*Vs), which dominates the 1.2 cruise margin.
V_stall = np.sqrt(2 * (togw_design * g) / (operating_atm.density() * S_w * CLmax_cruise))
opti.subject_to(airspeed >= loiter_stall_factor * V_stall)   # banked-loiter stall margin
opti.subject_to(airspeed >= crosswind_floor_ms)              # 15 kt crosswind headway floor
opti.subject_to(wing_airfoil.max_thickness() * chord_root >= 0.025)  # root must accomodate batteries (20mm)
# ROW-AWARE solar panel model (honest about taper): at each spanwise station the number of
# chordwise panel ROWS that fit is floor(local_chord / panel_pitch), so a full-chord inboard
# station carries 2 rows while a narrow tapered tip carries 1 or 0. We integrate the panel-able
# area across the half-span with a smooth row-counter (sum of soft steps), so taper is paid for
# in lost rows rather than by a rigid rectangular strip.
panel_pitch = 0.13            # m, 12.5 cm cell + gap
# Phase 6 -- packing density. The row-aware area already floors chordwise rows at panel_pitch;
# this factor is the REMAINING loss on top. Adam's reference tessellates floor(chord/pitch) x
# floor(span/pitch) with NO extra knockdown (effective ~1.0), so the old 0.85 was an extra
# conservatism not in the reference. Decomposed physically: spanwise cell fill 0.125/0.13 = 0.962
# x ~0.96 spar-rail + TE-curvature allowance = 0.92. Still below Adam's implicit 1.0, but ~8%
# denser than before -> more cells per wing area -> closes at shorter span.
panel_packing_eff = 0.92      # was 0.85; matches Adam's demonstrated density minus a spar/TE margin
_panel_rows_max = 4           # most rows the widest plausible chord could hold

def _rows_at(c):
    # Number of chordwise cell rows that PHYSICALLY fit in a band of width c. A row fits once the
    # band clears one cell WIDTH (solar_panel_side_length = 0.125 m), not a full pitch; each extra
    # row then needs +panel_pitch. This counts a band that just fits a cell as a full row (e.g. the
    # H-stab's 0.132 m band = 1 row, ~10 cells across span) instead of the old ~0.5 at the boundary.
    return sum(
        0.5 * (1 + np.tanh((c - (solar_panel_side_length + (r - 1) * panel_pitch)) * 120.0))
        for r in range(1, _panel_rows_max + 1)
    )

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
# H-stab chordwise rows: a SHARP count (a full row once the cell-able band clears one cell width),
# so the H-stab carries a clean 1-row strip of cells filling its ENTIRE span -- matching the
# viewer's floor(span/pitch) columns -- instead of the ~0.8-row under-count the wing's smooth model
# gave right at the boundary. hstab chord is derived (not a design lever), so the sharp step is safe.
_hstab_rows = sum(0.5 * (1 + np.tanh((hstab_panel_frac * hstab_chordlen - (solar_panel_side_length + (r - 1) * panel_pitch)) * 300.0)) for r in range(1, 3))
hstab_panel_area = _hstab_rows * panel_pitch * hstab_span
wing_panel_capacity = 2 * _panel_area_available * panel_packing_eff / panel_pitch**2   # cells the wing can hold
# H-stab holds a clean SINGLE-ROW strip: floor(span / pitch) cells per row (no 2-D packing derate,
# since a single row butts cell-to-cell along the span). A 1.5 m span at 0.13 m pitch => ~11 cells.
hstab_panel_capacity = _hstab_rows * hstab_span / panel_pitch                          # cells the hstab can hold
usable_panel_area = (2 * _panel_area_available + hstab_panel_area) * panel_packing_eff
opti.subject_to(solar_panel_n * panel_pitch**2 <= usable_panel_area)  # panels must physically fit (row-aware)
# Split the solved cell count: fill the wing first, overflow onto the hstab (reported numbers).
# Fill the wing first; overflow above wing capacity goes to the hstab, capped at what the hstab
# can physically hold (hstab_panel_capacity, which is ~0 when the hstab chord is too short for a
# panel row). Report-only, so hard min/max is fine. Viewer reads these exact numbers.
# Fill the H-stab FIRST (up to its capacity) so it is actively used across its full span, then put
# the remaining cells on the wing. (Was wing-first-overflow, which left the H-stab nearly bare.)
panels_on_hstab = np.minimum(solar_panel_n, hstab_panel_capacity)
panels_on_wing = solar_panel_n - panels_on_hstab
# Root chord must fit at least min_panel_rows inboard rows (plus a little for the aileron behind them).
opti.subject_to(chordlen >= min_panel_rows * panel_pitch + 0.03)
# Chordwise cell-row cap. A thin (<=2-row) wing packs cells end-to-end and is aerodynamically
# cleaner, BUT on the harder April day with only 3 packs at a 6 m span there isn't enough area for
# it to close -- the wing needs a fatter multi-row chord. So the cap is DISABLED here (set
# max_panel_rows small and re-enable the constraint only if you fly an easier day / accept >6 m span).
max_panel_rows = 3   # IPOPT converges here; thin needs DE (see notes)
opti.subject_to(chordlen <= (max_panel_rows + 1) * panel_pitch - 0.005)   # < 4 pitches => at most 3 rows
# Honest actual chordwise row count at the root, for the report/viewer (real value, not the floor).
solar_panel_n_rows = np.floor(chordlen / panel_pitch)
# Taper must leave the tip at least one panel-row wide (keeps the taper physical).
opti.subject_to(chord_tip >= panel_pitch)
opti.subject_to(hstab_chordlen >= 0.22)  # >=22cm so the 60% panel band (>=0.13m) fits a full row of cells
opti.subject_to(vstab_root_chord >= hstab_chordlen)  # enforce an actual taper (root >= tip) and avoid negative root chord

# Stability — longitudinal
# Minimum tail arm. The tail-VOLUME sizing (hstab_area = V_H*S_w*MAC/L_H) silently assumes a
# sensible arm; without this floor, min-mass shrinks the boom to a stub (~2 MAC) and pays for it
# with an oversized H-stab (cheap boom vs 1/L_H tail area) -- unrealistic (downwash interference,
# weak pitch damping, structural). Real sailplanes run ~3-4 MAC. Floor at 3 MAC for honest,
# slender proportions (Adam-like): smaller tail, and typically a smaller wing too.
TAIL_ARM_MIN_MAC = 3.0
opti.subject_to(L_H >= TAIL_ARM_MIN_MAC * MAC_w)
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
# Phase 4: pitch-acceleration authority -- full elevator at loiter speed must rotate the aircraft
# at >= 0.30 rad/s^2 (needs the Iyy inertia above). Checks maneuver capability, not just trim.
opti.subject_to(pitch_accel >= pitch_accel_min)

# Power (SOC window on DELIVERABLE energy)
soc_floor_wh = battery_capacity_effective * (1 - allowable_battery_depth_of_discharge)  # 20% SOC, in Wh
opti.subject_to(battery_states >= soc_floor_wh)
opti.subject_to(battery_states <= battery_capacity_effective)
# Two-day periodic gate (Phase 5): on a REPEATING design day, the cycle must not ratchet down --
# the SOC at the end of the 24 h march must be >= the start (start <= end). Combined with the hard
# 20% floor across every timestep, this is exactly Adam's "day-2 must close" test: a design that
# survives its first night from full but loses a little each day is rejected, because at steady
# state (start == end) it would eventually dip below the floor. One periodic cycle == steady state.
opti.subject_to(battery_states[0] <= battery_states[N-1])
# ENERGY MARGIN (rev7, emulating Adam Yang): deliverable energy above the 20% floor at the WORST
# point of the cycle. Instead of the hard max-min EPIGRAPH (energy_margin <= state - floor for
# every i), which makes IPOPT stall -- at a tight design the overnight SOC sits flat at the floor,
# so dozens of those constraints go active at once (a degenerate active set) -- we use a smooth
# SOFT-MIN of the SOC-above-floor trace. It has no auxiliary variable and no active-set degeneracy,
# just a differentiable objective, so the gradient solver converges even at the knife edge. The
# soft-min slightly UNDER-reports the true min (conservative); higher beta -> closer to the hard min.
_softmin_beta = 6.0
_soc_above = (battery_states - soc_floor_wh) / 100.0   # Wh/100, O(1) for numerical conditioning
# softmin(x) = -1/b * log(mean(exp(-b*x))). x >= 0 (hard floor holds), so exp(-b*x) in (0,1] -- safe.
energy_margin_wh = (-1.0 / _softmin_beta) * np.log(np.sum(np.exp(-_softmin_beta * _soc_above)) / N) * 100.0
# NOTE: no explicit daily-refill constraint. Minimizing span (hence solar area + battery mass)
# already drives the smallest pack that is fully cycled -- min SOC floors at 0.20 and the peak
# naturally reaches ~0.95, so the pack does not carry dead capacity. A fixed-index refill target
# only mis-pinned a morning point and inflated span with no physical benefit.
opti.subject_to(num_packs <= 3)   # user: limit battery to 3 packs
# Optional integer-pack pin for a discrete-reality solve (env FIXED_PACKS=2/3/4...).
if fixed_packs is not None:
    opti.subject_to(battery_capacity == float(fixed_packs) * pack_energy_Wh)
    print(f"[packs] pinned to {fixed_packs} packs = {float(fixed_packs)*pack_energy_Wh:.0f} Wh nameplate")



### SOLVE
# Objective: minimize mass (the design is feasibility-limited once the battery is bounded, so it
# just closes at ~0 margin; the lightest such design is the right target).
opti.minimize(-energy_margin_wh)   # MAXIMIZE energy margin (smooth soft-min -- see above)

try:
    sol = opti.solve(max_iter=5000)
    # WHOLE-PACK ROUND-UP: you cannot fly a fractional battery pack, so carry the NEXT integer
    # number of packs above what the continuous solve wants, then re-solve with that pinned so the
    # extra pack's mass (and the slack energy it buys) are properly accounted. (Skipped if the pack
    # count was already pinned via FIXED_PACKS.)
    if fixed_packs is None:
        _npacks_cont = float(sol.value(battery_capacity)) / pack_energy_Wh
        _n_int = int(onp.ceil(_npacks_cont - 1e-6))
        print(f"[packs] continuous {_npacks_cont:.3f} -> carrying {_n_int} whole packs; re-solving pinned")
        opti.set_initial(opti.x, sol.value(opti.x))                 # warm-start from the continuous solve
        opti.subject_to(battery_capacity == _n_int * pack_energy_Wh)
        sol = opti.solve(max_iter=5000)
    _bs = sol.value(battery_states); _cap = sol.value(battery_capacity_effective)
    print(f"[SOC] min={_bs.min()/_cap:.3f}  peak={_bs.max()/_cap:.3f} (i={int(_bs.argmax())})  dawn(i={i_dawn})={_bs[i_dawn]/_cap:.3f}  start={_bs[0]/_cap:.3f} end={_bs[-1]/_cap:.3f}")
    _margin_wh = float(_bs.min() - 0.20 * _cap)   # achieved margin above the 20% floor (Wh)
    print(f"[MARGIN] energy margin above 20% floor = {_margin_wh:.1f} Wh  (span={sol.value(wingspan):.3f} m, mass={sol.value(total_mass):.3f} kg)")
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
        "climb_rate_req": climb_rate_req,
        "climb_angle_at_loiter_deg": climb_angle_rad * 180 / np.pi,
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
        "Iyy_pitch": Iyy_pitch,           # lumped pitch inertia about CG (kg*m^2), Phase 4
        "pitch_accel_rad_s2": pitch_accel,  # full-elevator pitch authority (>= 0.30 floor)
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
        "fuselage_radius": fuselage_radius,
    },
    "Power": {
        "solar_panel_n": solar_panel_n,
        "panels_on_wing": panels_on_wing,
        "panels_on_hstab": panels_on_hstab,
        "wing_panel_capacity": wing_panel_capacity,
        "hstab_panel_capacity": hstab_panel_capacity,
        "n_mppt": n_mppt,
        "max_cells_per_string": max_cells_per_string,
        "energy_margin_wh": energy_margin_wh,
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
        "eta_prop_eff": eta_prop_eff,
        "prop_thrust_derate": prop_thrust_derate,
        "prop_power_inflate": prop_power_inflate,
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
    # Bake the airspeed thrust-lapse curve into the artifact (separation of concerns):
    # the prop physics stays here, on the design side, and the sim reads only this table.
    _thrust_curve = thrust_curve_from_soln(soln)
    if _thrust_curve is not None:
        soln.setdefault("Propulsion", {})["thrust_curve"] = _thrust_curve
    # Raw C_T(J)/C_P(J) prop map for the sim's real FGPropeller (higher fidelity).
    _prop_map = prop_coeff_table(soln)
    if _prop_map is not None:
        soln.setdefault("Propulsion", {})["prop_map"] = _prop_map
    write_json(run_dir / "soln.json", soln, indent=2)
    
    # Compute solved airplane once for all exports.
    airplane_sol = sol(airplane) if callable(sol) else airplane
    
    # Export an XFLR5 XML (synthesizes a single fin for the twin-fin geometry).
    try:
        export_xflr5_xml_from_soln(airplane_sol=airplane_sol, soln=soln, out_path=run_dir / "aircraft.xml")
    except Exception as e:
        print(f"[export_xflr5_xml_from_soln] skipped due to error: {e}")

    # Export the authoritative CAD/geometry mesh (twin fins + pods + booms) from the
    # SOLVED airplane, so the sim renders exactly this run's airframe.
    try:
        export_airplane_glb(airplane_sol, run_dir / "geometry.glb")
    except Exception as e:
        print(f"[export_airplane_glb] skipped due to error: {e}")

    # Draw airplane (needs the optional `pyvista` viewer; skip cleanly if absent).
    try:
        airplane_sol.draw()
    except Exception as e:
        print(f"[draw] skipped ({e})")