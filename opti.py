import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.library import power_solar, propulsion_electric, propulsion_propeller
from aerosandbox.atmosphere.atmosphere import Atmosphere
import numpy as onp
from pathlib import Path
from datetime import datetime, timezone

from lib.artifacts import process_raw_values, run_id_random, write_json
from lib.exports import export_xflr5_xml_from_soln, export_cadquery_step
from lib.models import apc_prop, hacker_motor_mass, hacker_motor_resistance



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
wing_airfoil = asb.Airfoil("s4110")
tail_airfoil = asb.Airfoil("naca0010")

# Main wing
polyhedral_angle = 3

# Tail sizing (via Tail Volume Coefficients)
V_H = 0.40 # Horizontal tail volume coefficient.
V_V = 0.015 # Vertical tail volume coefficient

# Static margin constraints (dimensionless, SM = (x_np - x_cg) / MAC)
SM_min = 0.05
SM_max = 0.5

# Stall / low-Re realism (AeroBuildup does not model stall)
CLmax_cruise = 1.25 # Obtained from S4110 airfoil data.
stall_speed_margin = 1.2 # V_cruise >= margin * V_stall

## Structural
structural_mass_markup = 1.2
boom_radius = 0.01 # radius of the boom in meters.
boom_spacing_frac = 0.25 # Fraction that the inboard wing spans.

## Propulsion
propeller_n = 2 # Number of motors / propellers.
advanced_ratio = 0.8 # ! This is a guess. Need to test to get actual values.

## Power
battery_voltage = 22.2
N = 180 # Number of discretization points.
time = onp.linspace(0, 24 * 60 * 60, N) # s (numeric; does not need to be symbolic)
dt = onp.diff(time)[0]  # s
solar_panel_side_length = 0.125 # m (12.5cm)
solar_panel_n_rows = 2
solar_encapsulation_eff_hit = 0.1 # Estimated 10% efficieincy loss from encapsulation.
solar_cell_efficiency = 0.243 * (1 - solar_encapsulation_eff_hit) # Efficiency of the solar cells.
energy_generation_margin = 1.05 # Losses in energy generation.
allowable_battery_depth_of_discharge = 0.85  # How much of the battery can you actually use?
pack_energy_Wh = 5 * 6 * 3.7  # 5 Ah, 6 cells, 3.7 V/cell

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
print(
    f"[solar_flux_profile] min={solar_flux_profile.min():.3f} W/m^2, "
    f"max={solar_flux_profile.max():.3f} W/m^2, "
    f"nan_count={onp.isnan(solar_flux_profile).sum()}"
)



### VARIABLES
## Performance
airspeed = opti.variable(init_guess=10, lower_bound=5, upper_bound=15, scale=5, category="airspeed")
togw_design = opti.variable(init_guess=8, lower_bound=1, upper_bound=togw_max, scale=1, category="togw_max")
power_out_max = opti.variable(init_guess=500, lower_bound=25*16, scale=100, category="power_out_max")

## Propulsion
# propeller_diameter = opti.variable(init_guess=0.5, lower_bound=0.1, upper_bound=0.8, scale=0.1, category="propeller_diameter")
propeller_diameter = opti.parameter(value=0.4)
motor_kv = opti.variable(init_guess=500, lower_bound=50, upper_bound=2000, scale=100, category="motor_kv")

## Avionics
solar_panel_n = opti.variable(init_guess=50, lower_bound=10, scale=10, category="solar_panel_n")
battery_capacity = opti.variable(init_guess=450, lower_bound=100, scale=100, category="battery_capacity") # initial battery energy in Wh
battery_states = opti.variable(n_vars=N, init_guess=500, scale=100, category="battery_states")

## Aerodynamics
# Main wing
wingspan = opti.variable(init_guess=6, lower_bound=2, upper_bound=8, scale=1, category="wingspan")
chordlen = opti.variable(init_guess=0.3, lower_bound=0.05, scale=0.1, category="chordlen")
struct_defined_aoa = opti.variable(init_guess=2, lower_bound=0, upper_bound=10, scale=1, category="struct_aoa")
cg_le_dist = 0.25 * chordlen # CG assumed at quarter-chord of main wing
# alpha_cruise = opti.variable(init_guess=3, lower_bound=-2, upper_bound=10, scale=5, category="alpha_cruise")

# Empennage (sized by tail volume coeffs)
hstab_AR = opti.variable(init_guess=6.0, lower_bound=2.0, upper_bound=20.0, scale=1, category="hstab_AR")
vstab_AR = 2.0
hstab_aoa = opti.variable(init_guess=-5, lower_bound=-10, upper_bound=5, scale=1, category="hstab_aoa")

# Structural
boom_length = opti.variable(init_guess=1, lower_bound=0.1, upper_bound=2, scale=1, category="boom_length")
boom_y = 0.5 * boom_spacing_frac * wingspan



### GEOMETRIES
main_wing = asb.Wing(
    name="Main Wing",
    symmetric=True,  # Should this wing be mirrored across the XZ plane?
    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            xyz_le=[
                0,
                0,
                0,
            ], # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            chord=chordlen,
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[0.00, boom_y, 0],
            chord=chordlen,
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[0.00, boom_y+wingspan*boom_spacing_frac/2, np.sind(polyhedral_angle)*(boom_y+wingspan*boom_spacing_frac/2)],
            chord=chordlen,
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
        asb.WingXSec(  # Tip
            xyz_le=[0.00, wingspan/2, np.sind(polyhedral_angle)*(wingspan*(1-boom_spacing_frac)/2)],
            chord=chordlen,
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
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

vert_stabilizer_L = asb.Wing(
    name="Vertical Stabilizer L",
    symmetric=False,
    xsecs=[
        asb.WingXSec(
            xyz_le=[0, 0, 0],
            chord=vstab_root_chord,
            twist=0,
            airfoil=tail_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[vstab_root_chord/4, 0, vstab_span],
            chord=hstab_chordlen,
            twist=0,
            airfoil=tail_airfoil,
        ),
    ],
).translate([boom_length, boom_y, 0])

vert_stabilizer_R = asb.Wing(
    name="Vertical Stabilizer R",
    symmetric=False,
    xsecs=[
        asb.WingXSec(
            xyz_le=[0, 0, 0],
            chord=vstab_root_chord,
            twist=0,
            airfoil=tail_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[vstab_root_chord/4, 0, vstab_span],
            chord=hstab_chordlen,
            twist=0,
            airfoil=tail_airfoil,
        ),
    ],
).translate([boom_length, -boom_y, 0])

# Booms
boom_L = asb.Fuselage(
    name="Boom L",
    xsecs=[
        asb.FuselageXSec(xyz_c=[boom_length * xi, boom_y, -0.02], radius=boom_radius)
        for xi in np.cosspace(0, 1, 20)
    ],
)
boom_R = asb.Fuselage(
    name="Boom R",
    xsecs=[
        asb.FuselageXSec(xyz_c=[boom_length * xi, -boom_y, -0.02], radius=boom_radius)
        for xi in np.cosspace(0, 1, 20)
    ],
)


left_pod = asb.Fuselage(  # left pod fuselage
    name="Left Fuse",
    xsecs=[
        asb.FuselageXSec(
            xyz_c=[0.5 * xi - 0.5, boom_y, -0.02],
            radius=0.4
            * asb.Airfoil("dae51").local_thickness(
                x_over_c=xi
            ),  # half a meter fuselage. Starting at LE and 0.5m forward
        )
        for xi in np.cosspace(0, 1, 30)
    ],
)

right_pod = asb.Fuselage(  # right pod fuselage
    name="Right Fuselage",
    xsecs=[
        asb.FuselageXSec(
            xyz_c=[0.5 * xi - 0.5, -boom_y, -0.02],
            radius=0.4
            * asb.Airfoil("dae51").local_thickness(
                x_over_c=xi
            ),  # half a meter fuselage. Starting at LE and 0.5m forward
        )
        for xi in np.cosspace(0, 1, 30)
    ],
)

# Model the airplane.
airplane = asb.Airplane(
    name="Rev 6",
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
        # alpha=alpha_cruise,  # deg (trim variable)
    ),
)
aero = vlm.run_with_stability_derivatives(
    alpha=True,
    beta=True,
    p=False,
    q=False,
    r=False
)



### STABILITY
static_margin = (aero["x_np"] - cg_le_dist) / main_wing.mean_aerodynamic_chord()



### PROPULSION
power_shaft_cruise = propulsion_propeller.propeller_shaft_power_from_thrust(
    # Guard against occasional negative drag during intermediate iterates
    # (can produce negative thrust/power and NaNs in downstream power-law models).
    thrust_force=aero["D"] / propeller_n,
    area_propulsive=np.pi / 4 * propeller_diameter ** 2,
    airspeed=airspeed,
    rho=operating_atm.density(),
    propeller_coefficient_of_performance=0.80  # calibrated to QProp output with Dongjoon
)

# Get cruise power.
rpm_cruise = 60 * airspeed / (advanced_ratio * propeller_diameter)
Q = power_shaft_cruise / (2 * np.pi * rpm_cruise / 60)
Kt = 60 / (2 * np.pi * motor_kv) 
i_cruise = Q / Kt
# The resistance scaling uses max_power**(-0.68); ensure it's strictly positive to avoid NaNs.
motor_resistance = hacker_motor_resistance(power_shaft_cruise * 2.5, motor_kv)
power_cruise = (power_shaft_cruise + (i_cruise**2) * motor_resistance) / 0.95 # 95% ESC efficiency.

# Propeller tip mach limit.
propeller_tip_mach = 0.7  # From Dongjoon, 4/30/20
propeller_rads_per_sec = propeller_tip_mach * operating_atm.speed_of_sound() / (propeller_diameter / 2)
propeller_rpm_max = propeller_rads_per_sec * 30 / np.pi

# Climb.
thrust_climb = togw_design * g * np.sin(min_climb_angle * np.pi / 180) + aero["D"]



### POWER
for i in range(N-1):
    solar_flux = solar_flux_profile[i]  # W / m^2 (numeric constant)
    solar_area = solar_panel_n * solar_panel_side_length**2 # m^2
    power_generated = solar_flux * solar_area * solar_cell_efficiency / energy_generation_margin
    power_used = (power_cruise*propeller_n + 8 + 1)  # 8W avionics, 1W for NavLights
    net_energy = (power_generated - power_used) * (dt / 3600)  # Wh

    battery_update = np.softmin(battery_states[i] + net_energy, battery_capacity, hardness=10)
    opti.subject_to(battery_states[i+1] == battery_update)



### Mass
# Power
mass_solar_cells = 0.0075 * solar_panel_n
mass_power_board = 0.050 * 2 # 75g estimate for each power board
mass_batteries = propulsion_electric.mass_battery_pack(
    battery_capacity_Wh=battery_capacity,
    battery_pack_cell_fraction=0.95
)
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
mass_avionics = mass_fc + mass_gps + mass_telemtry + mass_receiver + mass_power_board + mass_navlights + mass_pitot

# Actuators
mass_servos = .02 * 4 # 20g a servo

# Propulsion
# mass_motor_raw = hacker_motor_mass(
#     max_power=power_out_max,
#     kv=motor_kv
# )
mass_motor_raw = 0.275
mass_motors_mounted = mass_motor_raw * 1.1 * propeller_n # marign to account for motor mounts.
# mass_esc = propulsion_electric.mass_ESC(max_power=power_out_max)
mass_escs = 0.030 * propeller_n # Single ESC.
# mass_propellers = apc_prop(propeller_diameter)
mass_propellers = 0.055 * propeller_n # Single propeller.

# Structures
mass_main_wing =  main_wing.area("wetted")*1.2 * 0.70 # 700g per meter squared of planform area
mass_hstab = hor_stabilizer.area("wetted")*1.10 * 0.4 # 400g per meter squared of wetted area
mass_vstabs = (vert_stabilizer_L.area("wetted") + vert_stabilizer_R.area("wetted"))*1.10 * 0.4 # 400g per meter squared of wetted area
mass_booms = 2 * 0.09 * boom_length  # kg (90 g/m per boom)
mass_fuselages = 0.20 # rough estimate
mass_superstructures = 0.2 # Two superstructures.
mass_boom_vstab_interfaces = 0.1 # Two boom-vstab interfaces.
mass_vstab_hstab_interfaces = 0.06 # Two vstab-stab interfaces.

## Total
total_mass = (
    mass_solar_cells + 
    mass_power_board + 
    mass_batteries + 
    mass_wires + 
    mass_avionics + 
    mass_servos + 
    mass_motors_mounted + 
    mass_escs + 
    mass_propellers + 
    mass_main_wing + 
    mass_hstab + 
    mass_vstabs + 
    mass_booms + 
    mass_fuselages + 
    mass_superstructures + 
    mass_boom_vstab_interfaces + 
    mass_vstab_hstab_interfaces 
)



### CONSTRAINTS
# Mission
opti.subject_to(total_mass < togw_design)

# Aerodynamic
opti.subject_to(aero["L"] >= togw_design * g)
opti.subject_to(aero["Cm"] == 0)
opti.subject_to(aero["CL"] <= CLmax_cruise)

# Propulsion
opti.subject_to(rpm_cruise <= motor_kv * battery_voltage)  # Motor no-load speed limit
opti.subject_to(propeller_rpm_max >= rpm_cruise)
opti.subject_to(power_out_max >= power_cruise * propeller_n)
opti.subject_to(power_out_max >= thrust_climb * airspeed)

# Stall-speed margin (equivalent to constraining required CL below CLmax at cruise)
V_stall = np.sqrt(2 * (togw_design * g) / (operating_atm.density() * S_w * CLmax_cruise))
opti.subject_to(airspeed >= stall_speed_margin * V_stall)
opti.subject_to(chordlen >= solar_panel_n_rows * 0.13 + 0.05) # ! Justify this addition of 0.1
opti.subject_to(wing_airfoil.max_thickness() * chordlen >= 0.025)  # must accomodate batteries (20mm)
opti.subject_to(wingspan >= 0.13 * solar_panel_n / solar_panel_n_rows)  # Must be able to fit all of our solar panels 13cm each
opti.subject_to(hstab_chordlen >= 0.135) # 13.5cm to make room for panels, elevator, and manufacturing
opti.subject_to(vstab_root_chord >= hstab_chordlen)  # enforce an actual taper (root >= tip) and avoid negative root chord

# Stability
opti.subject_to(L_H >= 1e-3)
opti.subject_to(L_V >= 1e-3)
opti.subject_to(boom_length >= chordlen)
opti.subject_to(static_margin >= SM_min)
opti.subject_to(static_margin <= SM_max)

# Power
opti.subject_to(battery_states >= battery_capacity * (1-allowable_battery_depth_of_discharge))
opti.subject_to(battery_states <= battery_capacity)
opti.subject_to(battery_states[0] <= battery_states[N-1])
opti.subject_to(num_packs <= 8)



### SOLVE
opti.minimize(wingspan)

try:
    sol = opti.solve(max_iter=5000)
except Exception as e:
    print(f"[solve] failed: {e}")
    print("\n[DEBUG] Variable values at failure:")
    print(f"  airspeed: {opti.debug.value(airspeed):.3f} m/s")
    # print(f"  alpha_cruise: {opti.debug.value(alpha_cruise):.3f} deg")
    print(f"  togw_design: {opti.debug.value(togw_design):.3f} kg")
    print(f"  power_out_max: {opti.debug.value(power_out_max):.3f} W")
    print(f"  solar_panel_n: {opti.debug.value(solar_panel_n):.3f}")
    print(f"  battery_capacity: {opti.debug.value(battery_capacity):.3f} Wh")

    battery_states_vals = opti.debug.value(battery_states)
    print(f"  battery_states min: {float(onp.min(battery_states_vals)):.3f} Wh")
    print(f"  battery_states max: {float(onp.max(battery_states_vals)):.3f} Wh")
    print(f"  battery_states[0]: {float(battery_states_vals[0]):.3f} Wh")
    print(f"  battery_states[-1]: {float(battery_states_vals[-1]):.3f} Wh")

    print(f"  wingspan: {opti.debug.value(wingspan):.3f} m")
    print(f"  chordlen: {opti.debug.value(chordlen):.3f} m")
    print(f"  struct_defined_aoa: {opti.debug.value(struct_defined_aoa):.3f} deg")
    print(f"  hstab_AR: {opti.debug.value(hstab_AR):.3f}")
    print(f"  hstab_aoa: {opti.debug.value(hstab_aoa):.3f} deg")
    print(f"  boom_length: {opti.debug.value(boom_length):.3f} m")

    print(f"  num_packs: {opti.debug.value(num_packs):.3f}")
    print(f"  cg_le_dist: {opti.debug.value(cg_le_dist):.3f} m")
    print(f"  boom_y: {opti.debug.value(boom_y):.3f} m")
    print(f"  total_mass: {opti.debug.value(total_mass):.3f} kg")
    print(f"  static_margin: {opti.debug.value(static_margin):.3f}")
    print(f"  L_H: {opti.debug.value(L_H):.3f} m")
    print(f"  L_V: {opti.debug.value(L_V):.3f} m")
    print(f"  mass_wires: {opti.debug.value(mass_wires):.3f} kg")
    print(f"  power_cruise: {opti.debug.value(power_cruise):.3f} W")
    print(f"  power_shaft_cruise: {opti.debug.value(power_shaft_cruise):.3f} W")
    print(f"  aero Cm: {opti.debug.value(aero['Cm']):.6f}")
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
        # "alpha_cruise": alpha_cruise,
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
    },
    "Main Wing": {
        "wingspan": wingspan,
        "chordlen": chordlen,
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
    },
    "Power": {
        "solar_panel_n": solar_panel_n,
        "battery_capacity": battery_capacity,
        "battery_states": battery_states,
        "battery_voltage": battery_voltage,
        "solar_panel_side_length": solar_panel_side_length,
        "solar_panel_n_rows": solar_panel_n_rows,
        "solar_encapsulation_eff_hit": solar_encapsulation_eff_hit,
        "solar_cell_efficiency": solar_cell_efficiency,
        "energy_generation_margin": energy_generation_margin,
        "num_packs": num_packs,
    },
    "Propulsion": {
        "propeller_n": propeller_n,
        "propeller_diameter": propeller_diameter,
        "advanced_ratio": advanced_ratio,
        "motor_kv": motor_kv,
        "rpm_cruise": rpm_cruise,
        "i_cruise": i_cruise
    },
    "Masses": {
        "mass_solar_cells": mass_solar_cells,
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
    
    # # Export CAD as STEP file.
    # try:
    #     export_cadquery_step(airplane_sol=airplane_sol, out_path=run_dir / "airplane.step")
    # except Exception as e:
    #     print(f"[export_cadquery_step] skipped due to error: {e}")
    
    # # Save AeroSandbox airplane object.
    # try:
    #     airplane_sol.save(filename=str(run_dir / "airplane.aero"))
    # except Exception as e:
    #     print(f"[airplane.save] skipped due to error: {e}")
    
    # print(f"[artifacts] wrote: {run_dir / 'soln.json'}")

    # Draw airplane airplane_sol.draw()
    airplane_sol.draw() 