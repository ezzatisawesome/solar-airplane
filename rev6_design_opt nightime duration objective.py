import aerosandbox as asb
import aerosandbox.numpy as np
from lib.sun import calculate_sun_angles
from lib.aero import calculate_skin_friction
from datetime import datetime, timedelta


opti = asb.Opti()


### VARIABLES ###
#-- Aerodynamic variables --#
airspeed = opti.variable(init_guess=15, lower_bound=5, upper_bound=40)
wingspan = opti.variable(init_guess=3.5, lower_bound=2.5, upper_bound=6)
chordlen = opti.variable(init_guess=0.26)
cg_le_dist = opti.variable(init_guess=0, lower_bound=0)
struct_defined_aoa = opti.variable(init_guess=2, lower_bound=0, upper_bound=3)
hstab_span = opti.variable(init_guess=0.8, lower_bound=0.8, upper_bound=1.2)
hstab_chordlen = opti.variable(init_guess=0.2, lower_bound=0.15, upper_bound=0.4)
hstab_aoa = opti.variable(init_guess=-3, lower_bound=-8, upper_bound=-2)
boom_length = opti.variable(init_guess=2, lower_bound=1.5, upper_bound=4)

#-- Power --#
#nightime = opti.variable(init_guess=10, lower_bound=5)
battery_capacity = opti.variable(init_guess=700, lower_bound=111)
n_solar_panels = opti.variable(init_guess=20, lower_bound=10)


### CONSTANTS ###
#-- Airframe --#
vstab_span = 0.3
vstab_chordlen = 0.15
polyhedral_angle = 10

#-- Aerodynamic --#
wing_airfoil = asb.Airfoil("goe447")
tail_airfoil = asb.Airfoil("naca0010")

#-- Power --#


#-- Mission constants --#
mission_date = date_time = datetime(2025, 5, 15, 0)
lat = 37.398928
long = -122.152097


### PLANE GEOMETRY ###
main_wing = asb.Wing(
    name="Main Wing",
    symmetric=True,  # Should this wing be mirrored across the XZ plane?
    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            xyz_le=[
                0,
                0,
                0,
            ],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            chord=chordlen,
            twist=struct_defined_aoa,  # degrees
            airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
        ),
        asb.WingXSec(  # Mid
            xyz_le=[0.00, 0.5 * wingspan / 2, 0],
            chord=chordlen,
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
        asb.WingXSec(  # Tip
            xyz_le=[0.00, wingspan / 2, np.sin(10 * np.pi / 180) * 0.5 * wingspan / 2],
            chord=0.125,
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
    ],
)

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
).translate([boom_length, 0, 0])

vert_stabilizer = asb.Wing(
    name="Vertical Stabilizer",
    symmetric=False,
    xsecs=[
        asb.WingXSec(
            xyz_le=[0, 0, 0],
            chord=vstab_chordlen,
            twist=0,
            airfoil=tail_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[0.00, 0, vstab_span],
            chord=vstab_chordlen,
            twist=0,
            airfoil=tail_airfoil,
        ),
    ],
).translate([boom_length + hstab_chordlen, 0, 0])

main_fuselage = asb.Fuselage(  # main fuselage
    name="Fuselage",
    xsecs=[
        asb.FuselageXSec(
            xyz_c=[0.5 * xi, 0, 0],
            radius=0.6
            * asb.Airfoil("dae51").local_thickness(
                x_over_c=xi
            ),  # half a meter fuselage. Starting at LE and 0.5m forward
        )
        for xi in np.cosspace(0, 1, 30)
    ],
).translate([-0.5, 0, 0])

left_pod = asb.Fuselage(  # left pod fuselage
    name="Fuselage",
    xsecs=[
        asb.FuselageXSec(
            xyz_c=[0.2 * xi, 0.75, -0.02],
            radius=0.4
            * asb.Airfoil("dae51").local_thickness(
                x_over_c=xi
            ),  # half a meter fuselage. Starting at LE and 0.5m forward
        )
        for xi in np.cosspace(0, 1, 30)
    ],
)

right_pod = asb.Fuselage(  # right pod fuselage
    name="Fuselage",
    xsecs=[
        asb.FuselageXSec(
            xyz_c=[0.2 * xi, -0.75, -0.02],
            radius=0.4
            * asb.Airfoil("dae51").local_thickness(
                x_over_c=xi
            ),  # half a meter fuselage. Starting at LE and 0.5m forward
        )
        for xi in np.cosspace(0, 1, 30)
    ],
)

airplane = asb.Airplane(
    name="rev 6",
    xyz_ref=[cg_le_dist, 0, 0],  # CG location
    wings=[main_wing, hor_stabilizer, vert_stabilizer],
    fuselages=[main_fuselage, left_pod, right_pod],
)


### AIRFRAME SIMULATION ###
#-- Run Vortex Lattice Method --# 
vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=airspeed,  # m/s
    ),
)
aero = vlm.run_with_stability_derivatives()  # Returns a dictionary

# VLM does not calcualte parasitic drag, we must add this manually
CD0 = (
    calculate_skin_friction(chordlen, airspeed) * main_wing.area(type="wetted") / main_wing.area()
    + calculate_skin_friction(hstab_chordlen, airspeed)
    * hor_stabilizer.area(type="wetted")
    / main_wing.area()
    + calculate_skin_friction(vstab_chordlen, airspeed)
    * vert_stabilizer.area(type="wetted")
    / main_wing.area()
    + calculate_skin_friction(0.5, airspeed) * main_fuselage.area_wetted() / main_wing.area()
    + 2 * calculate_skin_friction(0.2, airspeed) * left_pod.area_wetted() / main_wing.area()
)

drag_parasite = 0.5 * 1.29 * airspeed**2 * main_wing.area() * CD0

aero["CD_tot"] = aero["CD"] + CD0
aero["D_tot"] = aero["D"] + drag_parasite
aero["power"] = (aero["D_tot"] * airspeed + 8) * 1.4  # 8w to run avionics 40% surplus for non ideal conditions and charging


### POWER ###
time_delta = 600  # seconds
N_steps = int(24 * 60 * 60 // time_delta)

battery_states = opti.variable(n_vars=(N_steps + 1))

opti.para

opti.subject_to(battery_states[0] == battery_capacity)
for i in range(N_steps):
    # Update current simulation time
    current_time = mission_date + timedelta(seconds=i * time_delta)
    
    # Power generation using time-updated sun angles
    elevation, azimuth = calculate_sun_angles(lat, long, current_time)
    panel_amps = max(0, 6 * np.sin(np.deg2rad(elevation)))  # constant numeric value at this step
    panel_wattage = panel_amps * 5.5
    
    # Net power (solar generation minus required power)
    power_net = panel_wattage * n_solar_panels - aero["power"]
    energy_net = power_net * time_delta / 3600.0  # convert from watts to Wh per time step

    opti.subject_to(battery_states[i + 1] == battery_states[i] + energy_net)
    opti.subject_to(battery_states[i + 1] > 0)


### WEIGHT ###
#-- Power mass --#
solar_cell_mass = 0.015 * n_solar_panels
n_batt_packs = battery_capacity / (3.7 * 6 * 5)
battery_mass = 0.450 * n_batt_packs

#-- Wing mass --#
foam_volume = main_wing.volume() + hor_stabilizer.volume() + vert_stabilizer.volume()
foam_mass = foam_volume * 30.0  # foam 30kg.m^2

#-- Structure mass --#
spar_mass = (wingspan / 2 + boom_length) * 0.09  # 90g/m carbon spar 22mm
fuselages_mass = 1.0  # 1kg for all fuselage pods

#-- TOTAL --#
weight = 9.81 * (solar_cell_mass + battery_mass + foam_mass + spar_mass + fuselages_mass)


### STABILITY ###
static_margin = (cg_le_dist - aero["x_np"]) / main_wing.mean_aerodynamic_chord()


### OPTIMIZE ###
opti.minimize(aero["power"])

#-- Aerodynamic --#
opti.subject_to(aero["L"] == weight) # Lift equal to weight
opti.subject_to(cg_le_dist <= 0.25 * chordlen)
opti.subject_to(wing_airfoil.max_thickness() * chordlen > 0.030)  # Must accomodate main spar (22mm)
opti.subject_to(wingspan > 0.13 * n_solar_panels)  # Must be able to fit all of our solar panels 13cm each
opti.subject_to(static_margin > 0.20)
opti.subject_to(static_margin < 0.35)

#-- Power --#
# opti.subject_to(sim_min_battery_cap > 0)

#-- SOLVE --#
sol = opti.solve()


### Print solutions ###
print(sol(airspeed))
print(sol(wingspan))
print(sol(chordlen))
print(sol(struct_defined_aoa))
print(sol(hstab_aoa))
print(sol(hstab_span))
print(sol(hstab_chordlen))
print(sol(solar_cell_mass), sol(n_solar_panels))
print(sol(battery_mass), sol(n_batt_packs))
print(sol(foam_mass))
print(sol(spar_mass))
print(sol(fuselages_mass))
print(sol(boom_length))
print(sol(static_margin))
print(sol(cg_le_dist))

for k, v in aero.items():
    print(f"{k.rjust(4)} : {sol(aero[k])}")

vlm = sol(vlm)
vlm.draw()

airplane = sol(airplane)
airplane.draw()