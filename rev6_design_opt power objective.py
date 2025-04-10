import aerosandbox as asb
import aerosandbox.numpy as np



def calculate_skin_friction(length):
    reynolds=1.29*airspeed*length/1.78e-5
    cf=0.455/(np.log10(reynolds)**2.58) #Prandtl-Schlichting’s approximation
    return cf

opti=asb.Opti(cache_filename="power_objective.json")

wing_airfoil = asb.Airfoil("sd7037")
tail_airfoil = asb.Airfoil("naca0010")

airspeed=opti.variable(init_guess=15,lower_bound=5,upper_bound=40)
wingspan=opti.variable(init_guess=3.5,lower_bound=2.5,upper_bound=10)
chordlen = opti.variable(init_guess=0.26)
boom_length = 1.5
hstab_span = 0.5
hstab_chordlen = 0.15
hstab_aoa = -2
vstab_span = 0.3
vstab_chordlen = 0.15
polyhedral_angle =10

main_wing = asb.Wing(
            name="Main Wing",
            symmetric=True,  # Should this wing be mirrored across the XZ plane?
            xsecs=[  # The wing's cross ("X") sections
                asb.WingXSec(  # Root
                    xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                    chord=chordlen,
                    twist=0,  # degrees
                    airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                ),
                asb.WingXSec(  # Mid
                    xyz_le=[0.00, 0.5*wingspan/2, 0],
                    chord=chordlen,
                    twist=0,
                    airfoil=wing_airfoil,
                ),
                asb.WingXSec(  # Tip
                    xyz_le=[0.00, wingspan/2, np.sin(10*np.pi/180)*0.5 * wingspan/2],
                    chord=0.125,
                    twist=0,
                    airfoil=wing_airfoil,
                ),
            ]
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
                    xyz_le=[0.0, hstab_span/2, 0],
                    chord=hstab_chordlen,
                    twist=hstab_aoa,
                    airfoil=tail_airfoil
                )
            ]
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
                    airfoil=tail_airfoil
                )
            ]
        ).translate([boom_length+hstab_chordlen, 0,0])

main_fuselage=asb.Fuselage( #main fuselage
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[0.5*xi, 0, 0],
                    radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi) #half a meter fuselage. Starting at LE and 0.5m forward
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        ).translate([-0.5,0,0])

left_pod = asb.Fuselage( #left pod fuselage
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[0.2*xi, 0.75, -0.02],
                    radius=0.4 * asb.Airfoil("dae51").local_thickness(x_over_c=xi) #half a meter fuselage. Starting at LE and 0.5m forward
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        )

right_pod = asb.Fuselage( #right pod fuselage
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[0.2*xi, -0.75, -0.02],
                    radius=0.4 * asb.Airfoil("dae51").local_thickness(x_over_c=xi) #half a meter fuselage. Starting at LE and 0.5m forward
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        )

airplane = asb.Airplane(
    name="rev 6",
    xyz_ref=[0.1*chordlen, 0, 0],  # CG location
    wings=[ main_wing, hor_stabilizer, vert_stabilizer],
    fuselages=[main_fuselage, left_pod, right_pod]
)


vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=airspeed,  # m/s 
    )
)


aero = vlm.run_with_stability_derivatives()  # Returns a dictionary
aerodyanmic_power=aero["D"]*airspeed

#VLM does not calcualte parasitic drag, we must add this manually
CD0 = (calculate_skin_friction(chordlen)*main_wing.area(type='wetted')/main_wing.area() +
    calculate_skin_friction(hstab_chordlen)*hor_stabilizer.area(type='wetted')/main_wing.area() +
    calculate_skin_friction(vstab_chordlen)*vert_stabilizer.area(type='wetted')/main_wing.area() +
    calculate_skin_friction(0.5)*main_fuselage.area_wetted()/main_wing.area()+
    2*calculate_skin_friction(0.2)*left_pod.area_wetted()/main_wing.area())

drag_parasite=0.5*1.29* airspeed**2 * main_wing.area() * CD0

aero["CD_tot"] = aero["CD"]+CD0
aero["D_tot"] = aero["D"] + drag_parasite
aero['power']=(aero['D_tot']*airspeed + 8 ) * 1.4 #8w to run avionics 40% surplus for non ideal conditions and charging



## POWER
from aerosandbox.library import power_solar

N = 100  # Number of discretization points
time = np.linspace(0, 24*60*60, N) # s
dt = np.diff(time)[0] # s
mission_date = 100
lat = 37.398928

battery_cap = opti.variable(init_guess=1500, lower_bound=100, scale=1000, category="power")  # initial battery energy in Wh
battery_states = opti.variable(n_vars = N, init_guess=500, scale=500, category="power")
opti.subject_to(battery_states[0] == battery_cap)

n_solar_panels = opti.variable(init_guess=40, lower_bound=10, scale=30, category="power")

for i in range(N-1):
    elevation = power_solar.solar_elevation_angle(lat, mission_date, time[i])
    panel_amps = np.softmax(0, 6 * np.sin(np.deg2rad(elevation)), hardness=10)
    panel_wattage = panel_amps * 0.55  # W per panel

    power_generated = panel_wattage * n_solar_panels
    power_used = aero["power"]
    net_energy = (power_generated - power_used) * (dt / 3600)  # Wh

    battery_update = np.softmin(battery_states[i] + net_energy, battery_cap, hardness=10)
    # The battery state integration constraint:
    opti.subject_to(battery_states[i+1] == battery_update)


opti.subject_to(battery_states > 0)


##Weight
solar_cell_weight = 0.015*n_solar_panels*9.81

current_draw = aero['power']/15 #4cell configuration is 15V
num_packs = battery_cap / (5 * 6 * 3.7)
battery_weight = num_packs * 0.450 * 9.81

foam_volume = main_wing.volume() + hor_stabilizer.volume() + vert_stabilizer.volume()
foam_weight = foam_volume * 30.0 * 9.81 #foam 30kg.m^2

spar_weights = (wingspan/2 + boom_length) * 0.09 * 9.81 #90g/m carbon spar 22mm

fuselages_weight = 1.0*9.81 #1kg for all fuselage pods

weight = solar_cell_weight + battery_weight + foam_weight+spar_weights+fuselages_weight

opti.minimize(aero['power'])

#Constraints
""""
Lift equal to weight
thickness of wing greater than spar
charge rate
battery capacity (duration)
aspect ratio
"""
opti.subject_to(aero["L"]==weight)
opti.subject_to(wing_airfoil.max_thickness()*chordlen>0.025) #must accomodate main spar (22mm)
opti.subject_to(wingspan>0.13*n_solar_panels)#Must be able to fit all of our solar panels 13cm each

sol=opti.solve()
opti.save_solution()
print(sol(airspeed))
print(sol(wingspan))
print(sol(chordlen))
print(sol(hstab_aoa))
print(sol(hstab_span))
print(sol(solar_cell_weight),sol(n_solar_panels))
print(sol(battery_weight),sol(num_packs))
print(sol(foam_weight))
print(sol(spar_weights))
print(sol(fuselages_weight))

for k, v in aero.items():
    print(f"{k.rjust(4)} : {sol(aero[k])}")

vlm=sol(vlm)
vlm.draw()

airplane=sol(airplane)
airplane.draw()