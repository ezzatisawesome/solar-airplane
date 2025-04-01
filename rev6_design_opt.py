import aerosandbox as asb
import aerosandbox.numpy as np



def calculate_skin_friction(length):
    reynolds=1.29*airspeed*length/1.78e-5
    cf=0.455/(np.log10(reynolds)**2.58) #Prandtl-Schlichting’s approximation
    return cf

opti=asb.Opti()

wing_airfoil = asb.Airfoil("sd7037")
tail_airfoil = asb.Airfoil("naca0010")
airspeed=opti.variable(init_guess=15,lower_bound=5,upper_bound=30)
wingspan=opti.variable(init_guess=3.4,lower_bound=3, upper_bound=4)
aoa=opti.variable(init_guess=3, lower_bound=-1,upper_bound=3)
chordlen = opti.variable(init_guess=0.26)
boom_length = 1.5
hstab_span = opti.variable(init_guess=0.5, lower_bound=0.4, upper_bound=0.8)
hstab_chordlen = 0.15
hstab_aoa = opti.variable(init_guess=-5, lower_bound=-5, upper_bound=-3)
vstab_span = 0.3
vstab_chordlen = 0.15
weight = 4*9.8
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
                    xyz_c=[0.2*xi, 0.75, 0],
                    radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi) #half a meter fuselage. Starting at LE and 0.5m forward
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        )

right_pod = asb.Fuselage( #right pod fuselage
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[0.2*xi, -0.75, 0],
                    radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi) #half a meter fuselage. Starting at LE and 0.5m forward
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        )
### Define the 3D geometry you want to analyze/optimize.
# Here, all distances are in meters and all angles are in degrees.
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
        alpha=aoa,  # degree
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
aero['power']=aero['D_tot']*airspeed + 8  #8w to run avionics


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
opti.subject_to(wing_airfoil.max_thickness()*chordlen>0.025) #must accomodate main spar (22nn)

sol=opti.solve()
print(sol(airspeed))
print(sol(wingspan))
print(sol(chordlen))
print(sol(hstab_aoa))
print(sol(hstab_span))

for k, v in aero.items():
    print(f"{k.rjust(4)} : {sol(aero[k])}")
