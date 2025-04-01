import aerosandbox as asb
import aerosandbox.numpy as np

wing_airfoil = asb.Airfoil("sd7037")
tail_airfoil = asb.Airfoil("naca0010")

opti=asb.Opti()

airspeed=opti.variable(init_guess=15,lower_bound=5,upper_bound=30)
wingspan=opti.variable(init_guess=3.4,lower_bound=3, upper_bound=4)
aoa=opti.variable(init_guess=3, lower_bound=-1,upper_bound=10)
chordlen = opti.variable(init_guess=0.26)
boom_length = 1.5
hstab_span = 0.5
hstab_chordlen = 0.15
hstab_aoa = -8
vstab_span = 0.3
vstab_chordlen = 0.15
weight = 4*9.8

### Define the 3D geometry you want to analyze/optimize.
# Here, all distances are in meters and all angles are in degrees.
airplane = asb.Airplane(
    name="rev 6",
    xyz_ref=[0.1*chordlen, 0, 0],  # CG location
    wings=[
        asb.Wing(
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
                    xyz_le=[0.00, wingspan/2, 0.1],
                    chord=0.125,
                    twist=0,
                    airfoil=wing_airfoil,
                ),
            ]
        ),
        asb.Wing(
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
        ).translate([boom_length, 0, 0]),
        asb.Wing(
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
    ],
    fuselages=[
        asb.Fuselage(
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[0.5*xi, 0, 0],
                    radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi) #half a meter fuselage. Starting at LE and 0.5m forward
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        )
    ]
)
vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=airspeed,  # m/s
        alpha=aoa,  # degree
    )
)


aero = vlm.run()  # Returns a dictionary
aerodyanmic_power=aero["D"]*airspeed


opti.minimize(aerodyanmic_power)

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
print(sol(aerodyanmic_power))
print(sol(wingspan))
print(sol(aoa))
print(sol(aero["D"]))
print(sol(chordlen))
#for k, v in aero.items():
#    print(f"{k.rjust(4)} : {v}")


