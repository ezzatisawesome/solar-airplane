import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.library.aerodynamics as lib_aero
from aerosandbox.library import mass_structural as lib_mass_struct
from aerosandbox.library import power_solar as lib_solar
from aerosandbox.library import propulsion_electric as lib_prop_elec
from aerosandbox.library import propulsion_propeller as lib_prop_prop
import aerosandbox.tools.units as u
import copy
from pathlib import Path


opti = asb.Opti(
    freeze_style='float',
    # variable_categories_to_freeze='all'
)

### CONSTANTS
# Power 
liion_21700_cell_m = 70 # grams
solar_panel_with_busbar_solder = 15 # grams

# Avionics
speedybee_f405_m = 55 # grams
gps_bn880_m = 12 # grams
telemtry_with_wire = 26 # grams
radio_receiver_with_wire = 18 # grams

# Propulsion
esc = 25 # ! grams, may change ESC
motor_410kv = 277 # grams

# Structures
carbon_22_per_meter = 90 # grams / m
carbon_20_per_meter = 80 # grams / m


# PARAMETERS
n_panels_spanwise_inboard = opti.parameter(12)  # number of panels per side
n_panels_chordwise_inboard = opti.parameter(1)
n_panels_spanwise_outboard = opti.parameter(12)  # number of panels per side
n_panels_chordwise_outboard = opti.parameter(1)



# MASS
design_mass_TOGW = opti.variable(
    init_guess=5,
    lower_bound=1e-3,
)





## Mass
# 