def apc_prop(diameter):
    """
    Estimate the mass of an APC Thin Electric propeller.

    Parameters
    ----------
    diameter_in : float
        Propeller diameter in meters.

    Returns
    -------
    mass_oz : float
        Estimated mass in ounces.
    mass_g : float
        Estimated mass in kg.
    """
    # Convert diameter form m to in
    diameter_in = diameter * 39.3701
    k = 0.0025340420925460466
    n = 2.4243828292426777
    m_oz = k * (diameter_in ** n)
    m_kg = m_oz * 0.0283495
    return m_kg


def apc_prop_carbon(diameter):
    """
    Estimate the mass of an APC Carbon Fiber (Electric Only) prop.

    Parameters
    ----------
    diameter_in : float
        Propeller diameter in meters.

    Returns
    -------
    mass_oz : float
        Estimated mass in ounces.
    mass_g : float
        Estimated mass in kg.
    """
    # Fitted coefficients from regression on catalog data
    k = 0.0028691378265985   # oz / in^n
    n = 2.255324275728351
    diameter_in = diameter * 39.3701
    mass_oz = k * (diameter_in ** n)
    mass_kg = mass_oz * 0.0283495
    return mass_kg


def hacker_motor_mass_kv(kv):
    """
    KV-only motor mass model fitted to your table.

    Parameters
    ----------
    kv_rpm_per_v : float or array-like
        Motor KV in rpm/V.

    Returns
    -------
    mass_kg : float or ndarray
        Predicted motor mass in kg.
    """
    mass_g = 1.146347e6 * (kv ** -1.349136)
    return mass_g / 1000.0  # g -> kg


def hacker_motor_mass(max_power, kv):
    """
    Predict BLDC motor mass (grams) from max power and KV.

    Parameters
    ----------
    max_power : float or array-like
        Maximum power in watts.
    kv_rpm : float or array-like
        Motor KV in rpm/V.

    Returns
    -------
    mass_g : float or ndarray
        Predicted motor mass in kg.
    """
    return (0.92 * (max_power ** 0.83) * (kv ** -0.47)) / 100


def hacker_motor_resistance(max_power, kv):
    """
    Predict BLDC motor internal resistance (ohms)
    using power- and KV-based scaling only.
    """
    return 0.060 * (max_power ** -0.68) * (kv ** 0.62)


def hacker_motor_i0(kv):
    """
    No-load (idle) current I0 in amps for a Hacker-class BLDC motor.

    FLAGGED: anchored to the A40-12L V4 kv410 datasheet point (I0 = 1.04 A @
    8.4 V) and scaled linearly with Kv (iron/windage losses rise with speed
    for a given stator). Replace with a per-motor measured value when known.

    Why it matters: at night cruise the propeller torque is tiny, so the
    torque-producing current Q/Kt is small and I0 becomes a *large* fraction
    of the motor current. Omitting it understates night bus watts — the one
    quantity that actually binds this design.
    """
    return 1.04 * (kv / 410.0)




### CONSTANTS
liion_21700_cell_mass = 70 # grams
solar_panel_with_busbar_solder_mass = 15 # grams
mppt_mass = 200 # grams

# Avionics
speedybee_f405_mass = 55 # grams
gps_bn880_mass = 12 # grams
telemtry_with_wire_mass = 26 # grams
radio_receiver_with_wire_mass = 18 # grams
wiring_mass = 50 # grams
servos_mass = 50 # grams

# Propulsion
esc_mass = 25 # ! grams, may change ESC
motor_410kv_mass = 277 # grams
propeller_mass = 52 # grams

# Structures
fuselage_mass = 100 # grams
motor_pods_mass = 60 # grams
carbon_22_per_meter_mass = 90 # grams / m
carbon_20_per_meter_mass = 80 # grams / m

# Wing
ngx_foam_density = 30.02 # kg / m^3