import aerosandbox.numpy as np
import matplotlib.pyplot as plt

import aerosandbox.numpy as np


# Function to calculate the declination (δ) of the Sun
def calculate_declination(day_of_year):
    return -23.45 * np.cos(np.radians((360 / 365) * (day_of_year + 10)))


# Function to calculate the local hour angle (γ)
def calculate_hour_angle(local_solar_time):
    return 15 * (local_solar_time - 12)


# Function to calculate the elevation (α) and azimuth (β) angles
def calculate_sun_angles(latitude, longitude, date_time):
    # Constants
    phi = latitude  # Latitude in degrees
    lambda_ = longitude  # Longitude in degrees
    day_of_year = date_time.timetuple().tm_yday  # Day of the year
    time = date_time.time()  # Time of day
    local_solar_time = date_time.hour + (date_time.minute / 60)  # Local time in hours

    # Declination (δ)
    declination = calculate_declination(day_of_year)

    # Local Hour Angle (γ)
    hour_angle = calculate_hour_angle(local_solar_time)

    # Elevation (α)
    alpha = np.degrees(
        np.asin(
            np.sin(np.radians(declination)) * np.sin(np.radians(phi))
            + np.cos(np.radians(declination))
            * np.cos(np.radians(phi))
            * np.cos(np.radians(hour_angle))
        )
    )

    # Azimuth (β)
    if hour_angle < 0:
        azimuth = np.degrees(
            np.acos(
                (
                    np.sin(np.radians(declination)) * np.cos(np.radians(phi))
                    - np.cos(np.radians(declination))
                    * np.sin(np.radians(phi))
                    * np.cos(np.radians(hour_angle))
                )
                / np.cos(np.radians(alpha))
            )
        )
    else:
        azimuth = 360 - np.degrees(
            np.acos(
                (
                    np.sin(np.radians(declination)) * np.cos(np.radians(phi))
                    - np.cos(np.radians(declination))
                    * np.sin(np.radians(phi))
                    * np.cos(np.radians(hour_angle))
                )
                / np.cos(np.radians(alpha))
            )
        )

    return alpha, azimuth


def propagate_eci(t, sma, ecc, inc, raan, aop, mu):
    """
    Propagates orbital position in ECI coordinates.

    Parameters:
    t (float): Time (in appropriate time units)
    sma (float): Semi-major axis in km
    ecc (float): Orbital eccentricity
    inc (float): Orbital inclination in radians
    raan (float): Right ascension in radians
    aop (float): Argument of periapsis in radians

    Returns:
    np.ndarray: Position vector in the ECI frame
    """
    # Calculate mean anomaly at time t
    T = 2 * np.pi * np.sqrt(sma**3 / mu)  # Orbital period in seconds
    mean_motion = 2 * np.pi / T  # Mean motion in rad/s
    mean_anomaly = mean_motion * t

    # Convert mean anomaly to eccentric anomaly using M2E function
    eccentric_anomaly = M2E(mean_anomaly, ecc, tol=1e-8)

    # Get perifocal coordinates (x, y) in the orbital plane
    x_peri, y_peri = OE2peri(eccentric_anomaly, sma, ecc)

    # Define the rotation matrix for argument of periapsis (R_aop)
    R_aop = np.array(
        [[np.cos(-aop), np.sin(-aop), 0], [-np.sin(-aop), np.cos(-aop), 0], [0, 0, 1]]
    )

    # Define the rotation matrix for inclination (R_inc)
    R_inc = np.array(
        [[1, 0, 0], [0, np.cos(-inc), np.sin(-inc)], [0, -np.sin(-inc), np.cos(-inc)]]
    )

    # Define the rotation matrix for right ascension (R_inc)
    R_raan = np.array(
        [
            [np.cos(-raan), np.sin(-raan), 0],
            [-np.sin(-raan), np.cos(-raan), 0],
            [0, 0, 1],
        ]
    )

    # Compute the perifocal position vector in the ECI frame
    r_peri = np.array([x_peri, y_peri, 0])

    # Perform the rotations
    r_eci = R_raan @ (R_inc @ (R_aop @ r_peri))

    return r_eci


def M2E(M, e, tol=1e-6):
    """
    Solves Kepler's equation for eccentric anomaly using Newton-Raphson iteration.

    Parameters:
    M (float): Mean anomaly [rad]
    e (float): Eccentricity of orbit
    tol (float): Tolerance for Newton-Raphson iteration

    Returns:
    float: Eccentric anomaly [rad]
    """
    M = M % (2 * np.pi)  # Ensure M is within [0, 2*pi]
    E = np.pi  # Initial guess
    d = float("inf")  # Initialize delta

    while abs(d) > tol:
        d = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= d

    return E


def OE2peri(E, a, e):
    """
    Converts orbital elements to perifocal coordinates.

    Parameters:
    E (float): Eccentric anomaly [rad]
    a (float): Semi-major axis of orbit [km]
    e (float): Eccentricity of orbit

    Returns:
    tuple: (x, y) perifocal coordinates [km]
    """
    x = a * (np.cos(E) - e)
    y = a * np.sqrt(1 - e**2) * np.sin(E)
    return x, y


def UT2GMST(cal, t):
    """
    Converts from UT1 time to GMST.

    Parameters:
    cal (list or tuple): Vector containing simulation start date as [YYYY, MM, DD]
                         where DD can be a decimal value accounting for the start time.
    t (float): Current time in seconds (from Simulink clock source)

    Returns:
    float: Current Greenwich Mean Sidereal Time [rad]
    """
    mjd = UT2MJD(cal, t / 3600)  # Convert seconds to hours for UT1_to_MJD
    gmst = MJD2GMST(mjd)  # Convert MJD to GMST
    return gmst


def UT2MJD(UT1, hr):
    """
    Converts UT1 date and time to Modified Julian Date (MJD).

    Parameters:
    UT1 (list or tuple): [MM, DD, YYYY]
    hr (float): Hour of the day (fractional hours allowed)

    Returns:
    float: Modified Julian Date
    """
    M = UT1[0]
    D = UT1[1] + hr / 24  # fractional day
    Y = UT1[2]

    # Adjust year and month
    if M <= 2:
        y = Y - 1
        m = M + 12
    else:
        y = Y
        m = M

    # Calculate B (no check for pre-1582)
    B = np.floor(y / 400) - np.floor(y / 100) + np.floor(y / 4)

    # MJD calculation
    mjd = (365 * y) - 679004 + np.floor(B) + np.floor(30.6001 * (m + 1)) + D

    return mjd


def MJD2GMST(mjd):
    """
    Converts UT1 time in MJD to GMST in radians.

    Parameters:
    mjd (float): Modified Julian Date

    Returns:
    float: Greenwich Mean Sidereal Time [rad]
    """
    mjd_shifted = mjd - 51544.5  # Jan 1, 2000 at 12:00 pm
    gmst = (4.894960891 + (1.0027379093 * (2 * np.pi) * mjd_shifted)) % (2 * np.pi)
    return gmst


def ECI2ECEF(r_eci, GMST):
    """
    Creates a rotation matrix from the Celestial Reference Frame (ECI)
    to the Terrestrial Reference Frame (ECEF).

    Parameters:
    GMST (float): Greenwich Mean Sidereal Time [rad]

    Returns:
    np.ndarray: 3x3 rotation matrix from ECI to ECEF
    """
    rot_eci_ecef = np.array(
        [[np.cos(GMST), np.sin(GMST), 0], [-np.sin(GMST), np.cos(GMST), 0], [0, 0, 1]]
    )
    return rot_eci_ecef @ r_eci


def ECEF2ENU(r_ecef, lat, long):
    E = np.array([-np.sin(lat), np.cos(lat), 0])
    N = np.array(
        [-np.sin(long) * np.cos(lat), -np.sin(long) * np.sin(lat), np.cos(long)]
    )
    U = np.array([np.cos(long) * np.cos(lat), np.cos(long) * np.sin(lat), np.sin(long)])
    R = np.vstack([E, N, U])

    r_station_XYZ = 6378 * np.array(
        [np.cos(long) * np.cos(lat), np.cos(long) * np.sin(lat), np.sin(long)]
    )
    r_satellite_XYZ = r_ecef

    enu = R @ (r_satellite_XYZ - r_station_XYZ)
    return enu


def ENU2SPHERICAL(r_enu):
    rE, rN, rU = r_enu

    azimuth = np.arctan2(rE, rN)
    elevation = np.arctan2(rU, np.sqrt(rE**2 + rN**2))

    # Ensure azimuth is in [0, 2π]
    if azimuth < 0:
        azimuth += 2 * np.pi

    return azimuth, elevation


if __name__ == "__main__":
    # Sun's orbital parameters wrt Earth
    mu = 1.327124e11
    sma = 149597870  # Semi-major axis in km
    ecc = 0.0167133  # Eccentricity
    inc = np.radians(0)  # Inclination in radians
    raan = np.radians(0)  # Right ascension in radians
    aop = np.radians(282.7685)  # Argument of perihelion in radians

    # Results
    r_eci = []
    r_ecef = []
    r_enu = []
    r_spherical = []

    # Timing
    day = 4
    month = 4
    year = 2025
    t0 = day * 24 * 60 * 60 + month * 31 * 24 * 60 * 60

    # Simulate
    hours = 1
    T = range(0, hours * 3600, 10)
    for t in T:
        gmst = UT2GMST([year, day, month], t)

        pos_eci = propagate_eci(t0 + t, sma, ecc, inc, raan, aop, mu)
        pos_ecef = ECI2ECEF(pos_eci, MJD2GMST(gmst))
        pos_enu = ECEF2ENU(pos_ecef, np.radians(37.399017), np.radians(-122.151495))
        pos_spherical = ENU2SPHERICAL(pos_enu)

        r_eci.append(pos_eci)
        r_ecef.append(pos_ecef)
        r_enu.append(pos_enu)
        r_spherical.append(pos_spherical)

    sun_eci = np.array(r_eci)
    sun_ecef = np.array(r_ecef)
    sun_enu = np.array(r_enu)
    sun_spherical = np.array(r_spherical)

    # --- Plotting ---
    plt.plot(T, sun_spherical[:, 1] * 180 / np.pi)
    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection='3d')

    # # Plot the Sun's ECEF path
    # ax.plot(sun_ecef[:, 0], sun_ecef[:, 1], sun_ecef[:, 2],
    #         c='orange', linewidth=0.6, label="Sun's trajectory (ECEF)")
    # ax.set_xlabel('X (km)')
    # ax.set_ylabel('Y (km)')
    # ax.set_zlabel('Z (km)')
    # ax.set_title('Sun\'s Position in ECEF Coordinates with Earth')
    # ax.legend()
    # ax.grid(True)

    plt.show()
