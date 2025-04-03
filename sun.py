import aerosandbox.numpy as np


# Sun's orbital parameters wrt Earth
sma = 149597870              # Semi-major axis in km
i = np.radians(23.4406)      # Inclination in radians
e = 0.0167133                # Eccentricity
aop = np.radians(282.7685)   # Argument of perihelion in radians
T = 365 * 24 * 60 * 60       # Orbital period in seconds
mean_motion = 2 * np.pi / T  # Mean motion in rad/s


def solve_kepler(M, e, tol=1e-6):
    E = M  # Initial guess
    for _ in range(10):  # Iterate using Newton's method
        E_new = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        if np.abs(E_new - E) < tol:
            break
        E = E_new
    return E

def propagate_eci(t):
    # Step 1: Compute mean anomaly
    mean_anomaly = mean_motion * t

    # Step 2: Get ecc anomaly
    eccentric_anomaly = solve_kepler(mean_anomaly, e)

    # Step 3: Get true anomaly
    true_anomaly = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(eccentric_anomaly / 2), np.sqrt(1 - e) * np.cos(eccentric_anomaly / 2))

    # Step 4: Get perifocal coordinates
    r = sma * (1 - e**2) / (1 + e * np.cos(true_anomaly))  # Orbital radius
    x_peri = r * np.cos(true_anomaly)
    y_peri = r * np.sin(true_anomaly)
    z_peri = 0

    # Step 5: Apply incliation and AoP rotations
    R_aop = np.array([[np.cos(aop), -np.sin(aop), 0],
                    [np.sin(aop), np.cos(aop), 0],
                    [0, 0, 1]])

    R_inc = np.array([[1, 0, 0],
                    [0, np.cos(i), -np.sin(i)],
                    [0, np.sin(i), np.cos(i)]])
    r_ecliptic = R_inc @ (R_aop @ np.array([x_peri, y_peri, z_peri]))

def eci2ecef(r):
    pass

def ecef2enu(r):
    pass

def ut12mjd(m, d, y, hr):
    d = d + hr/24

    if (m <= 2):
        y = y - 1
        m = m + 12
    
    b = np.floor(y/400) - np.floor(y/100) + np.floor(y/4)
    
    return (365*y) - 679004 + np.floor(b) + np.floor(30.6001 * (m+1)) + d

def mjd2gmst(t):
    return 280.4606 + 360.9856473 * t