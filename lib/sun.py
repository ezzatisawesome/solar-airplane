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