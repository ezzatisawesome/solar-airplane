from lib.sun import calculate_sun_angles
from datetime import datetime, timedelta
from aerosandbox import numpy as np

def estimate_solar_power(elevation_deg, max_current=6.0, max_voltage=0.55):
    if elevation_deg <= 0:
        return 0.0
    max_power = max_current * max_voltage
    angle_factor = np.sin(np.radians(elevation_deg))
    return max_power * angle_factor

# --- New: Total daily energy ---
def calculate_daily_energy(latitude, longitude, date, step_minutes=10):
    total_energy = 0.0
    time = datetime(date.year, date.month, date.day, 0, 0)
    end_time = time + timedelta(days=1)
    step = timedelta(minutes=step_minutes)

    previous_power = 0.0
    while time < end_time:
        elevation, _ = calculate_sun_angles(latitude, longitude, time)
        power = estimate_solar_power(elevation)
        # Trapezoidal integration
        energy = (power + previous_power) / 2 * (step_minutes / 60)
        total_energy += energy
        previous_power = power
        time += step

    return total_energy  # in watt-hours (Wh)

# --- Example usage ---
latitude = 38.6   # Saint Louis
longitude = -90.2
date = datetime(2023, 2, 20)

daily_energy = calculate_daily_energy(latitude, longitude, date)
print(f"Estimated Total Energy on {date.date()}: {daily_energy:.2f} Wh")