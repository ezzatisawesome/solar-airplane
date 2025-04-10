from aerosandbox.library import power_solar

mission_date = 100
lat = 37.398928

panel_wattage = 0.125**2 * power_solar.solar_flux(
    latitude=lat,
    day_of_year=mission_date,
    time=4*60*60,
    altitude=400,
    panel_azimuth_angle=0,
    panel_tilt_angle=0
) * 0.243*0.9


print(panel_wattage)