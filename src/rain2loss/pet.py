"""Potential evapotranspiration (PET) from air temperature (Oudin et al., 2005)."""
import numpy as np

SOLAR_CONSTANT = 0.0820   # MJ m-2 min-1 (FAO-56)
LATENT_HEAT = 2.45        # MJ kg-1, latent heat of vaporisation of water
WATER_DENSITY = 1000.0    # kg m-3


def extraterrestrial_radiation(doy, lat_deg):
    """Daily extraterrestrial radiation Re (MJ m-2 day-1), FAO-56 equations 21-25."""
    doy = np.asarray(doy, dtype=float)
    phi = np.deg2rad(lat_deg)                                       # latitude (rad)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)                  # inverse relative Earth-Sun distance
    delta = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)            # solar declination (rad)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1, 1))    # sunset hour angle (rad)
    return (24 * 60 / np.pi) * SOLAR_CONSTANT * dr * (
        ws * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(ws))


def oudin_pet(temp_c, doy, lat_deg):
    """PET (mm/day) from daily mean air temperature; zero when temp_c <= -5 degC."""
    temp_c = np.asarray(temp_c, dtype=float)
    re = extraterrestrial_radiation(doy, lat_deg)
    pet = re / (LATENT_HEAT * WATER_DENSITY) * (temp_c + 5) / 100 * 1000   # m/day -> mm/day
    return np.where(temp_c + 5 > 0, pet, 0.0)