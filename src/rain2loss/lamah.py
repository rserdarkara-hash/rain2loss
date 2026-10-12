"""Readers for the LamaH-CE dataset (daily resolution)."""
import pandas as pd
from pyproj import Transformer

from rain2loss.paths import LAMAH
from rain2loss.pet import oudin_pet

CZECH_REGIONS = ["OLM", "SBO", "SMO", "VYS", "ZLN"]


def read_lamah_ts(path):
    """Read a LamaH daily time series CSV into a date-indexed DataFrame."""
    df = pd.read_csv(path, sep=";")
    df.index = pd.to_datetime(dict(year=df.YYYY, month=df.MM, day=df.DD))
    return df.drop(columns=["YYYY", "MM", "DD"])


def load_attributes():
    """Catchment and gauge attributes, one row per gauge ID."""
    catchments = pd.read_csv(
        LAMAH / "A_basins_total_upstrm/1_attributes/Catchment_attributes.csv",
        sep=";", index_col="ID")
    gauges = pd.read_csv(
        LAMAH / "D_gauges/1_attributes/Gauge_attributes.csv",
        sep=";", index_col="ID")
    return pd.concat([catchments, gauges], axis=1)


def load_lamah(gauge_id, attrs=None):
    """Daily forcing plus observed discharge (m3/s and mm/day) for one gauge."""
    if attrs is None:
        attrs = load_attributes()
    forcing = read_lamah_ts(
        LAMAH / f"A_basins_total_upstrm/2_timeseries/daily/ID_{gauge_id}.csv")
    flow = read_lamah_ts(LAMAH / f"D_gauges/2_timeseries/daily/ID_{gauge_id}.csv")
    q = flow["qobs"].where(flow["qobs"] >= 0)    # -999 means missing
    area_km2 = attrs.loc[gauge_id, "area_gov"]
    out = forcing.copy()
    out["qobs_m3s"] = q
    out["qobs_mm"] = q / (area_km2 * 1e6) * 1000 * 86400
    return out


def gauge_latlon(gauge_id, attrs):
    """Gauge position in degrees. LamaH 'lon'/'lat' columns are EPSG:3035 metres, not degrees."""
    to_wgs84 = Transformer.from_crs("EPSG:3035", "EPSG:4326", always_xy=True)
    lon, lat = to_wgs84.transform(attrs.loc[gauge_id, "lon"], attrs.loc[gauge_id, "lat"])
    return lat, lon


def load_hbv_inputs(gauge_id, attrs=None):
    """Daily HBV forcing (prec, temp, pet) plus observed discharge qobs (mm/day) for one gauge."""
    if attrs is None:
        attrs = load_attributes()
    df = load_lamah(gauge_id, attrs)
    lat, _ = gauge_latlon(gauge_id, attrs)
    out = pd.DataFrame(index=df.index)
    out["prec"] = df["prec"]                      # mm/day, ERA5-Land catchment mean
    out["temp"] = df["2m_temp_mean"]              # degC
    out["pet"] = oudin_pet(out["temp"], out.index.dayofyear, lat)   # mm/day
    out["qobs"] = df["qobs_mm"]                   # mm/day
    return out


def load_cosero(gauge_id, attrs=None):
    """Discharge simulated by the COSERO benchmark model in LamaH-CE (m3/s -> mm/day)."""
    if attrs is None:
        attrs = load_attributes()
    df = read_lamah_ts(LAMAH / f"F_hydrol_model/2_timeseries/ID_{gauge_id}.csv")
    area_km2 = attrs.loc[gauge_id, "area_gov"]
    return df["Qsim_A"] / (area_km2 * 1e6) * 1000 * 86400