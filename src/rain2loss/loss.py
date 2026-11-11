"""Flood loss: hazard maps x exposure x vulnerability -> loss per return period, EP curve, AAL."""
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds

JRC_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-EFAS/flood_hazard/"
JRC_RETURN_PERIODS = [10, 20, 30, 40, 50, 75, 100, 200, 500]

# Huizinga et al. (2017), JRC global depth-damage functions: Europe, residential buildings
HUIZINGA_EUROPE_RESIDENTIAL = pd.DataFrame({
    "depth_m":      [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0],
    "damage_ratio": [0.0, 0.25, 0.40, 0.50, 0.60, 0.75, 0.85, 0.95, 1.00],
})
# Huizinga et al. (2017), sheet MaxDamage-Residential, Czech Republic, EUR per m2 (2010 prices)
MAX_DAMAGE_CZ_STRUCTURE = 391.5
MAX_DAMAGE_CZ_CONTENT = 195.7
MAX_DAMAGE_CZ_TOTAL = 587.2


def download_window(url, bounds, out_path):
    """Save the part of a remote cloud-optimised GeoTIFF that lies inside lon/lat bounds."""
    with rasterio.Env(GDAL_HTTP_MAX_RETRY=5, GDAL_HTTP_RETRY_DELAY=3):   # retry on network hiccups
        with rasterio.open("/vsicurl/" + url) as src:
            window = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
            data = src.read(1, window=window)
            profile = src.profile.copy()
            profile.update(width=data.shape[1], height=data.shape[0],
                           transform=src.window_transform(window), compress="deflate")
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data, 1)


def download_jrc_window(rp, bounds, out_path):
    """JRC v3 flood depth map of return period rp, cut to lon/lat bounds."""
    download_window(JRC_BASE + f"Europe_RP{rp}_filled_depth.tif", bounds, out_path)


def sample_depth(raster_path, points):
    """Flood depth (m) at point locations; 0 where the map is dry or has no data."""
    with rasterio.open(raster_path) as src:
        pts = points.to_crs(src.crs)
        coords = list(zip(pts.geometry.x, pts.geometry.y))
        depth = np.array([value[0] for value in src.sample(coords)], dtype=float)
        if src.nodata is not None:
            depth[depth == src.nodata] = np.nan
    depth = np.nan_to_num(depth, nan=0.0)
    return np.clip(depth, 0.0, None)


def damage_ratio(depth, curve=HUIZINGA_EUROPE_RESIDENTIAL):
    """Share of the building value lost at a given water depth (linear interpolation)."""
    return np.interp(depth, curve["depth_m"], curve["damage_ratio"])


def loss_table(buildings, depth_by_rp):
    """Total loss, number of flooded buildings and mean depth for every return period."""
    rows = []
    for rp, depth in depth_by_rp.items():
        loss = buildings["value_eur"].to_numpy() * damage_ratio(depth)
        rows.append({"rp": rp, "aep": 1 / rp, "flooded_buildings": int((depth > 0).sum()),
                     "mean_depth_m": depth[depth > 0].mean() if (depth > 0).any() else 0.0,
                     "loss_eur": loss.sum()})
    return pd.DataFrame(rows).sort_values("rp").reset_index(drop=True)


def average_annual_loss(rp, loss, protection_rp=None, tail="flat"):
    """AAL = area under the loss vs exceedance-probability curve (trapezoidal rule).

    protection_rp: floods more frequent than this cause no loss (simple flood defence).
    tail="flat" adds the loss of the rarest map for all rarer events; tail=None ignores them.
    """
    rp = np.asarray(rp, dtype=float)
    loss = np.asarray(loss, dtype=float).copy()
    order = np.argsort(rp)
    rp, loss = rp[order], loss[order]
    if protection_rp is not None:
        loss[rp <= protection_rp] = 0.0
    p = 1 / rp                                            # decreasing exceedance probability
    aal = np.sum((loss[:-1] + loss[1:]) / 2 * (p[:-1] - p[1:]))
    if tail == "flat":
        aal += loss[-1] * p[-1]                           # events rarer than the last map
    return aal