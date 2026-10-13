"""Raster and vector helpers: catchment polygons, gridded data to catchment means, reprojection."""
import geopandas as gpd
import numpy as np
import rasterio
import rioxarray  # noqa: F401  (registers the .rio accessor on xarray objects)
import xarray as xr
from exactextract import exact_extract
from rasterio.warp import Resampling, calculate_default_transform, reproject

from rain2loss.paths import LAMAH

CZ_BBOX = (12.0, 48.5, 19.0, 51.1)   # lon_min, lat_min, lon_max, lat_max around Czechia


def catchment_polygon(gauge_id, crs="EPSG:4326"):
    """LamaH-CE catchment (total upstream area) of one gauge as a GeoDataFrame."""
    basins = gpd.read_file(LAMAH / "A_basins_total_upstrm/3_shapefiles/Basins_A.shp")
    return basins.loc[basins["ID"] == gauge_id, ["ID", "geometry"]].to_crs(crs)


def open_eobs(path, variable="rr", bbox=CZ_BBOX):
    """Open an E-OBS NetCDF lazily, cut it to a lon/lat box and attach the CRS."""
    ds = xr.open_dataset(path, chunks="auto")          # lazy: data are read only when needed
    lon_min, lat_min, lon_max, lat_max = bbox
    grid = ds[variable].sel(longitude=slice(lon_min, lon_max), latitude=slice(lat_min, lat_max))
    grid = grid.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
    return grid.rio.write_crs("EPSG:4326")


def simple_mean(grid, polygon):
    """Unweighted mean of all cells touching the polygon (edge cells count fully)."""
    clipped = grid.rio.clip(polygon.geometry, polygon.crs, all_touched=True)
    return clipped.mean(dim=["latitude", "longitude"])


def cell_weights(grid, polygon):
    """Cells overlapping the polygon and their area weights (sum to 1), via exactextract."""
    template = grid.isel(time=0)
    cells = exact_extract(template, polygon.to_crs(grid.rio.crs),
                          ["coverage", "center_x", "center_y"], output="pandas")
    coverage = np.asarray(cells["coverage"].iloc[0])     # covered fraction of each cell (0-1)
    lon = np.asarray(cells["center_x"].iloc[0])
    lat = np.asarray(cells["center_y"].iloc[0])
    area = coverage * np.cos(np.deg2rad(lat))             # degree cells shrink towards the pole
    return lon, lat, area / area.sum()


def weighted_mean(grid, polygon):
    """Area-weighted catchment mean that counts partially covered edge cells correctly."""
    lon, lat, weights = cell_weights(grid, polygon)
    values = grid.sel(longitude=xr.DataArray(lon, dims="cell"),
                      latitude=xr.DataArray(lat, dims="cell"), method="nearest")
    return (values * xr.DataArray(weights, dims="cell")).sum("cell")


def reproject_raster(src_path, dst_path, dst_crs, resolution, resampling=Resampling.bilinear):
    """Reproject a single-band GeoTIFF to another CRS and resolution."""
    with rasterio.open(src_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, resolution=resolution)
        profile = src.profile.copy()
        profile.update(crs=dst_crs, transform=transform, width=width, height=height)
        with rasterio.open(dst_path, "w", **profile) as dst:
            reproject(source=rasterio.band(src, 1), destination=rasterio.band(dst, 1),
                      src_transform=src.transform, src_crs=src.crs,
                      dst_transform=transform, dst_crs=dst_crs, resampling=resampling)