import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from shapely.geometry import box

from rain2loss.gis import weighted_mean


@pytest.fixture
def grid():
    """2 x 2 cells of 0.1 degree, two days of rain."""
    lon = np.array([10.05, 10.15])
    lat = np.array([50.05, 50.15])
    time = pd.date_range("2000-01-01", periods=2, freq="D")
    values = np.array([[[1.0, 4.0], [9.0, 9.0]],
                       [[2.0, 8.0], [9.0, 9.0]]])          # dims: time, latitude, longitude
    da = xr.DataArray(values, coords={"time": time, "latitude": lat, "longitude": lon},
                      dims=("time", "latitude", "longitude"))
    return da.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude").rio.write_crs("EPSG:4326")


def test_partially_covered_cell_gets_partial_weight(grid):
    # polygon: the whole south-west cell plus half of the south-east cell
    polygon = gpd.GeoDataFrame(geometry=[box(10.0, 50.0, 10.15, 50.1)], crs="EPSG:4326")
    result = weighted_mean(grid, polygon).to_numpy()
    expected = [(2 * 1.0 + 1 * 4.0) / 3, (2 * 2.0 + 1 * 8.0) / 3]
    np.testing.assert_allclose(result, expected, rtol=1e-4)