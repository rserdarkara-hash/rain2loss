"""DEM-based catchment analysis with pysheds: conditioning, flow routing, delineation, HAND."""
import geopandas as gpd
import numpy as np
from pysheds.grid import Grid
from shapely.geometry import shape

# pysheds 0.5 still calls np.in1d, which numpy 2.4+ removed; np.isin does the same job here
if not hasattr(np, "in1d"):
    np.in1d = np.isin


def load_conditioned_dem(path):
    """Read a projected DEM and remove pits, depressions and flats so water can always flow."""
    grid = Grid.from_raster(path)
    dem = grid.read_raster(path)
    filled = grid.fill_depressions(grid.fill_pits(dem))
    conditioned = grid.resolve_flats(filled)
    return grid, dem, conditioned


def flow_direction_and_accumulation(grid, conditioned_dem):
    """D8 flow direction and the number of upstream cells draining through every cell."""
    fdir = grid.flowdir(conditioned_dem)
    acc = grid.accumulation(fdir)
    return fdir, acc


def cell_area_km2(grid):
    """Area of one grid cell in km2 (the DEM must be in a metric CRS such as EPSG:3035)."""
    return abs(grid.affine.a * grid.affine.e) / 1e6


def delineate_catchment(grid, fdir, acc, x, y, min_area_km2):
    """Snap the gauge to the nearest river cell larger than min_area_km2, then trace upstream."""
    river = acc > min_area_km2 / cell_area_km2(grid)
    x_snap, y_snap = grid.snap_to_mask(river, (x, y))
    catch = grid.catchment(x=x_snap, y=y_snap, fdir=fdir, xytype="coordinate")
    return catch, (x_snap, y_snap)


def catchment_polygon_from_mask(grid, catch):
    """Convert the boolean catchment raster into a single polygon (GeoDataFrame)."""
    shapes = [shape(geom) for geom, value in grid.polygonize(catch.astype(np.uint8)) if value == 1]
    polygon = max(shapes, key=lambda s: s.area)          # keep the main body, drop specks
    return gpd.GeoDataFrame(geometry=[polygon], crs=grid.crs.srs)


def stream_network(grid, fdir, acc, min_area_km2):
    """River lines for all cells draining more than min_area_km2."""
    streams = acc > min_area_km2 / cell_area_km2(grid)
    network = grid.extract_river_network(fdir, streams)
    return gpd.GeoDataFrame.from_features(network["features"], crs=grid.crs.srs)


def height_above_nearest_drainage(grid, fdir, conditioned_dem, acc, min_area_km2):
    """HAND: elevation of every cell above the stream cell it drains to (m)."""
    streams = acc > min_area_km2 / cell_area_km2(grid)
    return grid.compute_hand(fdir, conditioned_dem, streams)