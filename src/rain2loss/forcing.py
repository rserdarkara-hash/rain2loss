"""Catchment forcing from gridded products and simple bias correction (linear scaling)."""
import pandas as pd

from rain2loss.gis import open_eobs, weighted_mean


def eobs_catchment_series(path, polygon, variable="rr"):
    """Daily area-weighted catchment mean of one E-OBS variable as a pandas Series."""
    grid = open_eobs(path, variable=variable)
    return weighted_mean(grid, polygon).to_series().rename(variable)


def linear_scaling_factors(series, reference, kind="ratio"):
    """Monthly correction factors so that `series` matches the mean of `reference` (overlap period).

    kind="ratio" for precipitation (multiply), kind="difference" for temperature (add).
    """
    both = pd.DataFrame({"x": series, "ref": reference}).dropna()
    monthly = both.groupby(both.index.month).mean()
    if kind == "ratio":
        return monthly["ref"] / monthly["x"]
    return monthly["ref"] - monthly["x"]


def apply_linear_scaling(series, factors, kind="ratio"):
    """Apply monthly factors from linear_scaling_factors to a daily series."""
    month_factor = factors.reindex(series.index.month).to_numpy()
    return series * month_factor if kind == "ratio" else series + month_factor