"""Extreme value statistics for floods and a simple stochastic event set."""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gamma


def water_year(index, start_month=11):
    """Hydrological year of each date, labelled by the year in which it ends (CHMI: Nov-Oct)."""
    return index.year + (index.month >= start_month)


def annual_maxima(q, start_month=11, min_days=330):
    """Annual maximum series (AMS) with the date of each maximum; incomplete years are dropped."""
    q = q.dropna()
    wy = water_year(q.index, start_month)
    groups = q.groupby(wy)
    ams = pd.DataFrame({"peak": groups.max(), "date": groups.idxmax(), "n_days": groups.size()})
    return ams[ams["n_days"] >= min_days].drop(columns="n_days")


def lmoments(x):
    """First two sample L-moments and the L-skewness (l1, l2, t3), via probability-weighted moments."""
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    j = np.arange(1, n + 1)
    b0 = x.mean()
    b1 = np.sum((j - 1) / (n - 1) * x) / n
    b2 = np.sum((j - 1) * (j - 2) / ((n - 1) * (n - 2)) * x) / n
    l1, l2, l3 = b0, 2 * b1 - b0, 6 * b2 - 6 * b1 + b0
    return l1, l2, l3 / l2


def fit_gev_lmom(x):
    """GEV parameters from L-moments (Hosking 1990), in scipy's convention (c, loc, scale)."""
    l1, l2, t3 = lmoments(x)
    z = 2 / (3 + t3) - np.log(2) / np.log(3)
    c = 7.8590 * z + 2.9554 * z ** 2                 # Hosking's k, equal to scipy's c (= -xi)
    scale = l2 * c / ((1 - 2 ** (-c)) * gamma(1 + c))
    loc = l1 - scale * (1 - gamma(1 + c)) / c
    return c, loc, scale


def fit_gev_mle(x):
    """GEV parameters by maximum likelihood (scipy), as (c, loc, scale)."""
    return stats.genextreme.fit(np.asarray(x, dtype=float))


def return_level(params, T):
    """Discharge exceeded on average once every T years."""
    c, loc, scale = params
    return stats.genextreme.ppf(1 - 1 / np.asarray(T, dtype=float), c, loc, scale)


def return_period(params, value):
    """Return period (years) of a given discharge under a fitted GEV."""
    c, loc, scale = params
    return 1 / stats.genextreme.sf(value, c, loc, scale)


def empirical_return_periods(values):
    """Weibull plotting positions: T = (n + 1) / rank, rank 1 = largest value."""
    values = np.asarray(values, dtype=float)
    rank = stats.rankdata(-values, method="ordinal")
    return (len(values) + 1) / rank


def bootstrap_return_levels(x, T, fit=fit_gev_lmom, n_boot=1000, seed=42):
    """Return levels for many resampled series; percentiles of the result give confidence intervals."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    levels = np.empty((n_boot, len(T)))
    for i in range(n_boot):
        sample = rng.choice(x, size=len(x), replace=True)
        levels[i] = return_level(fit(sample), T)
    return pd.DataFrame(levels, columns=T)


def block_bootstrap_years(forcing, n_years=1000, start_month=11, seed=1):
    """Synthetic forcing made of whole hydrological years drawn at random from the record."""
    wy = water_year(forcing.index, start_month)
    days_per_year = pd.Series(wy).value_counts()
    complete = days_per_year[days_per_year >= 365].index.sort_values()
    rng = np.random.default_rng(seed)
    chosen = rng.choice(complete, size=n_years, replace=True)
    blocks = []
    for synthetic_year, source_year in enumerate(chosen, start=1):
        block = forcing[wy == source_year].copy()
        block["syn_year"] = synthetic_year
        block["source_year"] = source_year
        blocks.append(block)
    return pd.concat(blocks, ignore_index=True)


def scale_largest_storm(synthetic, window=5, sigma=0.2, seed=2):
    """Multiply the wettest `window`-day spell of every synthetic year by a random lognormal factor."""
    rng = np.random.default_rng(seed)
    out = synthetic.copy()
    factors = {}
    for year, block in out.groupby("syn_year"):
        rolling = block["prec"].rolling(window).sum()
        end = rolling.idxmax()                           # last day of the wettest spell
        factor = rng.lognormal(mean=0.0, sigma=sigma)
        out.loc[end - window + 1: end, "prec"] *= factor
        factors[year] = factor
    return out, pd.Series(factors, name="storm_factor")