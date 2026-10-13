"""Goodness-of-fit metrics for simulated vs observed discharge."""
import numpy as np


def _paired(sim, obs):
    """Return sim and obs as float arrays, keeping only days where both are present."""
    sim = np.asarray(sim, dtype=float)
    obs = np.asarray(obs, dtype=float)
    ok = ~np.isnan(sim) & ~np.isnan(obs)
    return sim[ok], obs[ok]


def nse(sim, obs):
    """Nash-Sutcliffe efficiency: 1 is perfect, 0 is as good as the observed mean."""
    sim, obs = _paired(sim, obs)
    return 1 - np.sum((sim - obs) ** 2) / np.sum((obs - obs.mean()) ** 2)


def kge(sim, obs, components=False):
    """Kling-Gupta efficiency (Gupta et al. 2009); optionally also r, alpha, beta."""
    sim, obs = _paired(sim, obs)
    r = np.corrcoef(sim, obs)[0, 1]          # timing and shape
    alpha = sim.std() / obs.std()            # variability ratio
    beta = sim.mean() / obs.mean()           # volume (bias) ratio
    value = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
    if components:
        return value, r, alpha, beta
    return value


def pbias(sim, obs):
    """Percent bias of the total volume; positive means the model overestimates."""
    sim, obs = _paired(sim, obs)
    return 100 * (sim.sum() - obs.sum()) / obs.sum()


def log_nse(sim, obs):
    """NSE of log-transformed flows: emphasises low flows."""
    sim, obs = _paired(sim, obs)
    eps = obs.mean() / 100                   # avoids log(0)
    return nse(np.log(sim + eps), np.log(obs + eps))


def evaluate(sim, obs):
    """All metrics in one dict, convenient for building tables."""
    k, r, alpha, beta = kge(sim, obs, components=True)
    return {"NSE": nse(sim, obs), "KGE": k, "r": r, "alpha": alpha, "beta": beta,
            "PBIAS": pbias(sim, obs), "logNSE": log_nse(sim, obs)}