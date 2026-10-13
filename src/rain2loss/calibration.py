"""Calibration of HBV with spotpy (SCE-UA, DDS) and a simple GLUE uncertainty analysis."""
import contextlib
import io

import numpy as np
import pandas as pd
import spotpy
from scipy.stats import qmc

from rain2loss.hbv import PARAMS, PARAM_NAMES, run_hbv
from rain2loss.metrics import kge


class HBVSetup:
    """spotpy setup: which parameters, how to simulate, what to compare, how to score."""

    def __init__(self, forcing, q_obs, eval_mask, zones=(0.0,), minimize=True):
        self.forcing = forcing                          # prec, temp, pet for the whole run
        self.q_obs = np.asarray(q_obs, dtype=float)
        self.mask = np.asarray(eval_mask, dtype=bool)   # True on the days that are scored
        self.zones = zones
        self.minimize = minimize                        # SCE-UA minimises, DDS maximises
        self.params = [spotpy.parameter.Uniform(name, low=lo, high=hi)
                       for name, (lo, hi, unit, meaning) in PARAMS.items()]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        params = dict(zip(PARAM_NAMES, np.asarray(vector, dtype=float)))
        q_sim = run_hbv(self.forcing, params, self.zones)["q_sim"].to_numpy()
        return q_sim[self.mask]

    def evaluation(self):
        return self.q_obs[self.mask]

    def objectivefunction(self, simulation, evaluation):
        score = kge(simulation, evaluation)
        return 1 - score if self.minimize else score


def best_params(results, minimize=True):
    """Best parameter set of a spotpy result array, as a dict."""
    i = np.nanargmin(results["like1"]) if minimize else np.nanargmax(results["like1"])
    return {name: float(results["par" + name][i]) for name in PARAM_NAMES}


def calibrate(algorithm, forcing, q_obs, eval_mask, zones=(0.0,), repetitions=5000, seed=42,
              quiet=True, **sample_kwargs):
    """Calibrate HBV with spotpy 'sceua' or 'dds'; returns (best parameter dict, all results)."""
    minimize = algorithm == "sceua"            # SCE-UA minimises 1 - KGE, DDS maximises KGE
    setup = HBVSetup(forcing, q_obs, eval_mask, zones, minimize=minimize)
    sampler_class = getattr(spotpy.algorithms, algorithm)
    sampler = sampler_class(setup, dbname="hbv_" + algorithm, dbformat="ram", save_sim=False)
    np.random.seed(seed)                       # reproducible random search
    output = contextlib.redirect_stdout(io.StringIO()) if quiet else contextlib.nullcontext()
    with output:                               # spotpy prints a lot; hide it when quiet=True
        sampler.sample(repetitions, **sample_kwargs)
    results = sampler.getdata()
    return best_params(results, minimize), results


def glue_sample(forcing, q_obs, eval_mask, zones=(0.0,), n=10000, seed=42):
    """Latin hypercube sample of parameter sets, each scored with KGE on the masked days."""
    lower = [PARAMS[name][0] for name in PARAM_NAMES]
    upper = [PARAMS[name][1] for name in PARAM_NAMES]
    unit_cube = qmc.LatinHypercube(d=len(PARAM_NAMES), seed=seed).random(n)
    samples = pd.DataFrame(qmc.scale(unit_cube, lower, upper), columns=PARAM_NAMES)
    obs = np.asarray(q_obs, dtype=float)[eval_mask]
    scores = []
    for params in samples.to_dict("records"):
        q_sim = run_hbv(forcing, params, zones)["q_sim"].to_numpy()[eval_mask]
        scores.append(kge(q_sim, obs))
    samples["KGE"] = scores
    return samples


def glue_bounds(forcing, behavioural, zones=(0.0,), lower_q=0.05, upper_q=0.95):
    """Daily prediction bounds (quantiles) from the behavioural parameter sets."""
    sims = np.array([run_hbv(forcing, params, zones)["q_sim"].to_numpy()
                     for params in behavioural[PARAM_NAMES].to_dict("records")])
    return pd.DataFrame({"lower": np.quantile(sims, lower_q, axis=0),
                         "median": np.quantile(sims, 0.5, axis=0),
                         "upper": np.quantile(sims, upper_q, axis=0)},
                        index=forcing.index)