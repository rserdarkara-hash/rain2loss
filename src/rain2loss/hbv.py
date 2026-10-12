"""HBV rainfall-runoff model, daily time step (HBV-light structure, Seibert & Vis 2012).

Snow is simulated separately in equal-area elevation zones; soil and response routines are lumped.
"""
import numpy as np
import pandas as pd
from numba import njit
from scipy.stats import norm

LAPSE_RATE = -0.006   # air temperature change with elevation (degC per m)

# name: (lower bound, upper bound, unit, meaning)
PARAMS = {
    "TT":     (-2.5, 2.5,   "degC",        "threshold temperature for snowfall and melt"),
    "CFMAX":  (0.5, 10.0,   "mm/degC/day", "degree-day factor: melt per degree above TT"),
    "SFCF":   (0.4, 1.4,    "-",           "snowfall correction factor"),
    "CFR":    (0.0, 0.1,    "-",           "refreezing coefficient"),
    "CWH":    (0.0, 0.2,    "-",           "liquid water holding capacity of the snowpack"),
    "FC":     (50.0, 500.0, "mm",          "maximum soil moisture storage (field capacity)"),
    "LP":     (0.3, 1.0,    "-",           "fraction of FC above which actual ET = PET"),
    "BETA":   (1.0, 6.0,    "-",           "shape of the soil moisture - recharge relation"),
    "PERC":   (0.0, 6.0,    "mm/day",      "maximum percolation from upper to lower zone"),
    "UZL":    (0.0, 100.0,  "mm",          "upper zone storage threshold for quick flow"),
    "K0":     (0.05, 0.9,   "1/day",       "recession coefficient, quick flow"),
    "K1":     (0.01, 0.5,   "1/day",       "recession coefficient, upper zone (interflow)"),
    "K2":     (0.001, 0.15, "1/day",       "recession coefficient, lower zone (baseflow)"),
    "MAXBAS": (1.0, 7.0,    "day",         "base length of the triangular routing function"),
}
PARAM_NAMES = list(PARAMS)

# A plausible first guess for a Central European mid-mountain catchment (not calibrated)
DEFAULT_PARAMS = {"TT": 0.0, "CFMAX": 3.5, "SFCF": 1.0, "CFR": 0.05, "CWH": 0.1,
                  "FC": 250.0, "LP": 0.7, "BETA": 2.0, "PERC": 1.5, "UZL": 20.0,
                  "K0": 0.2, "K1": 0.08, "K2": 0.02, "MAXBAS": 2.5}

STATE_NAMES = ["snow", "liquid", "sm", "suz", "slz"]
OUTPUT_NAMES = ["p_in", "melt", "recharge", "et_act", "q0", "q1", "q2", "q_gen"] + STATE_NAMES


def elevation_zones(elev_mean, elev_std, n_zones=3):
    """Elevation offsets (m) of equal-area zones from the mean, assuming normally distributed elevations."""
    return elev_std * norm.ppf((np.arange(n_zones) + 0.5) / n_zones)


@njit(cache=True)
def _hbv_loop(prec, temp, pet, p, dz, s0):
    """Day-by-day HBV water balance.

    p: 13 parameters (MAXBAS excluded), dz: elevation offsets of the snow zones (m),
    s0: initial storages snow, liquid, sm, suz, slz (mm).
    """
    tt, cfmax, sfcf, cfr, cwh = p[0], p[1], p[2], p[3], p[4]
    fc, lp, beta = p[5], p[6], p[7]
    perc_max, uzl, k0, k1, k2 = p[8], p[9], p[10], p[11], p[12]

    nz = dz.shape[0]
    snow = np.full(nz, s0[0])      # frozen water in the snowpack of each zone (mm)
    liquid = np.full(nz, s0[1])    # liquid water held in the snowpack of each zone (mm)
    sm, suz, slz = s0[2], s0[3], s0[4]

    n = prec.shape[0]
    out = np.zeros((n, 13))
    for t in range(n):
        # 1. Snow routine, separately in every elevation zone
        p_in = 0.0
        melt_total = 0.0
        to_soil = 0.0
        for z in range(nz):
            tz = temp[t] + LAPSE_RATE * dz[z]          # temperature of this zone
            if tz < tt:
                rain = 0.0
                snowfall = prec[t] * sfcf
            else:
                rain = prec[t]
                snowfall = 0.0
            snow[z] += snowfall
            melt = 0.0
            if tz > tt:
                melt = min(cfmax * (tz - tt), snow[z])
                snow[z] -= melt
                liquid[z] += melt
            else:
                refreeze = min(cfr * cfmax * (tt - tz), liquid[z])
                snow[z] += refreeze
                liquid[z] -= refreeze
            liquid[z] += rain
            release = max(liquid[z] - cwh * snow[z], 0.0)   # what the snowpack cannot hold
            liquid[z] -= release
            p_in += (rain + snowfall) / nz                 # zones have equal area
            melt_total += melt / nz
            to_soil += release / nz

        # 2. Soil routine: part of the input recharges the groundwater, the rest wets the soil
        recharge = to_soil * min(sm / fc, 1.0) ** beta
        sm += to_soil - recharge
        if sm > fc:
            recharge += sm - fc
            sm = fc
        et_act = min(pet[t] * min(sm / (lp * fc), 1.0), sm)
        sm -= et_act

        # 3. Response routine: upper zone (quick flow + interflow), lower zone (baseflow)
        suz += recharge
        perc = min(perc_max, suz)
        suz -= perc
        slz += perc
        q0 = k0 * max(suz - uzl, 0.0)
        suz -= q0
        q1 = k1 * suz
        suz -= q1
        q2 = k2 * slz
        slz -= q2

        out[t, 0] = p_in
        out[t, 1] = melt_total
        out[t, 2] = recharge
        out[t, 3] = et_act
        out[t, 4] = q0
        out[t, 5] = q1
        out[t, 6] = q2
        out[t, 7] = q0 + q1 + q2
        out[t, 8] = snow.mean()
        out[t, 9] = liquid.mean()
        out[t, 10] = sm
        out[t, 11] = suz
        out[t, 12] = slz
    return out


def triangular_weights(maxbas):
    """Daily weights of the HBV triangular routing function (they sum to 1)."""
    n = int(np.ceil(maxbas))
    edges = np.minimum(np.arange(n + 1, dtype=float), maxbas)   # day boundaries 0, 1, ..., n
    cdf = np.where(edges <= maxbas / 2,
                   2 * edges ** 2 / maxbas ** 2,
                   1 - 2 * (maxbas - edges) ** 2 / maxbas ** 2)
    return np.diff(cdf)


def route(q_gen, maxbas):
    """Spread the generated runoff over the following days with the triangular weights."""
    weights = triangular_weights(maxbas)
    return np.convolve(q_gen, weights)[: len(q_gen)]


def run_hbv(forcing, params, zones=(0.0,), states0=None):
    """Run HBV on a DataFrame with columns prec (mm/day), temp (degC) and pet (mm/day).

    params: dict {name: value}. zones: elevation offsets (m) of equal-area snow zones relative to
    the elevation of the temperature series; (0.0,) means one lumped zone.
    states0: initial storages in mm (snow, liquid, sm, suz, slz).
    Returns a DataFrame with fluxes, storages and the routed discharge q_sim (mm/day).
    """
    p = np.array([params[name] for name in PARAM_NAMES[:-1]], dtype=float)
    if states0 is None:
        states0 = [0.0, 0.0, 0.5 * params["FC"], 0.0, 0.0]   # half-full soil, other stores empty
    out = _hbv_loop(forcing["prec"].to_numpy(dtype=float),
                    forcing["temp"].to_numpy(dtype=float),
                    forcing["pet"].to_numpy(dtype=float),
                    p, np.asarray(zones, dtype=float), np.asarray(states0, dtype=float))
    sim = pd.DataFrame(out, index=forcing.index, columns=OUTPUT_NAMES)
    sim["q_sim"] = route(sim["q_gen"].to_numpy(), params["MAXBAS"])
    return sim


def water_balance_error(sim, states0):
    """Inputs - outputs - storage change over the whole run (mm); should be ~0."""
    storage_start = float(np.sum(states0))
    storage_end = sim[STATE_NAMES].iloc[-1].sum()
    return (sim["p_in"].sum() - sim["q_gen"].sum() - sim["et_act"].sum()
            - (storage_end - storage_start))