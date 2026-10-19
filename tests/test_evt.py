import numpy as np
import pandas as pd
import pytest
from scipy import stats

from rain2loss.evt import annual_maxima, fit_gev_lmom, return_level


def test_lmoment_fit_recovers_known_gev():
    true = (-0.1, 100.0, 30.0)                              # scipy convention (c, loc, scale)
    sample = stats.genextreme.rvs(*true, size=20000, random_state=1)
    c, loc, scale = fit_gev_lmom(sample)
    assert c == pytest.approx(true[0], abs=0.02)
    assert loc == pytest.approx(true[1], rel=0.02)
    assert scale == pytest.approx(true[2], rel=0.03)


def test_return_level_increases_with_T():
    levels = return_level((-0.1, 100.0, 30.0), [2, 10, 100])
    assert np.all(np.diff(levels) > 0)


def test_annual_maxima_water_year():
    days = pd.date_range("2000-11-01", "2002-10-31", freq="D")
    q = pd.Series(1.0, index=days)
    q["2001-03-15"] = 50.0                                  # belongs to water year 2001
    q["2001-11-20"] = 70.0                                  # November -> water year 2002
    ams = annual_maxima(q)
    assert ams.loc[2001, "peak"] == 50.0
    assert ams.loc[2002, "peak"] == 70.0