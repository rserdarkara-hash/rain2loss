import numpy as np
import pandas as pd
import pytest

from rain2loss.hbv import (DEFAULT_PARAMS, STATE_NAMES, run_hbv, triangular_weights,
                           water_balance_error)


@pytest.fixture
def forcing():
    """Three years of random but realistic daily weather."""
    rng = np.random.default_rng(0)
    n = 3 * 365
    days = pd.date_range("2000-01-01", periods=n, freq="D")
    temp = 8 - 10 * np.cos(2 * np.pi * days.dayofyear / 365) + rng.normal(0, 3, n)
    prec = rng.gamma(0.6, 6, n) * (rng.random(n) < 0.5)
    pet = np.clip(0.15 * (temp + 5), 0, None)
    return pd.DataFrame({"prec": prec, "temp": temp, "pet": pet}, index=days)


@pytest.mark.parametrize("maxbas", [1.0, 2.0, 2.5, 4.7, 7.0])
def test_routing_weights_sum_to_one(maxbas):
    assert triangular_weights(maxbas).sum() == pytest.approx(1.0)


@pytest.mark.parametrize("zones", [(0.0,), (-200.0, 0.0, 200.0)])
def test_water_balance_closes(forcing, zones):
    states0 = [10.0, 0.0, 100.0, 5.0, 20.0]
    sim = run_hbv(forcing, DEFAULT_PARAMS, zones, states0=states0)
    assert abs(water_balance_error(sim, states0)) < 1e-6


def test_storages_never_negative(forcing):
    sim = run_hbv(forcing, DEFAULT_PARAMS, (-200.0, 0.0, 200.0))
    assert (sim[STATE_NAMES] >= 0).all().all()
    assert (sim["q_sim"] >= 0).all()