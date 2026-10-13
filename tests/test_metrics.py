import numpy as np
import pytest

from rain2loss.metrics import kge, nse, pbias


def test_perfect_simulation():
    obs = np.array([1.0, 3.0, 2.0, 5.0])
    assert nse(obs, obs) == pytest.approx(1.0)
    assert kge(obs, obs) == pytest.approx(1.0)
    assert pbias(obs, obs) == pytest.approx(0.0)


def test_mean_flow_benchmark():
    obs = np.array([1.0, 3.0, 2.0, 5.0])
    mean_sim = np.full_like(obs, obs.mean())
    assert nse(mean_sim, obs) == pytest.approx(0.0)


def test_kge_components():
    obs = np.array([1.0, 3.0, 2.0, 5.0])
    value, r, alpha, beta = kge(2 * obs, obs, components=True)
    assert r == pytest.approx(1.0)
    assert alpha == pytest.approx(2.0)
    assert beta == pytest.approx(2.0)
    assert value == pytest.approx(1 - np.sqrt(2))


def test_pbias_sign():
    obs = np.array([1.0, 2.0, 3.0])
    assert pbias(obs * 1.1, obs) == pytest.approx(10.0)   # positive = overestimation


def test_missing_values_are_ignored():
    obs = np.array([1.0, np.nan, 3.0, 4.0])
    sim = np.array([1.0, 2.0, 3.0, np.nan])
    assert nse(sim, obs) == pytest.approx(1.0)