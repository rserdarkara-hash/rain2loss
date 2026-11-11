import pytest

from rain2loss.loss import average_annual_loss, damage_ratio


def test_damage_curve_end_points():
    assert damage_ratio(0.0) == 0.0
    assert damage_ratio(1.0) == pytest.approx(0.40)
    assert damage_ratio(10.0) == 1.0          # deeper than 6 m: total loss


def test_aal_by_hand():
    # losses 10 at T=10 and 20 at T=100: trapezoid (10 + 20) / 2 * (0.1 - 0.01) = 1.35
    assert average_annual_loss([10, 100], [10.0, 20.0], tail=None) == pytest.approx(1.35)
    # the flat tail adds 20 * 0.01 = 0.2; the order of the inputs does not matter
    assert average_annual_loss([100, 10], [20.0, 10.0]) == pytest.approx(1.55)


def test_protection_removes_frequent_losses():
    # loss at T=10 set to 0: (0 + 20) / 2 * 0.09 = 0.9
    assert average_annual_loss([10, 100], [10.0, 20.0], protection_rp=10, tail=None) == pytest.approx(0.9)