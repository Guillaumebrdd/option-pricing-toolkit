"""
Sanity tests for the Black-Scholes module.

Run with:  pytest tests/ -v
"""
import numpy as np
import pytest

from src import black_scholes as bs


# Standard test point (round numbers, ATM call)
S0, K, T, r, sigma = 100.0, 100.0, 1.0, 0.05, 0.20


def test_call_put_parity():
    """C - P = S0 * exp(-qT) - K * exp(-rT). Must hold exactly."""
    c = bs.call_price(S0, K, T, r, sigma)
    p = bs.put_price(S0, K, T, r, sigma)
    lhs = c - p
    rhs = S0 - K * np.exp(-r * T)
    assert abs(lhs - rhs) < 1e-10


def test_call_price_known_value():
    """Reference value from Hull Table 15.2 / standard texts."""
    c = bs.call_price(S0, K, T, r, sigma)
    assert abs(c - 10.4506) < 1e-3   # ~10.45


def test_delta_call_in_zero_one():
    """A call's delta must lie in (0, 1)."""
    d = bs.delta(S0, K, T, r, sigma, option_type="call")
    assert 0.0 < d < 1.0


def test_gamma_positive():
    """Gamma is positive for both calls and puts."""
    g = bs.gamma(S0, K, T, r, sigma)
    assert g > 0.0


def test_finite_difference_delta():
    """Analytical delta must match a centered finite difference to ~1e-5."""
    h = 1e-4
    fd_delta = (bs.call_price(S0 + h, K, T, r, sigma)
                - bs.call_price(S0 - h, K, T, r, sigma)) / (2 * h)
    an_delta = bs.delta(S0, K, T, r, sigma, option_type="call")
    assert abs(fd_delta - an_delta) < 1e-5


def test_finite_difference_vega():
    """Analytical vega must match a centered finite difference."""
    h = 1e-5
    fd_vega = (bs.call_price(S0, K, T, r, sigma + h)
               - bs.call_price(S0, K, T, r, sigma - h)) / (2 * h)
    an_vega = bs.vega(S0, K, T, r, sigma)
    assert abs(fd_vega - an_vega) < 1e-4


def test_vectorisation():
    """The functions must accept numpy arrays."""
    spots = np.array([80.0, 100.0, 120.0])
    prices = bs.call_price(spots, K, T, r, sigma)
    assert prices.shape == (3,)
    assert np.all(np.diff(prices) > 0)   # call price increases with S0


def test_invalid_inputs():
    """Negative or zero parameters should raise."""
    with pytest.raises(ValueError):
        bs.call_price(S0, K, T=-1.0, r=r, sigma=sigma)
    with pytest.raises(ValueError):
        bs.call_price(S0, K, T, r, sigma=0.0)
