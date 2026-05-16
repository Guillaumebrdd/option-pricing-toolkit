"""
Black-Scholes-Merton closed-form pricing and Greeks for European options.

All inputs are scalars or numpy arrays; the module is fully vectorized.

Notation
--------
S0    : spot price of the underlying
K     : strike price
T     : time to maturity (in years, T > 0)
r     : risk-free rate (continuously compounded, annualised)
sigma : volatility (annualised, sigma > 0)
q     : continuous dividend yield (default 0)

Reference: Hull (2018), Chapter 15. Joshi (2008), Chapter 6.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

ArrayLike = float | np.ndarray


# ---------------------------------------------------------------------------
# Core auxiliary quantities
# ---------------------------------------------------------------------------

def _d1_d2(
    S0: ArrayLike,
    K: ArrayLike,
    T: ArrayLike,
    r: ArrayLike,
    sigma: ArrayLike,
    q: ArrayLike = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the Black-Scholes auxiliary quantities d1 and d2."""
    S0, K, T, r, sigma, q = map(np.asarray, (S0, K, T, r, sigma, q))
    if np.any(T <= 0):
        raise ValueError("Time to maturity T must be strictly positive.")
    if np.any(sigma <= 0):
        raise ValueError("Volatility sigma must be strictly positive.")
    if np.any(S0 <= 0) or np.any(K <= 0):
        raise ValueError("Spot S0 and strike K must be strictly positive.")

    sqrtT = np.sqrt(T)
    d1 = (np.log(S0 / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * sqrtT)
    d2 = d1 - sigma * sqrtT
    return d1, d2


# ---------------------------------------------------------------------------
# Prices
# ---------------------------------------------------------------------------

def call_price(
    S0: ArrayLike, K: ArrayLike, T: ArrayLike,
    r: ArrayLike, sigma: ArrayLike, q: ArrayLike = 0.0,
) -> ArrayLike:
    """Price of a European call option under Black-Scholes-Merton."""
    d1, d2 = _d1_d2(S0, K, T, r, sigma, q)
    return np.asarray(S0) * np.exp(-np.asarray(q) * np.asarray(T)) * norm.cdf(d1) \
        - np.asarray(K) * np.exp(-np.asarray(r) * np.asarray(T)) * norm.cdf(d2)


def put_price(
    S0: ArrayLike, K: ArrayLike, T: ArrayLike,
    r: ArrayLike, sigma: ArrayLike, q: ArrayLike = 0.0,
) -> ArrayLike:
    """Price of a European put option under Black-Scholes-Merton."""
    d1, d2 = _d1_d2(S0, K, T, r, sigma, q)
    return np.asarray(K) * np.exp(-np.asarray(r) * np.asarray(T)) * norm.cdf(-d2) \
        - np.asarray(S0) * np.exp(-np.asarray(q) * np.asarray(T)) * norm.cdf(-d1)


# ---------------------------------------------------------------------------
# Greeks (analytical)
# ---------------------------------------------------------------------------
# NOTE: for a fully-fledged toolkit, also implement: rho, vanna, charm, vomma.
# Left as exercise — the patterns below should make the extension trivial.

def delta(
    S0: ArrayLike, K: ArrayLike, T: ArrayLike,
    r: ArrayLike, sigma: ArrayLike, q: ArrayLike = 0.0,
    option_type: str = "call",
) -> ArrayLike:
    """Delta = ∂Price / ∂S0."""
    d1, _ = _d1_d2(S0, K, T, r, sigma, q)
    disc_q = np.exp(-np.asarray(q) * np.asarray(T))
    if option_type == "call":
        return disc_q * norm.cdf(d1)
    elif option_type == "put":
        return disc_q * (norm.cdf(d1) - 1.0)
    raise ValueError("option_type must be 'call' or 'put'.")


def gamma(
    S0: ArrayLike, K: ArrayLike, T: ArrayLike,
    r: ArrayLike, sigma: ArrayLike, q: ArrayLike = 0.0,
) -> ArrayLike:
    """Gamma = ∂²Price / ∂S0² — identical for call and put."""
    d1, _ = _d1_d2(S0, K, T, r, sigma, q)
    S0, T, sigma, q = map(np.asarray, (S0, T, sigma, q))
    return np.exp(-q * T) * norm.pdf(d1) / (S0 * sigma * np.sqrt(T))


def vega(
    S0: ArrayLike, K: ArrayLike, T: ArrayLike,
    r: ArrayLike, sigma: ArrayLike, q: ArrayLike = 0.0,
) -> ArrayLike:
    """Vega = ∂Price / ∂σ — identical for call and put. Per 1.0 change in σ."""
    d1, _ = _d1_d2(S0, K, T, r, sigma, q)
    S0, T, q = map(np.asarray, (S0, T, q))
    return S0 * np.exp(-q * T) * norm.pdf(d1) * np.sqrt(T)


def theta(
    S0: ArrayLike, K: ArrayLike, T: ArrayLike,
    r: ArrayLike, sigma: ArrayLike, q: ArrayLike = 0.0,
    option_type: str = "call",
) -> ArrayLike:
    """Theta = ∂Price / ∂t. Annualised (divide by 365 for daily theta)."""
    d1, d2 = _d1_d2(S0, K, T, r, sigma, q)
    S0, K, T, r, sigma, q = map(np.asarray, (S0, K, T, r, sigma, q))
    term1 = -S0 * np.exp(-q * T) * norm.pdf(d1) * sigma / (2 * np.sqrt(T))
    if option_type == "call":
        return term1 - r * K * np.exp(-r * T) * norm.cdf(d2) \
            + q * S0 * np.exp(-q * T) * norm.cdf(d1)
    elif option_type == "put":
        return term1 + r * K * np.exp(-r * T) * norm.cdf(-d2) \
            - q * S0 * np.exp(-q * T) * norm.cdf(-d1)
    raise ValueError("option_type must be 'call' or 'put'.")
