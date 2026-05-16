# Option Pricing & Greeks Toolkit

A Python library implementing and cross-validating **analytical**, **lattice-based**, and **Monte Carlo** methods for European and American option pricing, with Greeks computation and implied-volatility surface construction from market data.

> **Status**: 🚧 Work in progress — see [Roadmap](#roadmap).

## Motivation

This project derives the Black-Scholes formula as the continuous-time limit of the Cox-Ross-Rubinstein binomial model, then implements three independent pricing methods (closed-form, lattice, Monte Carlo) and cross-validates them against each other. The goal is twofold: (i) build a rigorous reference toolkit for European, American, and path-dependent options, and (ii) calibrate a parametric volatility surface (SVI) to real market data and verify the resulting surface is arbitrage-free.

## Mathematical background

A short self-contained derivation is provided in [`docs/from_discrete_to_continuous.md`](docs/from_discrete_to_continuous.md). The key result is that under the risk-neutral measure $\mathbb{Q}$, with CRR parameters $u = e^{\sigma\sqrt{\Delta t}}$ and $d = 1/u$, the log-price $\log(S_N / S_0)$ converges in distribution to $\mathcal{N}\!\left((r - \sigma^2/2)T,\ \sigma^2 T\right)$ as $N \to \infty$, yielding the Black-Scholes formula by direct expectation:

$$C = S_0 \Phi(d_1) - K e^{-rT} \Phi(d_2), \quad d_{1,2} = \frac{\log(S_0/K) + (r \pm \sigma^2/2) T}{\sigma \sqrt{T}}.$$

## Methods implemented

| Method | European | American | Path-dependent | Greeks |
|---|---|---|---|---|
| Black-Scholes (closed-form) | ✅ | ❌ | ❌ | ✅ analytical |
| CRR / Leisen-Reimer binomial | ✅ | ✅ | ⚠️ limited | ✅ finite diff. |
| Monte Carlo (vanilla + variance reduction) | ✅ | ✅ via Longstaff-Schwartz | ✅ | ✅ pathwise / LRM |

## Project structure

```
option-pricing-toolkit/
├── README.md
├── requirements.txt
├── src/
│   ├── black_scholes.py       # closed-form pricing + Greeks
│   ├── binomial.py            # CRR + Leisen-Reimer
│   ├── monte_carlo.py         # MC + variance reduction + Longstaff-Schwartz
│   ├── implied_vol.py         # IV inversion (Brent, Newton)
│   ├── svi.py                 # SVI calibration + no-arbitrage checks
│   └── utils.py
├── tests/                     # pytest test suite
├── notebooks/                 # research and exploration
├── docs/                      # mathematical write-ups
└── results/figures/           # generated plots
```

## Installation

```bash
git clone https://github.com/<your-username>/option-pricing-toolkit.git
cd option-pricing-toolkit
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pytest                                                # run the test suite
```

## Roadmap

- [ ] **Phase 1** — Discrete-to-continuous derivation (write-up + reproducible code)
- [ ] **Phase 2** — Black-Scholes, CRR, Leisen-Reimer + cross-validation
- [ ] **Phase 3** — Monte Carlo, variance reduction, Longstaff-Schwartz, exotic options
- [ ] **Phase 4** — Implied volatility, SVI calibration on SPY options, no-arbitrage checks
- [ ] **Phase 5** — Heston (or Merton jumps) + Carr-Madan Fourier pricing
- [ ] **Phase 6** — Polish: tests ≥ 70% coverage, type hints, docstrings, showcase notebook

## Key results

*(To be filled in as phases complete — example placeholders below.)*

- Convergence of CRR to Black-Scholes: see `results/figures/crr_convergence.png`
- Monte Carlo variance reduction: factor ~30× speedup with antithetic + control variates
- SVI fit RMSE on SPY surface (2026-XX-XX snapshot): < X bps across all maturities

## References

- Hull, J. — *Options, Futures, and Other Derivatives* (10th ed.)
- Shreve, S. — *Stochastic Calculus for Finance II*
- Joshi, M. — *The Concepts and Practice of Mathematical Finance*
- Glasserman, P. — *Monte Carlo Methods in Financial Engineering*
- Gatheral, J. — *The Volatility Surface*
- Gatheral & Jacquier (2014) — *Arbitrage-free SVI volatility surfaces*

## License

MIT
