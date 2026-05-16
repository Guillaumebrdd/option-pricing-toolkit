# Project Roadmap — Actionable checklist

Each phase ends with a concrete deliverable. Check items off as you finish.

---

## Phase 0 — Python ecosystem onboarding  (3-5 days)

You already know C and theoretical CS — focus on the *libraries*, not the syntax.

- [ ] Install Python 3.11+ via `pyenv` or Anaconda
- [ ] Set up `venv`, VS Code with Python + Jupyter extensions
- [ ] `numpy`: read the official Absolute Beginner tutorial; do `numpy-100` exercises 1-40
- [ ] `scipy.stats`: practice with `norm.cdf`, `norm.pdf`, `norm.ppf`
- [ ] `scipy.optimize`: practice with `brentq` (root finding) and `minimize`
- [ ] `matplotlib`: practice `plot`, `subplot`, `plot_surface`
- [ ] `pandas`: read official "10 minutes to pandas" and try indexing/groupby
- [ ] **Deliverable:** notebook generating 100,000 normal samples, plotting histogram vs theoretical PDF, computing mean & variance with 95% CI

---

## Phase 1 — Discrete-to-continuous derivation  (week 1)

Leverage your measure-theory + martingale background. This is *your* phase.

- [ ] Read Shreve vol. II chapter 4 (skim — you know measure theory already)
- [ ] Read Hull chapters 13-15 for market intuition
- [ ] Rework on paper: CRR risk-neutral probability `q = (R - d)/(u - d)`
- [ ] Prove the convergence of `log(S_N / S_0)` to `N((r - σ²/2)T, σ²T)` under Q
- [ ] Derive the Black-Scholes formula from the log-normal expectation
- [ ] Complete the write-up `docs/from_discrete_to_continuous.md` (template provided)
- [ ] Implement the convergence numerically — first version of `binomial.py`
- [ ] **Deliverable:** 10-15 page Markdown write-up + reproducible convergence plot

---

## Phase 2 — Black-Scholes, CRR, Leisen-Reimer  (weeks 2-3)

- [ ] `src/black_scholes.py` — closed form prices (call/put) ✅ *starter provided*
- [ ] `src/black_scholes.py` — Greeks: Delta, Gamma, Vega, Theta ✅ *starter provided*
- [ ] Add Rho, Vanna, Charm, Vomma to `black_scholes.py`
- [ ] `src/binomial.py` — CRR pricer, vectorised
- [ ] `src/binomial.py` — American option support (early exercise check at each node)
- [ ] `src/binomial.py` — Leisen-Reimer parametrisation
- [ ] `notebooks/02_crr_vs_bs.ipynb` — convergence plot, oscillation visualisation
- [ ] `notebooks/03_american_premium.ipynb` — American put vs European put as a function of moneyness
- [ ] Tests: put-call parity, Greeks vs finite differences, CRR → BS within 1e-3 at N=5000
- [ ] **Deliverable:** working `bs`, `binomial` modules + 2 notebooks + passing test suite

---

## Phase 3 — Monte Carlo + exotics  (week 4)

- [ ] `src/monte_carlo.py` — vanilla MC for European options with 95% CI
- [ ] Antithetic variates implementation
- [ ] Control variates (use S_T itself as control)
- [ ] `notebooks/04_mc_variance_reduction.ipynb` — log-log plot of variance vs N for plain / antithetic / control
- [ ] Asian option (arithmetic average) pricing — only MC can do this
- [ ] Barrier option (knock-out call) pricing
- [ ] Longstaff-Schwartz for American options (regression on Laguerre polynomials)
- [ ] Compare LSM vs binomial for American put — should match
- [ ] **Deliverable:** `monte_carlo` module + 1 notebook + exotic options demo

---

## Phase 4 — Implied vol + SVI calibration  (week 5)

- [ ] `src/implied_vol.py` — Brent and Newton methods
- [ ] Compare convergence speed and robustness
- [ ] Fetch SPY option chain via `yfinance`
- [ ] Clean data: filter OTM only, remove zero-bid quotes, mid-price
- [ ] Plot smiles per maturity → see the shape
- [ ] `src/svi.py` — SVI parametrisation `w(k) = a + b(ρ(k-m) + √((k-m)² + σ²))`
- [ ] Calibrate via `scipy.optimize.minimize` (least squares on each maturity)
- [ ] Check static no-arbitrage: butterfly (g(k) ≥ 0) and calendar (w monotone in T)
- [ ] `notebooks/05_volatility_surface.ipynb` — 3D surface plot, observed vs fitted
- [ ] **Deliverable:** working SVI calibration on real SPY data + arbitrage report

---

## Phase 5 — Beyond Black-Scholes  (week 6, optional but recommended)

Pick **one** of these — depth over breadth:

### Option A: Heston stochastic volatility
- [ ] MC with Andersen QE scheme (handles variance ≥ 0)
- [ ] Semi-analytic pricing via Carr-Madan FFT
- [ ] Calibration to a few maturities — compare to SVI fit

### Option B: Merton jump-diffusion
- [ ] Analytical price as Poisson-weighted BS series
- [ ] MC implementation with compound Poisson jumps
- [ ] Show how the model naturally produces a smile

### Option C: Local volatility (Dupire)
- [ ] Compute local vol `σ_loc(K, T)` from the SVI surface
- [ ] MC pricing under local vol dynamics
- [ ] Verify the local-vol prices match the input vanilla prices

- [ ] **Deliverable:** 1 chosen model implemented + comparison notebook

---

## Phase 6 — Polish & publication  (week 6-7)

- [ ] Achieve ≥ 70% test coverage (`pytest --cov`)
- [ ] Add type hints to all public functions
- [ ] Add NumPy-style docstrings to all public functions
- [ ] Write `notebooks/00_overview.ipynb` — 10-cell project showcase (read first)
- [ ] Fill in the README "Key results" section with real numbers and figures
- [ ] Add 4-6 high-quality figures to `results/figures/`
- [ ] Final pass on README: fix typos, check all formulas render, verify install instructions
- [ ] Push to GitHub with a clean commit history
- [ ] Pin the repository to your GitHub profile
- [ ] **Deliverable:** publishable repo, ready to link in applications

---

## Recurring practices (do these throughout)

- **Commit often** — small, descriptive commits. Avoid one giant "final" commit.
- **Write tests as you code** — don't leave testing to the end.
- **Keep notebooks reproducible** — set a random seed, run top-to-bottom before committing.
- **AI assistants as tutors, not generators** — you must be able to re-explain every line in interviews.
- **Read in parallel** — Hull, Shreve, Joshi, Glasserman, Gatheral. 30 min/day adds up.
