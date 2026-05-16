# From Discrete to Continuous: A Derivation of Black-Scholes

## Goal

Derive the Black-Scholes call price formula as the continuous-time limit of the Cox-Ross-Rubinstein (CRR) binomial model, using only discrete-time martingale theory and the central limit theorem — without invoking Itô calculus.

The Itô / PDE derivation is given in [Appendix A](#appendix-a-the-itô-derivation) for completeness.

---

## 1. The CRR binomial model in discrete time

Fix a maturity $T > 0$ and partition $[0, T]$ into $N$ steps of size $\Delta t = T/N$. At each step, the underlying $S$ moves by a factor $u > 1$ (up) or $d < 1$ (down). We assume $0 < d < e^{r\Delta t} < u$ to rule out arbitrage. Denote $R = e^{r\Delta t}$ the per-step risk-free growth.

**Risk-neutral measure $\mathbb{Q}$.** There exists a unique probability $q \in (0,1)$ under which the discounted price $\widetilde{S}_n = R^{-n} S_n$ is a martingale:

$$\mathbb{E}^{\mathbb{Q}}[S_{n+1} \mid \mathcal{F}_n] = R \, S_n
\quad \Longleftrightarrow \quad
q u + (1 - q) d = R
\quad \Longleftrightarrow \quad
q = \frac{R - d}{u - d}.$$

**Pricing a European claim.** For a claim with payoff $h(S_N)$ at $T$, the no-arbitrage price at time 0 is

$$V_0 = R^{-N} \, \mathbb{E}^{\mathbb{Q}}[h(S_N)].$$

For a European call, $h(s) = (s - K)^+$.

---

## 2. CRR parametrisation

Among all choices of $(u, d)$ matching given moments, the CRR choice is

$$u = e^{\sigma \sqrt{\Delta t}}, \qquad d = e^{-\sigma \sqrt{\Delta t}} = 1/u.$$

This choice has the property that $\log u = -\log d = \sigma \sqrt{\Delta t}$, which makes the log-returns symmetric. Other parametrisations (Jarrow-Rudd, Tian, Leisen-Reimer) exist; CRR is the simplest.

---

## 3. The limit theorem

Let $X_k = \log(S_k / S_{k-1}) \in \{\log u, \log d\} = \{\sigma\sqrt{\Delta t}, -\sigma\sqrt{\Delta t}\}$. Under $\mathbb{Q}$, the $X_k$ are i.i.d. with

$$\mathbb{E}^{\mathbb{Q}}[X_k] = q \sigma\sqrt{\Delta t} - (1-q) \sigma\sqrt{\Delta t} = (2q - 1) \sigma\sqrt{\Delta t},$$

$$\operatorname{Var}^{\mathbb{Q}}[X_k] = \sigma^2 \Delta t - \big((2q-1)\sigma\sqrt{\Delta t}\big)^2 = \sigma^2 \Delta t \big(1 - (2q-1)^2\big) = 4 q (1-q) \sigma^2 \Delta t.$$

**Expanding $q$ as $\Delta t \to 0$.** Using $R = 1 + r\Delta t + O(\Delta t^2)$ and Taylor-expanding $u, d$:

$$q = \frac{R - d}{u - d}
= \frac{1 + r\Delta t - 1 + \sigma\sqrt{\Delta t} - \tfrac{1}{2}\sigma^2 \Delta t + O(\Delta t^{3/2})}{2\sigma\sqrt{\Delta t} + O(\Delta t^{3/2})}
= \frac{1}{2} + \frac{1}{2}\!\left(\frac{r - \tfrac{1}{2}\sigma^2}{\sigma}\right) \sqrt{\Delta t} + O(\Delta t).$$

So $2q - 1 = \left(\frac{r - \sigma^2/2}{\sigma}\right) \sqrt{\Delta t} + O(\Delta t)$, and:

$$N \cdot \mathbb{E}^{\mathbb{Q}}[X_k] = N (2q-1) \sigma \sqrt{\Delta t} \xrightarrow{N\to\infty} \left(r - \tfrac{\sigma^2}{2}\right) T,$$

$$N \cdot \operatorname{Var}^{\mathbb{Q}}[X_k] = 4 q (1-q) \sigma^2 N \Delta t \xrightarrow{N\to\infty} \sigma^2 T.$$

By the **Lindeberg-Lévy CLT** applied to the triangular array $\{X_k^{(N)}\}_{k=1,\dots,N}$:

$$\log\!\frac{S_N}{S_0} = \sum_{k=1}^N X_k \xrightarrow[N\to\infty]{\text{(d)}} \mathcal{N}\!\left((r - \tfrac{\sigma^2}{2}) T,\ \sigma^2 T\right) \quad \text{under } \mathbb{Q}.$$

---

## 4. The Black-Scholes formula

Under $\mathbb{Q}$, $S_T = S_0 \exp(Y)$ where $Y \sim \mathcal{N}(\mu, \sigma^2 T)$ with $\mu = (r - \sigma^2/2) T$. The discount factor is $R^{-N} = e^{-rT}$, so:

$$C_0 = e^{-rT} \, \mathbb{E}^{\mathbb{Q}}\!\left[(S_0 e^Y - K)^+\right].$$

This expectation is computed by splitting the indicator $\{S_0 e^Y > K\} = \{Y > \log(K/S_0)\}$ and applying the standard log-normal expectation identity. After standardising $Z = (Y - \mu)/(\sigma\sqrt{T}) \sim \mathcal{N}(0,1)$, one obtains:

$$\boxed{\ C_0 = S_0 \, \Phi(d_1) - K e^{-rT} \, \Phi(d_2), \qquad d_{1,2} = \frac{\log(S_0/K) + (r \pm \sigma^2/2) T}{\sigma \sqrt{T}}\ }$$

The put follows by put-call parity: $P_0 = C_0 - S_0 + K e^{-rT}$.

---

## 5. Numerical verification

The convergence $C_0^{\text{CRR}}(N) \to C_0^{\text{BS}}$ is implemented in `notebooks/01_crr_to_bs_convergence.ipynb`. The convergence is of order $O(1/N)$ but oscillatory; Leisen-Reimer's parametrisation removes the oscillation and accelerates convergence to $O(1/N^2)$.

---

## Appendix A: the Itô derivation

For reference, the standard route. Assume $S_t$ follows geometric Brownian motion under $\mathbb{Q}$:

$$dS_t = r S_t \, dt + \sigma S_t \, dW_t^{\mathbb{Q}}.$$

By Itô's lemma applied to $f(S, t) = \log S$:

$$d(\log S_t) = \left(r - \tfrac{\sigma^2}{2}\right) dt + \sigma \, dW_t^{\mathbb{Q}},$$

so $\log(S_T/S_0) \sim \mathcal{N}\!\big((r - \sigma^2/2) T,\ \sigma^2 T\big)$ under $\mathbb{Q}$, recovering the same law as the discrete limit. The call price is then computed identically.

Alternatively, one shows that any $C(S, t) \in C^{2,1}$ satisfying the **Black-Scholes PDE**

$$\frac{\partial C}{\partial t} + r S \frac{\partial C}{\partial S} + \tfrac{1}{2} \sigma^2 S^2 \frac{\partial^2 C}{\partial S^2} - r C = 0$$

with terminal condition $C(S, T) = (S - K)^+$ admits the same closed form, via the Feynman-Kac representation. Both routes coincide with the discrete-limit derivation given above.

---

## References

- Cox, J. C., Ross, S. A., Rubinstein, M. (1979). *Option pricing: a simplified approach*. JFE 7, 229-263.
- Shreve, S. (2004). *Stochastic Calculus for Finance, vols. I & II*. Springer. (Vol. I, ch. 2 for the discrete model; vol. II, ch. 4-5 for the continuous limit and Itô.)
- Hull, J. (2018). *Options, Futures, and Other Derivatives*, 10th ed., ch. 13-15.
- Joshi, M. S. (2008). *The Concepts and Practice of Mathematical Finance*, 2nd ed., ch. 3-6.
