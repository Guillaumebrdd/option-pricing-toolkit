# From Discrete to Continuous: A Derivation of Black-Scholes

## Goal

Derive the Black-Scholes formula as the continuous-time limit of the Cox-Ross-Rubinstein (CRR) binomial model, using only discrete-time martingale theory and the central limit theorem for triangular arrays — without invoking Itô calculus.

The Itô / PDE derivation is given in [Appendix A](#appendix-a-the-ito-derivation) as a parallel route.

---

## 1. Probabilistic framework of the CRR model
Fix a maturity $T > 0$ and partition $[0,T]$ into $n$ steps of length $\Delta t = T/n$. We work on a probability space $(\Omega, \mathcal{F}, \mathbb{P})$ carrying an i.i.d. sequence $(\xi_k)_{k=1,\ldots,n}$ taking values in $\{u, d\}$, with

$$u = e^{\sigma \sqrt{\Delta t}}, \qquad d = e^{-\sigma \sqrt{\Delta t}} = 1/u, \qquad R = e^{r \Delta t},$$

and we assume the no-arbitrage condition $0 < d < R < u$ (satisfied for $n$ large enough as soon as $\sigma > 0$). Under the historical measure $\mathbb{P}$, write $p = \mathbb{P}(\xi_k = u) \in (0,1)$.

The underlying price is defined recursively by

$$S_0 \text{ given}, \qquad S_{k+1} = S_k \cdot \xi_{k+1}, \qquad k = 0, \ldots, n-1,$$

and we equip $\Omega$ with the natural filtration

$$\mathcal{F}_{k} = \sigma(\xi_1, \ldots, \xi_k) = \sigma(S_1, \ldots, S_k), \qquad \mathcal{F}_{0} = \{\emptyset, \Omega\}.$$

In particular, $S\_{k}$ is $\mathcal{F}\_{k}$-measurable, while $\xi\_{k+1}$ is independent of $\mathcal{F}\_{k}$.

### Construction of the risk-neutral measure

We seek an equivalent probability measure $\mathbb{Q} \sim \mathbb{P}$ under which the discounted price $\widetilde{S}\_{k} = R^{-k} S\_{k}$ is a $(\mathcal{F}\_{k})$-martingale, and under which the $(\xi_k)$ remain i.i.d. with $\mathbb{Q}(\xi_k = u) = q \in (0,1)$.

The martingale condition

$$\mathbb{E}^{\mathbb{Q}}[\widetilde{S}_{k+1} \mid \mathcal{F}_{k}] = \widetilde{S}_{k}$$

rewrites, after multiplication by $R^{k+1}$, as

$$\mathbb{E}^{\mathbb{Q}}[S_{k+1} \mid \mathcal{F}_{k}] = R \cdot S_{k}.$$

Using the $\mathcal{F}\_{k}$-measurability of $S\_{k}$ (which factors out of the conditional expectation) and the independence of $\xi\_{k+1}$ from $\mathcal{F}\_{k}$ (which collapses the conditional expectation into a plain one):

$$\mathbb{E}^{\mathbb{Q}}[S_{k+1} \mid \mathcal{F}_{k}] = S_{k} \cdot \mathbb{E}^{\mathbb{Q}}[\xi_{k+1}] = S_{k} \cdot \big(q u + (1-q) d\big).$$

The martingale condition is therefore equivalent to $qu + (1-q)d = R$, which yields

$$q = \frac{R - d}{u - d}.$$

The arbitrage-free condition $d < R < u$ ensures $q \in (0,1)$, so $\mathbb{Q}$ is indeed an equivalent probability measure. Its uniqueness follows from the fact that the condition above determines $q$ uniquely: the CRR market is **complete**.

> **Remark.** Here we assume that $\mathbb{Q}$ exists due to the no-arbitrage assumption, as in discrete time. The continuous-time analogue of this result, the Fundamental Theorem of Asset Pricing, is significantly more delicate. The rigorous version (Delbaen-Schachermayer, 1994) establishes the equivalence between *No Free Lunch with Vanishing Risk* (NFLVR) and the existence of an equivalent local-martingale measure.

### Pricing a European claim

For a claim with payoff $h(S_n)$ at maturity, the no-arbitrage price at time $0$ is

$$V_0 = R^{-n} \, \mathbb{E}^{\mathbb{Q}}[h(S_n)] = e^{-rT} \, \mathbb{E}^{\mathbb{Q}}[h(S_n)],$$

by backward induction on the one-step pricing relation and the tower property of conditional expectations.

This is the formula we will study in the limit $n \to \infty$.

---

## 2. Centered reparametrisation

To study the asymptotic distribution of $S_n$ as $n \to \infty$, we face two technical obstacles. First, $S_n = S_0 \prod_{k=1}^n \xi_k$ is a *product* of i.i.d. random variables, whereas standard limit theorems are stated for *sums*. Second, taking logarithms gives

$$\log \frac{S_n}{S_0} = \sum_{k=1}^n \log \xi_k, \qquad \log \xi_k \in \{\sigma\sqrt{\Delta t},\, -\sigma\sqrt{\Delta t}\},$$

but the summands shrink to zero as $\Delta t \to 0$, so the classical CLT does not apply to them as stated.

The remedy is to separate the *amplitude* of each increment (which depends on $n$) from its *direction* (which does not). Define

$$Z_k = \frac{\log \xi_k}{\sigma \sqrt{\Delta t}} \in \{-1, +1\},$$

so that $\xi_k = e^{\sigma\sqrt{\Delta t}\,Z_k}$. As a deterministic transformation of i.i.d. variables, the sequence $(Z_k)$ remains i.i.d. under $\mathbb{Q}$, with

$$\mathbb{Q}(Z_k = +1) = q, \qquad \mathbb{Q}(Z_k = -1) = 1 - q.$$

Substituting $\sqrt{\Delta t} = \sqrt{T}/\sqrt{n}$, the terminal log-price rewrites as

$$\log \frac{S_n}{S_0} = \sigma \sqrt{\Delta t} \sum_{k=1}^n Z_k = \sigma \sqrt{T} \cdot \widetilde{Z}_n, \qquad \text{where } \widetilde{Z}_n := \frac{1}{\sqrt{n}} \sum_{k=1}^n Z_k.$$

The renormalised sum $\widetilde{Z}_n$ is now in the standard form to which the central limit theorem applies. Identifying its limiting distribution under $\mathbb{Q}$ — the object of the next two sections — will yield the limiting distribution of $\log(S_n/S_0)$, and hence of $S_n$, via the linear relation above.

---


## 3. Asymptotic moments of $Z_k$ under $\mathbb{Q}$

To apply a central limit theorem to $\widetilde{Z}_n$, we need the limiting behaviour of the first two moments of $Z_k$ under $\mathbb{Q}$ as $\Delta t \to 0$. The non-trivial point is that $q$ depends on $\Delta t$, so these moments depend on $n$ and must be computed asymptotically.

### Asymptotic expansion of $q$

We write $\varepsilon := \sqrt{\Delta t}$ and expand $u = e^{\sigma\varepsilon}$, $d = e^{-\sigma\varepsilon}$, $R = e^{r\varepsilon^2}$ to order $\varepsilon^2$:

$$u = 1 + \sigma\varepsilon + \tfrac{1}{2}\sigma^2\varepsilon^2 + O(\varepsilon^3), \quad d = 1 - \sigma\varepsilon + \tfrac{1}{2}\sigma^2\varepsilon^2 + O(\varepsilon^3), \quad R = 1 + r\varepsilon^2 + O(\varepsilon^4).$$

The numerator and denominator of $q = (R - d)/(u - d)$ then read

$$R - d = \sigma\varepsilon + \left(r - \tfrac{1}{2}\sigma^2\right)\varepsilon^2 + O(\varepsilon^3), \qquad u - d = 2\sigma\varepsilon + O(\varepsilon^3),$$

so that

$$q = \frac{1}{2} + \frac{1}{2}\cdot\frac{r - \sigma^2/2}{\sigma}\,\varepsilon + O(\varepsilon^2).$$

The key feature is that $q \to 1/2$ as $\Delta t \to 0$, but the *speed* of this convergence — the $\varepsilon$-coefficient — encodes the drift $(r - \sigma^2/2)$.

### Moments of $Z_k$

Since $Z_k \in \{-1, +1\}$ with $\mathbb{Q}(Z_k = 1) = q$:

$$\mathbb{E}^{\mathbb{Q}}[Z_k] = 2q - 1 = \frac{r - \sigma^2/2}{\sigma}\,\sqrt{\Delta t} + O(\Delta t),$$

$$\mathrm{Var}^{\mathbb{Q}}(Z_k) = 1 - (2q - 1)^2 = 4q(1-q) = 1 + O(\Delta t).$$

### Moments of $\widetilde{Z}_n$

By linearity of expectation and independence of the $(Z_k)$:

$$\mathbb{E}^{\mathbb{Q}}[\widetilde{Z}_n] = \sqrt{n}\cdot\mathbb{E}^{\mathbb{Q}}[Z_1] = \sqrt{n\,\Delta t}\cdot\frac{r - \sigma^2/2}{\sigma} + O(\sqrt{n}\,\Delta t) \xrightarrow[n \to \infty]{} \frac{\sqrt{T}\,(r - \sigma^2/2)}{\sigma},$$

$$\mathrm{Var}^{\mathbb{Q}}(\widetilde{Z}_n) = \mathrm{Var}^{\mathbb{Q}}(Z_1) \xrightarrow[n \to \infty]{} 1.$$

The renormalised sum $\widetilde{Z}_n$ has, in the limit, an explicit non-trivial mean and unit variance — precisely the moments of the Gaussian distribution we expect to recover via a central limit theorem.

---

## 4. Central limit theorem for a triangular array

We can now identify the limiting distribution of $\widetilde{Z}_n$. The technical subtlety is that the variables $(Z_k)$ are **independent but not identically distributed across $n$**: their common law under $\mathbb{Q}$ depends on $\Delta t = T/n$, hence on $n$. The classical Lindeberg-Lévy CLT (which assumes a fixed common distribution) does not apply.

The right framework is that of a **triangular array** $\{Z_k^{(n)}\}_{1 \leq k \leq n}$: for each $n$, we have a row of $n$ i.i.d. variables, but the rows themselves change with $n$. The appropriate result is the **Lindeberg-Feller central limit theorem** (see Billingsley, *Probability and Measure*, Theorem 27.2).

### Verifying the Lindeberg condition

Let $\mu_n := \mathbb{E}^{\mathbb{Q}}[Z_k^{(n)}]$ and $s_n^2 := \sum_{k=1}^n \mathrm{Var}^{\mathbb{Q}}(Z_k^{(n)}) = n \cdot \mathrm{Var}^{\mathbb{Q}}(Z_1^{(n)})$. The Lindeberg-Feller theorem requires that for every $\delta > 0$,

$$\frac{1}{s_n^2} \sum_{k=1}^n \mathbb{E}^{\mathbb{Q}}\left[ (Z_k^{(n)} - \mu_n)^2 \, \mathbf{1}_{|Z_k^{(n)} - \mu_n| > \delta s_n} \right] \xrightarrow[n \to \infty]{} 0.$$

In our case the variables are *uniformly bounded*: $|Z_k^{(n)} - \mu_n| \leq 2$ for all $n, k$. Meanwhile $s_n \sim \sqrt{n} \to \infty$, so for any fixed $\delta > 0$ the indicator $\mathbf{1}_{|Z_k^{(n)} - \mu_n| > \delta s_n}$ vanishes identically for $n$ large enough. The Lindeberg condition holds trivially.


### Conclusion of the CLT

Lindeberg-Feller then gives

$$\frac{1}{s_n}\sum_{k=1}^n \big(Z_k^{(n)} - \mu_n\big) \xrightarrow[n \to \infty]{(d)} \mathcal{N}(0, 1).$$

Since $s_n^2 = n \cdot \mathrm{Var}^{\mathbb{Q}}(Z_1^{(n)}) \to n$, dividing by $\sqrt{n}$ instead of $s_n$ changes nothing in the limit. Combining with the asymptotics of Section 3:

$$\widetilde{Z}_n \xrightarrow[n \to \infty]{(d)} \frac{\sqrt{T}\,(r - \sigma^2/2)}{\sigma} + \widetilde{Z}, \qquad \widetilde{Z} \sim \mathcal{N}(0, 1).$$

### Limiting distribution of the log-price


Recall that $\log(S_n / S_0) = \sigma\sqrt{T}\cdot \widetilde{Z}_n$. The linear transformation yields, under $\mathbb{Q}$:

$$\log\frac{S_n}{S_0} \xrightarrow[n \to \infty]{(d)} \mathcal{N}\left( \big(r - \tfrac{\sigma^2}{2}\big)\,T, \sigma^2 T \right).$$

In the continuous-time limit, the underlying is log-normal under $\mathbb{Q}$ with drift $(r - \sigma^2/2)T$ and variance $\sigma^2 T$ — the parameters of the Black-Scholes log-price distribution.


---

## 5. Convergence of the price: Put first, Call by parity

By the definition of convergence in distribution (the *Portmanteau theorem*), $\widetilde{Z}_n \xrightarrow{(d)} \widetilde{Z}$ implies $\mathbb{E}[f(\widetilde{Z}_n)] \to \mathbb{E}[f(\widetilde{Z})]$ **for every bounded continuous function** $f$.

### The Put: bounded payoff, direct convergence

The payoff of the European put, $(K - S_n)^+$, is bounded above by $K$. Writing $S_n = S_0 \exp(\sigma\sqrt{T} \, \widetilde{Z}_n)$, we have $(K - S_n)^+ = f(\widetilde{Z}_n)$ with $f(z) = (K - S_0 e^{\sigma \sqrt{T} z})^+$ continuous and bounded by $K$. Direct application:

$$P_0^{\text{CRR}}(n) = e^{-rT} \mathbb{E}^{\mathbb{Q}}\left[(K - S_n)^+\right] \xrightarrow[n \to \infty]{} e^{-rT} \mathbb{E}^{\mathbb{Q}}\left[(K - S_T^{\text{cont}})^+\right] =: P_0^{\text{BS}},$$

where $S_T^{\text{cont}} := S_0 \exp\big((r - \sigma^2/2)T + \sigma\sqrt{T}\, \widetilde{Z}\big)$ is the log-normal limit.

### The Call: unbounded payoff, parity workaround

The payoff $(S_n - K)^+$ is **not bounded** in $n$, so the Portmanteau theorem does not apply directly. Two routes:

1. **Strengthen the convergence** to convergence in $L^1$ via uniform integrability of $(S_n)$ under $\mathbb{Q}$, then apply Vitali's theorem.
2. **Use put-call parity**, valid at every $n$ by no-arbitrage:
   
$$C_0^{\text{CRR}}(n) - P_0^{\text{CRR}}(n) = S_0 - K e^{-rT}.$$

The right-hand side is independent of $n$. Passing to the limit on the left and using the Put convergence:

$$C_0^{\text{CRR}}(n) \xrightarrow[n \to \infty]{} P_0^{\text{BS}} + S_0 - K e^{-rT} =: C_0^{\text{BS}}.$$

The second route is more elegant and avoids any analytic argument.

---
## 6. The Black-Scholes formula

It remains to compute the limiting expectation explicitly. Set

$$Y := (r - \tfrac{\sigma^2}{2})T + \sigma\sqrt{T}\,\widetilde{Z}, \qquad \widetilde{Z} \sim \mathcal{N}(0, 1),$$

so that $S_T^{\text{cont}} = S_0 e^Y$ is log-normal with $Y \sim \mathcal{N}\big((r - \sigma^2/2)T,\, \sigma^2 T\big)$. The call price reads

$$C_0^{\text{BS}} = e^{-rT}\, \mathbb{E}^{\mathbb{Q}}\left[ (S_0 e^Y - K)^+ \right].$$

### Splitting on the exercise event

The payoff $(S_0 e^Y - K)^+$ vanishes on $\{S_0 e^Y \leq K\}$. We split the indicator $\mathbf{1}\_{S_0 e^Y > K}$, which is equivalent to $\mathbf{1}\_{Y > \log(K/S_0)}$, and write

$$\mathbb{E}^{\mathbb{Q}}\left[ (S_0 e^Y - K)^+ \right] = \underbrace{\mathbb{E}^{\mathbb{Q}}\left[ S_0 e^Y \,\mathbf{1}_{Y > \log(K/S_0)} \right]}_{=: A} - \underbrace{K \cdot \mathbb{Q}\left( Y > \log(K/S_0) \right)}_{=: B}.$$

The two terms $A$ and $B$ are now Gaussian integrals over a half-line, both computable in closed form.

### Standardisation

Let $\mu := (r - \sigma^2/2)T$ and $s := \sigma\sqrt{T}$, so $Y = \mu + s\widetilde{Z}$. The event $\{Y > \log(K/S_0)\}$ rewrites as $\{\widetilde{Z} > -d_2\}$ where

$$d_2 := \frac{\mu - \log(K/S_0)}{s} = \frac{\log(S_0/K) + (r - \sigma^2/2)T}{\sigma\sqrt{T}}.$$

By symmetry of the standard normal distribution, $\mathbb{Q}(\widetilde{Z} > -d_2) = \Phi(d_2)$, so

$$B = K \, \Phi(d_2).$$

### Computing $A$ by completing the square

We rewrite the integral defining $A$:

$$A = \int_{-d_2}^{+\infty} S_0 e^{\mu + s z} \cdot \frac{1}{\sqrt{2\pi}} e^{-z^2/2} \, dz = \frac{S_0\,e^{\mu}}{\sqrt{2\pi}} \int_{-d_2}^{+\infty} e^{s z - z^2/2}\, dz.$$

Completing the square: $sz - z^2/2 = s^2/2 - (z - s)^2/2$. Substituting and using $\mu + s^2/2 = rT$:

$$A = S_0\,e^{\mu + s^2/2} \cdot \frac{1}{\sqrt{2\pi}} \int_{-d_2}^{+\infty} e^{-(z-s)^2/2}\, dz = S_0\,e^{rT} \cdot \mathbb{Q}(\widetilde{Z}' > -d_2 - s),$$

where $\widetilde{Z}' := \widetilde{Z} - s \sim \mathcal{N}(0, 1)$. Setting $d_1 := d_2 + s = d_2 + \sigma\sqrt{T}$:

$$A = S_0\,e^{rT}\,\Phi(d_1).$$

### Final formula

Combining $A - B$, multiplying by the discount factor $e^{-rT}$, and using $d_1 = d_2 + \sigma\sqrt{T}$:

$$\boxed{ C_0^{\text{BS}} = S_0\,\Phi(d_1) - K e^{-rT}\,\Phi(d_2), \qquad d_{1,2} = \frac{\log(S_0/K) + \big(r \pm \tfrac{\sigma^2}{2}\big)T}{\sigma\sqrt{T}} }$$

The put follows by put-call parity:

$$P_0^{\text{BS}} = C_0^{\text{BS}} - S_0 + K e^{-rT} = K e^{-rT}\,\Phi(-d_2) - S_0\,\Phi(-d_1).$$

---

## 7. Numerical verification

The convergence $C_0^{\text{CRR}}(n) \to C_0^{\text{BS}}$ is implemented in [`notebooks/01_crr_to_bs_convergence.ipynb`](../notebooks/01_crr_to_bs_convergence.ipynb). The standard CRR parametrisation produces oscillatory convergence of order $O(1/n)$; the Leisen-Reimer parametrisation removes the oscillation and accelerates convergence to $O(1/n^2)$.

*(Figure to be added: `results/figures/crr_convergence.png`)*

---

## Appendix A: The Itô derivation

For reference, the standard route via continuous-time stochastic calculus.

Assume directly that under $\mathbb{Q}$, the underlying follows a geometric Brownian motion:

$$dS_t = r S_t \, dt + \sigma S_t \, dW_t^{\mathbb{Q}}.$$

By Itô's lemma applied to $f(s, t) = \log s$:

$$d(\log S_t) = \left(r - \tfrac{\sigma^2}{2}\right) dt + \sigma \, dW_t^{\mathbb{Q}},$$

so $\log(S_T / S_0) \sim \mathcal{N}\big((r - \sigma^2/2) T, \, \sigma^2 T\big)$ under $\mathbb{Q}$, **recovering exactly the distribution obtained as the discrete limit in Section 4**. The pricing formula follows identically.

Alternatively, $C(s, t) \in C^{2,1}$ satisfying the **Black-Scholes PDE**

$$\frac{\partial C}{\partial t} + r s \frac{\partial C}{\partial s} + \tfrac{1}{2} \sigma^2 s^2 \frac{\partial^2 C}{\partial s^2} - r C = 0$$

with terminal condition $C(s, T) = (s - K)^+$ admits the same closed form, via the Feynman-Kac representation. All three routes — discrete limit, risk-neutral expectation, PDE — coincide.

---

## References

- Cox, J. C., Ross, S. A., Rubinstein, M. (1979). *Option pricing: a simplified approach*. JFE 7, 229-263.
- Billingsley, P. (1995). *Probability and Measure*, 3rd ed. Wiley. (Theorem 27.2 for the Lindeberg-Feller CLT.)
- Shreve, S. (2004). *Stochastic Calculus for Finance, vols. I & II*. Springer.
- Hull, J. (2018). *Options, Futures, and Other Derivatives*, 10th ed., ch. 13-15.
- Joshi, M. S. (2008). *The Concepts and Practice of Mathematical Finance*, 2nd ed.
- Delbaen, F., Schachermayer, W. (1994). *A general version of the Fundamental Theorem of Asset Pricing*. Math. Ann. 300, 463-520.
