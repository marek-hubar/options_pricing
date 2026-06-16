### The calculator is live [here](https://marek-hubar.github.io/options_pricing/).

# Option Pricing Algorithms

This repository contains my C++ implementations of pricing European and (more importantly) American options.
All options are modelled using the Black-Scholes equation and the underlying price is assumed to be a geometric Brownian motion.

In particular for a the price $V = V(t, s)$ of an option the Black-Scholes formula is
```math
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial s^2} + (r - q) S \frac{\partial V}{\partial s} - r V = 0,
```
and for the underlying stock price $S = S(t)$ we assume:

```math
d S_t = (r - q) S_t d t + \sigma S_t d W_t 
```

where $r$ is the risk-free rate, $\sigma$ is the volatility, and $q$ is the continuous dividend yield.


## European Options
These are priced using an explicitly solved Black-Scholes equation.
The solution relies on transforming the Black-Scholes equation into the standard heat equation, solving this heat equation explicitly for the correspondingly transformed payoff kernel, and then transforming the solution back into the original coordinates.


## American Options
These cannot be priced using an explicit formula, so we have to calculate the solution numerically.
Our solution still relies on a suitable coordinate transformation to convert the Black-Scholes equation into a standard heat equation, which we then solve.
Namely, we invert and scale the time, and take the logarithm of "centered" space:
```math
\tau = \frac{1}{2} \sigma^2 (T - t),
```
```math
x = \log(\frac{s}{K}).
```


### Transformation to the Heat Equation

Assume the price transforms as $V(t, s) = K \cdot v(\tau, x)$ and $v(\tau, x) = e^{\alpha x + \beta\tau}\, u(\tau, x)$. 

First we change variables. Let $t^* = T - t$ be time to maturity, so $\frac{\partial}{\partial t} = -\frac{\partial}{\partial t^*}$. With $S = K e^x$ and $V(t^*, S) = K \cdot v(t^*, x)$, the partial derivatives become:

```math
S \frac{\partial V}{\partial S} = K \frac{\partial v}{\partial x}, \qquad
S^2 \frac{\partial^2 V}{\partial S^2} = K\left( \frac{\partial^2 v}{\partial x^2} - \frac{\partial v}{\partial x} \right)
```

Substituting into the Black-Scholes PDE and dividing by $K$:

```math
-\frac{\partial v}{\partial t^*} + \frac{1}{2}\sigma^2 \frac{\partial^2 v}{\partial x^2} + \left( r - q - \frac{1}{2}\sigma^2 \right) \frac{\partial v}{\partial x} - r v = 0
```

Now scale time by letting $\tau = \frac{1}{2}\sigma^2\, t^*$, so $\frac{\partial}{\partial t^*} = \frac{1}{2}\sigma^2 \frac{\partial}{\partial\tau}$. Dividing through by $\frac{1}{2}\sigma^2$ and defining the dimensionless constants:

```math
k_1 = \frac{2(r - q)}{\sigma^2}, \qquad k_2 = \frac{2r}{\sigma^2}
```

we obtain:

```math
-\frac{\partial v}{\partial \tau} + \frac{\partial^2 v}{\partial x^2} + (k_1 - 1)\frac{\partial v}{\partial x} - k_2\, v = 0
```

Now substitute $v = e^{\alpha x + \beta\tau}\, u$:

```math
\begin{aligned}
\frac{\partial v}{\partial \tau} &= e^{\alpha x + \beta\tau}(\beta u + \frac{\partial u}{\partial \tau}), \\[4pt]
\frac{\partial v}{\partial x} &= e^{\alpha x + \beta\tau}(\alpha u + \frac{\partial u}{\partial x}), \\[4pt]
\frac{\partial^2 v}{\partial x^2} &= e^{\alpha x + \beta\tau}(\alpha^2 u + 2\alpha \frac{\partial u}{\partial x} + \frac{\partial^2 u}{\partial x^2}).
\end{aligned}
```

After substituting and dividing by $e^{\alpha x + \beta\tau}$, the terms collect by derivative order:

```math
-\beta u - \frac{\partial u}{\partial \tau} + \alpha^2 u + 2\alpha \frac{\partial u}{\partial x} + \frac{\partial^2 u}{\partial x^2} + (k_1 - 1)\alpha u + (k_1 - 1)\frac{\partial u}{\partial x} - k_2\, u = 0
```

To recover the standard heat equation, we eliminate the $\frac{\partial u}{\partial x}$ term and the $u$ term by choosing $\alpha$ and $\beta$ appropriately:

```math
\begin{aligned}
\frac{\partial u}{\partial x} \text{ term:}&\quad 2\alpha + (k_1 - 1) = 0 \;\Longrightarrow\; \boxed{\alpha = -\frac{1}{2}(k_1 - 1)}, \\[8pt]
u \text{ term:}&\quad -\beta + \alpha^2 + (k_1 - 1)\alpha - k_2 = 0 \;\Longrightarrow\; \boxed{\beta = -\alpha^2 - k_2}.
\end{aligned}
```

With these choices, the $\frac{\partial u}{\partial x}$ and $u$ terms vanish, leaving the standard heat equation:

```math
\boxed{\frac{\partial u}{\partial \tau} = \frac{\partial^2 u}{\partial x^2}}
```

**Special case — no dividends ($q = 0$).** Then $k_1 = k_2 = k = \frac{2r}{\sigma^2}$, and:

```math
\alpha = -\frac{1}{2}(k - 1), \qquad
\beta = -\frac{1}{4}(k-1)^2 - k = -\frac{1}{4}(k+1)^2.
```

**Inverse transformation.** The option price is recovered by:

```math
V(T - t, S) = K \cdot v(\tau, x) = K \cdot e^{\alpha x + \beta\tau} \cdot u(\tau, x)
```

**Initial condition ($\tau = 0$).** The payoff at expiry gives $u(0, x)$:

```math
u(0, x) = e^{-\alpha x} \cdot \frac{\max(K - S, 0)}{K} = e^{-\alpha x} \cdot \max(1 - e^x, 0) \quad \text{(put)},
```

```math
u(0, x) = e^{-\alpha x} \cdot \frac{\max(S - K, 0)}{K} = e^{-\alpha x} \cdot \max(e^x - 1, 0) \quad \text{(call)}.
```


### Numerical methods
Finite difference method is used to calculate the values of $u$ on a discretized grid.
For the purposes of this project, we do not consider adaptive grid methods.
We focus on the specifics of time and space discretization in further sections.
For simpler derivations we assume a uniformity of the time grid.

We start with the transformed payoff $u(0, \cdot)$ at time $\tau = 0$ and evaluate this function at the spatial grid points.
For iteratively finding the values of $u$ on the grid points in the next time step, we employ the **Crank-Nicolson method**.
This is an implicit-explicit method, which has the advantage of being $O(N)$ per step ($N$ is the number of grid points), converging at $O(\Delta \tau^2)$, and being unconditionally stable.
It works by approximating the partial derivatives by finite differences.
Specifically, we can approximate
```math
\frac{\partial u(\tau_i, x_j)}{\partial \tau} \approx \frac{u(\tau_{i+1}, x_j) - u(\tau_i, x_j)}{\Delta\tau} \approx \frac{\partial u(\tau_{i+1}, x_j)}{\partial \tau}
```
and, denoting $h^L_j = x_j - x_{j-1}$, $h^R_j = x_{j+1} - x_j$, and
```math
\begin{aligned}
a_j &= \frac{2}{h^L_j (h^L_j + h^R_j)}$,\\
b_j &= \frac{-2}{h^L_j h^R_j}$,\\
c_j &= \frac{2}{h^R_j (h^L_j + h^R_j)}$,\\
\end{aligned}
```
```math
\begin{aligned}
\frac{\partial^2 u(\tau_i, x_j)}{\partial x^2} &\approx a_j u(\tau_i, x_{j-1}) + b_j u(\tau_i, x_{j}) + c_j u(\tau_i, x_{j+1}), \\
\frac{\partial^2 u(\tau_{i+1}, x_j)}{\partial x^2} &\approx a_j u(\tau_{i+1}, x_{j-1}) + b_j u(\tau_{i+1}, x_{j}) + c_j u(\tau_{i+1}, x_{j+1}).
\end{aligned}
```
An explicit method would equate the finite differences at time $i$, while an implicit method would equate the finite differences at time $i+1$.
This implicit-explicit method combines the two approaches to get the advantages of both:
```math
\frac{u^{i+1}_j - u^i_j}{\Delta\tau} = \frac{1}{2}\left(a_j u^{i+1}_{j-1} + b_j u^{i+1}_{j} + c_j u^{i+1}_{j+1}\right) + \frac{1}{2}\left(a_j u^i_{j-1} + b_j u^i_{j} + c_j u^i_{j+1}\right),
```
where we denote for all $k$ and $l$ the value of $u$ at time $\tau_k$ and space $x_l$ by $u^k_l$.
Denoting $\gamma = \frac{\Delta\tau}{2}$, this equation can be rewritten as:
```math
- \gamma a_j u^{i+1}_{j-1} + (1 - \gamma b_j) u^{i+1}_{j} - \gamma c_j u^{i+1}_{j+1} = \gamma a_j u^i_{j-1} + (1 + \gamma b_j) u^i_{j} + \gamma c_j u^i_{j+1}.
```
At time $i$ the RHS is known and we need to solve for the $u^{i+1} \in \mathbb{R}^N$.
To do this, we utilize the fact that on the LHS we multiply vector $u^{i+1}$ by a tridiagonal matrix.
Tridiagonal systems can be solved efficiently using the [Thomas matrix algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) in $O(N)$ time per step.


### Boundary Conditions

The Black-Scholes PDE is posed on the infinite domain $x \in (-\infty, \infty)$, but we must truncate to a finite computational window $[x_0, x_{N-1}]$. Truncation introduces artificial boundaries where we must supply conditions for the PDE.

We impose $\frac{\partial^2 V}{\partial S^2} = 0$ at the extreme stock prices. This is equivalent to assuming the option gamma vanishes far from the strike: for deep in-the-money options the price curve becomes locally linear (slope ±1), and for deep out-of-the-money options the price is essentially zero. In both cases $\frac{\partial^2 V}{\partial S^2} \approx 0$.

In the log-space coordinates, $\frac{\partial^2 V}{\partial S^2} = 0$ is not exactly $\frac{\partial^2 u}{\partial x^2} = 0$, but the cumulative error introduced by using a linearity condition on $u$ at a distance of $3\sigma\sqrt{T}$ from the center is negligible for practical option prices.

Discretely, the linearity condition means the boundary value must satisfy a linear extrapolation from its two interior neighbours. On a non-uniform grid the extrapolation weights depend on the local mesh ratios:

```math
c_L = \frac{x_1 - x_0}{x_2 - x_1}, \qquad
c_R = \frac{x_{N-1} - x_{N-2}}{x_{N-2} - x_{N-3}}
```

```math
u_0 = (1 + c_L) u_1 - c_L\, u_2, \qquad
u_{N-1} = (1 + c_R) u_{N-2} - c_R\, u_{N-3}
```

These relations are substituted directly into the tridiagonal system. For the first interior point ($j=1$), substituting $u_0$ into the LHS of the Crank-Nicolson scheme modifies the coefficients:

```math
b_1 \leftarrow b_1 + a_1 (1 + c_L), \qquad c_1 \leftarrow c_1 - a_1\, c_L, \qquad a_1 \leftarrow 0
```

An analogous modification is made for the last interior point ($j = N-2$). This keeps the system tridiagonal while enforcing the linearity condition at the boundaries.


### Early Exercise Constraint

European options can only be exercised at expiry $T$, whereas American options can be exercised at any $t \leq T$. The holder will exercise early whenever the continuation value (the price of holding the option) drops below the immediate exercise payoff. The price of an American option is therefore the *optimal stopping* value:

```math
V(t, S) = \sup_{\text{stopping times } \zeta \in [t,T]} \mathbb{E}^\mathbb{Q}\!\left[ e^{-r(\zeta - t)} \, \text{payoff}(S_\zeta) \;\middle|\; S_t = S \right]
```

In the discretized $u$-space, the early exercise condition is enforced after each Crank-Nicolson time step. The intrinsic payoff at expiry ($\tau = 0$) is pre-computed as:

```math
u_{\text{intrinsic}}(0, x_j) = e^{-\alpha x_j} \cdot \max(1 - e^{x_j}, 0) \quad \text{or} \quad e^{-\alpha x_j} \cdot \max(e^{x_j} - 1, 0)
```

At a later time $\tau > 0$, the intrinsic value gains an $e^{-\beta\tau}$ factor:

```math
u_{\text{intrinsic}}(\tau, x_j) = e^{-\beta\tau} \cdot u_{\text{intrinsic}}(0, x_j)
```

After advancing the PDE from $\tau_i$ to $\tau_{i+1}$, we overwrite every grid point with:

```math
u(\tau_{i+1}, x_j) \leftarrow \max\!\big( u(\tau_{i+1}, x_j),\; e^{-\beta \tau_{i+1}} \cdot u_{\text{intrinsic}}(0, x_j) \big)
```

This ensures the option price never falls below the immediate exercise value — the defining property of the American option.


### Spatial Grid Design

The spatial grid is constructed in log-price space $x = \log(S/K)$ rather than in $S$-space. This has three advantages:

1. **Well-conditioned coefficients.** The Black-Scholes coefficients in $x$ are constant after the transformation — the $a_j, b_j, c_j$ on a uniform $x$-grid depend only on the constant $\Delta x$ and $\Delta\tau$, making the Crank-Nicolson matrix stable.

2. **Natural centering.** For at-the-money options ($S = K$), the center of the grid $x = 0$ falls exactly at the current spot price, minimizing the interpolation error at the final pricing step.

3. **Symmetric coverage.** Stock prices are approximately log-normal, so a uniform grid in $x$ naturally covers the relevant price range with uniform resolution in standard deviations.

**Grid bounds: $\pm 3\sigma\sqrt{T}$.** The log-return over the option's lifetime follows a normal distribution with standard deviation $\sigma\sqrt{T}$. At $\pm 3\sigma\sqrt{T}$, the cumulative probability mass exceeds 99.7%, so the truncation error for practical option values is negligible. Narrower bounds would cut off the tails; wider bounds would waste nodes on prices far from the strike.

```math
x_{\min} = \log\!\left(\frac{S_0}{K}\right) - 3\sigma\sqrt{T}, \qquad
x_{\max} = \log\!\left(\frac{S_0}{K}\right) + 3\sigma\sqrt{T}
```

The grid is uniform: $x_j = x_{\min} + j \cdot \Delta x$ for $j = 0, \dots, N-1$, where $\Delta x = (x_{\max} - x_{\min}) / (N-1)$.


### Dividend Yield Extension

The continuous dividend yield $q$ enters the Black-Scholes PDE by replacing the drift $r$ with $r - q$ in the convection term:

```math
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + (r - q) S \frac{\partial V}{\partial S} - r V = 0
```

This propagates into the heat equation coefficients. With dividends, the dimensionless constants become:

```math
k_1 = \frac{2(r - q)}{\sigma^2}, \qquad k_2 = \frac{2r}{\sigma^2}
```

$\alpha$ depends on $k_1$ (drift with dividends), while $\beta$ depends on both $k_1$ and $k_2$. When $q = 0$, these reduce to the standard no-dividend case $k_1 = k_2 = k$. The discounting term always uses $r$ alone, since the risk-free rate governs the time value of money regardless of dividends.


### Implementation Notes

The repository contains two C++ pricers:

- **`EuropeanOptionPricer`** — Closed-form Black-Scholes solution. Computes both single prices and price surfaces over $(S, t)$ grids by direct evaluation of the Black-Scholes formula at each point.

- **`AmericanOptionPricer`** — Finite difference solver using Crank-Nicolson in the transformed $(x, \tau)$ coordinates. Both pricers are stateless: all grid buffers are passed explicitly from the caller.

The C++ code is compiled to WebAssembly using Emscripten and served as a static website on GitHub Pages. The web app overlays a Plotly 3D surface plot for visualization, with configurable spatial resolution, temporal resolution, option type, and market parameters (spot price, strike, volatility, risk-free rate, and dividend yield).

The `price_surface` in both classes populates three output buffers: stock prices $S$, times to maturity $t$, and a 2D array of option prices $V(S, t)$.

The Thomas solver implementation pre-computes the LU-like factorization of the tridiagonal matrix in its constructor (the `denom` and `c_prime` arrays). Each time step then requires only an $O(N)$ forward sweep and back substitution. This avoids recomputing the matrix inverse or factorization at every step.
