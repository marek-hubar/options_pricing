# Option Pricing Algorithms

This repository contains my C++ implementations of pricing European and (more importantly) American options.
All options are modelled using the Black-Scholes equation and the underlying price is assumed to have be a geometric Brownian motion.

In particular for a the price $V = V(t, s)$ of an option the Black-Scholes formula is
```math
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial s^2} + r S \frac{\partial V}{\partial s} - r V = 0,
```
and for the underlying stock price $S = S(t)$ we assume:

```math
d S_t = r S_t d t + \sigma S_t d W_t 
```

where $r$ is the risk-free rate, and $\sigma$ is the volatility.


## European Options
These are priced using an explicitly solved Black-Scholes equation.
The solution relies on transforming the Black-Scholes equation into the standard heat equation, solving this heat equation explicitly for the correspondingly transformed payoff kernel, and then transforming the solution back into the original coordinates.

## American Options
These cannot be priced using an explicit formula, so we have to calculate the solution numerically.
Our solution still relies on a suitable coordinate transformation to convert the black-scholes equation into a standard heat equation, which we then solve.
Namely, we invert and scale the time, and take the logarithm of "centered" space:
```math
\tau = \frac{1}{2} \sigma^2 (T - t),
```
```math
x = \log(\frac{s}{K}).
```
Further, we assume that the price $V$ can be expressed as 
```math
V(t, s) = e^{-\alpha x - \beta \tau} u(\tau, x),
```
where $\alpha$ and $\beta$ are some constants, which we will take a look at later.
Using these transformations we find the Black-Scholes equation for $V$ is equivalent to the standard heat equation
```math
\frac{\partial u}{\partial \tau} = \frac{\partial^2 u}{\partial x^2}
```
for $u$.

### Numerical methods
Finite difference method is used to calculate the values of $u$ on a discretized grid.
For the purposes of this project, we do not consider adaptive grid methods.
We focus on the specifics of time and space discretization in further sections.
For simpler derivations we assume a uniformity of the time grid.

We start with the transformed payoff $u(0, \cdot)$ at time $\tau = 0$ and evaluate this function at the spatial grid points.
For iteratively finding the values of $u$ on the grid points in the next time step, we employ the **Crank-Nicolson method**.
This is an implicit-explicit method, which has the advantage of being $O(N)$ per step ($N$ is the number of grid poitns), converging at $O(\Delta \tau)$, and being unconditionally stable.
It works by approximating the parital derivatives by finite differences.
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
where we denote for all $k$ and $l$ the value of $u$ at time $tau_k$ and space $x_l$ by $u^k_l$.
Denoting $\gamma = \frac{\Delta\tau}{2}$, this equation can be rewritten as:
```math
- \gamma a_j u^{i+1}_{j-1} + (1 - \gamma b_j) u^{i+1}_{j} - \gamma c_j u^{i+1}_{j+1} = \gamma a_j u^i_{j-1} + (1 + \gamma b_j) u^i_{j} + \gamma c_j u^i_{j+1}.
```
At time $i$ the RHS is known and we need to solve for the $u^{i+1} \in \mathbb{R}^N$.
To do this, we utilize the fact that the on the LHS we basically multiply vector $u^{i+1}$ by a tridiagonal matrix.
Tridiagonal systems can be solved efficiently using the [Thomas matrix algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) in $O(N)$ time.
