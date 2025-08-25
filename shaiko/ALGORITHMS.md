# Boost.Math::cubature — Algorithms

This document specifies the mathematics and core algorithms used in the proposed module. Citations point to standard references and software with compatible ideas/terminology.

---

## 1) Deterministic Adaptive Cubature (hypercubes)

### 1.1 Rules: Fully Symmetric Embedded Pairs (Genz–Malik / DCUHRE)

We employ **fully symmetric** cubature rules on the hypercube with **embedded pairs** for **a posteriori error estimation**. The canonical choice is degree‑7 with embedded degree‑5, and degree‑9 with embedded degree‑7; higher pairs (11/9, 13/11) are available in CUBA/Cuhre and the literature. The embedding means all degree‑k nodes are a subset of degree‑(k+2) nodes, enabling a cheap difference estimate. citeturn2view0turn12search4turn12search5

**Notes.** Construction of good fully symmetric rules and their properties are summarized by Espelid, Sørevik, Cools, Genz and others. We will include ready‑to‑use tables for d up to 15 (practical) and a generator for general d. citeturn12search0turn12search2turn12search8

### 1.2 Global Adaptive Subdivision (DCUHRE‑style)

We follow **DCUHRE**/**Cuhre**: maintain a max‑heap of subregions ordered by estimated error; repeatedly split the worst region along the axis maximizing an estimate of local variation; continue until the **global error estimate ≤ max(abs_tol, rel_tol·|I|)** or `max_eval` is reached. Error per region is the embedded‑rule difference, optionally safeguarded by a **spread** estimate. citeturn2view0turn3view0

**Axis selection.** Use component‑wise differences between symmetric node pairs to estimate directional variation; bisect the dimension with the largest weighted variation (as in Cuhre/DCUHRE). citeturn3view0

**Vector integrands.** Evaluate all components at once; region error is the norm (L∞ by default) of component‑wise embedded differences; the splitting criterion uses the component with the largest contribution.

**Pseudocode (simplified):**
```text
PQ ← { whole box R0 with estimate e0 from embedded pair }
I  ← 0; E ← 0
while E > tol and evals < max_eval:
  R ← pop_max(PQ)
  I ← I - I_R + I_R_refined           // keep global sum consistent
  E ← E - e_R + e_R_refined
  if e_R_refined > local_tol:         // local_tol ~ α·(tol * vol(R)/vol(total))
     (R1,R2) ← split(R, argmax_directional_variation(R))
     for Ri in {R1,R2}: push(PQ, estimate(Ri))
return (I, E)
```

### 1.3 Numerical details

- **Scaling to [0,1]^d** with Jacobian; affine invariance of rules respected.
- **Cancellation**: accumulate with Kahan or two‑sum for multiprecision friendliness.
- **Safeguards**: if embedded difference is spuriously small (near machine epsilon), fall back to **spread estimate** (max−min of partial sums) to avoid stopping too early, as in **Cuhre**. citeturn3view0

---

## 2) Sparse‑Grid Cubature (Smolyak)

### 2.1 Smolyak construction

Given nested 1‑D rules \( U^i \) of increasing accuracy (Clenshaw–Curtis or Gauss–Patterson), the **Smolyak operator** of level \( \ell \) in dimension \( d \) is
\[
A(\ell, d) = \sum_{\ell-d+1 \le |\boldsymbol{i}| \le \ell} (-1)^{\ell-|\boldsymbol{i}|} \binom{d-1}{\ell-|\boldsymbol{i}|}
\bigotimes_{j=1}^d U^{i_j},
\]
which yields a **sparse** combination of tensor rules with **nested abscissas**, drastically reducing points vs. full tensor products. citeturn1news0

### 2.2 1‑D bases & nestedness

- **Clenshaw–Curtis (CC)** rules are **nested** under 1→3→5→9→… growth (or other slow/exponential growth choices). We implement FFT/DCT‑based weight construction and cache weights; choose a slow growth rule suitable for sparse grids. citeturn10search0turn10search5
- **Gauss–Patterson (GP)** extend Gauss–Legendre with embedded nodes; usable but not perfectly nested across all orders—CC is the default for maximal reuse. citeturn10search7

### 2.3 Adaptivity

- **Level adaptivity**: increase \(\ell\) until a user budget or error heuristic is hit.
- **Dimension adaptivity** (optional): estimate **surpluses** per anisotropic index and refine where contribution is largest, following Gerstner–Griebel ideas. citeturn1news2

### 2.4 Error estimates

- Use **hierarchical surplus norms** and/or **difference between levels**; practical stop criteria mirror Tasmanian/SparseGrids literature (not a hard a priori bound).

---

## 3) (Randomized) Quasi‑Monte‑Carlo (QMC/RQMC)

### 3.1 Sobol’ points

We use `boost::random::sobol_engine` to generate base‑2 **Sobol’** low‑discrepancy points in \( [0,1)^d \). QMC integrates \(f\) by \( \hat{I}_N = \frac{1}{N}\sum f(u_i)\), often with superior error compared to IID MC. citeturn14view0

### 3.2 Owen scrambling (RQMC)

To obtain **unbiased** estimates and a **sample variance**, perform **Owen digital scrambling** of Sobol’ points and average over `k` replicates. The variance decays faster than IID MC for smooth integrands. We implement base‑2 per‑digit XOR scrambles (nested uniform scramble). citeturn1news2

**Error bar.** With replicates \(Y_r\), report \( \bar{Y} \pm t_{0.975,k-1} s/\sqrt{k} \).

### 3.3 Periodization / variance reduction

- **Tent (Baker) transform** \( \phi(u)=1-2|u-\tfrac12| \) improves QMC for non‑periodic integrands; enabled by default. citeturn2search3
- Optional: **shifts**, **digital shifts**, **Brownian bridge** (future).

### 3.4 Parallelism & determinism

We assign **disjoint index ranges** to threads to preserve order; with scrambling, each replicate uses an independent scramble seed. Reproducibility is achieved by fixed seeds and deterministic work partitioning.

---

## 4) Simplex Integration

Two options:
1. **Duffy transform** maps the unit cube to a simplex; good for endpoint/edge singularities. The Jacobian is handled analytically; robust and easy to implement. citeturn4search1
2. **Direct embedded rules + subdivision** as in **CUBPACK** (bisect longest edge, choose edge via smoothness indicators; use embedded pairs for error). We can start with Duffy + adaptive hypercube after mapping, then add native simplex rules in phase‑2. citeturn13search1

---

## 5) Domain & Improper Integral Transforms

- Infinite intervals: standard bijective maps \( (-\infty,\infty) \leftrightarrow (0,1) \) (e.g., `x = \tan(\pi(u-1/2))`, `x = u/(1-u)`), with weight/Jacobian; documented and selectable per algorithm. citeturn4search2
- Endpoint singularities: **tanh‑sinh** is 1‑D; for ND, we combine mild algebraic endpoint transforms per axis before cubature.
- Coordinate **shifts/rotations** are allowed via user functors for aligning features.

---

## 6) Vector‑Valued Integrands

Given \( f: \mathbb{R}^d \to \mathbb{R}^m \), evaluate at nodes once and accumulate **m** components with identical weights. Error estimates use L∞ (default) or user‑selected norms over components. The API provides overloads for `std::array<Real,M>`, `std::vector<Real>`, and a callback `(x, out*, m)` style (zero‑alloc).

---

## 7) Numerical Stability & Precision

- All sums use pairwise accumulation; Kahan/two‑sum toggled for large **m** or high precision.
- Works with `double`, `long double`, and `boost::multiprecision::cpp_dec_float_50/100`. citeturn6search1

---

## 8) Benchmarks & Canonical Tests

- **Genz test functions** across \(d=2..12\) with known/closed‑form references and parameterized difficulty. citeturn5search0
- **Gaussian/Keister‑type** integrals for sanity checks. citeturn5search1turn5search9
- Comparisons vs. **CUBA (Cuhre/Suave/Vegas)** and **cubature** (Steven G. Johnson) on vector/scalar integrands to ensure competitive performance/accuracy. citeturn3view0turn11search11
