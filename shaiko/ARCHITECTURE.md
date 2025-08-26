# Boost.Math::cubature — Architecture

## 0) Scope & Definitions of Done (DoD)

**Scope.** Header‑only, policy‑driven **N‑dimensional integration** over hyperrectangles and simplices:
- Deterministic **adaptive cubature** (Genz–Malik/DCUHRE‑style) with embedded error estimates.
- **Sparse‑grid (Smolyak)** cubature with nested 1‑D rules (Clenshaw–Curtis and Gauss–Patterson) and optional anisotropic levels.
- **Quasi‑Monte‑Carlo (QMC)** with Sobol’ points from Boost.Random and optional **Owen scrambling** (RQMC) for variance estimates; periodization transforms (tent/Baker).
- Domain transforms for **improper** limits and weak endpoint singularities; optional **Duffy** transform for simplices/edge singularities.
- Vector‑valued integrands (compute **all components with shared abscissas**).
- Concurrency via `boost::asio::thread_pool`; integrates with `boost::math::policies` and `Boost.Multiprecision`.
- Clean error policies, deterministic reproducibility, workspace reuse, cancellation via stop token.

**DoD (phase‑1):**
1. Hyperrectangles: adaptive cubature (Genz–Malik degree‑7/5, degree‑9/7 embedded), sparse grids (levels 1–7), QMC/RQMC (Sobol’, k‑replicates) with transforms.
2. Simplices: Duffy mapping + DCUHRE‑style subdivision with embedded rules.
3. Vector‑valued integrands and **n ≤ 32** dimensions supported (soft cap; higher possible with QMC).
4. Multiprecision (`cpp_dec_float_50`) support; error policies honored.
5. Deterministic parallel execution and cancellation; benchmarks + acceptance tests pass (see `TESTING_VALIDATION.md`).

Motivation and fit with existing Boost.Math quadrature & policies are documented in Boost docs and in literature. citeturn7view0turn8view0


> Docs sources
>
> - [Algorithms](ALGORITHMS.md)
> - [Architecture](ARCHITECTURE.md)
> - [API](API.md)
> - [Testing & Validation](TESTING_VALIDATION.md)

---

## 1) Package Layout (headers only)

```
boost/math/cubature/
  adaptive.hpp           // Genz–Malik/DCUHRE adaptive hypercube cubature
  sparse_grid.hpp        // Smolyak + nested 1-D rules (CC, GP)
  qmc.hpp                // Sobol’/RQMC driver + transforms
  simplex.hpp            // Simplex domains + Duffy + embedded rules
  transforms.hpp         // finite↔infinite, weak singular endpoints, Duffy, etc.
  regions.hpp            // region types: hypercube, simplex; subdivision utils
  concepts.hpp           // Integrand/Domain concepts (C++20 constraints)
  policies.hpp           // aliases/adapters to boost::math::policies
  workspace.hpp          // scratch memory, stacks, RNG/QMC state
  detail/
    gm_rules.hpp         // embedded Genz–Malik point/weight tables
    cc_rules.hpp         // Clenshaw–Curtis weights (cached / constexpr)
    gp_rules.hpp         // Gauss–Patterson weights
    sobol_owen.hpp       // Owen scrambling utilities for Sobol’ points
    adaptivity.hpp       // priority queues, error bookkeeping
    vectorize.hpp        // pack/unpack vector-valued integrands
```

**Header‑only.** Boost.Math is predominantly header‑only; we follow that pattern, with optional pre‑generated tables guarded by macros. citeturn16search1turn16search8

---

## 2) Public API (high level sketch)

```cpp
namespace boost::math::cubature {

// Domain types
template <class Real> struct hypercube { std::vector<Real> lower, upper; };
template <class Real> struct simplex  { std::vector<std::vector<Real>> vertices; /* d+1 vertices, each of length d */ };

// Vector-valued integrand support via traits:
template <class F, class Real>
concept ScalarIntegrand = /* f(std::span<const Real>) -> Real */;

template <class F, class Real>
concept VectorIntegrand = /* f(x, out*, ncomp) or returns std::array<Real,M> */;

// Policies
using default_policy = boost::math::policies::policy<>; // forwarded

// Adaptive cubature (hypercube)
template <class Real, class F, class Policy = default_policy>
auto integrate_adaptive(const F& f, const hypercube<Real>& box,
                        Real abs_tol, Real rel_tol,
                        std::size_t max_eval = 0,
                        Policy const& pol = Policy{}) -> result<Real>;

// Sparse-grid Smolyak (hypercube)
template <class Real, class F, class Policy = default_policy>
auto integrate_sparse_grid(const F& f, const hypercube<Real>& box,
                           unsigned level,
                           Policy const& pol = Policy{}) -> result<Real>;

// QMC / RQMC (hypercube); k replicates for variance/error estimate
template <class Real, class F, class Policy = default_policy, class Sobol = boost::random::sobol_engine<>>
auto integrate_qmc(const F& f, const hypercube<Real>& box,
                   std::size_t n_points_per_rep, unsigned replicates = 1,
                   bool tent_transform = true,
                   Policy const& pol = Policy{}) -> result<Real>;

// Simplex
template <class Real, class F, class Policy = default_policy>
auto integrate_simplex(const F& f, const simplex<Real>& tri_or_simplex,
                       Real abs_tol, Real rel_tol,
                       Policy const& pol = Policy{}) -> result<Real>;

// Vector-valued overloads return std::vector<Real> with same semantics.

// Parallelism
struct execution_options {
  std::size_t max_threads = std::thread::hardware_concurrency();
  std::uint64_t max_eval  = 0;
  bool deterministic      = true;  // fixed work splitting for determinism
};
template <class Real, class F, class Domain, class Policy = default_policy>
auto integrate(const F& f, const Domain& dom, /* algorithm tag */,
               execution_options const& ex = {}, Policy const& pol = Policy{});

} // namespace
```

- `integrate_adaptive` is modeled after DCUHRE/Cuhre with embedded error estimates from **Genz–Malik** rules. citeturn2view0turn3view0
- `integrate_sparse_grid` implements **Smolyak** with nested 1‑D rules (**Clenshaw–Curtis** default). citeturn1news0turn10search0
- `integrate_qmc` uses **Boost.Random Sobol’**; optional **Owen scrambling** for replicates. citeturn14view0turn1news2
- `integrate_simplex` supports triangles/tetrahedra/… via Duffy transform + embedded rules. citeturn4search1turn13search1

---

## 3) Concurrency & Determinism

- Threading via `boost::asio::thread_pool`; we post **region tasks** (adaptive) or **batch point blocks** (sparse‑grid/QMC). To keep **deterministic** results irrespective of thread count, we (a) use fixed traversal orders, (b) join partial sums in a reproducible tree, and (c) avoid data races via local accumulators. citeturn6search10
- Vector‑valued integrands are evaluated **once per node** and accumulated component‑wise.
- `execution_options::deterministic=true` fixes split sizes and Sobol’ index ranges per thread.

---

## 4) Memory & Workspace

- **`workspace`** holds scratch vectors (abscissas, function values), priority queues for adaptive subdivision, and QMC RNG state. Users can supply a workspace to amortize allocations across calls; otherwise a local workspace is constructed.
- Small‑dimensional coordinate buffers use SSO via `boost::container::small_vector` (implementation detail) to avoid heap traffic.
- For QMC, replicates share generated points when **tent transform** disabled; with Owen scrambling enabled, each replicate applies an independent digital scramble (different seed) to the same index set. citeturn1news2

---

## 5) Domains & Transforms

- **Hypercube** `[a_i,b_i]` — internally mapped to `[0,1]^d`, applying Jacobians.
- **Infinite limits** `(-∞,∞)`, `(0,∞)`, `(-∞,b]`, … use standard monotone maps with square‑integrable images, e.g., `x = tan(π(u−1/2))` or `x = u/(1−u)` with weight; sensible defaults per algorithm. citeturn4search2
- **Simplex** — Duffy map from `[0,1]^d` → simplex, good for **edge singularities**; embedded rules also available. citeturn4search1

For QMC, we expose **tent (Baker)** transform to improve non‑periodic integrands before Sobol’ sampling; toggled on by default. citeturn2search3

---

## 6) Policies & Error Semantics

- We forward to `boost::math::policies` for overflow/underflow/indeterminate handling and domain errors; user may pass a custom policy type to all integrators. citeturn0search1
- Each algorithm reports: estimate, **absolute error bound** (deterministic: difference of embedded rules or surplus norm; QMC: replicate‑wise standard error), evaluation counts, and a status flag (`ok`, `maxeval_reached`, `cancelled`).

---

## 7) GPU & Heterogeneous Notes (future‑proof)

- Boost.Math’s double‑exponential quadrature exposes **CUDA device** entry points; we adopt the same **free‑function on device** pattern should we later add device kernels (requires `__device__` integrands). Phase‑1 ships host‑only. citeturn8view0

---

## 8) Internal Modules (brief)

- `detail::gm_rules.hpp` — constexpr tables for degree‑7/5 and 9/7 fully symmetric rules; layout mirrors **Genz–Malik** and **Cuhre**. citeturn3view0turn12search5
- `detail::cc_rules.hpp` — FFT‑based construction of Clenshaw–Curtis weights with nested growth; cached by order. citeturn10search0
- `detail::sobol_owen.hpp` — per‑digit XOR scrambles (base‑2) à la **Owen**; integrator requests `k` replicates for variance estimates. citeturn1news2

---

## 9) Integration with the rest of Boost

- **Random**: use `boost::random::sobol_engine` for deterministically reproducible Sobol’ points (QMC/RQMC). citeturn14view0
- **Multiprecision**: tested with `cpp_dec_float_50/100`; algorithms avoid assuming IEEE‑754 and respect user precision. citeturn6search1
- **Asio**: `thread_pool` for parallel evaluation; minimal dependency surface. citeturn6search10

---

## 10) Example (sketch)

```cpp
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/qmc.hpp>
using Real = double;
auto f = [](std::span<const Real> x)->Real {
  return 1.0 / (1.0 - std::cos(x[0])*std::cos(x[1])*std::cos(x[2]));
};
boost::math::cubature::hypercube<Real> box{{0,0,0}, {M_PI,M_PI,M_PI}};
auto I1 = boost::math::cubature::integrate_adaptive(f, box, 1e-8, 1e-8);
auto [I2, se] = boost::math::cubature::integrate_qmc(f, box, /*N*/ 1<<16, /*rep*/ 8);
```

For a deeper API, see `API.md`.
