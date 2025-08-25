# Boost.Math::cubature — Technical Blueprint (engineering)

## 0) Principles

- **Single header per algorithm**, header‑only, no global state; reproducible by default.
- **Policy‑first**: every entry takes a `Policy` to integrate with `boost::math::policies`. citeturn0search1
- **Deterministic parallelism** using `boost::asio::thread_pool`. citeturn6search10
- **Dimension & precision agnostic**: no IEEE‑754 assumptions; support `cpp_dec_float_50`. citeturn6search1

---

## 1) Milestones & DoD

**M1 — Adaptive (hypercube)**  
- Implement degree‑7/5 embedded **Genz–Malik** rules (constexpr tables).  
- DCUHRE loop (global max‑heap, axis bisection, spread safeguard).  
- Vector‑valued integrands.  
- DoD: passes Genz suite up to d=10 at 1e‑8 on smooth cases; intrinsic error estimate monotone under refinement. citeturn2view0

**M2 — Sparse‑grid**  
- Clenshaw–Curtis 1‑D rules (nested slow‑growth), Smolyak assembly, surpluses.  
- DoD: matches reference results on smooth tensor integrands; points << tensor product. citeturn10search0

**M3 — QMC/RQMC**  
- Sobol’ generator (Boost.Random), tent transform, `k` replicates, Owen scramble.  
- DoD: empirical rate better than IID MC on smooth tests; CI on replicates sane. citeturn14view0turn1news2

**M4 — Simplex**  
- Duffy transform + adaptive hypercube; (optional) native embedded simplex rules.  
- DoD: triangles/tetrahedra validated on polynomial tests and Genz‑style functions. citeturn4search1

**M5 — Docs, examples, CI & benchmarks**  
- Complete docs, `examples/`, correctness + perf CI, comparison scripts vs. Cuba/cubature. citeturn3view0turn2view1

---

## 2) Data Structures

### 2.1 Region node (adaptive)
```cpp
template <class Real>
struct region {
  std::vector<Real> a, b;         // bounds (size d)
  Real estimate_coarse, estimate_fine;
  Real error;                      // |fine - coarse| (or safeguarded)
  unsigned split_dim;
  std::size_t evals;
};
```
Priority queue ordered by `error` (largest first). Regions store cached function values at mirrored nodes to avoid re‑evaluation after split when possible.

### 2.2 Sparse‑grid assembly
- Store **surplus coefficients** keyed by multi‑index; memoize 1‑D weights per order; combine by level with integer coefficients.
- Node de‑duplication by hashing coordinates (within ULP) to avoid double evaluation.

### 2.3 QMC state
- For each replicate: scramble seed, Sobol’ index range, per‑thread accumulators; optional **leapfrog**/block partitioning to keep cache‑friendly batches. citeturn14view0

---

## 3) Implementation Notes

### 3.1 Tables for embedded rules
- Adopt **Genz–Malik** node/weight sets (degree 7/5 and 9/7) that are used in **Cuhre**; place in `detail::gm_rules.hpp` as `constexpr` arrays with compile‑time dimension. citeturn12search5

### 3.2 Clenshaw–Curtis weights
- Build by DCT‑I following standard formulas; cache by order; ensure **nestedness** under slow growth (1, 3, 5, 9, …). citeturn10search0turn10search5

### 3.3 Owen scrambling
- Implement base‑2 **nested uniform scrambles** using random bitmasks per digit/coordinate; guarantee reproducibility via 64‑bit seeds. citeturn1news2

### 3.4 Duffy transform
- Provide `[0,1]^d → simplex` mapping and Jacobian; expose helpers `duffy_map(u)`/`duffy_jacobian(u)`. citeturn4search1

### 3.5 Improper maps
- Provide standard `finite_to_infinite(u)` maps with Jacobians selectable per axis; document conditioning trade‑offs. citeturn4search2

### 3.6 Precision & accumulation
- Use pairwise reductions; for multiprecision types, avoid fused‑multiply‑add assumptions; prefer `two_sum` helper.

### 3.7 Parallel evaluation
- Use `boost::asio::thread_pool` (`post`) with bounded task queue; deterministic **work stealing disabled** when `deterministic=true`. citeturn6search10

---

## 4) API Surface (condensed, see `API.md`)

- `integrate_adaptive`, `integrate_sparse_grid`, `integrate_qmc`, `integrate_simplex` with scalar & vector overloads.
- `execution_options` (threads, max_eval, deterministic).
- Result struct: value, error, evals, status.

---

## 5) Complexity & Costs

- **Adaptive**: \(O(P \log R)\) bookkeeping where \(P\) points per rule eval, \(R\) active regions; dominant cost is integrand evaluations.
- **Sparse grid**: \(O(N \log N)\) to assemble + evaluate (with hashing) where \(N\) is unique nodes.
- **QMC**: \(O(kN)\) function calls; small overhead for scrambling.

---

## 6) Acceptance Criteria (Phase‑1)

1. **Correctness**: on Genz suite (d ≤ 10), achieve requested tolerances for smooth instances; cross‑check vs **Cuhre** and **cubature** within 2× evaluations. citeturn2view1turn3view0
2. **Performance**: wall‑time within 1.5× of CUBA’s Cuhre on scalar integrands for d ≤ 8; QMC beats IID Monte‑Carlo of Boost.Math on smooth tests at equal budget. citeturn9view0
3. **Robustness**: deterministic across thread counts; cancellation support; passes multiprecision smoke tests. 

---

## 7) Risks & Mitigations

- **Rule availability for high d** → ship 7/5 and 9/7 up to practical d, fall back to tensor‑sparse (Smolyak) for higher d. citeturn12search0
- **Non‑smooth integrands** → prefer QMC/RQMC; expose transforms and user guidance. citeturn1news2
- **State bloat** (vector integrands) → zero‑alloc callback `(x, out*, m)` to avoid temporaries.

---

## 8) Example Usage (vector integrand)

```cpp
auto f = [](std::span<const double> x, double* out, std::size_t m) {
  out[0] = std::exp(-x[0]*x[0] - x[1]*x[1]);
  if (m > 1) out[1] = x[0] + x[1];
};
hypercube<double> box{{-5, -5}, {5, 5}};
auto [I, E, evals, status] = integrate_adaptive(f, box, 1e-9, 1e-9);
```

---

## 9) Files & Code Skeletons

```
include/boost/math/cubature/adaptive.hpp   // class adaptive_integrator<>
include/boost/math/cubature/sparse_grid.hpp
include/boost/math/cubature/qmc.hpp
include/boost/math/cubature/simplex.hpp
include/boost/math/cubature/transforms.hpp
include/boost/math/cubature/workspace.hpp
```

- Each header has a `detail::*` namespace for internals and a flat set of free functions for public API.

---

## 10) Future Extensions (Phase‑2)

- **Vegas+/Suave‑style** hybrids (adaptive importance + subregion sampling). citeturn11search8turn11search5
- **Dimension‑adaptive sparse grids** with anisotropic weights (Gerstner–Griebel). citeturn1news2
- **Device backends** following Boost.Math’s CUDA pattern for inner kernels. citeturn8view0
