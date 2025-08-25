# Boost.Math::cubature — Testing & Validation

## 1) Correctness testbeds

- **Genz test functions** (oscillatory, corner peak, Gaussian, product peak, continuous, discontinuous) across dimensions 2–12 with parameter sweeps. Include reference values (closed forms or high‑precision baselines). citeturn5search0
- **Gaussian/Keister** integrals for sanity (closed forms in 1‑D; higher‑D products / known references). citeturn5search1turn5search9
- **Simplex polynomials**: integrate exactly to machine epsilon using known formulas (degree ≤ embedded degree).
- **Vector‑valued**: compare each component vs. scalar runs on same node sets.

## 2) Cross‑library comparisons

- **CUBA/Cuhre, Suave, Vegas**; ensure comparable accuracy/time on smooth problems; check adaptive behavior (nregions, eval counts) and failure modes on discontinuities. citeturn3view0turn11search11
- **cubature** (Steven G. Johnson): compare to `hcubature_pcubature` on hyperrectangles; check vector‑valued interface parity. citeturn2view1

## 3) CI Matrix

- Compilers: GCC ≥ 11, Clang ≥ 14, MSVC ≥ 19.3.
- C++: 17 (min), 20 (preferred for concepts).
- OS: Linux (x64, ARM64), Windows, macOS.
- Types: `double`, `long double`, `cpp_dec_float_50` (smoke). citeturn6search1

## 4) Benchmarks

- Measure **evals/s**, **error vs. time**, **scaling with threads**, **determinism across thread counts**.
- Data: Genz suite; synthetic peaks; smooth tensor products; ill‑conditioned cases (document expected failures).

## 5) Acceptance thresholds

- Correctness: within requested tol on smooth tests; monotone error decrease under refinement.
- Performance: ≤1.5× Cuhre on representative tests (scalar), parity on vector‑valued; QMC beats IID MC (Boost naive MC) at equal budgets. citeturn9view0
