# Boost.Math::cubature **SimLab** — Hierarchical Prompt Roadmap (Pre‑Implementation V1)

_A phase‑structured, dependency‑aware plan that decomposes each capability into buildable sub‑prompts for a production‑grade N‑D integration suite in Boost.Math._

**Goal.** Stand up a deterministic, test‑driven development flow that exercises architecture → algorithms → API → tests → docs **before** full implementation is complete by using a prompt‑first approach. Each phase specifies concrete prompts to generate code, tests, and docs that align with Boost norms and our spec pack. The plan layers cleanly and is designed for iterative PRs against `boostorg/math`.

> Primary sources for scope, algorithms, API surface, and acceptance bars are the **Proposal Pack** and its constituent docs (Architecture, Algorithms, Technical Blueprint, API, Testing/Validation, References). fileciteturn0file0 fileciteturn0file6 fileciteturn0file5 fileciteturn0file4 fileciteturn0file3 fileciteturn0file2 fileciteturn0file1

---

## Legend

- **Status:** **(NOW)** implement immediately • **(SEQ)** queued next in sequence • **(LATER)** post‑V1 / optional stretch.
- **Each prompt** lists: **Context, Objective, Depends‑On, Prompts (A…Z), Deliverables, Acceptance, Risks/Notes, Out‑of‑Scope.**
- **Milestones** and a **global dependency graph** appear at the end.
- **Source of Truth.** This document is the sole planning artifact for the prompt‑driven build. Update statuses (NOW/SEQ/LATER), append “‑ done” to completed phase bullets, and adjust milestones/gates inline. Other tickets derive from this roadmap.

---

## Status & Next Actions (2025‑08‑25)

### Current Understanding (Specs + Targets)

- **Spec pack baselined.** Architecture, Algorithms, Technical Blueprint, API draft, and Testing/Validation are ready and cross‑consistent; references collected. fileciteturn0file6 fileciteturn0file5 fileciteturn0file4 fileciteturn0file3 fileciteturn0file2
- **Scope & DoD.** Phase‑1 DoD: adaptive (Genz–Malik/DCUHRE), sparse‑grid (Smolyak+CC), QMC/RQMC (Sobol+Owen), transforms (improper, Duffy), vector outputs, policies/multiprecision, deterministic parallelism, CI tests/benchmarks. fileciteturn0file6
- **Acceptance gates** defined in Testing/Validation + Blueprint (accuracy, monotone refinement, perf parity vs. Cuhre/cubature, determinism, CI matrix). fileciteturn0file2 fileciteturn0file4

### Immediate Next Actions

1) **P00 — Prompt Foundations & Determinism** (NOW): bootstrap repo scaffolding prompts, code‑style guards, and a reproducible “design PR” generator.
2) **P01 — Public API Shells + Concepts** (NOW): generate headers/guards/namespaces + placeholder bodies with TODOs matching `API.md`. fileciteturn0file3
3) **P02 — Embedded Rules (Genz–Malik 7/5, 9/7)** (NOW): derive constexpr tables + exactness tests; wire `detail::gm_rules.hpp`. fileciteturn0file5
4) **P03 — DCUHRE Controller (adaptive)** (SEQ): region queue, split policy, error aggregation; minimal demo on Genz functions. fileciteturn0file5
5) **P04 — Clenshaw–Curtis & Smolyak** (SEQ): nested CC, DCT‑I weights, Smolyak assembly + surplus tests. fileciteturn0file5

**Success Criteria (near‑term).** `integrate_adaptive` header compiles; rule tables pass polynomial exactness; a toy driver integrates Genz‑Gaussian at d≤3 to 1e‑10; CI green on GCC/Clang; docs stub build. fileciteturn0file2

---

## V1 Critical Path (Phased Prompts)

> Each **Pxx** block specifies concrete LLM prompts. Use them as PR‑sized vertical slices. Prefer **header‑only**, Boost policies, deterministic parallelism, and multiprecision readiness per spec. fileciteturn0file6 fileciteturn0file4

### **P00 — Prompt Foundations & Determinism** (NOW)

**Context.** We need reproducible, review‑ready generations: fixed style, error policies, tests, and doc stubs aligned with Boost.Math practice.  
**Objective.** Establish the “meta‑prompts” and scaffolding that all later prompts build upon.  
**Depends‑On.** Spec pack (Architecture/Blueprint/API/Tests). fileciteturn0file6 fileciteturn0file4 fileciteturn0file3 fileciteturn0file2

**Prompts**  
- **Prompt A — Repo/Headers Skeleton.**  
  _Context:_ Copy `Package Layout` and `Public API (high level sketch)` into system prompt; require header‑only, `detail::` internals, include guards, and `BOOST_MATH_CUBATURE_EXPERIMENTAL`. _Ask:_ “Generate empty headers exactly as listed with namespaces, doxygen stubs, and TODO blocks, no implementations; add unit test skeletons compiling a ‘hello integration’ stub.” fileciteturn0file6
- **Prompt B — Policies & Result Type._**  
  _Ask:_ “Create `policies.hpp` adapter aliases and `result<T>` with {value,error,evaluations,status}; integrate with `boost::math::policies::policy<>` and document semantics.” fileciteturn0file3
- **Prompt C — CI Matrix Seed.**  
  _Ask:_ “Emit a CMake + GitHub Actions matrix (GCC/Clang/MSVC; C++17/20; linux/macos/windows) with warnings‑as‑errors and Boost.Test bootstrap.” fileciteturn0file2

**Deliverables.** Header shells, minimal tests, CI bootstrap.  
**Acceptance.** All headers compile; CI green on skeleton; docs build; coding guidelines pass.  
**Risks/Notes.** Guard against accidental non‑header symbols.  
**Out‑of‑Scope.** Real math kernels.

---

### **P01 — Public API Shells + Concepts** (NOW)

**Context.** Lock visible surface early to minimize churn.  
**Objective.** Implement `concepts.hpp`, domain types, and function signatures from `API.md`, with doxygen and examples. fileciteturn0file3  
**Depends‑On.** P00.

**Prompts**  
- **Prompt A — Concepts & Domains.**  
  _Ask:_ “Implement `ScalarIntegrand`/`VectorIntegrand` concepts, `hypercube<T>` and `simplex<T>` domain types, with sanity checks & small examples.” fileciteturn0file3
- **Prompt B — Result & Status.**  
  _Ask:_ “Emit `result<Real>` and status codes; integrate into API signatures; show usage snippets from `API.md`.” fileciteturn0file3

**Deliverables.** `concepts.hpp`, domain structs, `result<>`.  
**Acceptance.** Examples compile & run; SFINAE/concepts behave; docs cross‑link.  
**Risks/Notes.** Be careful with `std::span` availability on older libstdc++; provide fallback typedefs.  
**Out‑of‑Scope.** Algorithm bodies.

---

### **P02 — Embedded Rules (Genz–Malik 7/5, 9/7) + Exactness Tests** (NOW)

**Context.** Deterministic adaptive needs embedded fully symmetric rules and intrinsic errors.  
**Objective.** Implement `detail::gm_rules.hpp` with constexpr tables and generators; add polynomial exactness tests up to degree. fileciteturn0file5  
**Depends‑On.** P01.

**Prompts**  
- **Prompt A — Tables & Layout.**  
  _Ask:_ “Emit constexpr nodes/weights for degree‑7/5 and 9/7 fully symmetric Genz–Malik pairs up to d=15; provide compile‑time accessors and unit tests for monomial exactness.” fileciteturn0file5
- **Prompt B — Generator (Optional).**  
  _Ask:_ “Provide a small offline generator (test‑only) that solves moment equations to reproduce tables; document derivation.” fileciteturn0file5

**Deliverables.** `detail/gm_rules.hpp`, tests.  
**Acceptance.** Exactness tests pass to machine epsilon; coverage on d=2..8 (quick) and d=15 (nightly).  
**Risks/Notes.** Code size from tables; gate larger degrees via macro.  
**Out‑of‑Scope.** Subdivision controller.

---

### **P03 — DCUHRE Adaptive Controller** (SEQ)

**Context.** Global max‑error refinement with axis bisection; normed errors for vector outputs.  
**Objective.** Implement `adaptive.hpp`: region PQ, embedded estimates, safeguarded error, axis‑selection, split, and global stop conditions; vector path and cancellation. fileciteturn0file5 fileciteturn0file6  
**Depends‑On.** P02.

**Prompts**  
- **Prompt A — Region & PQ.**  
  _Ask:_ “Implement `region{a,b,estimate_fine,estimate_coarse,error,split_dim,evals}` + PQ keyed by error; two‑sum accumulation; policy‑driven exceptions.” fileciteturn0file4
- **Prompt B — Axis Selection & Split.**  
  _Ask:_ “Compute directional variation from symmetric node pairs; split on argmax; bisect bounds; re‑use mirrored node evaluations when possible.” fileciteturn0file5
- **Prompt C — Global Loop.**  
  _Ask:_ “Write the DCUHRE loop with absolute/relative tol per `TESTING_VALIDATION.md`; expose `max_eval`, `max_regions`, `max_depth`.” fileciteturn0file2

**Deliverables.** `adaptive.hpp` initial implementation + tests on Genz functions.  
**Acceptance.** Meets requested tol on smooth Genz at d≤10; monotone error decrease; deterministic results across thread counts (single‑thread first). fileciteturn0file2  
**Risks/Notes.** Guard for spuriously low embedded differences → fallback spread estimate.  
**Out‑of‑Scope.** Parallelism (handled in P09).

---

### **P04 — Clenshaw–Curtis (CC) & Smolyak Sparse‑Grid** (SEQ)

**Context.** Medium‑dim smooth integrands need sparse grids with nested 1‑D rules.  
**Objective.** Implement `cc_rules.hpp` (DCT‑I weights), `sparse_grid.hpp` (Smolyak, fixed level), node de‑duplication, and basic surplus‑based error. fileciteturn0file5  
**Depends‑On.** P01.

**Prompts**  
- **Prompt A — CC weights with nested growth.**  
  _Ask:_ “Implement slow‑growth CC (1,3,5,9,…) weights via DCT‑I; cache by order; prove nestedness in tests.” fileciteturn0file5
- **Prompt B — Smolyak Assembly.**  
  _Ask:_ “Build \(A(\ell,d)\) with integer coefficients; hash nodes to dedupe (ULP‑aware); evaluate once; return value & cheap error via level difference.” fileciteturn0file5

**Deliverables.** `sparse_grid.hpp`, tests on tensor products.  
**Acceptance.** Matches reference integrals; point count ≪ tensor; deterministic node ordering. fileciteturn0file2  
**Risks/Notes.** GP rules optional; keep CC default.  
**Out‑of‑Scope.** Dimension‑adaptive surplus (P20).

---

### **P05 — QMC / RQMC (Sobol + Owen) + Tent Transform** (SEQ)

**Context.** High‑d or rough integrands benefit from (randomized) QMC.  
**Objective.** Implement `qmc.hpp`: Sobol engine (Boost.Random), replicates, Owen scrambling, tent transform, variance estimate; fixed index partitioning for determinism. fileciteturn0file6 fileciteturn0file5  
**Depends‑On.** P01.

**Prompts**  
- **Prompt A — Base QMC.**  
  _Ask:_ “Wrap `boost::random::sobol_engine`; map to [a,b]^d; implement tent transform; accumulate value; return estimate.” fileciteturn0file3
- **Prompt B — RQMC & Error Bars.**  
  _Ask:_ “Add Owen per‑digit scrambling (seeded); `k` replicates; report std error via replicate variance; deterministic block partition per thread.” fileciteturn0file5

**Deliverables.** `qmc.hpp` with tests on smooth vs IID MC.  
**Acceptance.** Convergence beats IID MC on smooth Genz; error bars sensible; deterministic across thread counts. fileciteturn0file2  
**Risks/Notes.** Scrambling cost vs benefit; allow k=1 fast path.  
**Out‑of‑Scope.** Brownian bridge.

---

### **P06 — Transforms: Improper Bounds & Simplex (Duffy)** (SEQ)

**Context.** Practical integrals require ∞ limits and simplex domains.  
**Objective.** Implement `transforms.hpp` mappings for infinite/semi‑infinite axes and the Duffy map; add helpers and Jacobians; demos and tests. fileciteturn0file5  
**Depends‑On.** P01.

**Prompts**  
- **Prompt A — Infinite maps.**  
  _Ask:_ “Provide standard bijections (e.g., `x = u/(1-u)`; `x = tan(pi*(u-0.5))`) with Jacobians; unit tests with known integrals.” fileciteturn0file5
- **Prompt B — Duffy.**  
  _Ask:_ “Implement `[0,1]^d → simplex` and Jacobian; add triangle/tetrahedra wrappers for `integrate_simplex` that call adaptive/QMC backends.” fileciteturn0file5

**Deliverables.** `transforms.hpp`, `simplex.hpp` scaffolding + tests.  
**Acceptance.** Polynomial integrals on simplex exact to tol; improper Gaussians match references. fileciteturn0file2  
**Risks/Notes.** Conditioning of `tan` mapping; document alternatives.  
**Out‑of‑Scope.** Native simplex embedded rules (P19).

---

### **P07 — Vector‑Valued Integrands (Zero‑Alloc Path)** (SEQ)

**Context.** Users often integrate m‑component outputs; we must evaluate once per node.  
**Objective.** Implement `(x, out*, m)` callback style; L∞/L2 norms; scatter/gather utilities; AoS↔SoA adapters. fileciteturn0file6  
**Depends‑On.** P01, P03/P04/P05 surfaces.

**Prompts**  
- **Prompt A — Callback contract.**  
  _Ask:_ “Add vector overloads to all entry points; zero allocations in hot path; comprehensive examples.” fileciteturn0file3
- **Prompt B — Norms & tests.**  
  _Ask:_ “Error aggregation via L∞ (default) and L2; unit tests comparing component‑wise vs scalar runs.” fileciteturn0file2

**Deliverables.** Vector path in all algorithms; docs/examples.  
**Acceptance.** No re‑evaluation per node; tests pass; perf within 3% vs scalar.  
**Risks/Notes.** Avoid temporary `std::vector` in inner loops.

---

### **P08 — Policies & Multiprecision Integration** (SEQ)

**Context.** Align with Boost policies and support `cpp_dec_float_50/100`.  
**Objective.** Thread policies through API; implement two‑sum/pairwise accumulation; multiprecision smoke tests. fileciteturn0file6  
**Depends‑On.** P01–P07.

**Prompts**  
- **Prompt A — Policy threading.**  
  _Ask:_ “Propagate `Policy` to error handling; add tests that flip policy behaviors.” fileciteturn0file3
- **Prompt B — Multiprecision.**  
  _Ask:_ “Compile & run with `cpp_dec_float_50`; add tests to ensure no FMA assumptions; document performance caveats.” fileciteturn0file2

**Deliverables.** Policy adherence; multiprecision tests.  
**Acceptance.** CI matrix green incl. multiprecision smoke; policies honored.

---

### **P09 — Deterministic Parallelism & Workspace** (SEQ)

**Context.** Scale with threads but preserve bitwise stability.  
**Objective.** `boost::asio::thread_pool` integration; deterministic split and reduction; reusable `workspace` to amortize allocations. fileciteturn0file6 fileciteturn0file4  
**Depends‑On.** P03/P04/P05.

**Prompts**  
- **Prompt A — Thread pool plumbing.**  
  _Ask:_ “Post region/batch tasks; bind index ranges per thread; reproducible tree reduction.” fileciteturn0file4
- **Prompt B — Workspace.**  
  _Ask:_ “Implement `workspace` with scratch buffers, PQ, Sobol/Owen state; API to reuse across calls.” fileciteturn0file6

**Deliverables.** Parallel paths; workspace type.  
**Acceptance.** Identical results across thread counts; speedup ≥ (0.7×threads) on compute‑bound demos.

---

### **P10 — Testing & Validation Pack (Genz/Keister/Simplex)** (SEQ)

**Context.** Lock correctness and stability.  
**Objective.** Implement tests/fixtures/bench harness as per `TESTING_VALIDATION.md`. fileciteturn0file2  
**Depends‑On.** P03–P06/P08/P09.

**Prompts**  
- **Prompt A — Genz Suite.**  
  _Ask:_ “Add oscillatory, corner/product peak, Gaussian, continuous/discontinuous tests 2–12D with references.” fileciteturn0file2
- **Prompt B — Gaussian/Keister & Simplex.**  
  _Ask:_ “Closed‑form Gaussians; Keister‐like; simplex polynomials exactness.” fileciteturn0file2
- **Prompt C — CI Matrix.**  
  _Ask:_ “Wire compilers/OS/types; nightly heavy runs.” fileciteturn0file2

**Deliverables.** Comprehensive tests + CI jobs.  
**Acceptance.** All gates green at target tols; determinism probes pass.

---

### **P11 — Benchmarks & Cross‑Library Comparisons** (SEQ)

**Context.** Demonstrate competitiveness and find regressions.  
**Objective.** Bench scripts comparing to CUBA/Cuhre and `cubature` (behavior only), and to Boost naive MC; plots + CSVs. fileciteturn0file2  
**Depends‑On.** P10.

**Prompts**  
- **Prompt A — Micro/Macro benches.**  
  _Ask:_ “Add eval/s, error vs time, scaling, determinism; emit JSON+CSV.” fileciteturn0file2
- **Prompt B — Comparators.**  
  _Ask:_ “Run external tools if available; otherwise use published results as reference; document methodology.” fileciteturn0file1

**Deliverables.** Bench harness + docs.  
**Acceptance.** ≤1.5× Cuhre on select tests; QMC beats IID MC at equal budget. fileciteturn0file2

---

### **P12 — Documentation (Quickbook) & Examples** (SEQ)

**Context.** Make the library discoverable and reviewable.  
**Objective.** Author a **Cubature** chapter (Overview, Adaptive, Sparse‑grid, QMC/RQMC, Transforms, Vector outputs, Parallelism, Examples) with runnable snippets. fileciteturn0file6  
**Depends‑On.** P01–P10.

**Prompts**  
- **Prompt A — Chapter Outline.**  
  _Ask:_ “Generate Quickbook stubs and example code mirroring `API.md`; include policy/multiprecision sections.” fileciteturn0file3
- **Prompt B — How‑to Guides.**  
  _Ask:_ “Write tasks: improper integrals, anisotropic sparse‑grid, vector‑valued, QMC with error bars.” fileciteturn0file6

**Deliverables.** Docs & examples.  
**Acceptance.** Docs build; examples compile & pass; linkcheck clean.

---

### **P13 — Review Hardening & Error Semantics** (SEQ)

**Context.** Align semantics (status codes, cancellation, errors) with Boost policies.  
**Objective.** Finalize error/status wording; add cancellation token; ensure policy hooks honored. fileciteturn0file6  
**Depends‑On.** P08–P10.

**Prompts**  
- **Prompt A — Semantics sweep.**  
  _Ask:_ “Consistency pass across algorithms; document behaviors; tests toggling policies.” fileciteturn0file3

**Deliverables.** Semantics doc + tests.  
**Acceptance.** Policy tests pass; reviewers sign‑off.

---

### **P14 — Packaging, CI Badges & Release Notes** (SEQ)

**Context.** Prepare for an experimental release branch.  
**Objective.** Add `BOOST_MATH_CUBATURE_EXPERIMENTAL` toggle, CI badges, and detailed release notes tying back to acceptance bars. fileciteturn0file4  
**Depends‑On.** P12–P13.

**Prompts**  
- **Prompt A — Release doc.**  
  _Ask:_ “Generate release notes mapping DoD to passing tests/benches; list known limitations.” fileciteturn0file2

**Deliverables.** Release notes + toggle.  
**Acceptance.** Branch builds; notes approved.

---

## Post‑V1 / Optional Stretch

### **P19 — Native Simplex Embedded Rules & Subdivision** (LATER)
- Add CUBPACK‑style embedded rules for triangles/tetrahedra; longest‑edge bisection; error indicators. fileciteturn0file5

### **P20 — Dimension‑Adaptive Sparse Grids** (LATER)
- Gerstner–Griebel surplus‑driven anisotropic refinement; per‑dimension weights and stop rules. fileciteturn0file5

### **P21 — Device Backends (CUDA/HIP) Pattern** (LATER)
- Follow Boost.Math DE GPU pattern for batched evaluations; device free‑function style, no host/device mixing in v1. fileciteturn0file6

---

## Milestones & Gating Tests

- **M0 (Week 1)** — **Scaffold**: P00–P01, API shells compile; CI green.
  - _Gate:_ headers & examples compile; docs build.
- **M1 (Week 2–3)** — **Rules + Adaptive**: P02–P03.
  - _Gate:_ Genz smooth tests pass to tol; monotone refinement; deterministic.
- **M2 (Week 3–4)** — **Sparse‑grid**: P04.
  - _Gate:_ tensor product checks & point‑reduction vs tensor.
- **M3 (Week 4–5)** — **QMC/RQMC**: P05.
  - _Gate:_ convergence beats IID MC; error bars stable.
- **M4 (Week 5–6)** — **Transforms & Vector**: P06–P07.
  - _Gate:_ simplex polynomials exact; vector path parity.
- **M5 (Week 6–7)** — **Policies/Parallel/Tests**: P08–P10.
  - _Gate:_ multiprecision smoke; deterministic across threads; CI matrix green.
- **M6 (Week 7–8)** — **Benches & Docs**: P11–P12.
  - _Gate:_ ≤1.5× Cuhre on benchmarks; docs/examples complete.
- **M7 (Week 9)** — **Harden & Release Notes**: P13–P14.

Gates tie directly to `TESTING_VALIDATION.md` and the Technical Blueprint acceptance bars. fileciteturn0file2 fileciteturn0file4

---

## Global Dependency Graph (top‑level)

```
P00 ──> P01 ──> P02 ──> P03 ───────────────────────┐
       │        │                                   ├─> P08 ──> P13 ──> P14
       │        └───────────────┐                   │
       │                        ├─> P04 ────────────┤
       │                        └─> P05 ────────────┤
       └─────────────────> P06 ─────────────────────┤
                              └─> P07 ──────────────┤
P03/P04/P05 ──> P09 ───────────┘                    ├─> P10 ──> P11 ──> P12
```

---

## Implementation Notes (Staff+)

- **Bitwise determinism**: fix traversal orders; disjoint Sobol index ranges per thread; tree reductions; two‑sum accumulation; workspace reuse. fileciteturn0file4
- **Error estimation**: embedded differences safeguarded by spread estimates; QMC uses replicate variance; sparse‑grid uses level deltas/surpluses. fileciteturn0file5
- **Vector path**: single integrand evaluation per node; L∞ default norm; configurable. fileciteturn0file6
- **Transforms**: ship safe defaults and document trade‑offs (e.g., tan map). fileciteturn0file5
- **Policies & multiprecision**: thread through all functions; avoid FMA assumptions; smoke tests for `cpp_dec_float_50`. fileciteturn0file2

---

## Prompt Style Guide (meta)

- **Grounding**: Always include short quotes or bullet references to the spec sections you are implementing (Architecture §1/§2, Algorithms §1/§2/§3, API signatures, Testing gates). fileciteturn0file6 fileciteturn0file5 fileciteturn0file3 fileciteturn0file2
- **Deliverables**: Headers first, tests next, docs last in each PR. No hidden global state. Header‑only.
- **Style**: Prefer free functions; `detail::` for internals; `constexpr` where safe; avoid dynamic allocation in hot paths.
- **Review hooks**: Ask for (a) monomial exactness prints for rules, (b) Genz function diffs, (c) determinism across threads (`N=1,2,4,8`).

---

## Run Manifest (Draft v0.1.0 — Code Gen Repro)

```
schema_version: 0.1.0
run_id: 2025-08-25_cubature_prompt_build
git_sha: <auto from repo>
compilers: [gcc-13, clang-17, msvc-19.3]
cpp_std: [17, 20]
generators:
  - name: p00_scaffold
    inputs: [ARCHITECTURE.md, API.md]
    seed: 12345
  - name: p02_rules
    inputs: [ALGORITHMS.md, TESTING_VALIDATION.md]
    seed: 202581
  - name: p03_adaptive
    inputs: [ALGORITHMS.md, TECHNICAL_BLUEPRINT.md]
    seed: 9981
targets:
  ci_matrix: true
  docs_build: true
artifacts:
  headers: true
  tests: true
  benches: true
  docs: true
```

Reproduce builds by fixing `run_id`, `git_sha`, and toolchain images in CI. fileciteturn0file2

---

## Open Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Rule table bloat | Code size concerns | Gate higher degrees by macro; optional generator; compress with `constexpr` symmetry | 
| Hidden nondeterminism | Flaky tests | Deterministic partition/reduction; explicit seeds; CI determinism probe | 
| Perf gaps vs Cuhre | Rejection during review | Micro‑optimize hot loops; cache nodes; batched vector path; provide benches | 
| Multiprecision regressions | CI red on Windows/MSVC | Add MSVC‑specific flags; reduce template instantiations; smoke tests | 
| Reviewer overload | Slows merge | PRs as vertical slices with tests/docs and clear acceptance criteria | 

---

## Out‑of‑Scope (V1)

- Vegas/Suave‑style importance sampling; native simplex embedded rules; GPU device kernels (see Post‑V1). fileciteturn0file4

---

_This roadmap is the authoritative plan for prompt‑driven implementation of Boost.Math::cubature. Keep it updated as phases land and acceptance tests evolve._
