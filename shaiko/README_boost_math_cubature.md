# Boost.Math::cubature — Proposal Pack · v0.9 (2025-08-25)

> **Goal.** Add a modern, production‑grade **N‑dimensional integration suite** to Boost.Math that complements existing 1‑D quadrature and naive Monte‑Carlo with **deterministic adaptive cubature**, **sparse‑grid Smolyak rules**, and **(randomized) quasi‑Monte‑Carlo** over hyperrectangles and simplices, plus robust mappings for improper/weakly singular integrals. The design is header‑only, policy‑driven, thread‑friendly, and integrates cleanly with `boost::math::policies`, `Boost.Multiprecision`, and `Boost.Random`.

**What’s inside**

- `ARCHITECTURE.md` — module layout, public API, policies, threading, workspace, GPU notes, Boost integration.
- `ALGORITHMS.md` — explicit math/algorithms (Genz–Malik/DCUHRE; Smolyak sparse grids; QMC/RQMC with Owen scrambling; simplex rules & Duffy; transforms).
- `TECHNICAL_BLUEPRINT.md` — step‑by‑step engineering plan, acceptance bars, milestones, perf/complexity, code skeletons, review plan.
- `API.md` — proposed user‑facing headers, function/class templates, examples, error semantics.
- `TESTING_VALIDATION.md` — test datasets, correctness & stability checks, benchmarks, CI matrix.
- `REFERENCES.md` — consolidated bibliography + anchor links to cited sources.

All documents are self‑contained and cross‑referenced with citations to Boost docs and peer‑reviewed sources.
