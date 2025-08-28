# Genz–Malik 7/5 Cubature (5D, [-1,1]^5) — **Exact Nodes & Weights**

**Rule degree:** 7 with embedded 5

**Abscissae (dimension‑independent):**

- λ₂² = 9/70 ⇒ λ₂ = 3/√70 ≈ 0.358568582800318059
- λ₃² = 9/10 (also used on diagonals) ⇒ λ₃ = 3/√10 ≈ 0.948683298050513768
- λ₅² = 9/19 ⇒ λ₅ = 3/√19 ≈ 0.688247201611685178

## Degree‑7 rule (group weights, apply to **every** node in the group)

- center: **W₀ = -242944/6561 ≈ -37.0285017528**
- axis λ₂: **W₁ = 31360/6561 ≈ 4.77975918305**  (count = 10)
- axis λ₃: **W₂ = -640/2187 ≈ -0.29263831733**  (count = 10)
- face‑diagonals λ₃: **W₃ = 6400/19683 ≈ 0.325153685922**  (count = 40)
- full‑diagonals λ₅: **W₄ = 6859/19683 ≈ 0.348473301834**  (count = 32)

> **Checks.** Node count = 93; constant function sum = W₀ + 2d·W₁ + 2d·W₂ + 2d(d−1)·W₃ + 2^d·W₄ = 32 (exactly 32).

## Embedded degree‑5 rule (no corners; different weights)

- center: **V₀ = -88672/729 ≈ -121.635116598**
- axis λ₂: **V₁ = 3920/243 ≈ 16.1316872428**  (count = 10)
- axis λ₃: **V₂ = -3760/729 ≈ -5.15775034294**  (count = 10)
- face‑diagonals λ₃: **V₃ = 800/729 ≈ 1.09739368999**  (count = 40)

> **Checks.** Node count = 61; constant function sum = V₀ + 2d·V₁ + 2d·V₂ + 2d(d−1)·V₃ = 32 (exactly 32).

## Download artifacts

- CSV (deg‑7 nodes): [genz_malik_5d_7_5_nodes_deg7.csv](sandbox:/mnt/data/genz_malik_5d_7_5_nodes_deg7.csv)
- CSV (deg‑5 nodes): [genz_malik_5d_7_5_nodes_deg5.csv](sandbox:/mnt/data/genz_malik_5d_7_5_nodes_deg5.csv)
- C++ header: [genz_malik_5d_7_5.hpp](sandbox:/mnt/data/genz_malik_5d_7_5.hpp)
