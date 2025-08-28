# Genz–Malik 7/5 Cubature (4D, [-1,1]^4) — **Exact Nodes & Weights**

**Rule degree:** 7 with embedded 5

**Abscissae (dimension‑independent):**

- λ₂² = 9/70 ⇒ λ₂ = 3/√70 ≈ 0.358568582800318059
- λ₃² = 9/10 (also used on diagonals) ⇒ λ₃ = 3/√10 ≈ 0.948683298050513768
- λ₅² = 9/19 ⇒ λ₅ = 3/√19 ≈ 0.688247201611685178

## Degree‑7 rule (group weights, apply to **every** node in the group)

- center: **W₀ = -92032/6561 ≈ -14.0271300107**
- axis λ₂: **W₁ = 15680/6561 ≈ 2.38987959153**  (count = 8)
- axis λ₃: **W₂ = 3520/19683 ≈ 0.178834527257**  (count = 8)
- face‑diagonals λ₃: **W₃ = 3200/19683 ≈ 0.162576842961**  (count = 24)
- full‑diagonals λ₅: **W₄ = 6859/19683 ≈ 0.348473301834**  (count = 16)

> **Checks.** Node count = 57; constant function sum = W₀ + 2d·W₁ + 2d·W₂ + 2d(d−1)·W₃ + 2^d·W₄ = 16 (exactly 16).

## Embedded degree‑5 rule (no corners; different weights)

- center: **V₀ = -12112/243 ≈ -49.8436213992**
- axis λ₂: **V₁ = 1960/243 ≈ 8.0658436214**  (count = 8)
- axis λ₃: **V₂ = -40/27 ≈ -1.48148148148**  (count = 8)
- face‑diagonals λ₃: **V₃ = 400/729 ≈ 0.548696844993**  (count = 24)

> **Checks.** Node count = 41; constant function sum = V₀ + 2d·V₁ + 2d·V₂ + 2d(d−1)·V₃ = 16 (exactly 16).

## Download artifacts

- CSV (deg‑7 nodes): [genz_malik_4d_7_5_nodes_deg7.csv](sandbox:/mnt/data/genz_malik_4d_7_5_nodes_deg7.csv)
- CSV (deg‑5 nodes): [genz_malik_4d_7_5_nodes_deg5.csv](sandbox:/mnt/data/genz_malik_4d_7_5_nodes_deg5.csv)
- C++ header: [genz_malik_4d_7_5.hpp](sandbox:/mnt/data/genz_malik_4d_7_5.hpp)
