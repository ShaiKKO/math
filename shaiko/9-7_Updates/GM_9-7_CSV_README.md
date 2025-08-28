# Genz–Malik 9/7 — Enumerated Nodes & Weights (CSV Pack, d=2..5)

This pack enumerates every evaluation node (fully symmetric orbits expanded) and the **per-node weights** for the **degree‑9 rule** with its **embedded degree‑7 rule** on the hypercube **[−1,1]^d**.  
To integrate on **[0,1]^d**, map nodes by `u = (x+1)/2` and **scale all weights by `2^−d`**.

> **Source of truth.** Constants and formulas are transcribed from **DCUHRE** (Algorithm 698) routines **D09HRE** (degree‑9 family) and **D07HRE** (degree‑7 family). See citations in the PR description.

## Files

- **2D**: [GM_9-7_d2_deg9.csv](GM_9-7_d2_deg9.csv) • [GM_9-7_d2_deg7.csv](GM_9-7_d2_deg7.csv)
- **3D**: [GM_9-7_d3_deg9.csv](GM_9-7_d3_deg9.csv) • [GM_9-7_d3_deg7.csv](GM_9-7_d3_deg7.csv)
- **4D**: [GM_9-7_d4_deg9.csv](GM_9-7_d4_deg9.csv) • [GM_9-7_d4_deg7.csv](GM_9-7_d4_deg7.csv)
- **5D**: [GM_9-7_d5_deg9.csv](GM_9-7_d5_deg9.csv) • [GM_9-7_d5_deg7.csv](GM_9-7_d5_deg7.csv)

Each CSV has columns:

```
group, group_label, node_index, x1, ... , xd, weight
```

- **group**: numeric group index mirroring DCUHRE column order.
- **group_label**: human‑readable orbit description.
- **node_index**: 1‑based index within its group.
- **weight**: per‑node weight on **[−1,1]^d**.

## Group Legend

**Degree‑9 (columns 1..9)**

1. `center` — (0,…,0)  
2. `axis_l1` — (±λ₁,0,…,0) (2d points)  
3. `axis_l2` — (±λ₂,0,…,0) (2d)  
4. `axis_l3` — (±λ₃,0,…,0) (2d)  
5. `axis_lp_null` — (±λₚ,0,…,0) (2d) **weight 0 in basic deg‑9**  
6. `pair_l1_l1` — (±λ₁,±λ₁,0,…,0) over unordered axis pairs (2d(d−1))  
7. `pair_l1_l2` — (±λ₁,±λ₂,0,…,0) over ordered axis pairs (4d(d−1))  
8. `triple_l1_l1_l1` — (±λ₁,±λ₁,±λ₁,0,…) over unordered triples (4/3 d(d−1)(d−2))  
9. `full_diag_l0` — (±λ₀,…,±λ₀) (2^d)

**Degree‑7 (columns 1..6)**

1. `center` — (0,…,0)  
2. `axis_l2` — (±λ₂,0,…,0) (2d)  
3. `axis_l1` — (±λ₁,0,…,0) (2d)  
4. `axis_lp_null` — (±λₚ,0,…,0) (2d) **weight 0 in basic deg‑7**  
5. `pair_l1_l1` — (±λ₁,±λ₁,0,…,0) (2d(d−1))  
6. `full_diag_l0` — (±λ₀,…,±λ₀) (2^d)

> λ₀² = 0.4707. The remaining radii (λ₁, λ₂, λ₃ for deg‑9; λ₁, λ₂ for deg‑7; and λₚ) are computed exactly as in **D09HRE/D07HRE** and then square‑rooted to obtain abscissae.

## Verification

We verified that **∑ weights = 2^d** (constant function over [−1,1]^d) and thus **= 1** after scaling by `2^−d` for [0,1]^d:

| d | deg‑9 sum | deg‑7 sum |
|---|-----------|-----------|
| 2 | 4.0 | 4.0 |
| 3 | 8.0 | 8.0 |
| 4 | 16.0 | 16.0 |
| 5 | 32.0 | 32.0 |

Floating‑point printing uses 18 decimals; recomputation yields the same totals within round‑off.

## Notes

- The `axis_lp_null` orbits are evaluated to support **null rules** used by DCUHRE’s error estimator; their weights are **0** in the **basic** 9 and 7 degree rules.
- To generate C++ `constexpr` arrays from these CSVs for Boost.Math, treat each row’s `weight` as the per‑node weight on [−1,1]^d (no extra division by orbit size).

