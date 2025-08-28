# Genz–Malik 2D **7/5** Cubature Rule (Exact Values)

**Domain:** standardized square \([-1,1]^2\).  
**Use:** degree‑7 rule with embedded degree‑5 rule for error estimation (the “7/5” pair).

---

## Abscissae (exact radicals)

Let
\[
\lambda_2=\sqrt{\tfrac{9}{70}}=\frac{3}{\sqrt{70}},\quad
\lambda_3=\sqrt{\tfrac{9}{10}}=\frac{3}{\sqrt{10}},\quad
\lambda_4=\lambda_3,\quad
\lambda_5=\sqrt{\tfrac{9}{19}}=\frac{3}{\sqrt{19}}.
\]

Node sets (2D):
- **Center**: \((0,0)\).
- **Axis (λ₂)**: \((\pm\lambda_2,0), (0,\pm\lambda_2)\) — 4 nodes.
- **Axis (λ₃)**: \((\pm\lambda_3,0), (0,\pm\lambda_3)\) — 4 nodes.
- **Diagonals (λ₄=λ₃)**: \((\pm\lambda_3,\pm\lambda_3)\) (all 4 sign combos) — 4 nodes.
- **Corners (λ₅)**: \((\pm\lambda_5,\pm\lambda_5)\) (all 4 sign combos) — 4 nodes.

Total nodes (degree‑7): **17**.  
Embedded degree‑5 rule uses all nodes **except** the λ₅ corners (so **13** nodes).

---

## Weights (per‑node, exact rationals, n=2)

### Degree‑7 rule
- \(w_1 = -\frac{1696}{2187}\) for \((0,0)\).
- \(w_2 = \frac{3920}{6561}\) for each axis‑λ₂ node.
- \(w_3 = \frac{1360}{6561}\) for each axis‑λ₃ node.
- \(w_4 = \frac{800}{19683}\) for each diagonal (λ₃,λ₃) node.
- \(w_5 = \frac{6859}{19683}\) for each corner (λ₅,λ₅) node.

**Check:** \( w_1 + 4w_2 + 4w_3 + 4w_4 + 4w_5 = 4 \) (area of \([-1,1]^2\)).

### Embedded degree‑5 rule (no λ₅ corners)
- \( \tilde w_1 = -\frac{3884}{729} \) for \((0,0)\).
- \( \tilde w_2 = \frac{490}{243} \) for each axis‑λ₂ node.
- \( \tilde w_3 = \frac{130}{729} \) for each axis‑λ₃ node.
- \( \tilde w_4 = \frac{100}{729} \) for each diagonal (λ₃,λ₃) node.

**Check:** \( \tilde w_1 + 4\tilde w_2 + 4\tilde w_3 + 4\tilde w_4 = 4 \).

---

## High‑precision numeric values (for sanity)

\(
\lambda_2 \approx 0.3585685828003180919906451539,\;
\lambda_3=\lambda_4 \approx 0.9486832980505137995996680633,\;
\lambda_5 \approx 0.6882472016116852977216287343.
\)

Degree‑7 per‑node weights:
- \(w_1 \approx -0.7754915409236397\)
- \(w_2 \approx 0.5974698978814205\)
- \(w_3 \approx 0.20728547477518672\)
- \(w_4 \approx 0.04064421074023269\)
- \(w_5 \approx 0.3484733018340700\)

Degree‑5 per‑node weights:
- \(\tilde w_1 \approx -5.327846364883402\)
- \(\tilde w_2 \approx 2.016460905349794\)
- \(\tilde w_3 \approx 0.17832647462277093\)
- \(\tilde w_4 \approx 0.13717421124828533\)

---

## Coordinate lists

### Degree‑7 (17 nodes)
- **Center**: \((0,0)\) → \(w_1\).
- **Axis λ₂**: \((\pm\lambda_2,0), (0,\pm\lambda_2)\) → \(w_2\) each.
- **Axis λ₃**: \((\pm\lambda_3,0), (0,\pm\lambda_3)\) → \(w_3\) each.
- **Diagonals λ₃**: \((\pm\lambda_3,\pm\lambda_3)\) → \(w_4\) each.
- **Corners λ₅**: \((\pm\lambda_5,\pm\lambda_5)\) → \(w_5\) each.

### Degree‑5 (13 nodes)
- As above **minus** the 4 λ₅ corner nodes; weights \( \tilde w_i\) as listed.

---

## Scaling to a general rectangle \([a_1,b_1]\times[a_2,b_2]\)

Let \(c_j=\tfrac{a_j+b_j}{2}\), \(d_j=\tfrac{b_j-a_j}{2}\).  
Map each node \(\mathbf{x}=(x_1,x_2)\in[-1,1]^2\) to \(\mathbf{y}=(c_1+d_1x_1,\; c_2+d_2x_2)\).  
Multiply every weight by \(d_1 d_2\).

---

## Sources

- A. C. Genz & A. A. Malik, **“Remarks on algorithm 006: An adaptive algorithm for numerical integration over an N‑dimensional rectangular region”**, *J. Comput. Appl. Math.*, 6(4), 1980. (explicit \( \lambda_i^2 \) and \( w_i \) formulas).  
- Verified against open-source implementations in SciPy and Stan Math (see comment block in the header for links).