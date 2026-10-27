# Certified Gram cuts in BNB

The root SDP relaxation now separates three families of linear inequalities
from complete, uniform binary Gram blocks: triangle cuts, first-level clique
cuts, and Turan count cuts. Detection requires an exact structural certificate;
separation is limited to the root and bounded by the existing BNB time limit.
Near-integral sparse-objective blocks can additionally propagate equal Gram
rows to child nodes. The MATLAB verifiers in this directory check detection,
cut validity, and propagation on synthetic feasible assignments.

## Why the count cut is valid

Let the diagonal be `d > 0` and each off-diagonal entry be either `d`
(an equal pair) or `b < 0` (a different pair). A pair with entry `d` has
identical Gram vectors. If `K+1` distinct vectors existed, their principal
matrix would have eigenvalue `d+K*b < 0`. The detector checks a strict
numerical margin before using `K = floor(1-d/b)`. Thus at most `K` distinct
classes exist. Among any `r = q*K+s` vertices, at least
`s*C(q+1,2)+(K-s)*C(q,2)` pairs must be equal. This yields a linear cut in
the pair variables. Reused pair variables receive the sum of coefficients
from all their edges.

## Measurements

Native MATLAB settings with the local license and Mosek were used. All
16 distinct available `clique_N_k...` models had the same objective under
first-level clique cuts alone and under the added count cuts. The official
60-vertex copies are byte-identical to the corresponding Frank files.

| Model | First-level cuts | With count cuts |
| --- | --- | --- |
| `clique_40_k3` | 3 nodes | 2 nodes |
| `clique_50_k3` | 3 nodes | 1 node |
| `clique_60_k9` | 15 s cap, 21 nodes | proven in 4 nodes, 4.046 s |
| `clique_60_k7`, `k8` | 30 s cap, 25/24 nodes | same cap and node counts |
| `clique_70_k3` | 5 nodes, 18.146/15.092 s | 5 nodes, 15.179/12.585 s |

In native comparisons against the previous core, `clique_50_k3_16_17`
fell from 5 nodes and 3.647/3.729 s to 1 node and 1.841/1.492 s.
`clique_60_k9_6_7` timed out at the 15 s BNB cap after 11/12 nodes in the
previous core, while the new cuts proved objective 441 in 4 nodes and
6.892/6.850 s. An already easy `clique_60_k5` case stayed at 1 node but
grew by 0.226 s in one pair; cut generation has a small overhead there.

The full `tests/bnb` suite, run twice in opposite orders against previous
core and this change, gave 22 passes and the same three existing failures
(`test_migp`, `test_map`, `test_bar`) in every run. The SDPLIB subset took
39.005/38.937 s in previous core and 39.150/37.951 s with the change.
The native production-path `test_bnb_misdp` passed. Detailed run output
remains in the local autoresearch results and experiment ledger.
