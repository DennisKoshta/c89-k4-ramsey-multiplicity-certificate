# C89 (K4 Ramsey multiplicity, upper bound): independent verification

Verifier: thread 93377e06. Search author: thread adac2d8f. Date: 2026-10-23.

I wrote `verify/c89_k4.py` from `FORMAT.md` alone: the frozen spec, git f27e15c, sha256 47e6e8d1…. I didn't read any of
the author's code (`problems/c89/*.py`) or logs. The §6 sign-refinement hint isn't used. The verifier assumes no
structure beyond §2.

## Verdicts

| run | result | wall time (1 core, nice 19) | peak RSS |
|---|---|---|---|
| `verify/test_c89_k4.py` (tiny brute force, mid-size dense cross-check, 17 red controls) | ALL OK | 24 s | n/a |
| Test vector, FORMAT §5: VSJ final3840 re-encoded with block = 20, claim U_VSJ | **PASS**, F(W) = U_VSJ exactly | 103 s | 236 MB |
| Record `c89_sign3_crude_0.03013888570.json.gz` (sha256 e37d811d…), `--beat U_VSJ` | **PASS, BEATS** | 371 s | 831 MB |
| Red controls on the record file (6 mutated copies) | all **FAIL** as required | 389 s (one full rescore) | n/a |

Record, computed rather than copied from the file:

- n = 7680, den = 65536, block = 40, C = 192, 1248 explicit blocks (13 per coarse row).
- S = 8307142402237973815060289107446603942297600
- F(W) = 1126728298914655736634696330762614467/37384537381450758925419573570384691200 ≈ 0.030138885695392054. This
  equals the claim exactly, and the claim is in lowest terms.
- **Margin U_VSJ − F(W) = 1049256017227412221891173997253/560768060721761383881293603555770368000 ≈ 1.871105e-9**,
  computed in exact rationals by the verifier.

Test vector: the source file `data/constructions/final3840.json.gz` (lukevs/k4-ramsey @2cf216b) decompresses to
JSON with sha256 05302cbc… (VSJ's identifier). All 3840 `block_weights` are equal, and the dense sha256 is 4093a165….
I read it as data with `python -I`, by `verify/vsj_testvector.py`, and executed none of its code. The score is
S = 519196432373018212667371525712277942005760, and F = U_VSJ exactly.

Red controls on the record file, each changing exactly one thing (`verify/record_controls.py`):

| control | verifier output |
|---|---|
| one entry of block (0,62) +1 | `sha256_dense mismatch` → FAIL |
| claim 2p/2q | `not in lowest terms (gcd 2)` → FAIL |
| one entry written as 65536.0 | `non-integer entry 65536.0` → FAIL |
| coarse (0,1) changed to den − v on one side | `coarse not symmetric at (0,1)` → FAIL |
| block (0,62) stored as [62, 0, transpose] | `bad indices 62, 0` → FAIL |
| claim (p+1)/q, reduced | scored in full; `claim − F = 1/q` → FAIL, DOES NOT BEAT |

Tiny controls in `test_c89_k4.py` also cover these cases:
- a boolean coarse entry, a float `den` (both as a Python float and as the literal `11.0`), and an entry above `den`;
- a null diagonal block, a missing block item, a duplicate block item, and a missing claim;
- a non-gzip file, and a malformed `--beat` (`1:2`, `1/0`), which is rejected before scoring.

## Method

M_ij = K_{g(i)g(j)} + N_ij. Here K is the C×C matrix of coarse constants, set to 0 on the null pairs H. N is the
explicit blocks, set to 0 elsewhere. Exactly one of the two terms is present for each pair. So expanding Π over the 6
K4 edges of (K + N) gives an exact sum over edge sets F ⊆ E(K4): the edges in F use N, the rest use K, and a vertex
with no F edge contributes `block`. Up to relabelling there are 11 classes:

empty ×1, edge ×6, path P3 ×12, two disjoint edges ×3, triangle ×4, star ×4, path P4 ×12, C4 ×3, triangle+pendant ×12,
K4−e ×6, K4 ×1.

Each class is a coarse contraction of K, the block sums, the block row sums and the block products, summed over
labelled walks in H. Coarse indices may repeat; a repeated index is a diagonal constant. The same routine runs on
den − M, and S = hom(K4, M) + hom(K4, den − M). This is not the author's centred expansion. It also doesn't use the
row-sum structure or the sign-refinement structure.

Exactness:
- Everything runs in int64 modulo the 6 largest primes below 2^25. Their product is 2^150.00, above the §3 bound
  2·n^4·den^6 = 2^148.63 for the record (2^144.6 for VSJ).
- Operands are reduced into [0,p) before each product.
- Every matmul splits its inner dimension so that k·(p−1)² < 2^63.
- The CRT runs in Python integers. Floats only print decimals.
- The script refuses `python -O`, and BLAS is pinned to 1 thread.

Cross-checks of the implementation:
1. 10 random tiny witnesses (n ≤ 12, including block = 1 and null-but-constant blocks) match the literal sum (★) over
   all n^4 quadruples in Python integers. Every one of the 11 classes is nonzero in at least one of them.
2. Three mid-size witnesses (n = 240 and 480, block 15 and 20, den 1000 and 65536) agree mod two primes with a dense
   route, Σ_i tr((diag(M_i)·M)^3), computed in exact float64 BLAS (n·p² < 2^53).
3. The VSJ test vector reproduces U_VSJ exactly.

## Spec notes (ambiguities I hit; none affects the verdict)

1. **A null entry whose block is actually constant.** §2 says a coarse entry is an int "if every entry … equals v;
   `null` otherwise", but the 5 must-reject rules don't forbid a null over a constant block. The verifier accepts it,
   and scores it like any other explicit block (the tiny tests include such blocks).
2. **Rule 4's symmetry is implied by rules 2–3.** The diagonal blocks are constant and (b,a) is defined as the
   transpose. So an "asymmetric block" can only show up as an asymmetric coarse entry or a block item with a > b.
   Both are rejected. The verifier still checks M == Mᵀ on the assembled matrix.
3. **Duplicate JSON keys.** The spec is silent. Python's default keeps the last one, so the verifier rejects them to
   be strict.
4. **`claim` syntax.** "p/q" is read as `[0-9]+/[0-9]+`, so signs and spaces are rejected and leading zeros are
   allowed. The `sha256_dense` hex string is compared case-insensitively. `schema` and `constant` must equal the
   values in the §2 table.
5. **The bound in §3, S ≤ 2·n^4·den^6, is loose by a factor of 2**, since Π w + Π(1−w) ≤ 1. It is used as stated,
   which only means more primes than necessary.

## Reproduce

```
python -I verify/test_c89_k4.py
python -I verify/vsj_testvector.py <lukevs/k4-ramsey@2cf216b checkout> vsj3840_b20.json.gz
python -I verify/c89_k4.py vsj3840_b20.json.gz                              # PASS (claim U_VSJ == F)
python -I verify/c89_k4.py c89_sign3_crude_0.03013888570.json.gz --beat U_VSJ   # PASS, BEATS
python -I verify/record_controls.py c89_sign3_crude_0.03013888570.json.gz OUTDIR   # then each ctrl_*: FAIL
```

Logs are in `/home/dennis/workspace/runs/record-hunt/c89_verify/`: `vsj_testvector.log`, `record.log` and
`controls.log`.
