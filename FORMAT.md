# C89 witness: format and exact scoring rule

Status: **FROZEN** on 2026-10-23, before any C89 claim; the freeze commit is the commit that introduced this line
(`git log -S'FROZEN on 2026-10-23' -- problems/c89/FORMAT.md`). Any later change needs a new commit, a dated note in
this file, and a message to the lead (5697e7f3). An independent verifier should be written from this file alone; it does
not need `k4eval.py` or any other code in this repository.

## 1. The constant and what a witness proves

C89 is the Ramsey multiplicity constant of K4 (teorth/optimizationproblems, constants/89a.md):

  C89 = lim_{n→∞} min k4 / C(n,4),

where the minimum is over all red/blue colourings of the edges of K_n, and k4 is the number of 4-sets whose six
edges all have the same colour.

An **equal-class step graphon** is a symmetric n×n matrix W = (w_ij) with entries in [0,1]. Each class i gets weight
1/n. The entry w_ij is the probability that an edge between a vertex of class i and a vertex of class j is red. The
diagonal entry w_ii is the probability for an edge between two distinct vertices of the same class.

Its **monochromatic K4 density** is

  F(W) = n^(-4) · Σ_{i1,i2,i3,i4 ∈ {0..n-1}} ( Π_{1≤s<t≤4} w_{i_s i_t} + Π_{1≤s<t≤4} (1 − w_{i_s i_t}) ).   (★)

The sum runs over **all n^4 ordered quadruples, including those with repeated indices**. A repeated index uses the
diagonal entry, e.g. i1 = i2 contributes the factor w_{i1 i1}. This is t(K4, W) + t(K4, 1 − W), and it is exactly
VSJ's definition (Van Seters–Jothimani, paper eq. (density)).

**Lemma (finite realization; VSJ Lemma "Finite realization").** For every n' ≥ 4, min k4 / C(n',4) ≤ F(W), hence
C89 ≤ F(W).

*Proof.* Give each vertex of K_{n'} an independent uniform class in {0..n−1}. Colour each edge red independently with
probability w_ij, where i and j are the classes of its ends. For four distinct vertices, the four classes are
independent and uniform (repetitions allowed), and the six edges are independent given the classes. So the
probability that they span a monochromatic K4 is exactly (★). By linearity, the expected proportion of monochromatic
4-sets is F(W), so some colouring does at least as well. ∎

So a witness proves **C89 ≤ F(W)**. Nothing else is claimed.

## 2. File format (schema `c89-step-graphon-v1`)

A witness file is gzip-compressed UTF-8 JSON. W = M/den, where M is a symmetric n×n integer matrix with entries in
[0, den]. The n classes are grouped, for compression only, into C = n/block consecutive **coarse groups** of `block`
classes each: class x lies in coarse group ⌊x/block⌋. M is stored by **coarse blocks** M[a·block.., b·block..]
(block×block).

| key | type | meaning |
|---|---|---|
| `schema` | string | `"c89-step-graphon-v1"` |
| `constant` | string | `"C_89"` |
| `n` | int > 0 | number of classes |
| `den` | int, 0 < den < 2^31 | common denominator (the bound keeps `sha256_dense`'s int32 encoding exact) |
| `block` | int > 0, divides n | coarse group size |
| `coarse` | C×C list | entry (a,b): an int v with 0 ≤ v ≤ den if every entry of coarse block (a,b) equals v; `null` otherwise |
| `blocks` | list of `[a, b, rows]` | one per `null` coarse entry with a < b; `rows` is a block×block list of ints, rows indexed by the class in group a, columns by the class in group b |
| `sha256_dense` | hex string | SHA-256 of the dense M as little-endian int32, row-major (n² values) |
| `claim` | string `"p/q"` | the exact F(W), in lowest terms, q > 0 (required in a record file) |

Any other key (e.g. `provenance`) is **not an input** and must be ignored.

Validity rules. A verifier must reject the file if any of these fails:

1. All of `n`, `den`, `block`, the coarse ints, and every entry of every `rows` are JSON integers. Floats and
   booleans are not allowed.
2. `coarse` is symmetric: (a,b) and (b,a) are both `null` or the same int. Every diagonal entry (a,a) is an int,
   i.e. diagonal coarse blocks are constant. That constant is also the fine diagonal w_ii for those classes.
3. `blocks` has exactly one item for each a < b with `coarse[a][b] = null`, and no other items. Block (b,a) is the
   transpose of block (a,b).
4. Every entry of the assembled M is in [0, den], and M is symmetric.
5. `sha256_dense` equals the SHA-256 of the assembled M (little-endian int32, row-major).

The format does not require constant row sums inside blocks, or any other structure. Such structure may be used to
compute (★) faster, but only after checking it.

## 3. Exact scoring rule

With integer entries, (★) becomes

  F(W) = S / (n^4 · den^6),  S = Σ_{i1..i4} ( Π_{s<t} M_{i_s i_t} + Π_{s<t} (den − M_{i_s i_t}) ),

summed over all n^4 ordered quadruples, repeated indices included, with M_ii from the matrix. S is a non-negative
integer with S ≤ 2·n^4·den^6.

The verifier computes F(W) **exactly**, as an integer S or a rational, and compares it with `claim` as rationals.

- Floating point may guide a computation, but it must not decide the verdict.
- Modular arithmetic with CRT is fine if the product of the moduli exceeds the bound on S, and every intermediate
  result is provably exact. For example, float64 matmuls of residues are exact only if k·(p−1)² < 2^53, where k is
  the inner length.

**Verdict PASS** iff the file is valid (§2), `claim` is present as the string "p/q" with p/q in lowest terms, and the
exact F(W) equals p/q.

## 4. The claim line

For a PASS file:

  **C89 ≤ claim = p/q ≈ 0.0301388… (10 significant digits, rounded up)**.

A claim **beats** the current record iff p/q < U_VSJ, where

  U_VSJ = 8450462766487926638466333426306607129 / 280384030360880691940646801777885184000 = 0.030138887566497220…

(VSJ 2026, lukevs/k4-ramsey @2cf216b, Theorem "Main bound"). The margin U_VSJ − p/q is reported from exact rationals.

## 5. Test vector

VSJ's witness itself, built by the verifier from their repository. We do not redistribute their files.

- Source: `data/constructions/final3840.json.gz` at lukevs/k4-ramsey @2cf216b, with n = 3840 and den = 65536. VSJ
  give the SHA-256 identifier 05302cbc635e939cc41f4ba019cdcba80199b0b83563100bac0b1a0d9fff1a29 for it.
- In the class order of that file, the dense SHA-256 of §2 is
  4093a1651998d485137e849795e68a3f9535b5c963fe2e3b776733f699b19a87.
- Encoded with block = 20, it must score exactly to U_VSJ.

## 6. Notes for an independent verifier (non-normative)

- The literal sum has n^4 ≈ 3.5·10^15 terms for n = 7680, so it is not feasible directly. Some valid routes:
  - Group the sum by coarse groups (a,b,c,d), with block-matrix products in exact or modular arithmetic.
  - Use the structure of a sign refinement. Our witnesses (n classes) have the form
    M'[2i+s, 2j+t] = M_ij + (−1)^(s+t)·A_ij for an N×N base table M with **N = n/2** (VSJ's, N = 3840), s, t ∈ {0,1},
    and an integer amplitude A with |A_ij| ≤ min(M_ij, den − M_ij). Then F(W') = F(W) + T3(A) + T4(A) (VSJ, Prop.
    "Sign refinement identity"), where W = M/den is the **N-class base** graphon and, with all sums over 0..N−1,
    - T3 = 4/N^4 · Σ_{ijk} a_ij a_jk a_ki (Σ_l w_il w_jl w_kl − q_il q_jl q_kl),
    - T4 = 3/N^4 · Σ_{ijkl} a_ij a_jk a_kl a_li (w_ik w_jl + q_ik q_jl),
    - where a = A/den and q = 1 − w. (The normaliser is the base size N, not the witness size n.)
  - The author's own route is different: a centred expansion around block row sums, with CRT over 10 primes below
    2^22.
- The witnesses we produce have constant row and column sums inside every coarse block, and constant diagonal
  coarse blocks equal to 0 (all within-class edges blue).
