# C89: Ramsey multiplicity constant of K₄, upper bound 0.03013888570

**Claim.** C₈₉ ≤ F = 1126728298914655736634696330762614467 / 37384537381450758925419573570384691200
= 0.030138885695392052467…, so C₈₉ ≤ 0.03013888570 (rounded up).

C₈₉ = lim_{n→∞} min k₄ / C(n,4). Here k₄ counts the 4-sets of a red/blue colouring of K_n whose six edges have the
same colour, and the minimum is over all colourings (constant 89 of
[teorth/optimizationproblems](https://github.com/teorth/optimizationproblems), page added in
[PR #255](https://github.com/teorth/optimizationproblems/pull/255)).

The previous bound is Van Seters–Jothimani [VSJ2026]:
U_VSJ = 8450462766487926638466333426306607129 / 280384030360880691940646801777885184000 = 0.030138887566497221…

The margin, computed in exact rationals, is

U_VSJ − F = 1049256017227412221891173997253 / 560768060721761383881293603555770368000 ≈ 1.8711051694·10⁻⁹
(relative 6.2·10⁻⁸).

The witness is an equal-weight step graphon on 7680 classes with entries in (1/65536)ℤ. F is its monochromatic K₄
density, computed exactly. By the finite-realization lemma (`FORMAT.md` §1, a first-moment argument), every such
graphon W gives C₈₉ ≤ F(W). The witness is a finite computation, so nothing here is asymptotic or conditional.

## What is new

[VSJ2026] obtain their bound from the 768-vertex Cayley graph of Parczyk–Pokutta–Spiegel–Szabó [PPSS2024]. They
apply two balanced sign refinements, ending with a 3840-class step graphon W = P/65536. A sign refinement doubles
every class i into (i, +) and (i, −) and sets

  W'((i,s),(j,t)) = w_ij + s·t·a_ij,  |a_ij| ≤ min(w_ij, 1 − w_ij).

VSJ prove (Prop. "Sign refinement identity") that F(W') − F(W) = T₃(a) + T₄(a). T₃ is a cubic form over triangles
and T₄ a quartic form over 4-cycles; the linear and quadratic terms cancel. So the first and second order
conditions at W say nothing about such a step.

We apply one more refinement to VSJ's final 3840-class table:

- The amplitudes are a_ij = ±min(w_ij, 1 − w_ij) on 124740 of the 214100 fractional pairs, and 0 elsewhere. Sign
  ascent sets a_ij = −sign(∂T₃/∂a_ij)·min(w_ij, 1 − w_ij), so a pair gets 0 where that derivative is exactly 0.
- The sign pattern is the best of eight runs of sign ascent on T₃, run over the 1152 triangles of fractional
  coarse blocks.
- T₃ ≈ −3.710·10⁻⁹ and T₄ ≈ +1.839·10⁻⁹, so F(W') − F(W) ≈ −1.871·10⁻⁹.
- The entries stay in (1/65536)ℤ.

The exact count above was then done by two separately written programs (see Verification).

## Contents

This repository contains **none of VSJ's table**. Their repository has no licence file, so you fetch their table
yourself; `rebuild/rebuild.py` then applies our amplitudes to it.

| path | what |
|---|---|
| `delta/c89_sign3_crude_0.03013888570.delta.json.gz` | our data only: the integer amplitudes A = 65536·a (1248 blocks of 20×20), the pinned identity of VSJ's base table, and the hashes and claim of the witness they produce |
| `rebuild/rebuild.py` | checks your copy of VSJ's `final3840.json.gz`, then writes our 7680-class witness `c89_sign3_crude_0.03013888570.json.gz`; the format is documented in the docstring |
| `rebuild/test_rebuild.py` | the real delta rebuilds; 9 tampered deltas are refused or fail the hash check |
| `FORMAT.md` | witness format (schema `c89-step-graphon-v1`) and exact scoring rule, frozen before any claim (SHA-256 `47e6e8d1…`) |
| `verifier/c89_k4.py` | independent exact verifier, written from `FORMAT.md` alone (see Verification) |
| `verifier/test_c89_k4.py` | its tests: brute force on tiny witnesses, a dense cross-check, rejection controls |
| `verifier/vsj_testvector.py` | builds FORMAT §5's test vector (VSJ's own table in our format) from your VSJ checkout |
| `verifier/record_controls.py` | writes six mutated copies of the witness; each must FAIL |
| `verifier/REVIEW.md` | the verifier's report |
| `expected/` | the outputs we obtained |
| `requirements.txt` | numpy (Python 3.12) |

`verifier/REVIEW.md` and `FORMAT.md` are reproduced verbatim. In them, `verify/…` is `verifier/…` here,
`problems/c89/FORMAT.md` is `FORMAT.md`, and "the lead" refers to our internal process.

## Reproduce

Use Python 3.12 with `pip install -r requirements.txt`. Each step is single-threaded.

```
sha256sum -c SHA256SUMS
git clone https://github.com/lukevs/k4-ramsey && git -C k4-ramsey checkout 2cf216ba885fbb1b2d1b2d4920efeaa6ba438842
VSJ=k4-ramsey/data/constructions/final3840.json.gz

# 1. rebuild the witness (about 10 s, 1.3 GB): all three hashes MATCH (expected/rebuild.txt)
python -I rebuild/rebuild.py $VSJ delta/c89_sign3_crude_0.03013888570.delta.json.gz .
#    file sha256 e37d811da4e8c261f30046c5453bb9c63fd9e9700989a33fccc948074c2b2a47
#    sha256_dense c7f0e7cb32258dcdf54003f8d214fcf16d401831f841d56366da8ece3611802e
#    (a different zlib may change the .gz bytes; the JSON sha256 and sha256_dense are what count)

# 2. exact score (about 6 min, 0.9 GB): claim == F(W) exactly, BEATS, VERDICT: PASS (expected/verifier_record.txt)
python -I verifier/c89_k4.py c89_sign3_crude_0.03013888570.json.gz --beat U_VSJ

# 3. test vector, FORMAT §5: VSJ's table scores exactly U_VSJ (about 2 min) (expected/verifier_testvector.txt)
python -I verifier/vsj_testvector.py k4-ramsey vsj3840_b20.json.gz
python -I verifier/c89_k4.py vsj3840_b20.json.gz

# 4. controls
python -I verifier/test_c89_k4.py
python -I rebuild/test_rebuild.py $VSJ delta/c89_sign3_crude_0.03013888570.delta.json.gz scratch/
python -I verifier/record_controls.py c89_sign3_crude_0.03013888570.json.gz controls/   # then each controls/ctrl_*: FAIL
```

## Verification

- **Independent verifier.** `verifier/c89_k4.py` was written by a separate Claude agent from `FORMAT.md` alone,
  without access to our search or scoring code.
  - Its route: expand hom(K₄, M) over the 64 edge subsets of K₄ (11 isomorphism classes) into coarse constants and
    explicit blocks, in int64 modulo six primes below 2²⁵, then CRT.
  - It uses neither the sign-refinement structure nor the row-sum structure.
  - It reproduces U_VSJ exactly on VSJ's table, and on our witness it gives exactly the claimed F, with
    S = 8307142402237973815060289107446603942297600 = F·7680⁴·65536⁶.
- **Our evaluator.** It uses a different route, a centred expansion around block row sums with CRT over ten primes
  below 2²², and gives the same exact fraction.

## Credits

- **The construction we refine:**
  - L. Van Seters and M. Jothimani [VSJ2026],
    *Twelve Clebsch graphs and a new upper bound for the Ramsey multiplicity of K₄*,
    [github.com/lukevs/k4-ramsey](https://github.com/lukevs/k4-ramsey) at commit `2cf216b`. This covers the
    3840-class table, the balanced sign refinement and its identity, and the exact-value conventions.
  - O. Parczyk, S. Pokutta, C. Spiegel and T. Szabó [PPSS2024], *New Ramsey multiplicity bounds and search
    heuristics*, Found. Comput. Math. (2024), [doi:10.1007/s10208-024-09675-6](https://doi.org/10.1007/s10208-024-09675-6).
    This is the 768-vertex Cayley graph that VSJ start from.

  Our contribution is the third sign refinement (the amplitudes in `delta/`), its exact evaluation, the format, and
  the checking.
- **The C89 page** of teorth/optimizationproblems, by JeetuSK0808 on GitHub, is [PR #255](https://github.com/teorth/optimizationproblems/pull/255).
- **AI use.** This bound was found and prepared with Claude (Anthropic, Claude Opus 5.5), directed by Dennis Koshta,
  who takes responsibility for it. That covers the amplitude search, the exact evaluator, the format, the rebuild
  tooling and the text.
  - The witness was checked by a separate Claude agent with its own verifier (`verifier/c89_k4.py`), written from
    the format specification without access to our code.
  - No human mathematician has checked the argument.

## Licence

MIT; see `LICENSE`. Copyright (c) 2026 Dennis Koshta. The licence covers the files in this repository. It does not
cover VSJ's table, which is not included and must be obtained from their repository.
