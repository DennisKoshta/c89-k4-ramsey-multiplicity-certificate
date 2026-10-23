#!/usr/bin/env python3
"""Self-tests and red controls for verify/c89_k4.py.  Run as a plain script: python -I verify/test_c89_k4.py

- green: random tiny witnesses (dense null pattern, so every one of the 11 edge-set classes occurs, with
  coarse indices repeating) scored against the literal sum (★) over all n^4 quadruples in Python integers;
- green: mid-size witnesses (n = 240, 480) against an independent dense route, S mod p via
  sum_i tr((diag(M_i) M)^3) in exact float64 BLAS (p < 2^20, n * p^2 < 2^53);
- red: invalid files and wrong claims must FAIL.
"""
import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import contextlib
import gzip
import hashlib
import io
import json
import os
import random
import sys
import tempfile
from fractions import Fraction
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c89_k4 as V  # noqa: E402

if not __debug__:
    sys.exit("refusing to run under python -O")


def encode(M, den, block, claim=None, force_null=()):
    """Witness dict for symmetric int matrix M (FORMAT §2). Blocks in force_null are stored explicitly."""
    n = M.shape[0]
    C = n // block
    coarse = [[None] * C for _ in range(C)]
    blocks = []
    for a in range(C):
        for b in range(a, C):
            B = M[a * block:(a + 1) * block, b * block:(b + 1) * block]
            if (B == B[0, 0]).all() and (a == b or (a, b) not in force_null):
                coarse[a][b] = coarse[b][a] = int(B[0, 0])
            else:
                assert a != b, "diagonal block not constant"
                blocks.append([a, b, B.tolist()])
    d = {"schema": V.SCHEMA, "constant": "C_89", "n": n, "den": den, "block": block, "coarse": coarse,
         "blocks": blocks, "sha256_dense": hashlib.sha256(M.astype("<i4").tobytes()).hexdigest()}
    if claim is not None:
        d["claim"] = f"{claim.numerator}/{claim.denominator}"
    return d


def brute_S(M, den):
    n = M.shape[0]
    L = M.tolist()
    S = 0
    for i, j, k, l in product(range(n), repeat=4):
        e = (L[i][j], L[i][k], L[i][l], L[j][k], L[j][l], L[k][l])
        r = q = 1
        for x in e:
            r *= x
            q *= den - x
        S += r + q
    return S


def dense_S_mod(M, den, p):
    def hom(A):
        A = A % p
        Af = A.astype(np.float64)
        t = 0
        for i in range(A.shape[0]):
            B = (A[i][:, None] * A) % p
            Bf = B.astype(np.float64)
            B2 = np.rint(Bf @ Bf).astype(np.int64) % p   # entries < n p^2 < 2^53: exact
            t = (t + int(((B2 * B.T) % p).sum() % p)) % p
        return t
    assert M.shape[0] * (p - 1) ** 2 < 2**53
    return (hom(M.astype(np.int64)) + hom(den - M.astype(np.int64))) % p


def random_witness(rng, C, block, den, pnull, pconst_in_null=0.15, extreme=False):
    n = C * block
    M = np.zeros((n, n), dtype=np.int64)
    force = set()
    vals = (lambda: rng.choice([0, den])) if extreme else (lambda: rng.randint(0, den))
    for a in range(C):
        for b in range(a, C):
            sl = (slice(a * block, (a + 1) * block), slice(b * block, (b + 1) * block))
            if a != b and rng.random() < pnull:
                if rng.random() < pconst_in_null:     # null entry whose block happens to be constant
                    B = np.full((block, block), vals())
                    force.add((a, b))
                else:
                    B = np.array([[rng.randint(0, den) for _ in range(block)] for _ in range(block)])
            else:
                B = np.full((block, block), vals())
            M[sl] = B
            M[sl[1], sl[0]] = B.T
    return M, force


def run_cli(d, extra=()):
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "w.json.gz")
        with gzip.open(path, "wt") as fh:
            json.dump(d, fh)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = V.main([path, *extra])
    return rc, out.getvalue()


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAILED'}] {name}" + (f"  ({detail})" if detail else ""))
    if not cond:
        check.bad += 1


check.bad = 0


def main():
    rng = random.Random(89)
    # ---- green: tiny, against (★) literally
    cases = [(4, 1, 7, 0.7), (3, 2, 5, 0.9), (4, 2, 11, 0.8), (3, 3, 4, 1.0), (5, 2, 9, 0.6), (4, 3, 3, 0.75),
             (2, 4, 6, 1.0), (6, 2, 13, 0.85), (1, 5, 8, 0.0), (5, 1, 1, 0.9)]
    for idx, (C, block, den, pnull) in enumerate(cases):
        M, force = random_witness(rng, C, block, den, pnull)
        S = brute_S(M, den)
        n = C * block
        F = Fraction(S, n**4 * den**6)
        d = encode(M, den, block, F, force)
        n_, den_, block_, Kc, nulls, _ = V.validate(json.loads(json.dumps(d)))
        S2, *_ = V.exact_S(n_, den_, block_, Kc, nulls, log=lambda *a: None)
        check(f"tiny brute force #{idx} C={C} block={block} den={den} null pairs={len(nulls) // 2}", S == S2,
              f"S={S}")
        if idx in (0, 2, 5):
            rc, out = run_cli(d)
            check(f"tiny #{idx} CLI PASS", rc == 0 and "VERDICT: PASS" in out)
    # every edge-set class contributes a nonzero term in some tiny case (so the brute-force checks exercise it)
    nonzero = set()
    for seed in range(6):
        M, force = random_witness(random.Random(seed), 7, 2, 9, 0.45 if seed % 2 else 0.9)
        _, _, _, Kc, nulls, _ = V.validate(json.loads(json.dumps(encode(M, 9, 2, None, force))))
        parts = {}
        V.hom_k4_mod(Kc, nulls, 2, 1048573, parts)
        nonzero |= {k for k, (_, v) in parts.items() if v}
    check("all 11 edge-set classes nonzero in tiny cases", len(nonzero) == 11, " ".join(sorted(nonzero)))

    # ---- green: mid-size, against the dense route mod p
    for (C, block, den, pnull, extreme) in [(12, 20, 65536, 0.35, True), (24, 20, 65536, 0.2, True),
                                             (16, 15, 1000, 0.4, False)]:
        M, force = random_witness(random.Random(C * block), C, block, den, pnull, extreme=extreme)
        n = C * block
        d = encode(M, den, block, None, force)
        _, _, _, Kc, nulls, _ = V.validate(json.loads(json.dumps(d)))
        S, *_ = V.exact_S(n, den, block, Kc, nulls, log=lambda *a: None)
        for p in (1048573, 1048571):
            check(f"mid-size n={n} block={block} den={den} null pairs={len(nulls) // 2} vs dense route mod {p}",
                  S % p == dense_S_mod(M, den, p))

    # ---- red controls on a tiny witness
    for seed in range(1000):   # want an explicit block, an off-diagonal constant and a 0/1 diagonal constant
        M, force = random_witness(random.Random(seed), 4, 2, 11, 0.5)
        g = encode(M, 11, 2)
        if (g["blocks"] and any(g["coarse"][a][b] is not None for a in range(4) for b in range(a + 1, 4))
                and any(g["coarse"][a][a] in (0, 1) for a in range(4))):
            break
    n = 8
    S = brute_S(M, 11)
    F = Fraction(S, n**4 * 11**6)
    good = encode(M, 11, 2, F, force)
    assert good["blocks"], "need at least one explicit block"

    def red(name, mutate, expect):
        d = json.loads(json.dumps(good))
        mutate(d)
        rc, out = run_cli(d)
        check(f"red: {name}", rc != 0 and "VERDICT: FAIL" in out and expect in out,
              next((l for l in out.splitlines() if l.startswith(("INVALID", "claim !="))), "?"))

    def perturb(d):
        r = d["blocks"][0][2]
        r[0][0] = r[0][0] + 1 if r[0][0] < 11 else r[0][0] - 1
    red("perturbed entry (sha mismatch)", perturb, "sha256_dense mismatch")

    def off_by_one(d):
        p, q = map(int, d["claim"].split("/"))
        w = Fraction(p + 1, q)   # off by one in p, then put back in lowest terms
        d["claim"] = f"{w.numerator}/{w.denominator}"
    red("wrong claim (p + 1, reduced)", off_by_one, "claim != F(W)")
    red("non-lowest-terms claim", lambda d: d.update(claim="/".join(str(2 * int(x)) for x in d["claim"].split("/"))),
        "not in lowest terms")

    def float_entry(d):
        d["blocks"][0][2][0][1] = float(d["blocks"][0][2][0][1])
    red("float entry in rows", float_entry, "non-integer entry")
    red("float den", lambda d: d.update(den=11.0), "den is not a JSON integer")

    def bool_coarse(d):
        a = next(i for i in range(4) if d["coarse"][i][i] in (0, 1))
        d["coarse"][a][a] = bool(d["coarse"][a][a])
    try:
        red("boolean coarse entry", bool_coarse, "neither int nor null")
    except StopIteration:
        print("[skip] boolean coarse entry: no 0/1 diagonal constant")

    def asym_coarse(d):
        a, b = next((a, b) for a in range(4) for b in range(a + 1, 4) if d["coarse"][a][b] is not None)
        d["coarse"][a][b] = (d["coarse"][a][b] + 1) % 12
    red("asymmetric coarse entry", asym_coarse, "coarse not symmetric")

    def lower_block(d):
        a, b, rows = d["blocks"][0]
        d["blocks"][0] = [b, a, [list(r) for r in zip(*rows)]]
    red("block stored as (b, a) with a < b", lower_block, "bad indices")

    def null_diag(d):
        d["coarse"][0][0] = None
        d["blocks"].append([0, 0, [[d["coarse"][1][1]] * 2] * 2])
    red("null (non-constant) diagonal block", null_diag, "diagonal coarse entry (0,0) is null")
    red("missing block item", lambda d: d["blocks"].pop(), "no blocks item")
    red("duplicate block item", lambda d: d["blocks"].append(d["blocks"][0]), "duplicate blocks item")
    red("entry above den", lambda d: d["blocks"][0][2][0].__setitem__(0, 12), "outside [0, den]")
    red("missing claim", lambda d: d.pop("claim"), "claim is not a string")

    # the literal JSON float spelling "65536.0" (not via Python float) and duplicate keys
    with tempfile.TemporaryDirectory() as td:
        txt = json.dumps(good).replace('"den": 11', '"den": 11.0')
        path = os.path.join(td, "w.json.gz")
        with gzip.open(path, "wt") as fh:
            fh.write(txt)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = V.main([path])
        check("red: literal 11.0 for den", rc != 0 and "not a JSON integer" in out.getvalue())

    # unreadable file and malformed --beat: clean FAIL / usage error, no traceback
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "w.json.gz")
        with open(path, "wb") as fh:
            fh.write(b"not gzip")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = V.main([path])
        check("red: non-gzip file", rc != 0 and "not a readable gzip file" in out.getvalue())
    for bad in ("1:2", "1/0"):
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                V.main(["/nonexistent", "--beat", bad])
            check(f"red: --beat {bad} rejected", False)
        except SystemExit as exc:
            check(f"red: --beat {bad} rejected before scoring", exc.code == 2)

    print("ALL OK" if check.bad == 0 else f"{check.bad} CHECK(S) FAILED")
    return 1 if check.bad else 0


if __name__ == "__main__":
    sys.exit(main())
