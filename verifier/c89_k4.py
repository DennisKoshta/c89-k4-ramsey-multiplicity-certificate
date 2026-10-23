#!/usr/bin/env python3
"""Independent verifier for C89 (K4 Ramsey multiplicity) step-graphon witnesses.

Written from problems/c89/FORMAT.md (frozen spec, schema c89-step-graphon-v1) alone.

    python -I verify/c89_k4.py WITNESS.json.gz [--beat P/Q]

Route (deliberately not the author's centred expansion).  Let g(i) be the coarse group of class i,
H the set of coarse pairs (a,b) whose block is stored explicitly (`null` in `coarse`), and write

    M_ij = K_{g(i)g(j)} + N_ij,

where K is the C x C matrix of coarse constants with K_ab = 0 on H, and N is M on the H blocks and 0
elsewhere.  For every pair exactly one of the two parts is present, so expanding the product over the
six edges of K4 gives

    hom(K4, M) = sum over edge sets F of K4 of  sum_{i1..i4} prod_{e in F} N_e * prod_{e not in F} K_e.

The 64 edge sets fall into 11 isomorphism classes; each class is summed once and weighted by its number
of labellings.  An edge in F must map to an H pair, and a vertex with no F edge contributes a factor
`block` for its fine index.  The same routine is run on M and on den - M (constant blocks den - K, H blocks
den - N), and S = hom(K4, M) + hom(K4, den - M).

Arithmetic: everything is done modulo primes p < 2^25 in int64.  Operands are reduced into [0, p) before
every product, so elementwise products are < 2^50.  Every matrix product goes through `mm`, which splits
the inner dimension into chunks of length k with k*(p-1)^2 < 2^63, so each partial dot product is exact.
S is reconstructed by CRT; the product of the primes exceeds 2*n^4*den^6 >= S >= 0, which pins down S.
Floats are used only to print approximate decimals.
"""
import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"  # one core; set before numpy is imported

import argparse
import gzip
import hashlib
import json
import math
import re
import sys
import time
from fractions import Fraction
from itertools import combinations

import numpy as np

if not __debug__:
    sys.exit("refusing to run under python -O")

SCHEMA = "c89-step-graphon-v1"
U_VSJ = Fraction(8450462766487926638466333426306607129, 280384030360880691940646801777885184000)
I64_LIMIT = 2**63 - 1


class Invalid(Exception):
    pass


# ----------------------------------------------------------------------------------------------------
# Validation (FORMAT §2)

def _is_int(x):
    return type(x) is int  # excludes bool, float


def _no_dup_keys(pairs):
    keys = [k for k, _ in pairs]
    if len(keys) != len(set(keys)):
        raise Invalid(f"duplicate JSON key among {sorted(keys)}")
    return dict(pairs)


def load(path):
    try:
        with gzip.open(path, "rb") as fh:
            raw = fh.read()
    except (OSError, EOFError) as exc:
        raise Invalid(f"not a readable gzip file: {exc}")
    try:
        return json.loads(raw.decode("utf-8"), object_pairs_hook=_no_dup_keys)
    except ValueError as exc:
        raise Invalid(f"not valid UTF-8 JSON: {exc}")


def validate(d):
    """Return (n, den, block, Kc, nulls, M) or raise Invalid.

    Kc: C x C int64 coarse constants (0 where null).  nulls: dict (a, b) -> block x block int64 for every
    ordered pair with a null coarse entry (both (a,b) and (b,a)).  M: dense int32 n x n."""
    if not isinstance(d, dict):
        raise Invalid("top level is not an object")
    for key in ("schema", "constant", "n", "den", "block", "coarse", "blocks", "sha256_dense"):
        if key not in d:
            raise Invalid(f"missing key {key!r}")
    if d["schema"] != SCHEMA:
        raise Invalid(f"schema {d['schema']!r}")
    if d["constant"] != "C_89":
        raise Invalid(f"constant {d['constant']!r}")
    n, den, block = d["n"], d["den"], d["block"]
    for name, v in (("n", n), ("den", den), ("block", block)):
        if not _is_int(v):
            raise Invalid(f"{name} is not a JSON integer: {v!r}")
    if n <= 0:
        raise Invalid("n <= 0")
    if not 0 < den < 2**31:
        raise Invalid("den out of (0, 2^31)")
    if block <= 0 or n % block:
        raise Invalid("block does not divide n")
    C = n // block
    coarse = d["coarse"]
    if not (isinstance(coarse, list) and len(coarse) == C and all(isinstance(r, list) and len(r) == C for r in coarse)):
        raise Invalid("coarse is not a C x C list")
    Kc = np.zeros((C, C), dtype=np.int64)
    isnull = np.zeros((C, C), dtype=bool)
    for a in range(C):
        for b in range(C):
            v = coarse[a][b]
            if v is None:
                isnull[a, b] = True
            elif _is_int(v):
                if not 0 <= v <= den:
                    raise Invalid(f"coarse[{a}][{b}] = {v} outside [0, den]")
                Kc[a, b] = v
            else:
                raise Invalid(f"coarse[{a}][{b}] is neither int nor null: {v!r}")
    for a in range(C):
        if isnull[a, a]:
            raise Invalid(f"diagonal coarse entry ({a},{a}) is null")
        for b in range(a + 1, C):
            if isnull[a, b] != isnull[b, a] or (not isnull[a, b] and coarse[a][b] != coarse[b][a]):
                raise Invalid(f"coarse not symmetric at ({a},{b})")
    blocks = d["blocks"]
    if not isinstance(blocks, list):
        raise Invalid("blocks is not a list")
    nulls = {}
    for item in blocks:
        if not (isinstance(item, list) and len(item) == 3):
            raise Invalid("blocks item is not [a, b, rows]")
        a, b, rows = item
        if not (_is_int(a) and _is_int(b) and 0 <= a < b < C):
            raise Invalid(f"blocks item has bad indices {a!r}, {b!r}")
        if not isnull[a, b]:
            raise Invalid(f"blocks item ({a},{b}) for a non-null coarse entry")
        if (a, b) in nulls:
            raise Invalid(f"duplicate blocks item ({a},{b})")
        if not (isinstance(rows, list) and len(rows) == block
                and all(isinstance(r, list) and len(r) == block for r in rows)):
            raise Invalid(f"blocks item ({a},{b}) is not block x block")
        for r in rows:
            for x in r:
                if not _is_int(x):
                    raise Invalid(f"blocks item ({a},{b}) has a non-integer entry {x!r}")
                if not 0 <= x <= den:
                    raise Invalid(f"blocks item ({a},{b}) has entry {x} outside [0, den]")
        B = np.array(rows, dtype=np.int64)
        nulls[(a, b)] = B
        nulls[(b, a)] = B.T.copy()
    for a in range(C):
        for b in range(a + 1, C):
            if isnull[a, b] and (a, b) not in nulls:
                raise Invalid(f"no blocks item for null coarse entry ({a},{b})")
    # rule 4/5: assemble the dense M
    M = np.empty((n, n), dtype=np.int32)
    for a in range(C):
        for b in range(C):
            sl = (slice(a * block, (a + 1) * block), slice(b * block, (b + 1) * block))
            M[sl] = nulls[(a, b)] if isnull[a, b] else Kc[a, b]
    if M.min() < 0 or M.max() > den:
        raise Invalid("M has entries outside [0, den]")
    if not np.array_equal(M, M.T):
        raise Invalid("M not symmetric")
    sha = hashlib.sha256(M.astype("<i4").tobytes()).hexdigest()
    if not isinstance(d["sha256_dense"], str) or sha != d["sha256_dense"].lower():
        raise Invalid(f"sha256_dense mismatch: file {d['sha256_dense']!r}, assembled {sha}")
    return n, den, block, Kc, nulls, M


def parse_claim(s):
    if not isinstance(s, str):
        raise Invalid("claim is not a string")
    m = re.fullmatch(r"([0-9]+)/([0-9]+)", s)
    if not m:
        raise Invalid(f"claim {s!r} is not 'p/q'")
    p, q = int(m.group(1)), int(m.group(2))
    if q == 0:
        raise Invalid("claim has q = 0")
    if math.gcd(p, q) != 1:
        raise Invalid(f"claim {s} is not in lowest terms (gcd {math.gcd(p, q)})")
    return Fraction(p, q)


# ----------------------------------------------------------------------------------------------------
# Exact modular arithmetic

def _is_prime(m):
    if m < 2:
        return False
    for f in range(2, math.isqrt(m) + 1):
        if m % f == 0:
            return False
    return True


def primes_below(limit, count):
    out, m = [], limit - 1
    while len(out) < count:
        if _is_prime(m):
            out.append(m)
        m -= 1
    return out


def mm(A, B, p):
    """Exact (A @ B) mod p for int64 arrays with entries in [0, p)."""
    k = A.shape[-1]
    step = I64_LIMIT // ((p - 1) ** 2)
    assert step >= 1
    if k <= step:
        return (A @ B) % p
    out = None
    for s in range(0, k, step):
        part = (A[..., s:s + step] @ B[..., s:s + step, :]) % p
        out = part if out is None else (out + part) % p
    return out


def hom_k4_mod(Kc, nulls, block, p, parts=None):
    """hom(K4, M) mod p, M given by coarse constants Kc (0 on null pairs) and the null blocks."""
    C = Kc.shape[0]
    b = block
    K = Kc % p
    nb = [sorted(c for (a2, c) in nulls if a2 == a) for a in range(C)]
    adj = {a: set(nb[a]) for a in range(C)}
    N = {key: v % p for key, v in nulls.items()}
    R = {key: v.sum(axis=1) % p for key, v in nulls.items()}  # r(a->c)[i] = sum_{j in c} N_ij, i in a
    S = np.zeros((C, C), dtype=np.int64)
    for (a, c), v in nulls.items():
        S[a, c] = int(v.sum()) % p
    RS = {a: np.stack([R[(a, c)] for c in nb[a]]) for a in range(C) if nb[a]}   # rows r(a->c), c in nb[a]
    bp = b % p
    tot = {}

    # T0: F empty.  b^4 * sum_a prod_6 K
    t = 0
    for a1 in range(C):
        Y = (K[a1][None, :] * K) % p                 # Y[a2, c] = K_{a1 c} K_{a2 c}
        Q = (mm(Y, K, p) * Y % p).sum(axis=1) % p    # Q[a2] = sum_{c,d} Y[a2,c] K_cd Y[a2,d]
        t = (t + int((K[a1] * Q % p).sum() % p)) % p
    tot["T0"] = (1, t * pow(bp, 4, p) % p)

    pairs = [(a, c) for a in range(C) for c in nb[a]]
    if pairs:
        PA = np.array(pairs, dtype=np.int64)
        Z = (K[PA[:, 0]] * K[PA[:, 1]]) % p          # z[c] = K_{a1 c} K_{a2 c}
        Sv = S[PA[:, 0], PA[:, 1]]
        qK = (mm(Z, K, p) * Z % p).sum(axis=1) % p
        qS = (mm(Z, S, p) * Z % p).sum(axis=1) % p
        t1 = int((Sv * qK % p).sum() % p)
        t2b = int((Sv * qS % p).sum() % p)
    else:
        t1 = t2b = 0
    # T1: F = {12}.  b^2 sum S12 K13 K14 K23 K24 K34
    tot["T1"] = (6, t1 * bp % p * bp % p)
    # T2b: F = {12, 34}.  sum S12 S34 K13 K14 K23 K24
    tot["T2b"] = (3, t2b)

    t2a = t3a = t3b = t4b = 0
    for a2 in range(C):
        L = nb[a2]
        if not L:
            continue
        RL = RS[a2]
        KL = K[L]
        c3 = mm((KL * K[a2][None, :]) % p, KL.T.copy(), p)          # c3[x,z] = sum_a4 K_{x a4} K_{a2 a4} K_{z a4}
        # T2a: F = {12, 23} with centre a2.  b * sum P(a1,a2,a3) K13 c3(a1,a2,a3)
        P = mm(RL, RL.T.copy(), p)
        t2a = (t2a + int((P * K[np.ix_(L, L)] % p * c3 % p).sum() % p)) % p
        # T3b: star centred at a2 (rename: centre = a2, leaves L).  sum_i R_x R_y R_z * K_xy K_xz K_yz
        RR = (RL[:, None, :] * RL[None, :, :]) % p                  # [x, y, i]
        T = mm(RR.reshape(-1, b), RL.T.copy(), p).reshape(len(L), len(L), len(L))
        KLL = K[np.ix_(L, L)]
        W = (KLL[:, :, None] * KLL[:, None, :]) % p * KLL[None, :, :] % p
        t3b = (t3b + int((T * W % p).sum() % p)) % p
        # T3a (triangle a1 a2 a3 with a2 = this vertex as vertex 1) and T4b (triangle + pendant at a3)
        for xi, a1 in enumerate(L):
            for zi, a3 in enumerate(L):
                if a3 not in adj[a1]:
                    continue
                X = mm(N[(a2, a1)], N[(a1, a3)], p)                 # (N_{a2a1} N_{a1a3})[i2, i3]
                tr = int((X * N[(a3, a2)].T % p).sum() % p)         # tr(N21 N13 N32)
                t3a = (t3a + tr * int(c3[xi, zi])) % p
                # T4b: triangle (a3; a2, a1) with pendant a4 at a3: sum_{i} diag_{i3}(N32 N21 N13) r(a3->a4)[i3] K_{a2 a4} K_{a1 a4}
                dg = (N[(a3, a2)] * X.T % p).sum(axis=1) % p        # dg[i3] = sum_{i2} N_{i3 i2} X[i2, i3]
                L4 = nb[a3]
                v = mm(RS[a3], dg[:, None], p)[:, 0]
                w = K[a2, L4] * K[a1, L4] % p
                t4b = (t4b + int((v * w % p).sum() % p)) % p
    tot["T2a"] = (12, t2a * bp % p)
    tot["T3a"] = (4, t3a * bp % p)
    tot["T3b"] = (4, t3b)
    tot["T4b"] = (12, t4b)

    # T3c: path 1-2-3-4, F = {12, 23, 34}.  sum r(a2->a1)^T N23 r(a3->a4) * K13 K14 K24
    t3c = 0
    for (a2, a3) in pairs:
        L1, L4 = nb[a2], nb[a3]
        G = mm(mm(RS[a2], N[(a2, a3)], p), RS[a3].T.copy(), p)               # G[x, y] = r(a2->x)^T N23 r(a3->y)
        Wt = (K[L1, a3][:, None] * K[np.ix_(L1, L4)]) % p * K[a2, L4][None, :] % p
        t3c = (t3c + int((G * Wt % p).sum() % p)) % p
    tot["T3c"] = (12, t3c)

    # T4a: 4-cycle 1-2-3-4-1.  sum tr(N12 N23 N34 N41) K13 K24
    # T5:  K4 minus {24}.  sum_{i1,i3} N13 (N12 N23)_{13} (N14 N43)_{13} K24
    t4a = t5 = 0
    for a1 in range(C):
        for a3 in range(C):
            cn = sorted(adj[a1] & adj[a3])
            if not cn:
                continue
            Xs = np.stack([mm(N[(a1, c)], N[(c, a3)], p) for c in cn])   # [c, i1, i3]
            Kcc = K[np.ix_(cn, cn)]
            flat = Xs.reshape(len(cn), -1)
            if K[a1, a3]:
                G = mm(flat, flat.T.copy(), p)                       # tr(A_c A_c'^T) = tr(N12 N23 N34 N41)
                t4a = (t4a + int(K[a1, a3]) * int((G * Kcc % p).sum() % p)) % p
            if a3 in adj[a1]:
                E = (flat * N[(a1, a3)].reshape(1, -1)) % p
                G = mm(E, flat.T.copy(), p)
                t5 = (t5 + int((G * Kcc % p).sum() % p)) % p
    tot["T4a"] = (3, t4a)
    tot["T5"] = (6, t5)

    # T6: all six edges in H.  sum_{i1,i2} N12 * x^T N34 y,  x = N13[i1]∘N23[i2],  y = N14[i1]∘N24[i2]
    t6 = 0
    for (a1, a2) in pairs:
        cn = sorted(adj[a1] & adj[a2])
        if len(cn) < 2:
            continue
        N12 = N[(a1, a2)].reshape(-1)
        Xc = {c: ((N[(a1, c)][:, None, :] * N[(a2, c)][None, :, :]) % p).reshape(b * b, b) for c in cn}
        for a3 in cn:
            for a4 in cn:
                if a4 not in adj[a3]:
                    continue
                v = (mm(Xc[a3], N[(a3, a4)], p) * Xc[a4] % p).sum(axis=1) % p
                t6 = (t6 + int((v * N12 % p).sum() % p)) % p
    tot["T6"] = (1, t6)

    if parts is not None:
        parts.update(tot)
    hom = 0
    for mult, val in tot.values():
        hom = (hom + mult * int(val)) % p
    return int(hom)


def exact_S(n, den, block, Kc, nulls, log=print):
    """S = sum over all n^4 quadruples of prod M + prod (den - M), computed by CRT."""
    bound = 2 * n**4 * den**6           # FORMAT §3: 0 <= S <= 2 n^4 den^6
    primes, prod = [], 1
    cand = primes_below(2**25, 64)
    for p in cand:
        primes.append(p)
        prod *= p
        if prod > bound:
            break
    assert prod > bound, "not enough primes"
    isnull = np.zeros_like(Kc, dtype=bool)
    for (a, c) in nulls:
        isnull[a, c] = True
    Kd = np.where(isnull, 0, den - Kc).astype(np.int64)
    Kc0 = np.where(isnull, 0, Kc).astype(np.int64)
    nulls_d = {key: den - v for key, v in nulls.items()}
    residues = []
    for p in primes:
        t = time.time()
        r = (hom_k4_mod(Kc0, nulls, block, p) + hom_k4_mod(Kd, nulls_d, block, p)) % p
        residues.append(r)
        log(f"  p = {p}: S mod p = {r}  ({time.time() - t:.1f} s)")
    # CRT (Garner-free, plain)
    S = 0
    for p, r in zip(primes, residues):
        m = prod // p
        S = (S + r * m * pow(m, -1, p)) % prod
    assert 0 <= S <= bound
    return S, primes, prod, bound


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("witness")
    ap.add_argument("--beat", help="record to compare with, as P/Q (e.g. U_VSJ)")
    ap.add_argument("--no-claim", action="store_true", help="score only; do not require a claim (test vector)")
    args = ap.parse_args(argv)
    rec = None
    if args.beat:
        m = re.fullmatch(r"\s*([0-9]+)\s*/\s*([0-9]+)\s*", args.beat)
        if args.beat == "U_VSJ":
            rec = U_VSJ
        elif m and int(m.group(2)) > 0:
            rec = Fraction(int(m.group(1)), int(m.group(2)))
        else:
            ap.error(f"--beat must be U_VSJ or P/Q with Q > 0, got {args.beat!r}")
    t0 = time.time()
    try:
        d = load(args.witness)
        n, den, block, Kc, nulls, M = validate(d)
        del M
        claim = None if args.no_claim and "claim" not in d else parse_claim(d.get("claim"))
    except Invalid as exc:
        print(f"INVALID: {exc}")
        print("VERDICT: FAIL")
        return 1
    print(f"valid: n = {n}, den = {den}, block = {block}, C = {n // block}, null blocks (a<b) = {len(nulls) // 2}"
          f"  ({time.time() - t0:.1f} s)")
    S, primes, prod, bound = exact_S(n, den, block, Kc, nulls)
    F = Fraction(S, n**4 * den**6)
    print(f"CRT over {len(primes)} primes < 2^25; product = 2^{math.log2(prod):.2f} > bound 2 n^4 den^6 = 2^{math.log2(bound):.2f}")
    print(f"S = {S}")
    print(f"F(W) = {F.numerator}/{F.denominator}")
    print(f"     ~ {F.numerator * 10**22 // F.denominator / 1e22:.20f}")
    ok = True
    if claim is not None:
        print(f"claim = {claim.numerator}/{claim.denominator}")
        if claim == F:
            print("claim == F(W) exactly")
        else:
            print(f"claim != F(W): claim - F = {claim - F}")
            ok = False
    if rec is not None:
        val = claim if claim is not None else F
        margin = rec - val
        print(f"record = {rec.numerator}/{rec.denominator}")
        print(f"margin = record - value = {margin.numerator}/{margin.denominator}")
        print(f"       ~ {float(margin):.6e}")
        print("BEATS" if (ok and val < rec) else "DOES NOT BEAT")
    print(f"VERDICT: {'PASS' if ok and claim is not None else ('SCORED (no claim)' if ok else 'FAIL')}"
          f"  ({time.time() - t0:.1f} s total)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
