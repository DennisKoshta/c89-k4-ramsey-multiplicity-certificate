#!/usr/bin/env python3
"""FORMAT §5 test vector: build VSJ's final3840 witness in schema c89-step-graphon-v1 (block = 20), with
claim U_VSJ, and write it out for verify/c89_k4.py.

    python -I verify/vsj_testvector.py VSJ_REPO_DIR OUT.json.gz

Reads data/constructions/final3840.json.gz (lukevs/k4-ramsey @2cf216b) as data only; nothing from that
repository is imported or executed.
"""
import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import gzip
import hashlib
import json
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c89_k4 import U_VSJ  # noqa: E402
from test_c89_k4 import encode  # noqa: E402

SRC_SHA = "05302cbc635e939cc41f4ba019cdcba80199b0b83563100bac0b1a0d9fff1a29"
DENSE_SHA = "4093a1651998d485137e849795e68a3f9535b5c963fe2e3b776733f699b19a87"


def main(repo, out):
    with gzip.open(os.path.join(repo, "data/constructions/final3840.json.gz"), "rb") as fh:
        raw = fh.read()
    sha = hashlib.sha256(raw).hexdigest()
    print(f"source JSON sha256 {sha} {'(matches VSJ id)' if sha == SRC_SHA else '(MISMATCH)'}")
    if sha != SRC_SHA:
        return 1
    d = json.loads(raw)
    print(f"source schema {d['schema']!r}, den {d['edge_probability_denominator']}")
    w = d["block_weights"]
    if not (len(w) == 3840 and all(type(x) is int and x == w[0] for x in w)):
        print("block_weights are not all equal: not an equal-class step graphon")
        return 1
    den = d["edge_probability_denominator"]
    rows = d["red_probability_numerators"]
    if not (type(den) is int and len(rows) == 3840 and all(len(r) == 3840 and all(type(x) is int for x in r) for r in rows)):
        print("unexpected shape or non-integer entries")
        return 1
    M = np.array(rows, dtype=np.int64)
    dense = hashlib.sha256(M.astype("<i4").tobytes()).hexdigest()
    print(f"dense sha256 {dense} {'(matches FORMAT §5)' if dense == DENSE_SHA else '(MISMATCH)'}")
    if dense != DENSE_SHA:
        return 1
    enc = encode(M, den, 20, U_VSJ)
    enc["provenance"] = "VSJ final3840 (lukevs/k4-ramsey @2cf216b), re-encoded by verify/vsj_testvector.py"
    with gzip.open(out, "wt") as fh:
        json.dump(enc, fh)
    print(f"wrote {out}: n = 3840, block = 20, null blocks (a<b) = {len(enc['blocks'])}, claim = U_VSJ")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
