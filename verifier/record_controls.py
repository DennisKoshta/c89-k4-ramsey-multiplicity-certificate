#!/usr/bin/env python3
"""Red controls on a real witness file: write mutated copies of it, for verify/c89_k4.py to reject.

    python -I verify/record_controls.py WITNESS.json.gz OUTDIR

Each copy changes exactly one thing.  `wrongclaim` (p + 1 over q, reduced) is valid in every other way, so
the verifier has to score it in full to reject it; the others fail validation or the claim parse.
"""
import gzip
import json
import os
import sys
from fractions import Fraction


def main(src, outdir):
    with gzip.open(src, "rt") as fh:
        base = fh.read()

    def write(name, d):
        path = os.path.join(outdir, f"ctrl_{name}.json.gz")
        with gzip.open(path, "wt") as fh:
            json.dump(d, fh)
        print(path)

    def fresh():
        return json.loads(base)

    d = fresh()
    r = d["blocks"][0][2]
    r[0][0] += 1 if r[0][0] < d["den"] else -1
    write("perturbed", d)

    d = fresh()
    p, q = map(int, d["claim"].split("/"))
    w = Fraction(p + 1, q)
    d["claim"] = f"{w.numerator}/{w.denominator}"
    write("wrongclaim", d)

    d = fresh()
    d["claim"] = f"{2 * p}/{2 * q}"
    write("nonlowest", d)

    d = fresh()
    d["blocks"][0][2][0][1] = float(d["blocks"][0][2][0][1])
    write("floatentry", d)

    d = fresh()
    C = len(d["coarse"])
    a, b = next((a, b) for a in range(C) for b in range(a + 1, C) if d["coarse"][a][b] is not None)
    d["coarse"][a][b] = d["den"] - d["coarse"][a][b]
    write("asymcoarse", d)

    d = fresh()
    a, b, rows = d["blocks"][0]
    d["blocks"][0] = [b, a, [list(x) for x in zip(*rows)]]   # stored as the (b, a) block
    write("lowerblock", d)
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
