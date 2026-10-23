"""Rebuild our C89 witness (schema c89-step-graphon-v1, FORMAT.md) from VSJ's published table plus our delta.

The delta (schema "c89-sign-delta-v1", gzipped JSON) holds only our own data, never VSJ's table:
  base         pinned identity of the base table: lukevs/k4-ramsey @ commit, path, SHA-256 of its decompressed JSON,
               n = 3840, den = 65536, and the SHA-256 of the dense table (little-endian int32, row-major)
  base_block   20: the amplitudes are stored by groups of 20 consecutive base classes
  amplitudes   [[a, b, rows], ...] for group pairs a < b with a nonzero amplitude block; rows is 20 x 20 ints.
               A is symmetric (block (b, a) is the transpose) and 0 everywhere else, including the diagonal groups.
  witness      what the rebuilt witness must be: block, sha256_dense, claim, extra keys (non-input, e.g. provenance),
               json_sha256 (decompressed JSON), file_sha256, file_name and gzip_mtime (for a byte-identical .gz)
The witness is the balanced sign refinement of VSJ's table P (VSJ paper, Prop. "sign refinement identity"):
  M'[2i+s, 2j+t] = P_ij + (-1)^(s+t) A_ij,   i, j < 3840,  s, t in {0, 1},   |A_ij| <= min(P_ij, den - P_ij),
so n' = 7680 classes with the same den.  Exit 0 iff the rebuilt witness has the recorded sha256_dense and JSON SHA-256.
Needs numpy only.  The VSJ file is read as data; nothing from that repository is run.
usage: python -I rebuild.py /path/to/k4-ramsey/data/constructions/final3840.json.gz DELTA.json.gz OUT_DIR"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import sys, io, gzip, json, hashlib
import numpy as np

BASE = {"repo": "https://github.com/lukevs/k4-ramsey", "commit": "2cf216ba885fbb1b2d1b2d4920efeaa6ba438842",
        "path": "data/constructions/final3840.json.gz",
        "json_sha256": "05302cbc635e939cc41f4ba019cdcba80199b0b83563100bac0b1a0d9fff1a29",
        "dense_sha256": "4093a1651998d485137e849795e68a3f9535b5c963fe2e3b776733f699b19a87",
        "n": 3840, "den": 65536}
BASE_BLOCK = 20


def dense_sha(M):
    return hashlib.sha256(np.ascontiguousarray(M, dtype="<i4").tobytes()).hexdigest()


def load_base(path):
    """VSJ's table as an int64 matrix, after checking it is exactly the pinned file."""
    raw = gzip.open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != BASE["json_sha256"]:
        sys.exit(f"base: decompressed JSON sha256 {sha} is not {BASE['json_sha256']} (VSJ final3840 @2cf216b)")
    d = json.loads(raw)
    w = d["block_weights"]
    if not (d["schema"] == "rational-step-graphon-v1" and d["edge_probability_denominator"] == BASE["den"]
            and len(w) == BASE["n"] and all(type(x) is int and x == w[0] for x in w)):
        sys.exit("base: unexpected schema, denominator or class weights")
    P = np.array(d["red_probability_numerators"], dtype=np.int64)
    if not (P.shape == (BASE["n"],) * 2 and (P == P.T).all() and P.min() >= 0 and P.max() <= BASE["den"]
            and dense_sha(P) == BASE["dense_sha256"]):
        sys.exit("base: table shape, symmetry, range or dense sha256 mismatch")
    return P


def amplitudes(delta, N):
    """Dense symmetric int64 amplitude matrix from the delta, with structural checks."""
    f = delta["base_block"]
    if not (delta["schema"] == "c89-sign-delta-v1" and delta["base"] == BASE and f == BASE_BLOCK):
        raise ValueError("delta: wrong schema or base")
    C = N // f
    A = np.zeros((C, C, f, f), np.int64)
    seen = set()
    for item in delta["amplitudes"]:
        if not (isinstance(item, list) and len(item) == 3):
            raise ValueError("delta: amplitude items must be [a, b, rows]")
        a, b, rows = item
        if not (type(a) is int and type(b) is int and 0 <= a < b < C) or (a, b) in seen:
            raise ValueError(f"delta: bad or duplicate group pair ({a}, {b})")
        if not (isinstance(rows, list) and len(rows) == f and all(isinstance(r, list) and len(r) == f and
                                                                   all(type(v) is int for v in r) for r in rows)):
            raise ValueError(f"delta: block ({a}, {b}) is not {f} x {f} ints")
        x = np.array(rows, dtype=np.int64)
        A[a, b], A[b, a] = x, x.T
        seen.add((a, b))
    return A.transpose(0, 2, 1, 3).reshape(N, N)


def encode(M, block, den, claim, extra):
    """The c89-step-graphon-v1 encoding (FORMAT.md section 2), keys in this fixed order."""
    n = M.shape[0]
    C = n // block
    B = M.reshape(C, block, C, block).transpose(0, 2, 1, 3)
    const = (B == B[:, :, :1, :1]).all(axis=(2, 3))
    assert const.diagonal().all(), "non-constant diagonal coarse block"
    coarse = [[int(B[a, b, 0, 0]) if const[a, b] else None for b in range(C)] for a in range(C)]
    blocks = [[a, b, B[a, b].astype(int).tolist()] for a in range(C) for b in range(a + 1, C) if not const[a, b]]
    d = {"schema": "c89-step-graphon-v1", "constant": "C_89", "n": n, "den": den, "block": block,
         "coarse": coarse, "blocks": blocks, "sha256_dense": dense_sha(M), "claim": claim}
    clash = set(extra) & set(d)
    if clash:
        raise ValueError(f"delta: extra keys {sorted(clash)} would overwrite witness fields")
    d.update(extra)
    return d


def main(base_path, delta_path, out_dir):
    P = load_base(base_path)
    N, den = BASE["n"], BASE["den"]
    print(f"base OK: {BASE['repo']} @{BASE['commit'][:7]} {BASE['path']}, JSON sha256 {BASE['json_sha256'][:16]}..., "
          f"dense sha256 {BASE['dense_sha256'][:16]}...", flush=True)
    delta = json.loads(gzip.open(delta_path, "rb").read())
    A = amplitudes(delta, N)
    if not (np.abs(A) <= np.minimum(P, den - P)).all():
        sys.exit("delta: an amplitude exceeds min(P, den - P)")
    wit = delta["witness"]
    sgn = np.array([[1, -1], [-1, 1]], np.int64)
    M = (P[:, None, :, None] + A[:, None, :, None] * sgn[None, :, None, :]).reshape(2 * N, 2 * N)
    del P, A
    assert M.min() >= 0 and M.max() <= den and (M == M.T).all()
    name = wit["file_name"]
    if not (isinstance(name, str) and name == os.path.basename(name) and name not in ("", ".", "..")
            and name.endswith(".json.gz")):
        raise ValueError("delta: file_name must be a plain *.json.gz file name")
    dense = dense_sha(M)                                     # from the rebuilt matrix itself, not from the delta
    d = encode(M, wit["block"], den, wit["claim"], wit["extra"])
    assert d["sha256_dense"] == dense
    data = json.dumps(d, separators=(",", ":")).encode()
    out = os.path.join(out_dir, name)
    # written as the record was (gzip.open(..., "wt", compresslevel=9) + json.dump): deflate output depends on how
    # the input arrives, so this streaming is what makes the .gz byte-identical, not only the JSON
    with open(out, "wb") as raw, gzip.GzipFile(filename=name, mode="wb", compresslevel=9,
                                               fileobj=raw, mtime=wit["gzip_mtime"]) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8") as fh:
        json.dump(d, fh, separators=(",", ":"))
    dense_ok = dense == wit["sha256_dense"]
    json_sha = hashlib.sha256(data).hexdigest()
    file_sha = hashlib.sha256(open(out, "rb").read()).hexdigest()
    print(f"wrote {out}: n {d['n']}, den {den}, block {d['block']}, {len(d['blocks'])} explicit blocks")
    print(f"sha256_dense   {dense}  {'MATCH' if dense_ok else 'MISMATCH'}")
    print(f"JSON sha256    {json_sha}  {'MATCH' if json_sha == wit['json_sha256'] else 'MISMATCH'}")
    print(f"file sha256    {file_sha}  {'MATCH' if file_sha == wit['file_sha256'] else 'differs (gzip/zlib build only; the JSON is what counts)'}")
    print(f"claim          {d['claim']}")
    return 0 if dense_ok and json_sha == wit["json_sha256"] else 1


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:4]))
