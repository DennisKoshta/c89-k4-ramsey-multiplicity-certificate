"""Controls for rebuild.py: the real delta rebuilds (exit 0); each tampered input is refused (nonzero exit / error).
usage: python -I test_rebuild.py /path/to/final3840.json.gz DELTA.json.gz SCRATCH_DIR"""
import os, sys, gzip, json, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rebuild


def run(base, delta_obj, scratch, name):
    p = os.path.join(scratch, f"delta_{name}.json.gz")
    with gzip.open(p, "wb") as fh:
        fh.write(json.dumps(delta_obj, separators=(",", ":")).encode())
    out = os.path.join(scratch, name)
    os.makedirs(out, exist_ok=True)
    try:
        rc = rebuild.main(base, p, out)
    except SystemExit as e:
        rc = f"refused: {e.code}"
    except (ValueError, AssertionError) as e:
        rc = f"refused: {type(e).__name__}: {e}"
    return rc


def main(base, delta_path, scratch):
    d = json.loads(gzip.open(delta_path, "rb").read())
    ok = True
    rc = run(base, d, scratch, "original")
    print(f"{'ok ' if rc == 0 else 'BAD'} original -> {rc}", flush=True)
    ok &= rc == 0

    def ctl(name, mut, want=None):
        """want=None: any refusal or nonzero exit; want=1: rebuilt, but the hash check must fail (exit 1)."""
        nonlocal ok
        e = copy.deepcopy(d)
        mut(e)
        rc = run(base, e, scratch, name.replace(" ", "_"))
        good = rc != 0 if want is None else rc == want
        ok &= good
        print(f"{'ok ' if good else 'BAD'} {name} -> {rc}", flush=True)

    def shrink(e):                                           # within the bound: only the hashes can catch it
        rows = e["amplitudes"][0][2]
        i, j = next((i, j) for i in range(len(rows)) for j in range(len(rows)) if rows[i][j] != 0)
        rows[i][j] -= 1 if rows[i][j] > 0 else -1

    ctl("amplitude moved 1 toward 0", shrink, want=1)
    ctl("extra overwrites sha256_dense",
        lambda e: e["witness"]["extra"].__setitem__("sha256_dense", e["witness"]["sha256_dense"]))
    ctl("file_name with a path", lambda e: e["witness"].__setitem__("file_name", "../escape.json.gz"))
    ctl("block dropped", lambda e: e["amplitudes"].pop(), want=1)
    ctl("duplicate block", lambda e: e["amplitudes"].append(copy.deepcopy(e["amplitudes"][0])))
    ctl("amplitude beyond bound", lambda e: e["amplitudes"][0][2][0].__setitem__(0, 70000))
    ctl("float amplitude", lambda e: e["amplitudes"][0][2][0].__setitem__(0, 1.0))
    ctl("other base commit", lambda e: e["base"].__setitem__("commit", "0" * 40))
    ctl("claim changed", lambda e: e["witness"].__setitem__("claim", "1/33"), want=1)
    print("ALL CONTROLS BEHAVE" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
