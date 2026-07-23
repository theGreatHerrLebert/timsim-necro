#!/usr/bin/env python3
"""Golden-run regression gate for the timsim v2 pipeline.

Reruns a small, frozen reference experiment (60 proteins) end to end through BOTH closed loops —
Thermo `.raw` DIA and Bruker `.d` DIA — searches each with DiaNN, scores against the render's own
answer key, and diffs the resulting metrics against a pinned baseline. The point is a NUMBER: every
change to a render/predictor/schema node lands against this gate, so a regression in detectable recall
or FDP is visible immediately rather than three features later.

Cheap because the DAG is content-addressed: the two loops share their whole feature-space (digest →
… → spectra), so it is computed once; only ccs (Bruker) + render + search + score differ.

Usage (via run.sh, which sets the venv + tool paths):
    python gate.py                 # run both loops, diff vs baseline.json, append to history.jsonl
    python gate.py --only thermo   # one instrument
    python gate.py --update-baseline   # rerun and REWRITE baseline.json (do this deliberately)

Big templates are NOT committed; point at them with env vars (defaults suit the dev box):
    GOLDEN_THERMO_TEMPLATE   a no-IMS Orbitrap/Astral .raw   (default /scratch/timsim-demo/orbi_hela_dia.raw)
    GOLDEN_BRUKER_REF        a reference DIA .d               (default .../MIDIA-250K-NOISELESS-FRAG.d)
    GOLDEN_WORKDIR           where the DAG materialises        (default /scratch/timsim-demo/golden-work)
"""
import argparse
import datetime
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLOW = HERE.parent / "flow" / "timsim_flow.py"
CONFIG = HERE / "config"
BASELINE = HERE / "baseline.json"
HISTORY = HERE / "history.jsonl"

# A metric moves by less than this (fraction) → "no change"; recall drop / FDP rise beyond it → regression.
EPS = 0.003  # 0.3 percentage points

THERMO_TEMPLATE = os.environ.get("GOLDEN_THERMO_TEMPLATE", "/scratch/timsim-demo/orbi_hela_dia.raw")
BRUKER_REF = os.environ.get("GOLDEN_BRUKER_REF",
                            "/scratch/timsim-demo/MIDIA-250K-NOISELESS-FRAG/MIDIA-250K-NOISELESS-FRAG.d")
WORKDIR = Path(os.environ.get("GOLDEN_WORKDIR", "/scratch/timsim-demo/golden-work"))


def render_proteome_spec() -> Path:
    """Materialise the proteome spec from its template, pointing at the frozen FASTA (portable — no
    absolute path baked into the committed config)."""
    tmpl = (CONFIG / "proteome.toml.tmpl").read_text()
    spec = WORKDIR / "proteome.toml"
    spec.write_text(tmpl.replace("__FASTA__", str(CONFIG / "tiny.fasta")))
    return spec


def run_pipeline(kind: str, proteome_spec: Path) -> dict:
    """Run one instrument's closed loop (render → search → score) and return its parsed metrics.json."""
    assert kind in ("thermo", "bruker")
    measurement = (["--thermo-template", THERMO_TEMPLATE] if kind == "thermo"
                   else ["--bruker-reference", BRUKER_REF])
    score_dir = "score" if kind == "thermo" else "score_bruker"
    cmd = [
        sys.executable, str(FLOW),
        "--outdir", str(WORKDIR),
        "--proteome-spec", str(proteome_spec),
        "--mods", str(CONFIG / "mods_basic.toml"),
        "--design-spec", str(CONFIG / "tiny_design.toml"),
        *measurement,
        "--frag-model", "",
        "--search-fasta", str(CONFIG / "tiny.fasta"),
        "--samples", "A_R1",
    ]
    print(f"  [{kind}] running closed loop …", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        tail = "\n".join((r.stdout + r.stderr).splitlines()[-25:])
        raise SystemExit(f"[{kind}] pipeline FAILED (exit {r.returncode}):\n{tail}")
    hits = glob.glob(str(WORKDIR / score_dir / "*" / "metrics.json"))
    if not hits:
        raise SystemExit(f"[{kind}] no metrics.json under {WORKDIR/score_dir} — did the score node run?")
    return json.loads(Path(sorted(hits)[0]).read_text())


def summarize(m: dict) -> dict:
    """The gate's headline numbers — the strictest (detectable) recall, FDP, ID counts, decile ladder."""
    det = m["hierarchy"][-1]
    return {
        "n_ids": m["diann_ids"],
        "correct": m["correct"],
        "false": m["false"],
        "fdp": round(m["fdp"], 5),
        "detectable_recall": round(det["recall"], 5),
        "detectable_denom": det["size"],
        "deciles": [round(d["recall"], 4) for d in m["recall_by_abundance_decile"]],
    }


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def fmt_delta(cur: float, base: float, higher_is_better: bool) -> str:
    d = cur - base
    if abs(d) < EPS:
        return "[no change]"
    good = (d > 0) if higher_is_better else (d < 0)
    arrow = f"{d*100:+.1f}pp"
    return f"[improved {arrow}]" if good else f"[REGRESSION {arrow}]"


def report(kind: str, cur: dict, base: dict | None) -> bool:
    """Print one instrument's line, return True if it REGRESSED vs baseline."""
    label = {"thermo": "Thermo DIA", "bruker": "Bruker DIA"}[kind]
    rc = cur["detectable_recall"]
    fdp = cur["fdp"]
    if base is None:
        print(f"  {label:11s} recall(det) {rc*100:5.1f}%   FDP {fdp*100:.2f}%   n_ids {cur['n_ids']:,}   [no baseline]")
        return False
    br, bf = base["detectable_recall"], base["fdp"]
    rflag = fmt_delta(rc, br, higher_is_better=True)
    fflag = fmt_delta(fdp, bf, higher_is_better=False)
    print(f"  {label:11s} recall(det) {br*100:5.1f}% -> {rc*100:5.1f}% {rflag:22s}"
          f"  FDP {bf*100:.2f}% -> {fdp*100:.2f}% {fflag}")
    return "REGRESSION" in rflag or "REGRESSION" in fflag


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=["thermo", "bruker"], help="run just one instrument")
    ap.add_argument("--update-baseline", action="store_true", help="rewrite baseline.json from this run")
    a = ap.parse_args()

    WORKDIR.mkdir(parents=True, exist_ok=True)
    spec = render_proteome_spec()
    kinds = [a.only] if a.only else ["thermo", "bruker"]

    results = {k: summarize(run_pipeline(k, spec)) for k in kinds}

    baseline = json.loads(BASELINE.read_text()) if BASELINE.exists() else {}
    print()
    regressed = False
    for k in kinds:
        regressed |= report(k, results[k], baseline.get(k))

    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    commit = git_commit()
    with HISTORY.open("a") as f:
        f.write(json.dumps({"ts": stamp, "commit": commit, **results}) + "\n")
    print(f"\n  logged → {HISTORY.name}  (commit {commit}, {stamp})")

    if a.update_baseline:
        merged = {**baseline, **results}  # keep the other instrument if --only was used
        BASELINE.write_text(json.dumps(merged, indent=2) + "\n")
        print(f"  baseline UPDATED → {BASELINE.name}")
    elif regressed:
        print("\n  ✗ REGRESSION vs baseline — investigate before committing.")
        sys.exit(1)
    elif baseline:
        print("\n  ✓ no regression vs baseline.")


if __name__ == "__main__":
    main()
