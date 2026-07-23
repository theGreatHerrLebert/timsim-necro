#!/usr/bin/env python3
"""Golden-run regression gate for the timsim v2 pipeline — two axes against one fixed oracle.

The render's answer key is a fixed ground truth. Holding one side still and varying the other measures
two different things, so the gate has two modes:

  SIM AXIS  (default)   fix the tool (DiaNN), vary the render/predictors/schema  → simulation realism
  TOOL AXIS (--tool-axis)  fix a FROZEN rendered dataset, vary the search engine → tool regression + benchmark

Sim axis reruns the whole pipeline end to end (both Thermo `.raw` and Bruker `.d` DIA loops), from a CLEAN
work dir — necroflow content-addresses on inputs/config/command, not the binary, so only a clean run
reflects current render/predictor code. The two loops still share their feature space within that run.

Tool axis runs only search→score against a dataset you froze with `--freeze` (a render you trust), so it
is cheap (~40s, no re-simulation) and isolates the SOFTWARE: DiaNN-vs-Sage-vs-FragPipe, or your own
tool's version bump, on identical ground truth.

Usage (via run.sh, which sets venv + tool paths):
    python gate.py                       # sim axis: both loops, diff baseline.sim_axis, log history
    python gate.py --only bruker         # sim axis, one instrument
    python gate.py --update-baseline     # rerun and REWRITE the relevant baseline section
    python gate.py --freeze              # snapshot the current renders as the frozen tool-axis dataset
    python gate.py --tool-axis --tool diann          # tool axis on the frozen dataset (both instruments)
    python gate.py --tool-axis --tool diann --only bruker

Env (defaults suit the dev box):
    GOLDEN_THERMO_TEMPLATE  no-IMS Orbitrap/Astral .raw   GOLDEN_BRUKER_REF  reference DIA .d
    GOLDEN_WORKDIR          sim-axis DAG scratch          GOLDEN_DATASET     frozen tool-axis dataset
    TIMSIM_DIANN            DiaNN binary                  DOTNET_ROOT        .NET (Thermo .raw only)
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLOW = HERE.parent / "flow" / "timsim_flow.py"
CONFIG = HERE / "config"
BASELINE = HERE / "baseline.json"
HISTORY = HERE / "history.jsonl"
DATASET_MANIFEST = HERE / "dataset.json"

EPS = 0.003  # 0.3 pp: smaller move = "no change"; recall drop / FDP rise beyond it = regression

THERMO_TEMPLATE = os.environ.get("GOLDEN_THERMO_TEMPLATE", "/scratch/timsim-demo/orbi_hela_dia.raw")
BRUKER_REF = os.environ.get("GOLDEN_BRUKER_REF",
                            "/scratch/timsim-demo/MIDIA-250K-NOISELESS-FRAG/MIDIA-250K-NOISELESS-FRAG.d")
WORKDIR = Path(os.environ.get("GOLDEN_WORKDIR", "/scratch/timsim-demo/golden-work"))
DATASET = Path(os.environ.get("GOLDEN_DATASET", "/scratch/timsim-demo/golden-dataset"))
TOOLWORK = Path(os.environ.get("GOLDEN_TOOLWORK", "/scratch/timsim-demo/golden-toolwork"))
DIANN = os.environ.get("TIMSIM_DIANN", "/home/administrator/dia-nn/diann-2.5.0/diann-linux")
DOTNET = os.environ.get("DOTNET_ROOT", os.path.expanduser("~/.dotnet"))

# DiaNN library-free params — MUST match the DAG's search / search_bruker nodes so a tool-axis DiaNN run
# on a frozen render reproduces the sim-axis number (self-consistency), and so tool comparisons are fair.
DIANN_ARGS = [
    "--fasta-search", "--predictor", "--gen-spec-lib", "--threads", "16",
    "--met-excision", "--cut", "K*,R*", "--missed-cleavages", "2",
    "--min-pep-len", "7", "--max-pep-len", "30", "--var-mods", "1", "--unimod35",
]
QVALUE = 0.01
LABEL = {"thermo": "Thermo DIA", "bruker": "Bruker DIA"}


# ── shared helpers ────────────────────────────────────────────────────────────

def render_proteome_spec(dest: Path) -> Path:
    """Materialise the proteome spec from its template → the frozen FASTA (portable: no absolute path
    baked into the committed config)."""
    dest.mkdir(parents=True, exist_ok=True)
    spec = dest / "proteome.toml"
    spec.write_text((CONFIG / "proteome.toml.tmpl").read_text().replace("__FASTA__", str(CONFIG / "tiny.fasta")))
    return spec


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()[:16]


def summarize(m: dict) -> dict:
    det = m["hierarchy"][-1]
    return {
        "n_ids": m["diann_ids"], "correct": m["correct"], "false": m["false"],
        "fdp": round(m["fdp"], 5),
        "detectable_recall": round(det["recall"], 5), "detectable_denom": det["size"],
        "deciles": [round(d["recall"], 4) for d in m["recall_by_abundance_decile"]],
    }


def score(report: Path, truth: Path, peptides: Path, out: Path) -> dict:
    """Run the instrument-agnostic timsim_eval scorer → parsed metrics."""
    cmd = [sys.executable, "-m", "timsim_eval.v2_thermo_eval",
           "--report", str(report), "--truth", str(truth), "--peptides", str(peptides),
           "--fdr", str(QVALUE), "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("score FAILED:\n" + "\n".join((r.stdout + r.stderr).splitlines()[-20:]))
    return json.loads(out.read_text())


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
    return f"[improved {d*100:+.1f}pp]" if good else f"[REGRESSION {d*100:+.1f}pp]"


def report_line(label: str, cur: dict, base: dict | None) -> bool:
    rc, fdp = cur["detectable_recall"], cur["fdp"]
    if base is None:
        print(f"  {label:11s} recall(det) {rc*100:5.1f}%   FDP {fdp*100:.2f}%   n_ids {cur['n_ids']:,}   [no baseline]")
        return False
    rflag = fmt_delta(rc, base["detectable_recall"], True)
    fflag = fmt_delta(fdp, base["fdp"], False)
    print(f"  {label:11s} recall(det) {base['detectable_recall']*100:5.1f}% -> {rc*100:5.1f}% {rflag:20s}"
          f"  FDP {base['fdp']*100:.2f}% -> {fdp*100:.2f}% {fflag}")
    return "REGRESSION" in rflag or "REGRESSION" in fflag


def load_baseline() -> dict:
    return json.loads(BASELINE.read_text()) if BASELINE.exists() else {}


def log_history(kind: str, results: dict) -> str:
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    commit = git_commit()
    with HISTORY.open("a") as f:
        f.write(json.dumps({"ts": stamp, "commit": commit, "axis": kind, **results}) + "\n")
    print(f"\n  logged → {HISTORY.name}  (commit {commit}, {stamp})")
    return commit


# ── SIM AXIS: vary the render, fix DiaNN ──────────────────────────────────────

def run_pipeline(kind: str, spec: Path) -> dict:
    measurement = (["--thermo-template", THERMO_TEMPLATE] if kind == "thermo"
                   else ["--bruker-reference", BRUKER_REF])
    score_dir = "score" if kind == "thermo" else "score_bruker"
    cmd = [sys.executable, str(FLOW), "--outdir", str(WORKDIR),
           "--proteome-spec", str(spec), "--mods", str(CONFIG / "mods_basic.toml"),
           "--design-spec", str(CONFIG / "tiny_design.toml"), *measurement,
           "--frag-model", "", "--search-fasta", str(CONFIG / "tiny.fasta"), "--samples", "A_R1"]
    print(f"  [{kind}] running closed loop …", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"[{kind}] pipeline FAILED (exit {r.returncode}):\n"
                         + "\n".join((r.stdout + r.stderr).splitlines()[-25:]))
    hits = glob.glob(str(WORKDIR / score_dir / "*" / "metrics.json"))
    if not hits:
        raise SystemExit(f"[{kind}] no metrics.json under {WORKDIR/score_dir}")
    return json.loads(Path(sorted(hits)[0]).read_text())


def sim_axis(kinds: list[str], update: bool, clean: bool) -> None:
    if clean:
        shutil.rmtree(WORKDIR, ignore_errors=True)  # necroflow caches on inputs, NOT the binary — start clean
    spec = render_proteome_spec(WORKDIR)
    results = {k: summarize(run_pipeline(k, spec)) for k in kinds}

    base = load_baseline()
    sim_base = base.get("sim_axis", {})
    print()
    regressed = any(report_line(LABEL[k], results[k], sim_base.get(k)) for k in kinds)
    log_history("sim", results)

    if update:
        base.setdefault("sim_axis", {}).update(results)
        BASELINE.write_text(json.dumps(base, indent=2) + "\n")
        print(f"  baseline.sim_axis UPDATED → {BASELINE.name}")
    elif regressed:
        sys.exit("\n  ✗ REGRESSION vs baseline — investigate before committing.")
    elif sim_base:
        print("\n  ✓ no regression vs baseline.")


# ── FREEZE: snapshot a trustworthy render as the tool-axis fixture ─────────────

def freeze() -> None:
    """Copy the current sim-axis renders + truths + shared peptides to the pinned GOLDEN_DATASET and
    record a small committed manifest (paths + content hashes). The tool axis reads THIS, never the
    work dir, so it stays fixed while render code moves."""
    peps = sorted(glob.glob(str(WORKDIR / "digest" / "*" / "peptides.parquet")))
    if not peps:
        raise SystemExit(f"no peptides in {WORKDIR} — run the sim axis first, then --freeze")
    peptides_src = Path(peps[0])
    sources = {
        "thermo": (glob.glob(str(WORKDIR / "render_thermo" / "*" / "data.raw")),
                   glob.glob(str(WORKDIR / "render_thermo" / "*" / "truth.parquet")), False),
        "bruker": (glob.glob(str(WORKDIR / "render" / "*" / "data.d")),
                   glob.glob(str(WORKDIR / "render" / "*" / "truth.parquet")), True),
    }
    manifest = {"frozen_at": datetime.datetime.now().isoformat(timespec="seconds"),
                "render_commit": git_commit(), "instruments": {}}
    for inst, (render_hits, truth_hits, is_dir) in sources.items():
        if not render_hits or not truth_hits:
            print(f"  [{inst}] no render in work dir — skipping")
            continue
        dest = DATASET / inst
        shutil.rmtree(dest, ignore_errors=True)
        dest.mkdir(parents=True, exist_ok=True)
        render_src = Path(render_hits[0])
        render_dst = dest / render_src.name
        (shutil.copytree if is_dir else shutil.copy2)(render_src, render_dst)
        shutil.copy2(truth_hits[0], dest / "truth.parquet")
        shutil.copy2(peptides_src, dest / "peptides.parquet")
        manifest["instruments"][inst] = {
            "render": str(render_dst), "is_dir": is_dir,
            "truth": str(dest / "truth.parquet"), "peptides": str(dest / "peptides.parquet"),
            "fasta": str(CONFIG / "tiny.fasta"),
            "truth_sha": sha256(dest / "truth.parquet"), "peptides_sha": sha256(dest / "peptides.parquet"),
        }
        print(f"  [{inst}] frozen → {dest}  (truth {manifest['instruments'][inst]['truth_sha']})")
    DATASET_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"  manifest → {DATASET_MANIFEST.name}  (render_commit {manifest['render_commit']})")


# ── TOOL AXIS: fix the frozen dataset, vary the search engine ─────────────────

def run_diann(inst: str, ds: dict, outdir: Path) -> Path:
    """DiaNN library-free on the frozen render. Bruker `.d` is native; a Thermo `.raw` needs .NET."""
    outdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    if not ds["is_dir"]:  # Thermo .raw → .NET runtime
        env["DOTNET_ROOT"] = DOTNET
        env["PATH"] = f"{DOTNET}:{env['PATH']}"
    cmd = [DIANN, "--f", ds["render"], "--fasta", ds["fasta"],
           "--out", str(outdir / "report.parquet"), "--qvalue", str(QVALUE), *DIANN_ARGS]
    log = outdir / "diann.log"
    with log.open("w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=lf, env=env)
    if r.returncode != 0:
        raise SystemExit(f"[{inst}] DiaNN FAILED (exit {r.returncode}) — see {log}")
    return outdir / "report.parquet"


TOOLS = {"diann": run_diann}  # sage / fragpipe slot in here via timsim_eval.{sage,fragpipe}_executor


def tool_axis(tool: str, insts: list[str], update: bool) -> None:
    if not DATASET_MANIFEST.exists():
        raise SystemExit("no frozen dataset — run `gate.py --freeze` after a sim-axis run first.")
    if tool not in TOOLS:
        raise SystemExit(f"tool '{tool}' not wired yet (have: {', '.join(TOOLS)}). "
                         f"sage/fragpipe use timsim_eval.{tool}_executor — add a runner to TOOLS.")
    manifest = json.loads(DATASET_MANIFEST.read_text())["instruments"]
    runner = TOOLS[tool]
    results = {}
    for inst in insts:
        if inst not in manifest:
            print(f"  [{inst}] not in frozen dataset — skipping")
            continue
        ds = manifest[inst]
        print(f"  [{inst}] {tool} on frozen render …", flush=True)
        wd = TOOLWORK / inst / tool
        report = runner(inst, ds, wd)
        m = score(report, Path(ds["truth"]), Path(ds["peptides"]), wd / "metrics.json")
        results[inst] = summarize(m)

    base = load_baseline()
    tool_base = base.get("tool_axis", {})
    print()
    regressed = False
    for inst in results:
        b = tool_base.get(inst, {}).get(tool)
        regressed |= report_line(f"{LABEL[inst].split()[0]}/{tool}", results[inst], b)
    log_history(f"tool:{tool}", results)

    if update:
        ta = base.setdefault("tool_axis", {})
        for inst, summ in results.items():
            ta.setdefault(inst, {})[tool] = summ
        BASELINE.write_text(json.dumps(base, indent=2) + "\n")
        print(f"  baseline.tool_axis[*][{tool}] UPDATED → {BASELINE.name}")
    elif regressed:
        sys.exit("\n  ✗ REGRESSION vs baseline — investigate before committing.")
    elif any(tool_base.get(i, {}).get(tool) for i in results):
        print("\n  ✓ no regression vs baseline.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=["thermo", "bruker"], help="one instrument")
    ap.add_argument("--update-baseline", action="store_true", help="rewrite the relevant baseline section")
    ap.add_argument("--freeze", action="store_true", help="snapshot current renders as the tool-axis dataset")
    ap.add_argument("--tool-axis", action="store_true", help="run the tool axis on the frozen dataset")
    ap.add_argument("--tool", default="diann", help="search engine for the tool axis (default: diann)")
    ap.add_argument("--no-clean", action="store_true", help="sim axis: reuse the work dir (faster; NOT a true reading)")
    a = ap.parse_args()

    insts = [a.only] if a.only else ["thermo", "bruker"]
    if a.freeze:
        freeze()
    elif a.tool_axis:
        tool_axis(a.tool, insts, a.update_baseline)
    else:
        sim_axis(insts, a.update_baseline, clean=not a.no_clean)


if __name__ == "__main__":
    main()
