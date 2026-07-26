#!/usr/bin/env python3
"""P1.4 — HeLa complexity ramp: empirical FDP vs peptide density, factorial over noise conditions.

For each (density level = --max-peptides, noise condition), run the necroflow Bruker-DIA pipeline (render →
DiaNN → score), then aggregate each run's `metrics.json` + `truth.parquet` into the empirical-FDP-vs-density
curve. See RAMP.md.

Key design (codex-reviewed):
- Empirical FDP (a run's realized false-discovery proportion vs the answer key), plotted against the NOMINAL
  DiaNN q-value — the calibration gap IS the "FDR inflation".
- BOTH views: fdp_raw = (false + background_subtracted)/(diann_ids + background_subtracted) is the PRIMARY
  end-to-end result; fdp_sub = false/diann_ids (the reported, background-adjusted) is the secondary
  attribution. (A2's noise-only control subtracts real blank IDs; showing only that hides part of the
  inflation A2 causes.)
- Factorial: noiseless / A1 / A1+A2 — the A1-only vs A1+A2 gap attributes the chemical-background effect.
- Clean sweep: --max-peptides is a seeded-shuffle prefix (NESTED across levels for a fixed seed); the design's
  per-protein load keeps per-peptide abundance stable, so only density varies.
- Density axis = actual rendered precursor count + detectable count + a co-elution proxy (precursors/RT-bin).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

# noise condition -> extra flow flags
CONDITIONS = {
    "noiseless": [],
    "A1": ["--noise-mz-ppm", "6.5", "--noise-frag-ppm", "6.5"],
    "A1A2": ["--noise-mz-ppm", "6.5", "--noise-frag-ppm", "6.5", "--noise-real-data"],
}


def run_flow(outdir, level, cond_flags, a) -> Path:
    """Run one pipeline (a density level under one noise condition) and return its outdir."""
    lvl_dir = Path(outdir) / f"n{level}_{a._cond_name}"
    log = Path(f"{lvl_dir}.log")
    cmd = [
        a.flow_python, str(Path(a.flow) / "timsim_flow.py"),
        "--bruker-reference", a.reference_d,
        "--proteome-spec", a.proteome_spec, "--mods", a.mods, "--design-spec", a.design_spec,
        "--search-fasta", a.search_fasta, "--samples", "A_R1",
        "--max-peptides", str(level), "--noise-seed", str(a.seed),
        "--outdir", str(lvl_dir),
        *cond_flags,
    ]
    env = dict(os.environ, PATH=f"{a.predict_bin}:{os.environ['PATH']}", TIMSIM_BIN=a.timsim_bin)
    print(f"  [{a._cond_name} n={level}] running flow -> {lvl_dir}", flush=True)
    with open(log, "w") as fh:
        subprocess.run(cmd, cwd=str(Path(a.flow) / "configs"), env=env, stdout=fh, stderr=fh, check=True)
    return lvl_dir


def read_metrics(lvl_dir) -> dict:
    m = sorted(glob.glob(f"{lvl_dir}/score*/**/metrics.json", recursive=True))
    if not m:
        raise FileNotFoundError(f"no score metrics.json under {lvl_dir}")
    d = json.load(open(m[-1]))
    diann, false = d["diann_ids"], d["false"]
    bg = d.get("background_subtracted") or 0
    fdp_sub = false / max(1, diann)
    fdp_raw = (false + bg) / max(1, diann + bg)  # un-subtract the real blank IDs
    # detectable recall = the strictest hierarchy level
    rec = d["hierarchy"][-1]["recall"] if d.get("hierarchy") else None
    return {"diann_ids": diann, "false": false, "background_subtracted": bg,
            "fdp_raw": fdp_raw, "fdp_sub": fdp_sub, "recall": rec, "q": d.get("q_threshold", 0.01)}


def read_density(lvl_dir) -> dict:
    t = sorted(glob.glob(f"{lvl_dir}/render/**/truth.parquet", recursive=True))
    if not t:
        return {}
    df = pq.read_table(t[-1]).to_pandas()
    present = df[df["abundance"] > 0]
    det = present[present["in_any_window"] & present["has_ms2"] & (present["abundance"] > 1e-3)]
    # co-elution proxy: mean precursors per 1-second RT bin (detectable set)
    coel = None
    if len(det):
        bins = (det["rt_seconds"] // 1).astype(int)
        coel = float(bins.value_counts().mean())
    return {"precursors": int(len(present)), "detectable": int(len(det)), "coelution_per_s": coel}


def svg_curve(rows, path):
    """Self-contained SVG: empirical FDP (raw, solid) + subtracted (dashed) vs detectable density, per
    condition, with the nominal-q line. Recall on a secondary implied scale is omitted for clarity."""
    W, H, ml, mr, mt, mb = 720, 420, 70, 160, 40, 60
    xs = [r["detectable"] for r in rows if r.get("detectable")]
    ys = [r["fdp_raw"] for r in rows] + [r["q"] for r in rows]
    if not xs:
        Path(path).write_text("<svg xmlns='http://www.w3.org/2000/svg'/>"); return
    xmin, xmax = min(xs), max(xs)
    ymax = max(max(ys) * 1.15, 0.02)
    def X(v): return ml + (v - xmin) / max(1, xmax - xmin) * (W - ml - mr)
    def Y(v): return H - mb - v / ymax * (H - mt - mb)
    colors = {"noiseless": "#6b7280", "A1": "#2563eb", "A1A2": "#dc2626"}
    p = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' font-family='sans-serif' font-size='12'>"]
    p.append(f"<rect width='{W}' height='{H}' fill='white'/>")
    p.append(f"<text x='{W/2}' y='20' text-anchor='middle' font-size='14' font-weight='bold'>Empirical FDP vs density (HeLa ramp)</text>")
    # axes
    p.append(f"<line x1='{ml}' y1='{H-mb}' x2='{W-mr}' y2='{H-mb}' stroke='#333'/><line x1='{ml}' y1='{mt}' x2='{ml}' y2='{H-mb}' stroke='#333'/>")
    p.append(f"<text x='{(ml+W-mr)/2}' y='{H-20}' text-anchor='middle'>detectable precursors</text>")
    p.append(f"<text transform='translate(20,{(mt+H-mb)/2}) rotate(-90)' text-anchor='middle'>empirical FDP</text>")
    # nominal q line
    q = rows[0]["q"]
    p.append(f"<line x1='{ml}' y1='{Y(q):.1f}' x2='{W-mr}' y2='{Y(q):.1f}' stroke='#999' stroke-dasharray='4 3'/>")
    p.append(f"<text x='{W-mr+4}' y='{Y(q)+4:.1f}' fill='#999'>nominal q={q*100:.0f}%</text>")
    # y ticks
    for k in range(5):
        v = ymax * k / 4
        p.append(f"<text x='{ml-8}' y='{Y(v)+4:.1f}' text-anchor='end' fill='#666'>{v*100:.1f}%</text>")
    conds = sorted({r["condition"] for r in rows})
    for ci, cond in enumerate(conds):
        pts = sorted([r for r in rows if r["condition"] == cond and r.get("detectable")], key=lambda r: r["detectable"])
        col = colors.get(cond, "#000")
        for style, key in (("", "fdp_raw"), ("stroke-dasharray='5 4'", "fdp_sub")):
            d = " ".join(f"{'M' if i==0 else 'L'}{X(r['detectable']):.1f},{Y(r[key]):.1f}" for i, r in enumerate(pts))
            p.append(f"<path d='{d}' fill='none' stroke='{col}' stroke-width='2' {style}/>")
        for r in pts:
            p.append(f"<circle cx='{X(r['detectable']):.1f}' cy='{Y(r['fdp_raw']):.1f}' r='3' fill='{col}'/>")
        p.append(f"<text x='{W-mr+4}' y='{mt+18+ci*16}' fill='{col}'>{cond}</text>")
    p.append(f"<text x='{W-mr+4}' y='{mt+18+len(conds)*16+10}' fill='#333'>solid=raw FDP</text>")
    p.append(f"<text x='{W-mr+4}' y='{mt+18+len(conds)*16+26}' fill='#333'>dashed=bg-subtracted</text>")
    p.append("</svg>")
    Path(path).write_text("\n".join(p))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="ramp", description="HeLa complexity ramp — empirical FDP vs density")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--levels", type=int, nargs="+", default=[2000, 8000, 20000, 40000],
                    help="--max-peptides density levels (nested for a fixed seed)")
    ap.add_argument("--conditions", nargs="+", default=["noiseless", "A1", "A1A2"], choices=list(CONDITIONS))
    ap.add_argument("--seed", type=int, default=41)
    ap.add_argument("--reference-d", required=True)
    ap.add_argument("--search-fasta", required=True)
    ap.add_argument("--proteome-spec", default="hela_proteome.toml")
    ap.add_argument("--mods", default="mods_basic.toml")
    ap.add_argument("--design-spec", default="design_hela.toml")
    ap.add_argument("--flow", default=str(Path(__file__).resolve().parent.parent / "flow"))
    ap.add_argument("--flow-python", default=sys.executable)
    ap.add_argument("--predict-bin", required=True, help="dir with timsim-ccs/rt/fragments on PATH")
    ap.add_argument("--timsim-bin", default="/scratch/timsim-demo/timsim-cli/target/release")
    ap.add_argument("--skip-run", action="store_true", help="only re-aggregate existing level dirs")
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for cond in a.conditions:
        a._cond_name = cond
        for level in a.levels:
            lvl_dir = a.outdir / f"n{level}_{cond}"
            if not a.skip_run:
                try:
                    run_flow(a.outdir, level, CONDITIONS[cond], a)
                except subprocess.CalledProcessError as e:
                    print(f"  FAILED [{cond} n={level}]: see {lvl_dir}.log", flush=True)
                    continue
            try:
                row = {"condition": cond, "level": level, **read_density(lvl_dir), **read_metrics(lvl_dir)}
                rows.append(row)
                rec = "—" if row["recall"] is None else f"{row['recall']*100:.1f}%"
                print(f"  [{cond} n={level}] detectable={row.get('detectable')} "
                      f"FDP_raw={row['fdp_raw']*100:.2f}% FDP_sub={row['fdp_sub']*100:.2f}% recall={rec}", flush=True)
            except FileNotFoundError as e:
                print(f"  no metrics for [{cond} n={level}]: {e}", flush=True)

    (a.outdir / "ramp.json").write_text(json.dumps(rows, indent=2))
    # markdown table
    lines = ["| condition | max_peptides | precursors | detectable | coelution/s | FDP_raw | FDP_sub | recall |",
             "|---|--:|--:|--:|--:|--:|--:|--:|"]
    for r in sorted(rows, key=lambda r: (r["condition"], r["level"])):
        rec = "—" if r["recall"] is None else f"{r['recall']*100:.1f}%"
        lines.append(f"| {r['condition']} | {r['level']:,} | {r.get('precursors',0):,} | {r.get('detectable',0):,} | "
                     f"{r.get('coelution_per_s') or 0:.1f} | {r['fdp_raw']*100:.2f}% | {r['fdp_sub']*100:.2f}% | {rec} |")
    (a.outdir / "ramp.md").write_text("\n".join(lines) + "\n")
    svg_curve(rows, a.outdir / "ramp.svg")
    print(f"\n  -> {a.outdir}/ramp.json  ramp.md  ramp.svg  ({len(rows)} points)")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
