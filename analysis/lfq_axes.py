#!/usr/bin/env python3
"""Axes of variation in real DIA data — measured from the LFQ Benchmark Gen Beta DIA-NN reports.

WHY THIS EXISTS
---------------
The v2 realism audit left one blocker: the render is correct at the bottom and middle of the
intensity distribution and wrong in the upper tail, which is the signature of SIGNAL SPREADING — too
much of one peptide's signal piled into too few bins. Measuring that was assumed to need a partner
to characterise real data.

It does not. DIA-NN already reports, per precursor per run: `FWHM` (chromatographic peak width in
minutes), `RT.Start`/`RT.Stop` (elution extent), `IM` (ion mobility), and several intensity columns.
The LFQ Benchmark Gen Beta supplement ships those reports for FIVE current instruments. So the RT
axis of the spreading question is answerable from a 3.4 GB extract instead of a 2.9 TB download.

WHAT IT REPORTS
---------------
Per (instrument, method, gradient), over target precursors at Q.Value <= --qvalue:
  * ID density      — precursors and protein groups per run
  * RT spreading    — FWHM quantiles, and elution extent RT.Stop-RT.Start
  * intensity shape — p10/p50/p90/p99 of Ms1.Area and Precursor.Quantity, and p99/p50 as a
                      SCALE-INVARIANT tail statistic (the one our render fails)
  * mobility        — IM range, where the instrument has it

DENOMINATOR WARNING
-------------------
These runs are 5/11/15/30-min gradients at 50 ng / 250 pg, and our own cohorts are typically a
~40-min gradient at 200 ng. Peak width scales with gradient length, so DO NOT compare our sigma against a
5-min run without matching. The `--gradient-min` column is printed for exactly that reason: compare
like with like, or state the mismatch.

Usage
    lfq_axes.py --root <dir with datasets_diann_output/> [--qvalue 0.01] [--out summary.tsv]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

# Columns we actually need. Reading a subset keeps a 100 MB report inside a sane memory budget and
# makes the whole sweep runnable on the dev box rather than only on a cluster.
COLS = ["Run", "Precursor.Id", "Protein.Group", "Genes", "Q.Value", "Decoy",
        "RT", "RT.Start", "RT.Stop", "FWHM", "IM", "Ms1.Area", "Precursor.Quantity"]


def parse_name(path: Path) -> dict:
    """Instrument / method / gradient / load out of the directory and file name."""
    inst = path.parent.name
    stem = path.stem
    # Tolerate `5min`, `5-min`, `5_min`, `0.25ng`, `250PG` — the LFQ set uses one convention, but a
    # silent "unknown" here turns into a NaN gradient that groups wrongly downstream.
    grad = (re.search(r"(\d+(?:\.\d+)?)[-_ ]?min", stem, re.I) or [None, None])[1]
    load = (re.search(r"(\d+(?:\.\d+)?(?:pg|ng|ug))", stem, re.I) or [None, None])[1]
    # Longest-first, so "diagonalPASEF" is not swallowed by "diaPASEF" and the bare "DIA" fallback
    # only fires when no specific method matched. Several Astral files (e.g. "DIANN_Optimized_DI")
    # do not carry a method token at all — they are plain Astral DIA, so "DI"/"DIANN" must not be
    # read as a method, and the instrument directory settles it.
    method = None
    for m in ("diagonalPASEF", "diaPASEF", "ZTScanDIA", "ZenoSWATH", "PASEF"):
        if m.lower() in stem.lower():
            method = m
            break
    if method is None:
        method = "DIA" if "Astral" in inst or "DIA" in stem.upper() else "?"
    return {"instrument": inst, "method": method,
            "gradient_min": float(grad) if grad else np.nan, "load": load or "?", "file": path.name}


def q(a: np.ndarray, p: float) -> float:
    return float(np.percentile(a, p)) if len(a) else np.nan


def read_report(path: Path) -> pd.DataFrame:
    """DIA-NN writes .parquet or .tsv depending on version and flags; accept both."""
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t", low_memory=False)
    try:
        return pd.read_parquet(path, columns=list(COLS))
    except (ValueError, KeyError):         # a report LACKS a column -> read everything instead
        return pd.read_parquet(path)
    # Anything else (corrupt file, unreadable parquet) propagates: swallowing it here would retry a
    # potentially enormous all-column read and then surface as a confusing missing-column error.


def summarise(path: Path, qvalue: float, label: str | None = None) -> dict | None:
    df = read_report(path)
    meta = parse_name(path)
    if label:
        # An explicit label beats the filename grammar, which only fits the LFQ-benchmark naming.
        meta["instrument"] = label
        meta["method"] = meta.get("method", "?")
    if "Decoy" in df:
        df = df[df["Decoy"] == 0]
    if "Q.Value" in df:
        df = df[df["Q.Value"] <= qvalue]
    if df.empty:
        return None
    nruns = df["Run"].nunique()
    out = dict(meta)
    out["runs"] = nruns
    out["prec_per_run"] = round(len(df) / nruns)
    out["pg_per_run"] = round(df.groupby("Run")["Protein.Group"].nunique().mean()) if "Protein.Group" in df else np.nan

    # --- RT spreading: the axis the simulator gets wrong -------------------------------------
    if "FWHM" in df:
        f = df["FWHM"].to_numpy(float); f = f[np.isfinite(f) & (f > 0)]
        # DIA-NN reports FWHM in MINUTES; seconds is the unit our --sigma-seconds is in.
        out["fwhm_s_p10"], out["fwhm_s_p50"], out["fwhm_s_p90"] = (round(q(f, p) * 60, 2) for p in (10, 50, 90))
        # A Gaussian's sigma is FWHM / (2*sqrt(2*ln2)); this is the number to compare with
        # `timsim-render --sigma-seconds` (default 3.0).
        out["sigma_s_equiv"] = round(q(f, 50) * 60 / 2.3548, 2)
    if {"RT.Start", "RT.Stop"} <= set(df.columns):
        e = (df["RT.Stop"] - df["RT.Start"]).to_numpy(float); e = e[np.isfinite(e) & (e > 0)]
        out["elution_s_p50"] = round(q(e, 50) * 60, 2)

    # --- intensity shape: p99/p50 is invariant under any pure rescale, so it is comparable
    #     across instruments AND across our own calibration sweeps ------------------------------
    for col, tag in (("Ms1.Area", "ms1"), ("Precursor.Quantity", "quant")):
        if col not in df:
            continue
        v = df[col].to_numpy(float); v = v[np.isfinite(v) & (v > 0)]
        if not len(v):
            continue
        out[f"{tag}_p50"] = f"{q(v, 50):.3e}"
        out[f"{tag}_p99_over_p50"] = round(q(v, 99) / q(v, 50), 1)
        out[f"{tag}_max_over_p50"] = round(float(v.max()) / q(v, 50), 1)

    if "IM" in df:
        im = df["IM"].to_numpy(float); im = im[np.isfinite(im) & (im > 0)]
        if len(im):
            out["im_p10"], out["im_p90"] = round(q(im, 10), 4), round(q(im, 90), 4)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, help="directory to sweep for *report.parquet")
    ap.add_argument("--report", type=Path, action="append", default=[],
                    help="a single DIA-NN report (.parquet or .tsv) to score — repeatable. Use this "
                         "to put OUR OWN render's report on exactly the same footing as the real "
                         "instruments, which is the only way the comparison means anything.")
    ap.add_argument("--label", help="instrument/method label for --report (default: from the path)")
    ap.add_argument("--qvalue", type=float, default=0.01)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()

    # Resolve before de-duplicating: the same file reachable as `--report x.parquet` and via
    # `--root .` is ONE report, and scoring it twice would put a duplicate row into the summary and
    # quietly bias any aggregate taken over it. Track explicitness separately from path identity —
    # `p in a.report` compares Path VALUES, so a differently-spelled path to the same file would
    # both duplicate the row and mislabel it.
    explicit = {p.resolve() for p in a.report}
    seen, reports = set(), []
    for p in list(a.report) + (sorted(a.root.rglob("*report.parquet")) + sorted(a.root.rglob("*report.tsv")) if a.root else []):
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        reports.append(rp)
    if not reports:
        raise SystemExit("nothing to score: pass --root and/or --report")
    if a.label and len(explicit) != 1:
        raise SystemExit(f"--label names one report, but {len(explicit)} were passed with --report")
    print(f"[lfq-axes] {len(reports)} reports, target precursors at Q.Value <= {a.qvalue}\n", flush=True)
    rows = []
    for p in reports:
        try:
            r = summarise(p, a.qvalue, a.label if p in explicit else None)
        except Exception as e:                       # one bad report must not kill the sweep
            print(f"  !! {p.name}: {type(e).__name__}: {e}", flush=True)
            continue
        if r:
            rows.append(r)
            print(f"  ok {r['instrument']:<24}{r['method']:<16}{r['gradient_min']:>5.0f}min "
                  f"{r['load']:>6}  prec/run={r['prec_per_run']:>7}  "
                  f"FWHM_p50={r.get('fwhm_s_p50', float('nan')):>6}s  "
                  f"sigma~{r.get('sigma_s_equiv', float('nan')):>5}s", flush=True)
    df = pd.DataFrame(rows)
    if a.out:
        df.to_csv(a.out, sep="\t", index=False)
        print(f"\n[lfq-axes] -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
