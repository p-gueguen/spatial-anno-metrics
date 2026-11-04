"""Compare the scCube-world verdicts with the toy world's.

    python benchmarks/worlds/sccube/compare.py --toy3 <toy sweep run with --seeds 3>

Reads ``benchmarks/synthetic/results/verdicts.csv`` (toy, 5 seeds), the optional toy run at 3 seeds (same seed
count as the scCube runs, to tell a world effect from a seed-count effect) and ``results/<pattern>/verdicts.csv``.
Writes ``results/verdict_counts.csv`` (verdicts per source over the cells that carry a claim) and
``results/verdict_changes.csv`` (every scenario x metric cell whose verdict differs from the toy's in at least one
pattern; ``toy3_same`` says whether the toy at 3 seeds keeps the 5-seed verdict, so a change there is not a
seed-count artefact).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
TOY = HERE.parents[1] / "synthetic" / "results" / "verdicts.csv"
PATTERNS = ["stripes", "domains", "fine", "mixed"]
CLAIM = ["detects", "BLIND", "FALSE-ALARM", "WRONG-WAY", "erratic", "stable", "drifts-up"]


def load(p):
    d = pd.read_csv(p)
    return d.set_index(["scenario", "metric"])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy3", default=None)
    ap.add_argument("--results", default=str(HERE / "results"))
    a = ap.parse_args(argv)
    res = Path(a.results)
    src = {"toy (5 seeds)": load(TOY)}
    if a.toy3:
        src["toy (3 seeds)"] = load(Path(a.toy3) / "verdicts.csv")
    for p in PATTERNS:
        f = res / p / "verdicts.csv"
        if f.exists():
            src[f"sccube:{p}"] = load(f)
    counts = []
    for name, d in src.items():
        c = d[d["expect"].isin(["down", "flat"])]["verdict"].value_counts().reindex(CLAIM, fill_value=0)
        counts.append({"source": name, **c.to_dict(), "n_claim_cells": int(c.sum())})
    pd.DataFrame(counts).to_csv(res / "verdict_counts.csv", index=False)
    print(pd.DataFrame(counts).to_string(index=False))
    toy = src["toy (5 seeds)"]
    rows = []
    for key in toy.index:
        r = {"scenario": key[0], "metric": key[1], "expect": toy.loc[key, "expect"],
             "toy": toy.loc[key, "verdict"], "toy_v0": toy.loc[key, "v0"], "toy_vmax": toy.loc[key, "v_max"]}
        if "toy (3 seeds)" in src and key in src["toy (3 seeds)"].index:
            r["toy3_same"] = src["toy (3 seeds)"].loc[key, "verdict"] == r["toy"]
        changed = False
        for p in PATTERNS:
            d = src.get(f"sccube:{p}")
            if d is None or key not in d.index:
                continue
            r[p] = d.loc[key, "verdict"]
            r[f"{p}_v0"], r[f"{p}_vmax"] = d.loc[key, "v0"], d.loc[key, "v_max"]
            changed |= r[p] != r["toy"]
        if changed:
            rows.append(r)
    ch = pd.DataFrame(rows)
    ch.to_csv(res / "verdict_changes.csv", index=False, float_format="%.4g")
    claim = ch[ch["expect"].isin(["down", "flat"])]
    print(f"{len(ch)} cells change in >= 1 pattern; {len(claim)} of them carry a claim")
    for p in PATTERNS:
        if p in ch:
            n = (claim[p] != claim["toy"]).sum()
            robust = ((claim[p] != claim["toy"]) & claim.get("toy3_same", True)).sum()
            print(f"  {p}: {n} claim cells differ from toy-5 ({robust} where toy-3 keeps the toy-5 verdict)")


if __name__ == "__main__":
    main()
