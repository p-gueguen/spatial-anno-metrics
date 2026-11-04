"""Compare the scDesign3 world's verdicts and ISC regime calls with the toy world.

    python compare.py [--toy DIR]

``--toy``: a toy-world sweep run with the same number of seeds (default ``results/toy_3seeds``, made with
``python benchmarks/synthetic/sweep.py --seeds 3 --workers 16 --out <dir>`` at the same commit). The
committed 5-seed toy results (benchmarks/synthetic/results) are compared too, to show how many cells flip
with the seed count alone. Writes ``results/verdict_changes.csv``, ``results/verdict_counts.csv`` and
``results/isc_regimes.csv`` (the paper's fixed-threshold quadrant call, ``quadrant`` in the check files).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
TOY5 = HERE.parents[1] / "synthetic" / "results"
CLAIM = ["detects", "BLIND", "FALSE-ALARM", "WRONG-WAY", "erratic", "stable"]


def counts(v):
    c = v[v["expect"].isin(["down", "flat"])]["verdict"].value_counts()
    return c.reindex(CLAIM).fillna(0).astype(int)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--toy", default=str(RES / "toy_3seeds"))
    a = p.parse_args(argv)
    toy3 = pd.read_csv(Path(a.toy) / "verdicts.csv")
    toy5 = pd.read_csv(TOY5 / "verdicts.csv")
    scd = pd.read_csv(RES / "bands" / "verdicts.csv")
    keep = ["scenario", "metric", "expect", "verdict", "v0", "v_max", "threshold", "first_level"]
    m = (scd[keep].merge(toy3[keep], on=["scenario", "metric"], suffixes=("_scd3", "_toy3"), how="outer")
         .merge(toy5[["scenario", "metric", "verdict"]].rename(columns={"verdict": "verdict_toy5"}),
                on=["scenario", "metric"], how="left"))
    m["changed"] = m["verdict_scd3"] != m["verdict_toy3"]
    m["toy_seed_flip"] = m["verdict_toy3"] != m["verdict_toy5"]
    m.sort_values(["changed", "scenario", "metric"], ascending=[False, True, True]).to_csv(
        RES / "verdict_changes.csv", index=False, float_format="%.4g")
    pd.DataFrame({"toy, 5 seeds (committed)": counts(toy5), "toy, 3 seeds": counts(toy3),
                  "scDesign3 bands, 3 seeds": counts(scd)}).to_csv(RES / "verdict_counts.csv")
    print(pd.DataFrame({"toy5": counts(toy5), "toy3": counts(toy3), "scd3": counts(scd)}))
    print(f"cells: {len(m)}, changed toy3 -> scd3: {int(m['changed'].sum())}, toy3 vs toy5: "
          f"{int(m['toy_seed_flip'].sum())}")
    regs = []
    for name, path in (("toy", Path(a.toy)), ("scDesign3 bands", RES / "bands"),
                       ("scDesign3 dirichlet", RES / "dirichlet")):
        f = path / "check_isc_regimes.csv"
        if f.exists():
            d = pd.read_csv(f)
            # the paper's quadrant: ``quadrant`` since the calibrated ``regime`` (older files: ``regime``);
            # the calibrated call per world is in ../results/cross_world_isc_regimes.csv
            if "quadrant" in d:
                d = d.drop(columns="regime").rename(columns={"quadrant": "regime"})
            regs.append(d.assign(world=name))
    r = pd.concat(regs)
    # the labels a perturbation touched: merged pairs ("x+y"), split halves ("x_a" / "x_b"); else every label
    t = r["type"].astype(str)
    r["labels"] = "all"
    r.loc[t.str.contains(r"\+", regex=True), "labels"] = "merged"
    r.loc[t.str.endswith(("_a", "_b")), "labels"] = "split halves"
    r = r[(r["labels"] != "all") | ~r["case"].str.startswith(("merge", "over_split"))]
    tab = (r.groupby(["world", "case", "labels", "regime"]).size().unstack(fill_value=0)
           .join(r.groupby(["world", "case", "labels"])[["iscg", "iscl"]].mean().round(3)))
    tab.to_csv(RES / "isc_regimes.csv")
    with pd.option_context("display.width", 200):
        print(tab)


if __name__ == "__main__":
    main()
