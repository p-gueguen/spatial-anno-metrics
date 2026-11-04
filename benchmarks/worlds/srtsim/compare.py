"""Which metric verdicts change between the toy world and the SRTsim world.

    python benchmarks/worlds/srtsim/compare.py [--toy benchmarks/synthetic/results/verdicts.csv]
                                               [--toy3 <toy sweep rerun with --seeds 3>/verdicts.csv]

Writes ``verdict_changes.csv`` (every (scenario, metric) whose verdict differs, with both worlds' baseline and
max-strength values) and ``verdict_summary.md`` (tallies over the cells that carry a claim, and the changes).
``--toy3``: the toy sweep rerun with the same number of seeds as this world; a change that also appears
between the toy at 5 and at 3 seeds is flagged ``seed_count_sensitive`` (the move threshold is 3 x the seed SD
at strength 0, so the seed count alone can flip a borderline cell).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GOOD = {"detects", "stable"}
CLAIM = ["detects", "BLIND", "FALSE-ALARM", "WRONG-WAY", "erratic", "stable"]


def tally(v):
    c = v[v["expect"].isin(["down", "flat"])]["verdict"].value_counts()
    return {k: int(c.get(k, 0)) for k in CLAIM}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--toy", default=str(ROOT / "benchmarks" / "synthetic" / "results" / "verdicts.csv"))
    p.add_argument("--toy3", default=None)
    p.add_argument("--world", default=str(HERE / "results" / "verdicts.csv"))
    p.add_argument("--out", default=str(HERE / "results"))
    a = p.parse_args(argv)
    toy, w = pd.read_csv(a.toy), pd.read_csv(a.world)
    key = ["scenario", "metric"]
    m = toy.merge(w, on=key, suffixes=("_toy", "_srtsim"), how="outer")
    if a.toy3:
        t3 = pd.read_csv(a.toy3)[key + ["verdict"]].rename(columns={"verdict": "verdict_toy3"})
        m = m.merge(t3, on=key, how="left")
        m["seed_count_sensitive"] = m["verdict_toy3"] != m["verdict_toy"]
    ch = m[m["verdict_toy"] != m["verdict_srtsim"]].copy()
    ch["direction"] = [("better" if b in GOOD and x not in GOOD else "worse" if x in GOOD and b not in GOOD
                        else "other") for x, b in zip(ch["verdict_toy"], ch["verdict_srtsim"])]
    cols = key + ["expect_toy", "verdict_toy", "verdict_srtsim", "direction", "v0_toy", "v_max_toy",
                  "v0_srtsim", "v_max_srtsim", "threshold_toy", "threshold_srtsim", "first_level_toy",
                  "first_level_srtsim"] + (["verdict_toy3", "seed_count_sensitive"] if a.toy3 else [])
    ch = ch[cols].sort_values(["direction", "scenario", "metric"])
    out = Path(a.out)
    ch.to_csv(out / "verdict_changes.csv", index=False, float_format="%.4g")
    tt, tw = tally(toy), tally(w)
    lines = ["| verdict (cells with a claim) | toy (5 seeds) | SRTsim (3 seeds) |", "|---|---|---|"]
    lines += [f"| {k} | {tt[k]} | {tw[k]} |" for k in CLAIM]
    if a.toy3:
        t3 = tally(pd.read_csv(a.toy3))
        lines[0] += " toy (3 seeds) |"
        lines[1] += "---|"
        lines[2:] = [f"{ln} {t3[k]} |" for ln, k in zip(lines[2:], CLAIM)]
    lines += ["", f"{len(ch)} of {len(m)} (scenario, metric) cells change verdict "
                  f"({(ch['direction'] == 'better').sum()} to detects/stable, "
                  f"{(ch['direction'] == 'worse').sum()} away from it, {(ch['direction'] == 'other').sum()} other)"
                  + (f"; {int(ch['seed_count_sensitive'].sum())} of them also differ between the toy at 5 and "
                     "at 3 seeds." if a.toy3 else "."), ""]
    lines += ["| scenario | metric | toy | SRTsim | toy v0 → vmax | SRTsim v0 → vmax |" + (" seeds? |" if a.toy3 else ""),
              "|---|---|---|---|---|---|" + ("---|" if a.toy3 else "")]
    for r in ch.itertuples():
        lines.append(f"| {r.scenario} | `{r.metric}` | {r.verdict_toy} | {r.verdict_srtsim} | "
                     f"{r.v0_toy:.3g} → {r.v_max_toy:.3g} | {r.v0_srtsim:.3g} → {r.v_max_srtsim:.3g} |"
                     + ((" yes |" if r.seed_count_sensitive else " |") if a.toy3 else ""))
    (out / "verdict_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
