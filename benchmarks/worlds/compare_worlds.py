"""Put the committed sweep verdicts of every world side by side (benchmarks/worlds/README.md, cross-world table).

    python benchmarks/worlds/compare_worlds.py [--out benchmarks/worlds/results]

Reads only committed results (no sweep is run):

* toy: ``benchmarks/synthetic/results`` (5 seeds);
* scDesign3: ``scdesign3/results/bands`` (3 seeds);
* SRTsim: ``srtsim/results`` (3 seeds, a different crop set per seed) and, when present,
  ``srtsim/results/fixed_crops_0`` (the crops of seed 0 for every seed, only the count draw varies);
* scCube: ``sccube/results/{stripes,domains,fine,mixed}`` (3 seeds each).

Writes ``cross_world_verdicts.csv`` (one row per scenario x metric, a verdict and a clean value per world),
``cross_world_matrix.md`` (metric x scenario, claim cells only), ``cross_world_counts.csv`` (verdict counts per
world) and ``cross_world_baselines.csv`` (clean values of every metric per world).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCES = {
    "toy": HERE.parent / "synthetic" / "results",
    "scdesign3": HERE / "scdesign3" / "results" / "bands",
    "srtsim": HERE / "srtsim" / "results",
    "srtsim_fixed": HERE / "srtsim" / "results" / "fixed_crops_0",
    "sccube_stripes": HERE / "sccube" / "results" / "stripes",
    "sccube_domains": HERE / "sccube" / "results" / "domains",
    "sccube_fine": HERE / "sccube" / "results" / "fine",
    "sccube_mixed": HERE / "sccube" / "results" / "mixed",
}
SCCUBE = ["sccube_stripes", "sccube_domains", "sccube_fine", "sccube_mixed"]
REALISTIC = ["scdesign3", "srtsim", *SCCUBE]
SYMBOL = {"detects": "✅", "BLIND": "❌", "WRONG-WAY": "⛔", "erratic": "〰️", "stable": "✔", "FALSE-ALARM": "⚠️",
          "drifts-up": "↗", "worse": "↓", "better": "↑", "flat": "·", "n/a": "∅"}
SCENARIOS = ["label_noise", "sibling_confusion", "type_swap", "patch_error", "over_split", "merge", "spillover",
             "depth_loss", "ambient", "batch", "coord_shuffle", "coexpression", "rare_type", "infiltration",
             "panel_growth"]


def load() -> pd.DataFrame:
    frames = []
    for world, d in SOURCES.items():
        f = d / "verdicts.csv"
        if not f.exists():
            continue
        v = pd.read_csv(f)[["scenario", "metric", "expect", "verdict", "v0", "v_max", "threshold"]]
        frames.append(v.assign(world=world))
    return pd.concat(frames, ignore_index=True)


def sccube_summary(row: pd.Series) -> str:
    """One symbol when the four patterns agree, else the counts, e.g. ``✅3 ❌1``."""
    vs = [row[w] for w in SCCUBE if isinstance(row.get(w), str)]
    if not vs:
        return ""
    counts = pd.Series(vs).value_counts()
    if len(counts) == 1:
        return SYMBOL.get(counts.index[0], counts.index[0])
    return " ".join(f"{SYMBOL.get(k, k)}{n}" for k, n in counts.items())


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=str(HERE / "results"))
    out = Path(p.parse_args(argv).out)
    out.mkdir(parents=True, exist_ok=True)
    df = load()
    worlds = [w for w in SOURCES if w in set(df.world)]
    wide = df.pivot_table(index=["scenario", "metric"], columns="world", values="verdict", aggfunc="first")
    v0 = df.pivot_table(index=["scenario", "metric"], columns="world", values="v0", aggfunc="first")
    expect = df[df.world == "toy"].set_index(["scenario", "metric"])["expect"]
    tab = wide[worlds].join(v0[worlds].add_prefix("v0_")).join(expect).reset_index()
    tab["sccube"] = tab.apply(sccube_summary, axis=1)
    real = [w for w in REALISTIC if w in worlds]
    same = tab[real].nunique(axis=1, dropna=False) == 1
    tab["all_realistic_agree"] = same
    tab["realistic_vs_toy"] = [
        ("same as toy" if r[real[0]] == r["toy"] else "all realistic differ from toy") if s
        else ("toy verdict kept in no realistic world" if all(r[w] != r["toy"] for w in real) else "mixed")
        for s, (_, r) in zip(same, tab.iterrows())]
    tab.to_csv(out / "cross_world_verdicts.csv", index=False)

    claim = tab[tab.expect != "-"]
    cols = ["toy", "scdesign3", "srtsim"]
    lines = ["| metric | " + " | ".join(SCENARIOS) + " |", "|---|" + "---|" * len(SCENARIOS)]
    for m in sorted(claim.metric.unique()):
        cells = []
        for s in SCENARIOS:
            r = claim[(claim.metric == m) & (claim.scenario == s)]
            if r.empty:
                cells.append("")
                continue
            r = r.iloc[0]
            cells.append(" / ".join(SYMBOL.get(r[c], str(r[c])) if isinstance(r[c], str) else "-" for c in cols)
                         + " / " + (r["sccube"] or "-"))
        lines.append(f"| `{m}` | " + " | ".join(cells) + " |")
    (out / "cross_world_matrix.md").write_text(
        "Each cell: toy / scDesign3 (bands) / SRTsim / scCube (one symbol if the 4 patterns agree, else counts). "
        "Claim cells only.\n\n" + "\n".join(lines) + "\n")

    counts = (claim.melt(id_vars=["scenario", "metric"], value_vars=worlds, var_name="world", value_name="verdict")
              .dropna().groupby(["world", "verdict"]).size().unstack(fill_value=0).reindex(worlds))
    counts.to_csv(out / "cross_world_counts.csv")
    base = (df[df.scenario == "label_noise"].pivot_table(index="metric", columns="world", values="v0",
                                                         aggfunc="first")[worlds])
    base.to_csv(out / "cross_world_baselines.csv", float_format="%.4g")
    print(counts.to_string())
    print(tab.groupby("realistic_vs_toy").size().to_string())


if __name__ == "__main__":
    main()
