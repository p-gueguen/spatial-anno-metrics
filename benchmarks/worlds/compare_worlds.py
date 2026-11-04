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
world), ``cross_world_baselines.csv`` (clean values of every metric per world) and
``cross_world_isc_regimes.csv`` (ISC ``regime`` and ``quadrant`` calls per world and case, from each world's
``check_isc_regimes.csv``; the scDesign3 Dirichlet variant is included there).
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
    isc = isc_regimes()
    if not isc.empty:
        isc.to_csv(out / "cross_world_isc_regimes.csv", index=False)
        print(isc.to_string(index=False))


ISC_SOURCES = {**SOURCES, "scdesign3_dirichlet": HERE / "scdesign3" / "results" / "dirichlet"}
ISC_CASES = {"clean": "clean", "over_split": "over_split 1 (Epi halves)", "merge": "merge",
             "untouched": "", "label_noise": "label_noise 20%", "batch": "batch sd 1.0"}


def _calls(s: pd.Series) -> str:
    return ", ".join(f"{k} {n}" for k, n in s.value_counts().items())


def isc_regimes() -> pd.DataFrame:
    """ISC regime calls per world from each ``check_isc_regimes.csv`` that records both calls.

    ``regime`` is the calibrated call, ``quadrant`` the paper's fixed 0.74 / 0.45 guides. For ``over_split`` only
    the two halves are counted, for ``merge`` only the merged labels (``merge 1`` and ``merge 3`` together),
    for ``untouched`` the other labels of those two cases; otherwise every scored label. ``wrong_regime`` /
    ``wrong_quadrant`` count definite wrong actions, graded as in ``benchmarks/isc_regimes_real.py``:
    ``merge`` / ``remove`` / ``split`` on the clean and batch sections and for untouched labels (``inconclusive``
    is not counted), ``merge`` / ``split`` under label noise (a degraded label may read ``inconclusive`` or
    ``remove``), anything but ``merge`` for the halves; none for merged labels (not identifiable from ISC).
    """
    rows = []
    for world, d in ISC_SOURCES.items():
        f = d / "check_isc_regimes.csv"
        if not f.exists():
            continue
        r = pd.read_csv(f)
        if "quadrant" not in r:
            continue
        t = r["type"].astype(str)
        made = t.str.endswith(("_a", "_b")) | t.str.contains("+", regex=False)
        for key, case in ISC_CASES.items():
            sel = r["case"].str.startswith(case) if key == "merge" else r["case"] == case
            if key == "over_split":
                sel &= t.str.endswith(("_a", "_b"))
            elif key == "merge":
                sel &= t.str.contains("+", regex=False)
            elif key == "untouched":   # the labels an over_split or merge left as they were
                sel = r["case"].str.startswith(("over_split", "merge")) & ~made
            x = r[sel]
            if x.empty:
                continue
            row = dict(world=world, case=key, seeds=x["seed"].nunique(), n=len(x),
                       regime=_calls(x["regime"]), quadrant=_calls(x["quadrant"]),
                       iscg=round(x["iscg"].mean(), 3), iscl=round(x["iscl"].mean(), 3))
            for call in ("regime", "quadrant"):
                c = x[call]
                if key in ("clean", "batch", "untouched"):
                    bad = c.isin(["merge", "remove", "split"])
                elif key == "label_noise":
                    bad = c.isin(["merge", "split"])
                elif key == "over_split":
                    bad = c != "merge"
                else:
                    bad = pd.Series(False, index=c.index)
                row[f"wrong_{call}"] = int(bad.sum())
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
