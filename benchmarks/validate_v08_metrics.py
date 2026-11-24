"""Validation benchmark for v0.8.0 OpenProblems-inspired spatial and reference metrics.

Runs:
1. Synthetic 15-scenario toy sweep via ``benchmarks/synthetic/sweep.py`` comparing
   the three new metrics side-by-side against existing battery metrics.
2. Real-data perturbations across:
   - 10x Xenium breast crop (``tests/goldens/signalqc_spatialqm/``)
   - Yerly et al. 2022 human BCC cohort (``tests/goldens/isc_yerly2022/``,
     split into ref patients vs query patients)
   - MERFISH mouse brain 4-section cohort
     (public MERFISH deposit, Animal 1, Bregma 0.01 / 0.06 / 0.11 / 0.16;
     location set with ``SAM_MERFISH_BASE`` and ``SAM_MERFISH_PREFIX``)

Outputs under ``benchmarks/v08_validation/``:
- ``toy_verdicts.csv``
- ``real_perturbation.csv``
- ``summary.md``
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sam_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/sam_numba")
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks" / "synthetic"))
sys.path.insert(0, str(ROOT / "benchmarks"))

import real_perturbation as rp  # noqa: E402
import spatial_anno_metrics as sam  # noqa: E402
import sweep  # noqa: E402

OUT_DIR = ROOT / "benchmarks" / "v08_validation"
MERFISH_BASE = Path(os.environ.get("SAM_MERFISH_BASE", "data/merfish_bregma"))
MERFISH_PREFIX = os.environ.get("SAM_MERFISH_PREFIX", "merfish_Animal1_Bregma_")
MERFISH_SECTIONS = ["0.01", "0.06", "0.11", "0.16"]

COMPARE_METRICS = [
    "rel_pw_ct",
    "resid_moran.explained",
    "resid_moran.residual",
    "spatial_org.corr",
    "spatial_org.offdiag",
    "sc.mean_coherence",
    "sc.pas",
    "sc.chaos",
    "ns.n_suspicious",
    "nmp",
    "coex.by_type",
    "isc.integrated",
    "ext.macro_f1",
]

V08_ORIENTATION: dict[str, int] = {
    "rel_pw_ct": 1,
    "resid_moran.explained": 1,
    "resid_moran.residual": -1,
    "resid_moran.raw": 1,
    "spatial_org.corr": 1,
    "spatial_org.spearman": 1,
    "spatial_org.offdiag": 1,
    "sc.mean_coherence": 1,
    "sc.pas": -1,
    "sc.chaos": -1,
    "nc.mean": -1,
    "nmp.salas": 1,
    "coex.by_type": 1,
    "isc.integrated": 1,
    "ext.macro_f1": 1,
}


def _ad_rel_pw(a: ad.AnnData, ref: ad.AnnData) -> dict[str, float]:
    res = sam.relative_pairwise_ct_expression(
        a, ref, label_key="label", ref_label_key="cell_type", layer="counts"
    )
    return {"rel_pw_ct": res.get("relative_pairwise_ct_expression", np.nan)}


def _ad_resid_moran(a: ad.AnnData, ref: ad.AnnData) -> dict[str, float]:
    b = a.copy()
    for k in list(b.obsp):
        del b.obsp[k]
    res = sam.residual_spatial_autocorrelation(b, label_key="label", layer="counts", graph=None)
    return {
        "resid_moran.explained": res.get("explained_fraction", np.nan),
        "resid_moran.residual": res.get("residual_median_moran", np.nan),
    }


def _ad_spatial_org(a: ad.AnnData, ref: ad.AnnData) -> dict[str, float]:
    res = sam.spatial_organization_consistency(a, label_key="label", sample_key="sample")
    return {
        "spatial_org.corr": res.get("mean_corr", np.nan),
        "spatial_org.offdiag": res.get("offdiag_mean_corr", np.nan),
    }


def run_synthetic_part(seeds: int = 3, workers: int = 8) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run 15-scenario toy sweep with the three new v0.8 metrics included."""
    tmp_out = OUT_DIR / "sweep_tmp"
    extra_adapters = {
        "rel_pairwise_ct": _ad_rel_pw,
        "residual_spatial_autocorrelation": _ad_resid_moran,
        "spatial_org_consistency": _ad_spatial_org,
    }
    extra_metrics = {
        "rel_pw_ct": (1, {sweep.L, sweep.PUR, sweep.REF}),
        "resid_moran.explained": (1, {sweep.L, sweep.SP}),
        "resid_moran.residual": (-1, {sweep.L, sweep.SP}),
        "spatial_org.corr": (1, {sweep.L, sweep.SP}),
        "spatial_org.offdiag": (1, {sweep.L, sweep.SP}),
    }
    df, sc = sweep.main(
        argv=["--seeds", str(seeds), "--workers", str(workers), "--out", str(tmp_out)],
        extra_adapters=extra_adapters,
        extra_metrics=extra_metrics,
    )
    sub = sc[sc["metric"].isin(COMPARE_METRICS)].copy()
    # sweep.DATA_EXPECT leaves spatial_org.* as 'flat' under coord_shuffle; mark true semantic expectation.
    mask_coord = sub["metric"].isin(["spatial_org.corr", "spatial_org.offdiag"]) & (
        sub["scenario"] == "coord_shuffle"
    )
    sub["openproblems_expect"] = sub["expect"]
    sub.loc[mask_coord, "openproblems_expect"] = "down (spatial)"
    sub.to_csv(OUT_DIR / "toy_verdicts.csv", index=False)
    err_df = pd.read_csv(tmp_out / "errors.csv")
    return df, sub, err_df


def _score_real(df: pd.DataFrame, expectations: dict[tuple[str, str], set[str]]) -> pd.DataFrame:
    """Compute per-(dataset, scenario, metric) monotonicity and threshold verdicts."""
    rows = []
    for (ds, scen, metric), g in df.groupby(["dataset", "scenario", "metric"], sort=False):
        levels = sorted(g["level"].unique())
        o = V08_ORIENTATION[metric]
        m = g.groupby("level")["value"].agg(["mean", "std"]).reindex(levels)
        v0 = float(m["mean"].iloc[0])
        v_end = float(m["mean"].iloc[-1])
        bad = -o * (m["mean"] - v0)
        thr = max(3.0 * float(np.nan_to_num(m["std"].iloc[0])), 0.02 * max(abs(v0), 1e-3), 1e-6)
        ok = ~bad.isna()
        rho = (
            float(spearmanr(np.arange(ok.sum()), bad[ok].to_numpy())[0])
            if ok.sum() > 2 and np.nanstd(bad[ok]) > 0
            else (1.0 if ok.sum() == 2 and bad[ok].iloc[-1] > thr else 0.0)
        )
        moved = float(np.nanmax(np.abs(bad.to_numpy()))) > thr if ok.any() else False
        expected_set = expectations.get((ds, scen), set())
        if scen == "nuisance_shuffle":
            verdict = "stable" if not moved else "FALSE-ALARM"
            exp = "flat"
        elif metric in expected_set:
            exp = "down"
            if moved and float(bad.iloc[-1]) > thr and rho >= 0.8:
                verdict = "detects"
            elif moved and float(bad.iloc[-1]) < -thr:
                verdict = "WRONG-WAY"
            else:
                verdict = "weak/flat"
        else:
            exp = "-"
            if not moved:
                verdict = "flat"
            elif float(bad.iloc[-1]) < -thr:
                verdict = "improves (inverted)"
            else:
                verdict = "moved"
        first_lv = next((lv for lv, b in zip(levels, bad) if pd.notna(b) and b > thr), None)
        rows.append(
            {
                "dataset": ds,
                "scenario": scen,
                "metric": metric,
                "expect": exp,
                "v0": v0,
                "v_end": v_end,
                "delta": v_end - v0,
                "spearman_rho": rho,
                "threshold": thr,
                "first_level": first_lv,
                "verdict": verdict,
            }
        )
    return pd.DataFrame(rows)


def _perturb_xenium_ext(base: ad.AnnData, scenario: str, level: float, seed: int) -> ad.AnnData:
    if scenario in {"spillover", "coord_shuffle", "nuisance_shuffle"}:
        return rp.perturb_xenium(base, scenario, level, seed)
    a = base.copy()
    rng = np.random.default_rng(seed)
    lab = a.obs["truth"].astype(str).to_numpy().copy()
    types = sorted(set(lab))
    if scenario == "label_noise":
        n_flip = int(round(level * len(lab)))
        if n_flip > 0:
            idx = rng.choice(len(lab), n_flip, replace=False)
            for i in idx:
                lab[i] = rng.choice([t for t in types if t != lab[i]])
    elif scenario == "merge":
        pairs = [("Epithelial", "Fibroblast"), ("Tcell", "Myeloid")]
        for t1, t2 in pairs[: int(level)]:
            lab[(lab == t1) | (lab == t2)] = f"{t1}+{t2}"
    elif scenario == "type_swap":
        pairs = [("Epithelial", "Tcell"), ("Fibroblast", "Myeloid")]
        for t1, t2 in pairs[: int(level)]:
            m1, m2 = lab == t1, lab == t2
            lab[m1], lab[m2] = t2, t1
    a.obs["label"] = lab
    return a


def _f(res: dict, key: str) -> float:
    val = res.get(key)
    return np.nan if val is None else float(val)


def _eval_xenium_v08(a: ad.AnnData) -> dict[str, float]:
    b = a.copy()
    for k in list(b.obsp):
        del b.obsp[k]
    rm = sam.residual_spatial_autocorrelation(b, label_key="label", layer="counts", graph=None)
    b2 = a.copy()
    for k in list(b2.obsp):
        del b2.obsp[k]
    sc = sam.spatial_coherence(b2, "label", n_neighs=10)
    b3 = a.copy()
    for k in list(b3.obsp):
        del b3.obsp[k]
    nc = sam.neighbor_contamination(b3, rp.XENIUM_MARKERS, "label", layer="counts", n_perm=19)
    return {
        "resid_moran.residual": _f(rm, "residual_median_moran"),
        "resid_moran.explained": _f(rm, "explained_fraction"),
        "resid_moran.raw": _f(rm, "raw_median_moran"),
        "sc.mean_coherence": _f(sc, "mean_coherence"),
        "sc.pas": _f(sc, "pas"),
        "sc.chaos": _f(sc, "chaos"),
        "nc.mean": _f(nc, "mean_contamination"),
    }


def _load_yerly_split() -> tuple[ad.AnnData, ad.AnnData]:
    full = rp.load_yerly()
    ref_pats = {"BCC1", "BCC2", "BCC3"}
    is_ref = full.obs["patient_bcc"].isin(ref_pats).to_numpy()
    ref = full[is_ref].copy()
    ref.obs["cell_type"] = ref.obs["truth"].astype(str)
    qry = full[~is_ref].copy()
    return ref, qry


def _eval_yerly_v08(q: ad.AnnData, ref: ad.AnnData) -> dict[str, float]:
    rel = sam.relative_pairwise_ct_expression(
        q, ref, label_key="label", ref_label_key="cell_type", layer="counts"
    )
    nmp_s = sam.nmp(q, ref, "label", "cell_type", score="salas", n_permutations=10, layer="counts")
    coex = sam.coexpression_similarity(
        q, ref, label_key="label", ref_label_key="cell_type", layer="counts"
    )
    isc = sam.inter_sample_consistency(q, "label", "patient_bcc", layer="counts")
    ext = sam.external_scores(y_true=q.obs["truth"], y_pred=q.obs["label"])
    return {
        "rel_pw_ct": _f(rel, "relative_pairwise_ct_expression"),
        "nmp.salas": _f(nmp_s, "nmp"),
        "coex.by_type": _f(coex, "coexpression_similarity"),
        "isc.integrated": _f(isc, "integrated"),
        "ext.macro_f1": _f(ext, "macro_f1"),
    }


def _load_merfish_cohort() -> tuple[ad.AnnData, list[tuple[str, str]]]:
    parts = []
    for sec in MERFISH_SECTIONS:
        q_path = MERFISH_BASE / "datasets" / f"{MERFISH_PREFIX}{sec}" / "query.h5ad"
        s_path = MERFISH_BASE / "solutions" / f"{MERFISH_PREFIX}{sec}" / "solution.parquet"
        a = ad.read_h5ad(q_path)
        sol = pd.read_parquet(s_path).set_index("cell_id").loc[a.obs_names]
        a.obs["truth"] = sol["label"].astype(str).to_numpy()
        a.obs["label"] = a.obs["truth"].copy()
        a.obs["sample"] = sol["section"].astype(str).to_numpy()
        a.obs_names = [f"{sec}__{cid}" for cid in a.obs_names]
        parts.append(a)
    cohort = ad.concat(parts, join="inner", merge="same")
    # Pick top 4 abundant brain cell types -> 2 deterministic biological merge pairs.
    counts = cohort.obs["truth"].value_counts()
    top4 = list(counts.index[:4])
    merge_pairs = [(top4[0], top4[1]), (top4[2], top4[3])]
    return cohort, merge_pairs


def _perturb_merfish(
    base: ad.AnnData, scenario: str, level: float, seed: int, merge_pairs: list[tuple[str, str]]
) -> ad.AnnData:
    a = base.copy()
    rng = np.random.default_rng(seed)
    lab = a.obs["truth"].astype(str).to_numpy().copy()
    types = sorted(set(lab))
    if scenario == "label_noise":
        n_flip = int(round(level * len(lab)))
        if n_flip > 0:
            idx = rng.choice(len(lab), n_flip, replace=False)
            for i in idx:
                lab[i] = rng.choice([t for t in types if t != lab[i]])
    elif scenario == "merge":
        for t1, t2 in merge_pairs[: int(level)]:
            lab[(lab == t1) | (lab == t2)] = f"{t1}+{t2}"
    elif scenario == "coord_shuffle":
        if level > 0:
            xy = a.obsm["spatial"].copy()
            samples = a.obs["sample"].astype(str).to_numpy()
            for s in pd.unique(samples):
                idx = np.flatnonzero(samples == s)
                n_shuf = int(round(level * len(idx)))
                if n_shuf > 1:
                    pick = rng.choice(idx, n_shuf, replace=False)
                    xy[pick] = xy[rng.permutation(pick)]
            a.obsm["spatial"] = xy
    elif scenario == "nuisance_shuffle":
        if level > 0:
            order = rng.permutation(a.n_obs)
            a = a[order].copy()
            lab = a.obs["label"].astype(str).to_numpy()
    a.obs["label"] = lab
    return a


def _eval_merfish_v08(a: ad.AnnData) -> dict[str, float]:
    org = sam.spatial_organization_consistency(a, label_key="label", sample_key="sample")
    # Also run single-section residual Moran vs spatial_coherence (deterministic across nuisance_shuffle).
    sec0 = sorted(a.obs["sample"].astype(str).unique())[0]
    s0 = a[a.obs["sample"].astype(str) == sec0].copy()
    for k in list(s0.obsp):
        del s0.obsp[k]
    rm = sam.residual_spatial_autocorrelation(s0, label_key="label", layer=None, graph=None, exclude_controls=False)
    s0_sc = a[a.obs["sample"].astype(str) == sec0].copy()
    for k in list(s0_sc.obsp):
        del s0_sc.obsp[k]
    sc = sam.spatial_coherence(s0_sc, "label", n_neighs=6)
    return {
        "spatial_org.corr": _f(org, "mean_corr"),
        "spatial_org.spearman": _f(org, "mean_spearman"),
        "spatial_org.offdiag": _f(org, "offdiag_mean_corr"),
        "resid_moran.residual": _f(rm, "residual_median_moran"),
        "resid_moran.explained": _f(rm, "explained_fraction"),
        "sc.mean_coherence": _f(sc, "mean_coherence"),
        "sc.pas": _f(sc, "pas"),
        "sc.chaos": _f(sc, "chaos"),
    }


def run_real_part(seeds: int = 3) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Run real-data perturbation suites on Xenium, Yerly BCC split, and MERFISH Bregma cohort."""
    rows = []

    xenium = rp.load_xenium()
    xen_scenarios = {
        "label_noise": [0.0, 0.05, 0.1, 0.25, 0.5],
        "merge": [0, 1, 2],
        "type_swap": [0, 1, 2],
        "spillover": [0.0, 0.1, 0.2, 0.4],
        "coord_shuffle": [0.0, 0.25, 0.5, 1.0],
        "nuisance_shuffle": [0, 1, 2],
    }
    for scen, levels in xen_scenarios.items():
        for lv in levels:
            for seed in range(seeds):
                a = _perturb_xenium_ext(xenium, scen, lv, seed)
                for m, v in _eval_xenium_v08(a).items():
                    rows.append(("xenium_breast", scen, lv, seed, m, v))

    ref_y, qry_y = _load_yerly_split()
    yerly_scenarios = {
        "label_noise": [0.0, 0.05, 0.1, 0.25, 0.5],
        "merge": [0, 1, 2, 3],
        "over_split": [0, 1, 2, 3],
        "type_swap": [0, 1, 2, 3],
        "spillover": [0.0, 0.1, 0.2, 0.4],
        "nuisance_shuffle": [0, 1, 2],
    }
    for scen, levels in yerly_scenarios.items():
        for lv in levels:
            for seed in range(seeds):
                q = rp.perturb_yerly(qry_y, scen, lv, seed)
                for m, v in _eval_yerly_v08(q, ref_y).items():
                    rows.append(("yerly2022_bcc_split", scen, lv, seed, m, v))

    merfish, merge_pairs = _load_merfish_cohort()
    merfish_scenarios = {
        "label_noise": [0.0, 0.05, 0.1, 0.25, 0.5],
        "merge": [0, 1, 2],
        "coord_shuffle": [0.0, 0.5, 1.0],
        "nuisance_shuffle": [0, 1, 2],
    }
    for scen, levels in merfish_scenarios.items():
        for lv in levels:
            for seed in range(seeds):
                m_a = _perturb_merfish(merfish, scen, lv, seed, merge_pairs)
                for m, v in _eval_merfish_v08(m_a).items():
                    rows.append(("merfish_bregma4", scen, lv, seed, m, v))

    val_df = pd.DataFrame(rows, columns=["dataset", "scenario", "level", "seed", "metric", "value"])
    expectations = {
        ("xenium_breast", "label_noise"): {
            "resid_moran.residual", "resid_moran.explained", "sc.mean_coherence", "sc.pas", "sc.chaos"
        },
        ("xenium_breast", "merge"): {"resid_moran.residual", "resid_moran.explained"},
        ("xenium_breast", "type_swap"): {"resid_moran.residual", "resid_moran.explained"},
        ("xenium_breast", "spillover"): {"nc.mean", "resid_moran.residual"},
        ("xenium_breast", "coord_shuffle"): {
            "resid_moran.explained", "resid_moran.raw", "sc.mean_coherence", "sc.pas", "sc.chaos"
        },
        ("yerly2022_bcc_split", "label_noise"): {
            "rel_pw_ct", "nmp.salas", "coex.by_type", "isc.integrated", "ext.macro_f1"
        },
        ("yerly2022_bcc_split", "merge"): {"rel_pw_ct", "coex.by_type", "ext.macro_f1"},
        ("yerly2022_bcc_split", "over_split"): {"isc.integrated", "ext.macro_f1"},
        ("yerly2022_bcc_split", "type_swap"): {
            "rel_pw_ct", "nmp.salas", "coex.by_type", "ext.macro_f1"
        },
        ("yerly2022_bcc_split", "spillover"): {"rel_pw_ct", "nmp.salas", "coex.by_type"},
        ("merfish_bregma4", "label_noise"): {
            "spatial_org.corr",
            "spatial_org.offdiag",
            "resid_moran.residual",
            "resid_moran.explained",
            "sc.mean_coherence",
            "sc.pas",
            "sc.chaos",
        },
        ("merfish_bregma4", "merge"): {"resid_moran.residual", "resid_moran.explained"},
        ("merfish_bregma4", "coord_shuffle"): {
            "spatial_org.corr",
            "spatial_org.offdiag",
            "resid_moran.explained",
            "sc.mean_coherence",
            "sc.pas",
            "sc.chaos",
        },
    }
    verdicts = _score_real(val_df, expectations)
    verdicts.to_csv(OUT_DIR / "real_perturbation.csv", index=False)
    val_df.to_csv(OUT_DIR / "real_values.csv", index=False)

    meta = {
        "xenium_shape": f"{xenium.n_obs} cells x {xenium.n_vars} features",
        "yerly_ref": f"{ref_y.n_obs} cells (BCC1-3), {ref_y.obs['cell_type'].nunique()} types",
        "yerly_qry": f"{qry_y.n_obs} cells (BCC4-7), {qry_y.obs['truth'].nunique()} types",
        "merfish": (
            f"{merfish.n_obs} cells x {merfish.n_vars} genes across "
            f"{len(MERFISH_SECTIONS)} sections; merge_pairs={merge_pairs}"
        ),
    }
    return val_df, verdicts, meta


def _fmt_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, (float, np.floating)):
                cells.append("nan" if pd.isna(v) else f"{v:.4f}")
            else:
                cells.append("nan" if pd.isna(v) else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_summary(
    toy_df: pd.DataFrame,
    toy_sub: pd.DataFrame,
    toy_errs: pd.DataFrame,
    real_verdicts: pd.DataFrame,
    meta: dict[str, str],
) -> None:
    """Write human-readable validation report to ``benchmarks/v08_validation/summary.md``."""
    v08_only = toy_sub[
        toy_sub["metric"].isin(
            [
                "rel_pw_ct",
                "resid_moran.explained",
                "resid_moran.residual",
                "spatial_org.corr",
                "spatial_org.offdiag",
            ]
        )
    ].copy()
    piv = v08_only.pivot(index="metric", columns="scenario", values="verdict").fillna("")
    # Highlight merge & infiltration side-by-side vs existing spatial metrics
    comp = toy_sub[
        toy_sub["scenario"].isin(["merge", "infiltration", "sibling_confusion", "type_swap", "panel_growth"])
    ][["scenario", "metric", "expect", "trend", "verdict", "v0", "v_max", "rho"]].sort_values(
        ["scenario", "metric"]
    )

    real_show = real_verdicts[
        [
            "dataset",
            "scenario",
            "metric",
            "expect",
            "verdict",
            "v0",
            "v_end",
            "spearman_rho",
            "first_level",
        ]
    ].sort_values(["dataset", "scenario", "metric"])

    err_count = len(toy_errs)
    new_err = toy_errs[
        toy_errs["adapter"].isin(
            ["rel_pairwise_ct", "residual_spatial_autocorrelation", "spatial_org_consistency"]
        )
    ]

    # Measure merge & panel_growth nuance explicitly from toy_df
    pg = toy_df[toy_df["metric"] == "rel_pw_ct"]
    pg_lvl = pg[pg["scenario"] == "panel_growth"].groupby("level")["value"].mean().to_dict()
    pg_str = ", ".join(f"{k}: {v:.4f}" for k, v in sorted(pg_lvl.items()))

    lines = [
        "# spatial-anno-metrics v0.8.0 Validation Summary",
        "",
        "Generated by `benchmarks/validate_v08_metrics.py` (3 seeds across synthetic toy sweep and real cohorts).",
        "",
        "## Datasets Tested",
        "- **Synthetic toy sweep**: 15 scenarios x 3 seeds via `benchmarks/synthetic/sweep.py`",
        f"- **Xenium breast crop**: `{meta['xenium_shape']}` (`tests/goldens/signalqc_spatialqm`)",
        f"- **Yerly et al. 2022 BCC (patient split)**: ref `{meta['yerly_ref']}`, query `{meta['yerly_qry']}`",
        f"- **MERFISH mouse brain cohort**: `{meta['merfish']}`",
        f"- **Sweep runtime errors**: {err_count} total ({len(new_err)} in new v0.8.0 adapters; includes edge cases check)",
        "",
        "## 1. Synthetic 15-Scenario Toy Verdicts (v0.8.0 metrics)",
        "",
        _fmt_table(piv.reset_index()),
        "",
        "Note: `sweep.py` expects non-`spatial_coherence`/`signal_qc` metrics to remain `flat` on `coord_shuffle`, so a drop in `spatial_org.*` or `resid_moran.explained` under coordinate scrambling registers as `FALSE-ALARM` in the raw `sweep` table even though spatial drop is the desired behavior.",
        "",
        "## 2. Key Regimes: Merge, Infiltration, Sibling Confusion, Type Swap, Panel Growth",
        "",
        _fmt_table(comp),
        "",
        f"Observed `panel_growth` mean `rel_pw_ct` by extra background genes: `{pg_str}`.",
        "Notice the known methodological caveat from unweighted per-gene L1 normalization (`nf_g`): adding pure-noise/uninformative genes dilutes or perturbs the score unless filtered to marker/expressed panels.",
        "",
        "## 3. Real-Data Perturbation Verdicts",
        "",
        _fmt_table(real_show),
        "",
    ]
    (OUT_DIR / "summary.md").write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Running synthetic 15-scenario toy sweep...")
    toy_df, toy_sub, toy_errs = run_synthetic_part(seeds=3, workers=8)
    print("Running real-data perturbations (Xenium, Yerly BCC split, MERFISH Bregma 4-section)...")
    _, real_verdicts, meta = run_real_part(seeds=3)
    write_summary(toy_df, toy_sub, toy_errs, real_verdicts, meta)
    print(f"Done. Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
