"""Run every public metric over the synthetic scenarios and say which ones behave.

    python benchmarks/synthetic/sweep.py [--seeds 3] [--workers 24] [--out benchmarks/synthetic/results]

For each scenario x strength x seed it simulates a section (``sim.py``), runs every metric, and scores
each metric's trajectory against what it SHOULD do: move with the damage (labels wrong / data degraded),
or stay flat (negative controls: nothing is wrong). Writes ``values.csv`` (raw values, with the number of
UserWarnings the metric's adapter raised in ``n_user_warnings``), ``verdicts.csv`` / ``matrix.md``
(verdicts), ``edge_cases.md``, ``check_*.csv`` (ISC regimes, conformal coverage, panel resolvability,
deconvolution), ``errors.csv`` (every adapter error, edge cases included, with its full traceback) and
``warnings.csv`` (every warning message, counted). Extra metrics can be plugged in with
``main(extra_adapters=...)`` without touching this file.
"""
from __future__ import annotations

import argparse
import os
import sys
import traceback
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):  # one thread per worker
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import sim  # noqa: E402

import spatial_anno_metrics as sam  # noqa: E402

# --------------------------------------------------------------------------- #
# metric adapters: fn(adata, ref) -> {name: value}
# METRICS[name] = (orientation, families): +1 higher = better, -1 lower = better, 0 descriptive
# --------------------------------------------------------------------------- #
METRICS: dict[str, tuple[int, set]] = {}
ADAPTERS: dict = {}


def adapter(name, metrics):
    def deco(fn):
        ADAPTERS[name] = fn
        for m, (o, fam) in metrics.items():
            METRICS[m] = (o, set(fam))
        return fn
    return deco


def _fresh_graph(a):
    for k in list(a.obsp.keys()):
        del a.obsp[k]


L, SP, PUR, SIG, TR, REF = "label", "spatial", "purity", "signal", "truth", "ref"


@adapter("internal_validity", {
    "iv.silhouette": (1, [L]), "iv.neighborhood_purity": (1, [L]), "iv.orbital_medoid": (1, [L]),
    "iv.ward_propmatch": (1, [L]), "iv.avg_similarity": (1, [L]), "iv.davies_bouldin": (-1, [L]),
    "iv.calinski_harabasz": (1, [L]), "iv.integrated": (1, [L])})
def _iv(a, ref):
    r = sam.internal_validity(a, "label", embedding="X_pca")
    return {f"iv.{k}": r.get(k, np.nan) for k in ["silhouette", "neighborhood_purity", "orbital_medoid", "ward_propmatch",
                                       "avg_similarity", "davies_bouldin", "calinski_harabasz", "integrated"]}


@adapter("inter_sample_consistency", {"isc.iscg": (1, [L, "isc"]), "isc.iscl": (1, [L, "isc"]),
                                      "isc.integrated": (1, [L, "isc"])})
def _isc(a, ref):
    r = sam.inter_sample_consistency(a, "label", "sample", layer="counts")
    return {f"isc.{k}": r.get(k, np.nan) for k in ["iscg", "iscl", "integrated"]}


@adapter("marker_program_fidelity", {"fid.mean_auc": (1, [L, "marker"]), "fid.mean_cohens_d": (1, [L, "marker"])})
def _fid(a, ref):
    r = sam.marker_program_fidelity(a, "label", a.uns["marker_sets"])
    return {"fid.mean_auc": r.get("mean_auc"), "fid.mean_cohens_d": r.get("mean_cohens_d")}


@adapter("marker_gene_overlap", {"mgo.mean_overlap": (1, [L, "marker"])})
def _mgo(a, ref):
    return {"mgo.mean_overlap": sam.marker_gene_overlap(a, "label", a.uns["marker_sets"], n_top=10).get("mean_overlap")}


@adapter("purity", {"crisp_purity": (1, [PUR]), "mecr": (-1, [PUR]), "pmp": (1, [L, PUR])})
def _pur(a, ref):
    ms = a.uns["marker_sets"]
    return {"crisp_purity": sam.crisp_purity(a, ms).get("crisp_purity"), "mecr": sam.mecr(a, ms).get("mecr"),
            "pmp": sam.pmp(a, ms, "label").get("mean_pmp")}


@adapter("nmp", {"nmp": (-1, [L, PUR, REF])})
def _nmp(a, ref):
    return {"nmp": sam.nmp(a, ref, "label", "cell_type").get("mean_nmp")}


@adapter("ncp", {"ncp": (1, [PUR, REF])})
def _ncp(a, ref):
    return {"ncp": sam.ncp(a, ref).get("ncp", np.nan)}


@adapter("neighbor_contamination", {"neighbor_contamination": (-1, [L, PUR, SP])})
def _nc(a, ref):
    _fresh_graph(a)
    r = sam.neighbor_contamination(a, a.uns["marker_sets"], "label")
    return {"neighbor_contamination": r.get("mean_contamination")}


@adapter("spatial_coherence", {"sc.mean_coherence": (1, [L, SP]), "sc.pas": (-1, [L, SP])})
def _sc(a, ref):
    _fresh_graph(a)
    r = sam.spatial_coherence(a, "label")
    return {"sc.mean_coherence": r.get("mean_coherence"), "sc.pas": r.get("pas")}


@adapter("neighborhood_sanity", {"ns.n_suspicious": (-1, [L, SP])})
def _ns(a, ref):
    _fresh_graph(a)
    r = sam.neighborhood_sanity(a, "label")
    return {"ns.n_suspicious": len(r["suspicious"]) if r["status"] == "ok" else np.nan}


@adapter("per_cell_confidence", {"conf.mean_margin": (1, ["marker"]), "conf.mean_entropy": (-1, ["marker"])})
def _conf(a, ref):
    r = sam.per_cell_confidence(a, a.uns["marker_sets"])
    return {"conf.mean_margin": r.get("mean_margin"), "conf.mean_entropy": r.get("mean_entropy")}


def marker_argmax(a):
    """The label_stability annotator: argmax of mean log-normalised marker expression."""
    X = a.layers["counts"] if "counts" in a.layers else a.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    Ln = np.log1p(X / lib * 1e4)
    col = {g: i for i, g in enumerate(a.var_names)}
    ms = a.uns["marker_sets"]
    S = np.column_stack([Ln[:, [col[g] for g in ms[t] if g in col]].mean(1) for t in ms])
    return np.array(list(ms))[S.argmax(1)]


@adapter("label_stability", {"stab.flip_rate": (-1, ["marker"])})
def _stab(a, ref):
    return {"stab.flip_rate": sam.label_stability(a, marker_argmax, reps=3)["mean_flip_rate"]}


@adapter("signal_qc", {"sig.moran_gap": (1, [SIG, SP]), "sig.frac_above_ctrl": (1, [SIG, SP]),
                       "sig.snr": (1, [SIG]), "sig.zero_fraction": (-1, [SIG]),
                       "sig.detection_entropy": (0, [SIG]), "sig.tx_per_area": (1, [SIG])})
def _sig(a, ref):
    _fresh_graph(a)
    m = sam.moran_signal(a)
    return {"sig.moran_gap": m.get("separation", np.nan),
            "sig.frac_above_ctrl": m.get("frac_genes_above_control_p95", np.nan),
            "sig.snr": sam.signal_to_noise(a).get("median_mean_signal_ratio"),
            "sig.zero_fraction": sam.sparsity(a).get("zero_fraction"),
            "sig.detection_entropy": sam.detection_entropy(a).get("median_detection_entropy"),
            "sig.tx_per_area": sam.tx_per_area(a).get("median_tx_per_area")}


@adapter("external", {"ext.macro_f1": (1, [L, TR]), "ext.ari": (1, [L, TR]), "ext.ecs": (1, [L, TR]),
                      "ext.hier_acc": (1, [L, TR]), "ext.comp_l1": (-1, [L, TR]), "ext.avg_bio": (1, [L, TR])})
def _ext(a, ref):
    p, t = a.obs["label"], a.obs["truth"]          # raw columns: a NaN label must reach the label check
    e = sam.external_scores(y_true=t, y_pred=p)
    return {"ext.macro_f1": e["macro_f1"], "ext.ari": e["ari"], "ext.ecs": e["ecs"],
            "ext.hier_acc": sam.hierarchical_accuracy(y_true=t, y_pred=p,
                                                      hierarchy=a.uns["hierarchy"])["hierarchical_accuracy"],
            "ext.comp_l1": sam.composition_accuracy(y_true=t, y_pred=p)["l1"],
            "ext.avg_bio": sam.avg_bio(a, "label", "truth").get("avg_bio")}


def run_adapters(a, ref, adapters):
    """Run every adapter on ``a``. Returns ``(values {adapter: {metric: value}}, n_user_warnings {adapter: n},
    warnings [(adapter, category, message)], errors {adapter: (one-line message, full traceback)})``. Warnings
    are recorded, not shown or ignored: every occurrence counts."""
    vals, nwarn, warned, errs = {}, {}, [], {}
    for name, fn in adapters.items():
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                vals[name] = {k: (np.nan if v is None else float(v)) for k, v in fn(a, ref).items()}
            except Exception as e:  # noqa: BLE001 - a crash is a finding, record it
                errs[name] = (f"{type(e).__name__}: {str(e).splitlines()[0][:120] if str(e) else ''}",
                              traceback.format_exc())
        nwarn[name] = sum(issubclass(w.category, UserWarning) for w in caught)
        warned += [(name, w.category.__name__, str(w.message).splitlines()[0][:200] if str(w.message) else "")
                   for w in caught]
    return vals, nwarn, warned, errs


def _flat(vals):
    return {m: v for metrics in vals.values() for m, v in metrics.items()}


def _reference_for(a):
    w = a.uns["world"]
    return sim.reference(int(w["seed"]), coexpr=float(w["coexpr"]), extra_bg=int(w["extra_bg"]))


def _task(args):
    scen, level, seed, adapters = args
    fn = sim.SCENARIOS[scen][0]
    a = fn(level, seed)
    vals, nwarn, warned, errs = run_adapters(a, _reference_for(a), adapters)
    rows = [(scen, level, seed, m, v, nwarn[ad]) for ad, metrics in vals.items() for m, v in metrics.items()]
    return (rows, [(scen, level, seed, k, one, tb) for k, (one, tb) in errs.items()],
            [("scenarios", *w) for w in warned])


# --------------------------------------------------------------------------- #
# what each metric SHOULD do in each scenario
# --------------------------------------------------------------------------- #
DATA_EXPECT = {   # data-side scenarios: families (or metric names) that should degrade; the rest - no claim
    "spillover": {PUR},
    "depth_loss": {"sig.tx_per_area", "sig.zero_fraction"},
    "ambient": {"sig.snr", "sig.moran_gap", "sig.frac_above_ctrl"},
    "batch": {"isc"},
    "coord_shuffle": {"sc.mean_coherence", "sc.pas", "ns.n_suspicious", "sig.moran_gap", "sig.frac_above_ctrl"},
}


def expectation(metric, scen, metrics=None):
    o, fam = (metrics or METRICS)[metric]
    kind = sim.SCENARIOS[scen][2]
    if o == 0:
        return None
    if kind == "label":
        return "down" if L in fam else "flat"
    if kind == "null":
        return "flat" if (L in fam or PUR in fam) else None
    want = DATA_EXPECT.get(scen, set())
    if metric in want or fam & want:
        return "down"
    if scen == "coord_shuffle":
        return "flat"                        # nothing but the coordinates changed
    return None


def score(df, extra_metrics=None):
    """Per (metric, scenario): oriented change, monotonicity, min detectable level, verdict. A move counts
    when it exceeds 3 x the seed SD at strength 0 AND 5 % of the metric's range over the whole sweep."""
    from scipy.stats import spearmanr
    metrics = {**METRICS, **(extra_metrics or {})}
    span = df.groupby("metric")["value"].agg(lambda v: np.nanmax(v) - np.nanmin(v) if v.notna().any() else 0)
    rows = []
    for (scen, metric), g in df.groupby(["scenario", "metric"]):
        if metric not in metrics:
            continue
        o = metrics[metric][0] or 1
        levels = sim.SCENARIOS[scen][1]
        m = g.groupby("level")["value"].agg(["mean", "std"]).reindex(levels)
        if m["mean"].isna().all():
            rows.append(dict(scenario=scen, metric=metric, verdict="n/a", detail="all NaN"))
            continue
        v0, sd0 = m["mean"].iloc[0], np.nan_to_num(m["std"].iloc[0])
        bad = -o * (m["mean"] - v0)                          # > 0 = metric reports worse quality
        thr = max(3 * sd0, 0.05 * span[metric], 1e-9)
        ok = ~m["mean"].isna()
        rho = (spearmanr(np.arange(ok.sum()), bad[ok].to_numpy())[0]
               if ok.sum() > 2 and np.nanstd(bad[ok]) > 0 else 0.0)
        worst = bad.iloc[-1] if ok.iloc[-1] else np.nan
        moved = np.nanmax(np.abs(bad.to_numpy())) > thr
        if not moved:
            trend = "flat"
        elif worst > thr and rho >= 0.8:
            trend = "worse"
        elif worst < -thr and rho <= -0.8:
            trend = "better"
        else:
            trend = "erratic"
        first = next((lv for lv, b in zip(levels, bad) if b > thr), None)
        exp = expectation(metric, scen, metrics)
        verdict = {("down", "worse"): "detects", ("down", "flat"): "BLIND", ("down", "better"): "WRONG-WAY",
                   ("down", "erratic"): "erratic", ("flat", "flat"): "stable",
                   ("flat", "worse"): "FALSE-ALARM", ("flat", "better"): "drifts-up",
                   ("flat", "erratic"): "FALSE-ALARM"}.get((exp, trend), trend)
        rows.append(dict(scenario=scen, metric=metric, expect=exp or "-", trend=trend, verdict=verdict,
                         v0=v0, v_max=m["mean"].iloc[-1], rho=rho, threshold=thr,
                         first_level=first if trend == "worse" else None))
    return pd.DataFrame(rows)


SYMBOL = {"detects": "✅", "BLIND": "❌", "WRONG-WAY": "⛔", "erratic": "〰️", "stable": "✔", "FALSE-ALARM": "⚠️",
          "drifts-up": "↗", "worse": "↓", "better": "↑", "flat": "·", "n/a": "∅"}


def matrix_md(sc):
    scen = list(sim.SCENARIOS)
    cell = sc.assign(c=[SYMBOL.get(v, v) + (f" {fl:g}" if v == "detects" and fl is not None and not pd.isna(fl)
                                            else "") for v, fl in zip(sc.verdict, sc.first_level)])
    M = cell.pivot(index="metric", columns="scenario", values="c").reindex(columns=scen).fillna("")
    head = "| metric | " + " | ".join(f"{s}<br>*{sim.SCENARIOS[s][2]}*" for s in scen) + " |"
    lines = [head, "|" + "---|" * (len(scen) + 1)]
    lines += [f"| `{m}` | " + " | ".join(M.loc[m]) + " |" for m in M.index]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# dedicated checks
# --------------------------------------------------------------------------- #
def isc_regimes(seeds):
    rows = []
    cases = [("clean", sim.label_noise, 0), ("label_noise 20%", sim.label_noise, 0.2),
             ("over_split 1 (Epi halves)", sim.over_split, 1), ("merge 1 (Fibro+Epi)", sim.merge, 1),
             ("merge 3 (incl. T_CD4+T_CD8)", sim.merge, 3), ("batch sd 1.0", sim.batch, 1.0)]
    for name, fn, s in cases:
        for seed in range(seeds):
            r = sam.inter_sample_consistency(fn(s, seed), "label", "sample", layer="counts")
            for t, reg in r.get("regime", {}).items():
                rows.append(dict(case=name, seed=seed, type=t, regime=reg,
                                 iscg=r["per_type_iscg"][t], iscl=r["per_type_iscl"][t]))
    return pd.DataFrame(rows)


def _lognorm(X):
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    return np.log1p(X / lib * 1e4)


def conformal_check(seeds):
    """Coverage of 90 % conformal sets: calibrated on the reference (not exchangeable with the section)
    vs on a labelled 30 % of the section; marginal vs Mondrian vs clustered Mondrian on the rare type; query
    depth shift; calibration error (ECE) of the classifier on the section. ``held_out`` = "B": B is left out
    of the reference, so the section's B cells are a truly novel type - do empty sets and conformal novelty
    p-values (``nov_*``: PCA reconstruction error, BH 10 %) single them out, or also flag known cells that
    only lost depth? ``nov_*`` use the reference with each cell thinned to a depth drawn from the section's
    per-cell depths; ``nov_flag_known_unmatched`` the reference as it is (400 counts per cell vs 150)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    def thin(X, depths, rng):                      # each reference cell thinned to a depth drawn from the section's
        X = np.asarray(X.todense() if hasattr(X, "todense") else X).astype(int)
        target = rng.choice(depths, len(X))[:, None]
        return rng.binomial(X, np.clip(target / np.maximum(X.sum(1, keepdims=True), 1), 0, 1)).astype(float)

    rows = []
    for seed in range(seeds):
        for held_out in ["", "B"]:
            for depth_keep in [1.0, 0.5, 0.25]:
                for rare in [0.15, 0.02]:
                    if (depth_keep != 1.0 or held_out) and rare != 0.15:
                        continue
                    a = sim.rare_type(rare, seed) if depth_keep == 1.0 else sim.depth_loss(1 - depth_keep, seed)
                    ref = sim.reference(seed)
                    genes = list(ref.var_names)
                    kr = ref.obs["cell_type"].to_numpy() != held_out
                    Xr = ref[:, genes].layers["counts"][kr]
                    clf = LogisticRegression(max_iter=500, C=0.3).fit(_lognorm(Xr), ref.obs["cell_type"][kr])
                    cls = clf.classes_
                    Xq = a[:, genes].layers["counts"]
                    Pq = clf.predict_proba(_lognorm(Xq))
                    yq = a.obs["truth"].to_numpy()
                    novel = yq == held_out
                    rng = np.random.default_rng(seed)
                    cal = rng.random(a.n_obs) < 0.3
                    ref_cal = sim.reference(seed + 100)
                    kc = ref_cal.obs["cell_type"].to_numpy() != held_out
                    Pr = clf.predict_proba(_lognorm(ref_cal[:, genes].layers["counts"][kc]))
                    ce = sam.conformal.calibration_error(Pq[~novel], cls, yq[~novel])
                    Xm = _lognorm(thin(Xr, np.asarray(Xq.sum(1)).ravel(), np.random.default_rng(seed + 200)))
                    rn = sam.conformal.conformal_novelty(_lognorm(Xq), Xm, method="pca", seed=seed)
                    nv, pv = rn["novel"], rn["pvalues"]
                    nv_raw = sam.conformal.conformal_novelty(_lognorm(Xq), _lognorm(Xr), method="pca", seed=seed)["novel"]
                    nov = dict(ece=ce["ece"], classwise_ece=ce["classwise_ece"],
                               nov_auc=roc_auc_score(novel, -pv) if held_out else np.nan,
                               nov_flag_known=nv[~novel].mean(), nov_flag_novel=nv[novel].mean() if held_out else np.nan,
                               nov_fdp=(~novel[nv]).mean() if held_out and nv.any() else np.nan,
                               nov_flag_known_unmatched=nv_raw[~novel].mean())
                    cals = {"reference": (Pr, ref_cal.obs["cell_type"].to_numpy()[kc], Pq, yq)}
                    if not held_out:                       # a labelled section subset cannot hold an unknown type
                        cals["section 30%"] = (Pq[cal], yq[cal], Pq[~cal], yq[~cal])
                    for cal_name, (Pc, yc, Pqq, yqq) in cals.items():
                        for cc in (False, True, "clustered"):
                            r = sam.conformal_prediction_sets(Pqq, cls, Pc, yc, alpha=0.1, class_conditional=cc,
                                                              y_query=yqq)
                            known = yqq != held_out
                            row = dict(seed=seed, query_depth=depth_keep, T_CD8_share=rare, calibration=cal_name,
                                       mondrian=cc,
                                       coverage=np.mean([yqq[i] in r["sets"][i] for i in np.flatnonzero(known)]),
                                       T_CD8_coverage=r["per_type_coverage"].get("T_CD8", np.nan),
                                       mean_set_size=r["summary"]["mean_set_size"], held_out=held_out,
                                       covgap=np.mean([abs(v - 0.9) for t, v in r["per_type_coverage"].items()
                                                       if t != held_out]),
                                       pct_empty_known=np.mean(r["set_size"][known] == 0),
                                       pct_empty_novel=np.mean(r["set_size"][~known] == 0) if held_out else np.nan)
                            rows.append({**row, **nov})
    return pd.DataFrame(rows)


def panel_check(seeds):
    """panel_resolvability vs panel content (fraction of private markers kept) x target depth."""
    rows = []
    ms = sim.marker_sets()
    private = [g for t, gs in ms.items() for g in gs if g not in sim.T_LINEAGE]
    for seed in range(seeds):
        ref = sim.reference(seed, n_per_type=200)
        rng = np.random.default_rng(seed)
        for keep in [1.0, 0.5, 0.25, 0.0]:
            drop = set(rng.choice(private, int(round((1 - keep) * len(private))), replace=False))
            panel = [g for g in ref.var_names if g not in drop]
            for depth in [200, 50, 20, 5]:
                r = sam.panel_resolvability(ref, "cell_type", panel, target_depth=depth, seed=seed)
                cd8 = r["per_type"].get("T_CD8", {})
                rows.append(dict(seed=seed, markers_kept=keep, target_depth=depth, macro_f1=r.get("macro_f1"),
                                 frac_resolvable=r.get("frac_resolvable"), T_CD8_f1=cd8.get("f1"),
                                 T_CD8_confused_with=cd8.get("confused_with")))
    return pd.DataFrame(rows)


def deconv_check(seeds):
    """deconvolution_metrics on 25-cell pseudo-spots: prediction = (1-s) truth + s Dirichlet noise."""
    rows = []
    for seed in range(seeds):
        a = sim.simulate(seed)
        xy = a.obsm["spatial"]
        spot = (xy[:, 0] // 5).astype(int) * 1000 + (xy[:, 1] // 5).astype(int)
        T = pd.crosstab(spot, a.obs["truth"].to_numpy(), normalize="index")
        rng = np.random.default_rng(seed)
        noise = pd.DataFrame(rng.dirichlet(np.ones(T.shape[1]), T.shape[0]), index=T.index, columns=T.columns)
        for s in [0, 0.1, 0.25, 0.5, 1.0]:
            r = sam.deconvolution_metrics(true_prop=T, pred_prop=(1 - s) * T + s * noise)
            rows.append(dict(seed=seed, noise=s, r2=r["r2"], jsd=r["jsd"], rmse=r["rmse"]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# edge cases: does each metric return a sane value, NaN, or crash?
# --------------------------------------------------------------------------- #
def _edge(name):
    a = sim.simulate(0, n_cells=600)
    ms = {k: list(v) for k, v in a.uns["marker_sets"].items()}
    lab = a.obs["label"].astype(str).to_numpy().copy()
    if name == "one_label":
        lab[:] = "Epi"
    elif name == "singleton_type":
        b = np.flatnonzero(lab == "B")
        lab[b[1:]] = "Myeloid"
    elif name == "two_types":
        lab = np.where(np.isin(lab, ["Fibro", "Epi"]), "Stroma", "Immune")
    elif name == "empty_marker_set":
        ms["B"] = []
    elif name == "markers_off_panel":
        ms = {k: [g + "_absent" for g in v] for k, v in ms.items()}
    elif name == "zero_count_cells":
        X = np.asarray(a.layers["counts"].todense())
        X[::10] = 0
        a.layers["counts"] = X
        a = sim.finish(a)
    elif name == "unused_category":
        a.obs["label"] = pd.Categorical(lab, categories=sorted(set(lab)) + ["Ghost"])
    elif name == "integer_labels":
        code = {t: i for i, t in enumerate(sim.TYPES)}
        a.obs["label"] = [code[x] for x in lab]
        a.obs["truth"] = [code[x] for x in a.obs["truth"]]
        a.uns["hierarchy"] = {str(code[k]): v for k, v in sim.HIERARCHY.items()}
        ms = {str(code[k]): v for k, v in ms.items()}
    elif name == "nan_labels":
        a.obs["label"] = pd.Series(lab, index=a.obs_names).where(np.arange(a.n_obs) % 20 != 0)
    elif name == "no_control_probes":
        a = a[:, ~a.var["control"].to_numpy()].copy()
    elif name == "tiny_30_cells":
        a = a[np.random.default_rng(0).choice(a.n_obs, 30, replace=False)].copy()
        lab = a.obs["label"].astype(str).to_numpy()
    elif name == "duplicate_coords":
        a.obsm["spatial"] = np.round(a.obsm["spatial"] / 8) * 8
    elif name == "all_zero_marker_gene":
        X = np.asarray(a.layers["counts"].todense())
        X[:, list(a.var_names).index("B_m0")] = 0
        a.layers["counts"] = X
        a = sim.finish(a)
    elif name == "dense_matrix":
        a.layers["counts"] = np.asarray(a.layers["counts"].todense())
        a.X = a.layers["counts"].copy()
    elif name == "unicode_labels":
        ren = {t: f"{t} (cellule é+)" for t in sim.TYPES}
        lab = np.array([ren[x] for x in lab])
        a.obs["truth"] = [ren[x] for x in a.obs["truth"]]
        a.uns["hierarchy"] = {ren[k]: v for k, v in sim.HIERARCHY.items()}
        ms = {ren[k]: v for k, v in ms.items()}
    elif name == "extra_marker_type":
        ms["Ghost"] = ["NegCtrl0", "NegCtrl1"]       # an absent type's markers: on panel, ambient level only
    elif name == "no_spatial_coords":
        del a.obsm["spatial"]
    elif name == "lognorm_X_no_counts":
        a.X = _lognorm(a.layers["counts"])
        del a.layers["counts"]
    if name not in ("unused_category", "integer_labels", "nan_labels"):
        a.obs["label"] = lab
    a.uns["marker_sets"] = ms
    return a


EDGES = ["one_label", "singleton_type", "two_types", "empty_marker_set", "markers_off_panel", "zero_count_cells",
         "unused_category", "integer_labels", "nan_labels", "no_control_probes", "tiny_30_cells",
         "duplicate_coords", "all_zero_marker_gene", "dense_matrix", "unicode_labels", "extra_marker_type",
         "no_spatial_coords", "lognorm_X_no_counts"]


def _edge_task(args):
    name, adapters = args
    a = _edge(name)
    vals, _, warned, errs = run_adapters(a, sim.reference(0), adapters)
    return name, vals, errs, [(f"edge:{name}", *w) for w in warned]


def edge_md(results, base):
    """Per edge case x adapter: ok / NaN / error / differs-from-clean (for cases that should not change values)."""
    same_as_clean = {"unused_category", "dense_matrix", "unicode_labels", "extra_marker_type", "integer_labels",
                     "lognorm_X_no_counts"}   # here ≠ = values silently computed on log data
    names = list(OWNER)
    lines = ["| edge case | " + " | ".join(names) + " |", "|" + "---|" * (len(names) + 1)]
    for name, vals, errs in results:                 # vals: {metric: value}; errs: {adapter: one-line message}
        cells = []
        for ad_name in names:
            if ad_name in errs:
                cells.append(f"💥 {errs[ad_name].split(':')[0]}")
                continue
            ms = [m for m in vals if _owner(m) == ad_name]
            nan = [m for m in ms if pd.isna(vals[m])]
            diff = [m for m in ms if name in same_as_clean and not pd.isna(vals[m]) and not pd.isna(base.get(m))
                    and abs(vals[m] - base[m]) > 1e-6 * max(1, abs(base[m]))]
            cells.append("≠ " + ",".join(x.split(".")[-1] for x in diff) if diff
                         else ("NaN " + ",".join(x.split(".")[-1] for x in nan) if nan else "ok"))
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    errors = [f"- **{n}** / `{k}`: {e}" for n, _, errs in results for k, e in errs.items()]
    return "\n".join(lines) + "\n\nErrors:\n\n" + ("\n".join(errors) or "none")


def _owner(metric):
    for ad_name, fam in OWNER.items():
        if metric in fam:
            return ad_name
    return None


OWNER: dict[str, set] = {}   # adapter -> the metrics it returned on the clean section (empty if it failed)


# --------------------------------------------------------------------------- #
def main(argv=None, extra_adapters=None, extra_metrics=None):
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--workers", type=int, default=24)
    p.add_argument("--out", default=str(Path(__file__).parent / "results"))
    p.add_argument("--only", nargs="*", default=None, help="subset of scenarios")
    args = p.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adapters = {**ADAPTERS, **(extra_adapters or {})}
    scen = args.only or list(sim.SCENARIOS)
    tasks = [(s, lv, seed, adapters) for s in scen for lv in sim.SCENARIOS[s][1] for seed in range(args.seeds)]
    rows, errs, warned = [], [], Counter()
    with ProcessPoolExecutor(args.workers) as ex:
        for r, e, w in ex.map(_task, tasks):
            rows += r
            errs += e
            warned.update(w)
        edges = list(ex.map(_edge_task, [(n, adapters) for n in EDGES]))
    df = pd.DataFrame(rows, columns=["scenario", "level", "seed", "metric", "value", "n_user_warnings"])
    df.to_csv(out / "values.csv", index=False)
    sc = score(df, extra_metrics)
    sc.to_csv(out / "verdicts.csv", index=False)
    (out / "matrix.md").write_text(matrix_md(sc) + "\n")
    base_vals, _, base_warned, _ = run_adapters(sim.simulate(0, n_cells=600), sim.reference(0), adapters)
    OWNER.update({name: set(base_vals.get(name, {})) for name in adapters})
    warned.update(("edge:base", *w) for w in base_warned)
    for name, vals, e, w in edges:
        errs += [("edge", name, 0, k, one, tb) for k, (one, tb) in e.items()]
        warned.update(w)
    pd.DataFrame(errs, columns=["scenario", "level", "seed", "adapter", "error", "traceback"]).to_csv(
        out / "errors.csv", index=False)
    pd.DataFrame([(*k, n) for k, n in sorted(warned.items())],
                 columns=["source", "adapter", "category", "message", "count"]).to_csv(out / "warnings.csv", index=False)
    edge_rows = [(name, _flat(vals), {k: one for k, (one, _) in e.items()}) for name, vals, e, _ in edges]
    (out / "edge_cases.md").write_text(edge_md(edge_rows, _flat(base_vals)) + "\n")
    if extra_adapters is None:
        checks = {"isc_regimes": isc_regimes(args.seeds), "conformal": conformal_check(args.seeds),
                  "panel_resolvability": panel_check(args.seeds), "deconvolution": deconv_check(args.seeds)}
        for k, v in checks.items():
            v.to_csv(out / f"check_{k}.csv", index=False)
    return df, sc


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
