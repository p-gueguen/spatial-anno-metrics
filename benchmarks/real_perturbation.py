"""Real-data perturbation benchmark on the committed Yerly 2022 BCC and Xenium breast datasets.

    python benchmarks/real_perturbation.py [--seeds 3] [--out benchmarks/real_perturbation/results]

Evaluates the metric battery under controlled perturbations of two committed real datasets:
1. **Yerly et al. 2022 human BCC scRNA-seq cohort** (``tests/goldens/isc_yerly2022/``, 2,620 cells x 2,000
   HVGs across 7 patients and 9 cell types): ``label_noise`` (T1), ``over_split`` (T2), ``merge`` (T3),
   ``type_swap`` (T4), ``spillover`` (T5), ``depth_loss`` (T6), and ``nuisance_shuffle`` (T8 negative control).
2. **10x Xenium human breast cancer crop** (``tests/goldens/signalqc_spatialqm/``, 2,000 cells x 541 features:
   313 genes + 228 negative controls, with cell centroids and areas): ``ambient`` (T5), ``depth_loss`` (T6),
   ``coord_shuffle`` (T7), and ``nuisance_shuffle`` (T8 negative control).

Writes ``values.csv``, ``verdicts.csv``, and ``matrix.md`` scored by monotonicity (Spearman rho), minimum
detectable perturbation level, and constant-fit on negative controls (catalog section 7, roadmap #24).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sam_mpl")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import spatial_anno_metrics as sam  # noqa: E402

GOLDENS = ROOT / "tests" / "goldens"
OBO = ROOT / "tests" / "data" / "cl_tiny.obo"

YERLY_MARKERS: dict[str, list[str]] = {
    "Tcell": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD8A"],
    "Bcell": ["CD37", "LTB", "P2RY10"],
    "PlasmaCell": ["JCHAIN", "XBP1", "CD27"],
    "MoMacDC": ["LYZ", "C1QA", "C1QB", "AIF1", "TYROBP", "CD68"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "ACKR1", "PLVAP"],
    "Normal_Kerat": ["KRT14", "KRT5", "KRT15", "DSP", "SFN"],
    "Tumor_KC_Basal": ["KRT17", "TACSTD2", "PTCH1"],
    "Melanocyte": ["MLANA", "PMEL", "TYR", "DCT", "MITF"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}

YERLY_CL: dict[str, str] = {
    "Tcell": "CL:0000084",
    "Bcell": "CL:0000236",
    "MoMacDC": "CL:0000235",
    "Normal_Kerat": "CL:0000312",
    "Fibroblast": "CL:0000057",
}

XENIUM_MARKERS: dict[str, list[str]] = {
    "Epithelial": ["EPCAM", "KRT7", "KRT8", "ERBB2", "CDH1"],
    "Tcell": ["CD3D", "CD3E", "CD8A", "IL7R", "TRAC"],
    "Bcell": ["MS4A1", "CD79A", "BANK1", "MZB1"],
    "Myeloid": ["CD68", "CD163", "LYZ", "C1QA", "ITGAX"],
    "Fibroblast": ["LUM", "POSTN"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR"],
}


def _compute_pca(X: sp.csr_matrix, n_components: int = 20) -> np.ndarray:
    """Compute log1p(1e4 * counts / total) PCA embedding from a CSR count matrix."""
    from sklearn.decomposition import PCA

    tot = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(1e4, tot, out=np.zeros_like(tot), where=tot > 0)
    L = X.copy()
    L.data *= np.repeat(scale, np.diff(L.indptr))
    np.log1p(L.data, out=L.data)
    return PCA(n_components=n_components, svd_solver="arpack", random_state=0).fit_transform(L)


def load_yerly() -> ad.AnnData:
    """Load the Yerly et al. 2022 BCC cohort with PCA embedding and raw counts."""
    gdir = GOLDENS / "isc_yerly2022"
    with np.load(gdir / "counts.npz", allow_pickle=False) as z:
        X = sp.csr_matrix((z["data"].astype(np.float64), z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        genes = z["genes"].astype(str)
    obs = pd.read_csv(gdir / "obs.csv", dtype=str).set_index("cell_id")
    a = ad.AnnData(X=X.copy(), obs=obs)
    a.var_names = genes
    a.layers["counts"] = X.copy()
    a.obsm["X_pca"] = _compute_pca(X)
    a.obs["truth"] = a.obs["annotation"].to_numpy()
    a.obs["label"] = a.obs["annotation"].to_numpy()
    return a


def load_xenium() -> ad.AnnData:
    """Load the 10x Xenium human breast cancer crop with spatial coordinates and marker-derived labels."""
    gdir = GOLDENS / "signalqc_spatialqm"
    with np.load(gdir / "counts.npz", allow_pickle=False) as c:
        X = sp.csr_matrix((c["data"].astype(np.float64), c["indices"], c["indptr"]), shape=tuple(c["shape"]))
        names = c["feature_name"].astype(str)
        ftypes = c["feature_type"].astype(str)
    cells = pd.read_csv(gdir / "cells.csv")
    a = ad.AnnData(X=X.copy(), var=pd.DataFrame({"feature_type": ftypes}, index=names))
    a.layers["counts"] = X.copy()
    a.obs_names = cells["cell_id"].astype(str).to_numpy()
    a.obs["cell_area"] = cells["cell_area"].to_numpy()
    a.obs["total_counts"] = cells["total_counts"].to_numpy()
    a.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()
    col = {g: i for i, g in enumerate(names)}
    scores = np.column_stack([
        np.asarray(X[:, [col[g] for g in gs]].mean(axis=1)).ravel() for gs in XENIUM_MARKERS.values()
    ])
    types = np.array(list(XENIUM_MARKERS))
    a.obs["truth"] = types[scores.argmax(axis=1)]
    a.obs["label"] = a.obs["truth"].copy()
    return a


def perturb_yerly(base: ad.AnnData, scenario: str, level: float, seed: int) -> ad.AnnData:
    """Apply perturbation ``scenario`` at ``level`` to a copy of the Yerly 2022 cohort."""
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
    elif scenario == "over_split":
        order = ["Tcell", "Tumor_KC_Basal", "Fibroblast", "Endothelial"]
        for t in order[:int(level)]:
            m = np.flatnonzero(lab == t)
            half = rng.choice(m, len(m) // 2, replace=False)
            lab[half] = f"{t}_split2"
    elif scenario == "merge":
        pairs = [("Tcell", "Bcell"), ("Fibroblast", "Endothelial"), ("Normal_Kerat", "Tumor_KC_Basal")]
        for t1, t2 in pairs[:int(level)]:
            lab[(lab == t1) | (lab == t2)] = f"{t1}+{t2}"
    elif scenario == "type_swap":
        pairs = [("Tcell", "Normal_Kerat"), ("Fibroblast", "MoMacDC"), ("Endothelial", "Tumor_KC_Basal")]
        for t1, t2 in pairs[:int(level)]:
            m1, m2 = lab == t1, lab == t2
            lab[m1], lab[m2] = t2, t1
    elif scenario == "spillover":
        if level > 0:
            X = a.layers["counts"].toarray()
            perm = rng.permutation(len(lab))
            mixed = rng.binomial(X.astype(int), 1.0 - level) + rng.binomial(X[perm].astype(int), level)
            a.layers["counts"] = sp.csr_matrix(mixed.astype(np.float64))
            a.X = a.layers["counts"].copy()
            a.obsm["X_pca"] = _compute_pca(a.layers["counts"])
    elif scenario == "depth_loss":
        if level > 0:
            X = a.layers["counts"].tocsr().copy()
            X.data = rng.binomial(np.rint(X.data).astype(int), 1.0 - level).astype(np.float64)
            X.eliminate_zeros()
            a.layers["counts"] = X
            a.X = X.copy()
            a.obsm["X_pca"] = _compute_pca(a.layers["counts"])
    elif scenario == "nuisance_shuffle":
        if level > 0:
            order = rng.permutation(a.n_obs)
            a = a[order].copy()
            lab = a.obs["label"].astype(str).to_numpy()
    a.obs["label"] = lab
    return a


def perturb_xenium(base: ad.AnnData, scenario: str, level: float, seed: int) -> ad.AnnData:
    """Apply spatial/signal perturbation ``scenario`` at ``level`` to a copy of the Xenium crop."""
    from scipy.spatial import cKDTree

    a = base.copy()
    rng = np.random.default_rng(seed)
    if scenario == "ambient":
        if level > 0:
            X = a.layers["counts"].toarray()
            noise = rng.poisson(level, size=X.shape).astype(np.float64)
            a.layers["counts"] = sp.csr_matrix(X + noise)
            a.X = a.layers["counts"].copy()
            a.obs["total_counts"] = np.asarray(a.layers["counts"].sum(axis=1)).ravel()
    elif scenario == "spillover":
        if level > 0:
            X = a.layers["counts"].toarray()
            _, nn = cKDTree(a.obsm["spatial"]).query(a.obsm["spatial"], k=2)
            nbr = nn[:, 1]
            mixed = rng.binomial(X.astype(int), 1.0 - level) + rng.binomial(X[nbr].astype(int), level)
            a.layers["counts"] = sp.csr_matrix(mixed.astype(np.float64))
            a.X = a.layers["counts"].copy()
            a.obs["total_counts"] = np.asarray(a.layers["counts"].sum(axis=1)).ravel()
    elif scenario == "depth_loss":
        if level > 0:
            X = a.layers["counts"].tocsr().copy()
            X.data = rng.binomial(np.rint(X.data).astype(int), 1.0 - level).astype(np.float64)
            X.eliminate_zeros()
            a.layers["counts"] = X
            a.X = X.copy()
            a.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    elif scenario == "coord_shuffle":
        if level > 0:
            xy = a.obsm["spatial"].copy()
            n_shuf = int(round(level * a.n_obs))
            idx = rng.choice(a.n_obs, n_shuf, replace=False)
            xy[idx] = xy[rng.permutation(idx)]
            a.obsm["spatial"] = xy
    elif scenario == "nuisance_shuffle":
        if level > 0:
            order = rng.permutation(a.n_obs)
            a = a[order].copy()
    return a


def eval_yerly(a: ad.AnnData, ref: ad.AnnData) -> dict[str, float]:
    """Run cohort, internal, ontology, fidelity, purity, and external metrics on Yerly 2022."""
    iv = sam.internal_validity(a, "label", embedding="X_pca", subsample=1000)
    isc = sam.inter_sample_consistency(a, "label", "patient_bcc", layer="counts")
    fid = sam.marker_program_fidelity(a, "label", YERLY_MARKERS, layer="counts")
    pmp = sam.pmp(a, YERLY_MARKERS, "label", layer="counts")
    crisp = sam.crisp_purity(a, YERLY_MARKERS, layer="counts", n_null=0)
    mecr = sam.mecr(a, YERLY_MARKERS, layer="counts", n_null=0)
    nmp_s = sam.nmp(a, ref, "label", "truth", score="salas", n_permutations=10, layer="counts")
    ncp_m = sam.ncp(a, ref, score="moldia", layer="counts")
    onto = sam.ontology_consistency(a, "label", OBO, YERLY_CL, embedding="X_pca", n_perm=99)
    onto_rwr = sam.ontology_consistency(a, "label", OBO, YERLY_CL, embedding="X_pca", distance="rwr", n_perm=99)
    ext = sam.external_scores(y_true=a.obs["truth"], y_pred=a.obs["label"])
    return {
        "iv.integrated": iv.get("integrated", np.nan),
        "isc.iscg": isc.get("iscg", np.nan),
        "isc.iscl": isc.get("iscl", np.nan),
        "isc.integrated": isc.get("integrated", np.nan),
        "fid.mean_auc": fid.get("mean_auc", np.nan),
        "pmp": pmp.get("mean_pmp", np.nan),
        "crisp_purity": crisp.get("crisp_purity", np.nan),
        "mecr": mecr.get("mecr", np.nan),
        "nmp.salas": nmp_s.get("nmp", np.nan),
        "ncp.moldia": ncp_m.get("ncp", np.nan),
        "onto.rho": onto.get("rho", np.nan),
        "onto.rwr_rho": onto_rwr.get("rho", np.nan),
        "ext.macro_f1": ext.get("macro_f1", np.nan),
        "ext.ari": ext.get("ari", np.nan),
        "ext.ecs": ext.get("ecs", np.nan),
    }


def eval_xenium(a: ad.AnnData) -> dict[str, float]:
    """Run signal-QC, spatial coherence, neighbor contamination, and purity metrics on Xenium crop."""
    for k in list(a.obsp):
        del a.obsp[k]
    m = sam.moran_signal(a, method="spatialqm", layer="counts")
    snr = sam.signal_to_noise(a, layer="counts")
    sp_res = sam.sparsity(a, layer="counts")
    tpa = sam.tx_per_area(a, layer="counts")
    sc = sam.spatial_coherence(a, "label", n_neighs=10)
    nc = sam.neighbor_contamination(a, XENIUM_MARKERS, "label", layer="counts", n_perm=49)
    crisp = sam.crisp_purity(a, XENIUM_MARKERS, layer="counts", n_null=0)
    return {
        "sig.moran_gene": m.get("median_gene_moran", np.nan),
        "sig.snr_pooled": snr.get("pooled_mean_signal_ratio", np.nan),
        "sig.zero_fraction": sp_res.get("zero_fraction", np.nan),
        "sig.tx_per_area": tpa.get("mean_tx_per_area", np.nan),
        "sc.mean_coherence": sc.get("mean_coherence", np.nan),
        "sc.pas": sc.get("pas", np.nan),
        "nc.mean": nc.get("mean_contamination", np.nan),
        "crisp_purity": crisp.get("crisp_purity", np.nan),
    }


ORIENTATION: dict[str, int] = {
    "iv.integrated": 1, "isc.iscg": 1, "isc.iscl": 1, "isc.integrated": 1, "fid.mean_auc": 1,
    "pmp": 1, "crisp_purity": 1, "mecr": -1, "nmp.salas": 1, "ncp.moldia": 1,
    "onto.rho": 1, "onto.rwr_rho": 1, "ext.macro_f1": 1, "ext.ari": 1, "ext.ecs": 1,
    "sig.moran_gene": 1, "sig.snr_pooled": 1, "sig.zero_fraction": -1, "sig.tx_per_area": 1,
    "sc.mean_coherence": 1, "sc.pas": -1, "nc.mean": -1,
}

YERLY_SCENARIOS: dict[str, tuple[list[float], set[str]]] = {
    "label_noise": ([0.0, 0.05, 0.1, 0.25, 0.5],
                    {"iv.integrated", "isc.iscg", "isc.iscl", "isc.integrated", "fid.mean_auc", "pmp",
                     "nmp.salas", "onto.rho", "onto.rwr_rho", "ext.macro_f1", "ext.ari", "ext.ecs"}),
    "over_split": ([0, 1, 2, 3], {"isc.iscl", "isc.integrated", "ext.macro_f1", "ext.ari", "ext.ecs"}),
    "merge": ([0, 1, 2, 3], {"fid.mean_auc", "pmp", "ext.macro_f1", "ext.ari", "ext.ecs"}),
    "type_swap": ([0, 1, 2, 3],
                  {"fid.mean_auc", "pmp", "nmp.salas", "onto.rho", "onto.rwr_rho", "ext.macro_f1"}),
    "spillover": ([0.0, 0.05, 0.1, 0.2, 0.4],
                  {"crisp_purity", "mecr", "pmp", "nmp.salas", "ncp.moldia"}),
    "depth_loss": ([0.0, 0.25, 0.5, 0.75], {"iv.integrated", "fid.mean_auc"}),
    "nuisance_shuffle": ([0, 1, 2], set()),
}

XENIUM_SCENARIOS: dict[str, tuple[list[float], set[str]]] = {
    "ambient": ([0.0, 0.05, 0.15, 0.3, 0.6],
                {"sig.moran_gene", "sig.snr_pooled", "crisp_purity"}),
    "spillover": ([0.0, 0.05, 0.1, 0.2, 0.4],
                  {"nc.mean", "crisp_purity"}),
    "depth_loss": ([0.0, 0.25, 0.5, 0.75, 0.9], {"sig.zero_fraction", "sig.tx_per_area"}),
    "coord_shuffle": ([0.0, 0.1, 0.25, 0.5, 1.0],
                      {"sig.moran_gene", "sc.mean_coherence", "sc.pas"}),
    "nuisance_shuffle": ([0, 1, 2], set()),
}


def run_benchmark(seeds: int = 3, out_dir: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the real-data perturbation suite and return ``(values_df, verdicts_df)``."""
    yerly = load_yerly()
    xenium = load_xenium()
    rows = []
    for scen, (levels, _) in YERLY_SCENARIOS.items():
        for lv in levels:
            for seed in range(seeds):
                a = perturb_yerly(yerly, scen, lv, seed)
                for m, v in eval_yerly(a, yerly).items():
                    rows.append(("yerly2022_bcc", scen, lv, seed, m, float(v) if v is not None else np.nan))
    for scen, (levels, _) in XENIUM_SCENARIOS.items():
        for lv in levels:
            for seed in range(seeds):
                a = perturb_xenium(xenium, scen, lv, seed)
                for m, v in eval_xenium(a).items():
                    rows.append(("xenium_breast", scen, lv, seed, m, float(v) if v is not None else np.nan))
    df = pd.DataFrame(rows, columns=["dataset", "scenario", "level", "seed", "metric", "value"])

    verdicts = []
    for (ds, scen, metric), g in df.groupby(["dataset", "scenario", "metric"]):
        spec = YERLY_SCENARIOS if ds == "yerly2022_bcc" else XENIUM_SCENARIOS
        levels, expected_set = spec[scen]
        o = ORIENTATION[metric]
        m = g.groupby("level")["value"].agg(["mean", "std"]).reindex(levels)
        v0 = float(m["mean"].iloc[0])
        v_end = float(m["mean"].iloc[-1])
        bad = -o * (m["mean"] - v0)
        thr = max(3.0 * float(np.nan_to_num(m["std"].iloc[0])), 0.02 * max(abs(v0), 1e-3), 1e-6)
        ok = ~bad.isna()
        rho = (float(spearmanr(np.arange(ok.sum()), bad[ok].to_numpy())[0])
               if ok.sum() > 2 and np.nanstd(bad[ok]) > 0 else 0.0)
        moved = float(np.nanmax(np.abs(bad.to_numpy()))) > thr
        if scen == "nuisance_shuffle":
            verdict = "stable" if not moved else "FALSE-ALARM"
        elif metric in expected_set:
            verdict = "detects" if (moved and float(bad.iloc[-1]) > thr and rho >= 0.8) else "weak/flat"
        else:
            verdict = "flat" if not moved else "moved"
        first_lv = next((lv for lv, b in zip(levels, bad) if b > thr), None)
        verdicts.append(dict(
            dataset=ds, scenario=scen, metric=metric,
            expect="down" if metric in expected_set else ("flat" if scen == "nuisance_shuffle" else "-"),
            verdict=verdict, v0=v0, v_max=v_end, rho=rho, first_level=first_lv,
        ))
    vdf = pd.DataFrame(verdicts)
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_dir / "values.csv", index=False)
        vdf.to_csv(out_dir / "verdicts.csv", index=False)
        lines = [
            "# Real-Data Perturbation Benchmark Summary",
            "",
            "Perturbations applied to committed real datasets (`tests/goldens/isc_yerly2022/` and "
            "`tests/goldens/signalqc_spatialqm/`).",
            "",
            "| dataset | scenario | metric | expect | verdict | v0 | v_max | rho | first_level |",
            "|---|---|---|---|---|---:|---:|---:|---:|",
        ]
        for r in vdf.itertuples(index=False):
            if r.expect in ("down", "flat"):
                fl = f"{r.first_level:g}" if r.first_level is not None and not pd.isna(r.first_level) else "-"
                lines.append(
                    f"| `{r.dataset}` | `{r.scenario}` | `{r.metric}` | {r.expect} | **{r.verdict}** | "
                    f"{r.v0:.4f} | {r.v_max:.4f} | {r.rho:.2f} | {fl} |"
                )
        (out_dir / "matrix.md").write_text("\n".join(lines) + "\n")
    return df, vdf


def main(argv=None) -> None:
    """CLI entry point for the real-data perturbation benchmark."""
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--out", default=str(Path(__file__).parent / "real_perturbation" / "results"))
    args = p.parse_args(argv)
    run_benchmark(seeds=args.seeds, out_dir=Path(args.out))


if __name__ == "__main__":
    main()
