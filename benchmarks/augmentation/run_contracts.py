"""CLI runner to verify metric behavioral contracts across local datasets and report Top-5 ROI status."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sam_mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/sam_numba")
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks" / "augmentation"))

import datasets as ds_reg  # noqa: E402
import spatial_anno_metrics as sam  # noqa: E402


def _make_evaluator(cohort: dict) -> dict:
    label_key = cohort["label_key"]
    sample_key = cohort["sample_key"]
    layer = cohort["layer"]
    markers = cohort["marker_sets"]

    def _eval(a, ref):
        out: dict[str, float] = {}
        if "truth" in a.obs:
            out["ext.macro_f1"] = float(sam.external_scores(y_true=a.obs["truth"], y_pred=a.obs[label_key])["macro_f1"])
        if "X_pca" in a.obsm:
            iv = sam.internal_validity(a, label_key, embedding="X_pca", subsample=1000)
            out["iv.integrated"] = np.nan if iv.get("integrated") is None else float(iv["integrated"])
        if markers:
            fid = sam.marker_program_fidelity(a, label_key, markers, layer=layer)
            out["fid.mean_auc"] = np.nan if fid.get("mean_auc") is None else float(fid["mean_auc"])
            pmp = sam.pmp(a, markers, label_key, layer=layer)
            out["pmp"] = np.nan if pmp.get("mean_pmp") is None else float(pmp["mean_pmp"])
        if ref is not None:
            nmp = sam.nmp(a, ref, label_key, "cell_type", n_permutations=10, layer=layer)
            out["nmp"] = np.nan if nmp.get("nmp") is None else float(nmp["nmp"])
            rel = sam.relative_pairwise_ct_expression(
                a, ref, label_key=label_key, ref_label_key="cell_type", layer=layer
            )
            val = rel.get("relative_pairwise_ct_expression")
            out["rel_pw_ct"] = np.nan if val is None else float(val)
        if sample_key and sample_key in a.obs and a.obs[sample_key].nunique() >= 2:
            isc = sam.inter_sample_consistency(a, label_key, sample_key, layer=layer)
            if isc.get("status") == "ok":
                out["isc.integrated"] = float(isc["integrated"])
        if "spatial" in a.obsm:
            target = a
            if sample_key and sample_key in a.obs:
                sec0 = sorted(a.obs[sample_key].astype(str).unique())[0]
                target = a[a.obs[sample_key].astype(str) == sec0].copy()
            rm = sam.residual_spatial_autocorrelation(
                target,
                label_key=label_key,
                layer=layer,
                exclude_controls=(layer is not None),
            )
            out["resid_moran.explained"] = (
                np.nan if rm.get("explained_fraction") is None else float(rm["explained_fraction"])
            )
            out["resid_moran.residual"] = (
                np.nan if rm.get("residual_median_moran") is None else float(rm["residual_median_moran"])
            )
            if sample_key and sample_key in a.obs and a.obs[sample_key].nunique() >= 2:
                org = sam.spatial_organization_consistency(a, label_key=label_key, sample_key=sample_key)
                out["spatial_org.corr"] = np.nan if org.get("mean_corr") is None else float(org["mean_corr"])
                out["spatial_org.offdiag"] = (
                    np.nan if org.get("offdiag_mean_corr") is None else float(org["offdiag_mean_corr"])
                )
        return out

    return {"battery": _eval}


def _md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        cells = [f"{row[c]:.4f}" if isinstance(row[c], (float, np.floating)) and pd.notna(row[c]) else str(row[c]) for c in cols]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    """Execute behavioral contract checks on all local cohorts and write report files."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=ROOT / "benchmarks" / "augmentation" / "results")
    ap.add_argument("--seeds", type=int, default=2)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    cohorts = ds_reg.load_local_cohorts()
    all_summaries = []
    cohort_rows = []

    for cohort in cohorts:
        ev = _make_evaluator(cohort)
        res = sam.verify_metric_contracts(
            cohort["adata"],
            ev,
            reference=cohort["reference"],
            seeds=tuple(range(args.seeds)),
            label_key=cohort["label_key"],
            sample_key=cohort["sample_key"],
            layer=cohort["layer"],
            perturbation_kwargs=cohort["perturbation_kwargs"],
        )
        s = res["summary"].copy()
        s.insert(0, "dataset", cohort["name"])
        all_summaries.append(s)
        cohort_rows.append(
            {
                "dataset": cohort["name"],
                "status": res["status"],
                "n_checks": res["n_checks"],
                "n_passed": res["n_passed"],
                "n_failed": res["n_failed"],
                "n_caveats": res["n_caveats"],
                "n_skipped": res["n_skipped"],
            }
        )

    report_df = pd.concat(all_summaries, ignore_index=True)
    report_df.to_csv(args.out / "contract_report.csv", index=False)
    cohort_df = pd.DataFrame(cohort_rows)
    roi_df = ds_reg.roi_inventory(args.cache_dir)

    md = [
        "# Metric Augmentation & Behavioral Contract Report",
        "",
        "## Cohort Summary",
        "",
        _md_table(cohort_df),
        "",
        "## Top 5 Highest-ROI OpenProblems Dataset Picks",
        "",
        _md_table(roi_df),
        "",
        "## Per-Contract Verdicts",
        "",
        _md_table(report_df[["dataset", "metric", "perturbation", "role", "verdict", "v0", "v_end", "rho", "first_level"]]),
        "",
    ]
    (args.out / "summary.md").write_text("\n".join(md))
    print(cohort_df.to_string(index=False))


if __name__ == "__main__":
    main()
