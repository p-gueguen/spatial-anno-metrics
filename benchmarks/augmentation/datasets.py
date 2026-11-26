"""Top-5 highest-ROI OpenProblems dataset registry and local dataset loader for contract validation."""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MERFISH_BASE = Path(os.environ.get("SAM_MERFISH_BASE", "data/merfish_bregma"))
MERFISH_PREFIX = os.environ.get("SAM_MERFISH_PREFIX", "merfish_Animal1_Bregma_")

TOP5_OPENPROBLEMS_ROI: list[dict[str, Any]] = [
    {
        "rank": 1,
        "dataset_id": "zenodo_spatial/merfish/human_cortex_{2,3,4}",
        "source": "Zenodo 12785822 (OpenProblems task_ist_preprocessing / datasets)",
        "scale": "3-4k cells x ~3,550 genes per section (cluster_L1..L3 hierarchy)",
        "primary_targets": ["imputation benchmark (clean same-cell panel holdout, ~300 -> ~3.5k genes)", "annotation benchmark (3-level human cortex hierarchy)"],
        "contract_metrics": ["residual_spatial_autocorrelation", "spatial_organization_consistency", "hierarchical_accuracy"],
        "cache_glob": "merfish_human_cortex*.h5ad",
    },
    {
        "rank": 2,
        "dataset_id": "zenodo_spatial_slidetags/slidetags/{human_cortex,human_skin_melanoma,human_tonsil,mouse_embryo}",
        "source": "Zenodo 12785822 (OpenProblems datasets)",
        "scale": "4-6k nuclei x 15-19k WTA genes with 2D spatial coords (cell_type)",
        "primary_targets": ["spillover benchmark (zero-lateral-spillover nuclear WTA base)", "imputation benchmark (clean WTA holdouts)", "label-swap check (controlled spillover injection)"],
        "contract_metrics": ["residual_spatial_autocorrelation", "relative_pairwise_ct_expression", "neighbor_contamination"],
        "cache_glob": "*slidetags*.h5ad",
    },
    {
        "rank": 3,
        "dataset_id": "cellxgene_census/{hlca,tabula_sapiens,dkd,immune_cell_atlas}",
        "source": "CZI CELLxGENE Census (native CL: cell_type_ontology_term_id)",
        "scale": "Atlas-scale multi-donor references with standardized CL terms",
        "primary_targets": ["annotation benchmark (native CL terms, no manual crosswalk tables)", "spatial-anno-metrics ontology_consistency"],
        "contract_metrics": ["ontology_consistency", "inter_sample_consistency", "relative_pairwise_ct_expression"],
        "cache_glob": "*census*.h5ad",
    },
    {
        "rank": 4,
        "dataset_id": "task_ist_preprocessing/allen_brain_cell_atlas_merfish_combined",
        "source": "Brain Image Library + Yao scRNA-seq (OpenProblems task_ist_preprocessing)",
        "scale": "4 whole mouse brain MERFISH sections x ~1,100 genes + matched Yao scRNA-seq",
        "primary_targets": ["annotation benchmark (non-cancer multi-section hierarchical MERFISH)", "imputation benchmark", "spatial_organization_consistency"],
        "contract_metrics": ["spatial_organization_consistency", "residual_spatial_autocorrelation", "relative_pairwise_ct_expression"],
        "cache_glob": "*allen_brain*merfish*.h5ad",
    },
    {
        "rank": 5,
        "dataset_id": "tenx_visium/post_xenium/{human_colon_cancer,mouse_brain}",
        "source": "10x Genomics post-Xenium Visium paired with Xenium (OpenProblems task_ist_preprocessing)",
        "scale": "Same-section post-Xenium WTA paired with Xenium panel on colon cancer and mouse brain",
        "primary_targets": ["imputation benchmark (cross-section WTA)", "same-section panel separability (real vs simulated)"],
        "contract_metrics": ["relative_pairwise_ct_expression", "residual_spatial_autocorrelation", "moran_signal"],
        "cache_glob": "*post_xenium*.h5ad",
    },
]


def roi_inventory(cache_dir: Path | None = None) -> pd.DataFrame:
    """Return inventory table for the Top-5 OpenProblems ROI datasets, checking an optional local cache."""
    rows = []
    for item in TOP5_OPENPROBLEMS_ROI:
        hits = sorted(cache_dir.glob(item["cache_glob"])) if cache_dir and cache_dir.is_dir() else []
        rows.append(
            {
                "rank": item["rank"],
                "dataset_id": item["dataset_id"],
                "targets": "; ".join(item["primary_targets"]),
                "contract_metrics": ", ".join(item["contract_metrics"]),
                "local_status": f"cached ({len(hits)} files)" if hits else "pending_download",
            }
        )
    return pd.DataFrame(rows)


def load_local_cohorts() -> list[dict[str, Any]]:
    """Load committed real goldens (Xenium breast, Yerly BCC split) plus local MERFISH Bregma if mounted."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import real_perturbation as rp

    cohorts: list[dict[str, Any]] = []

    xenium = rp.load_xenium()
    cohorts.append(
        {
            "name": "xenium_breast_golden",
            "adata": xenium,
            "reference": None,
            "label_key": "label",
            "sample_key": None,
            "layer": "counts",
            "marker_sets": rp.XENIUM_MARKERS,
            "perturbation_kwargs": {"merge": {"pairs": [("Epithelial", "Fibroblast"), ("Tcell", "Myeloid")]}},
        }
    )

    full = rp.load_yerly()
    is_ref = full.obs["patient_bcc"].isin({"BCC1", "BCC2", "BCC3"}).to_numpy()
    ref_y = full[is_ref].copy()
    ref_y.obs["cell_type"] = ref_y.obs["truth"].astype(str)
    qry_y = full[~is_ref].copy()
    cohorts.append(
        {
            "name": "yerly2022_bcc_split",
            "adata": qry_y,
            "reference": ref_y,
            "label_key": "label",
            "sample_key": "patient_bcc",
            "layer": "counts",
            "marker_sets": rp.YERLY_MARKERS,
            "perturbation_kwargs": {"merge": {"pairs": [("Tcell", "Bcell"), ("Fibroblast", "Endothelial")]}},
        }
    )

    secs = ["0.01", "0.06", "0.11", "0.16"]
    if all((MERFISH_BASE / "datasets" / f"{MERFISH_PREFIX}{s}" / "query.h5ad").is_file() for s in secs):
        parts = []
        for s in secs:
            q_p = MERFISH_BASE / "datasets" / f"{MERFISH_PREFIX}{s}" / "query.h5ad"
            s_p = MERFISH_BASE / "solutions" / f"{MERFISH_PREFIX}{s}" / "solution.parquet"
            a = ad.read_h5ad(q_p)
            sol = pd.read_parquet(s_p).set_index("cell_id").loc[a.obs_names]
            a.obs["truth"] = sol["label"].astype(str).to_numpy()
            a.obs["label"] = a.obs["truth"].copy()
            a.obs["sample"] = sol["section"].astype(str).to_numpy()
            a.obs_names = [f"{s}__{cid}" for cid in a.obs_names]
            parts.append(a)
        merfish = ad.concat(parts, join="inner", merge="same")
        counts = merfish.obs["truth"].value_counts()
        top4 = list(counts.index[:4])
        cohorts.append(
            {
                "name": "merfish_bregma4_allen_proxy",
                "adata": merfish,
                "reference": None,
                "label_key": "label",
                "sample_key": "sample",
                "layer": None,
                "marker_sets": None,
                "perturbation_kwargs": {"merge": {"pairs": [(top4[0], top4[1]), (top4[2], top4[3])]}},
            }
        )
    return cohorts
