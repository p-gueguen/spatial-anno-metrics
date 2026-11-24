# Augmentation & Behavioral Contract Verification Harness

Generic `AnnData` perturbations (`sam.perturb_anndata`) and declarative behavioral contract checks (`sam.verify_metric_contracts`) verify that every metric degrades monotonically on the failure modes it is built to catch and remains invariant on negative controls.

## Quick Run

```bash
PYTHONPATH=src python benchmarks/augmentation/run_contracts.py
```

Outputs are written to `benchmarks/augmentation/results/contract_report.csv` and `benchmarks/augmentation/results/summary.md`.

## Top-5 Highest-ROI OpenProblems Dataset Picks

Tracked in [`datasets.py`](datasets.py) (`TOP5_OPENPROBLEMS_ROI`):

1. `zenodo_spatial/merfish/human_cortex_{2,3,4}` (Zenodo 12785822): 3-4k cells x ~3,550 genes with `cluster_L1..L3` hierarchy, for imputation (clean same-cell panel holdout, ~300 -> ~3.5k genes) and 3-level hierarchical annotation.
2. `zenodo_spatial_slidetags/slidetags/{human_cortex,human_skin_melanoma,human_tonsil,mouse_embryo}` (Zenodo 12785822): zero-lateral-spillover nuclear WTA base on real 2D coordinates, for spillover, imputation and label-swap checks.
3. `cellxgene_census/{hlca,tabula_sapiens,dkd,immune_cell_atlas}`: native `CL:` ontology terms, so annotation scoring and `ontology_consistency` need no manual crosswalk tables.
4. `task_ist_preprocessing/allen_brain_cell_atlas_merfish_combined`: 4 whole mouse brain MERFISH sections (~1,100 genes + Yao scRNA-seq) for multi-section `spatial_organization_consistency` and hierarchical annotation.
5. `tenx_visium/post_xenium/{human_colon_cancer,mouse_brain}` paired with Xenium: same-section post-Xenium WTA, for cross-section imputation and same-section panel separability.
