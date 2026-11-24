# Metric Augmentation & Behavioral Contract Report

## Cohort Summary

| dataset | status | n_checks | n_passed | n_failed | n_caveats | n_skipped |
|---|---|---|---|---|---|---|
| xenium_breast_golden | ok | 24 | 19 | 0 | 5 | 0 |
| yerly2022_bcc_split | ok | 28 | 17 | 0 | 3 | 8 |
| merfish_bregma4_allen_proxy | ok | 23 | 18 | 0 | 5 | 0 |

## Top 5 Highest-ROI OpenProblems Dataset Picks

| rank | dataset_id | targets | contract_metrics | local_status |
|---|---|---|---|---|
| 1 | zenodo_spatial/merfish/human_cortex_{2,3,4} | imputation benchmark (clean same-cell panel holdout, ~300 -> ~3.5k genes); annotation benchmark (3-level human cortex hierarchy) | residual_spatial_autocorrelation, spatial_organization_consistency, hierarchical_accuracy | pending_download |
| 2 | zenodo_spatial_slidetags/slidetags/{human_cortex,human_skin_melanoma,human_tonsil,mouse_embryo} | spillover benchmark (zero-lateral-spillover nuclear WTA base); imputation benchmark (clean WTA holdouts); label-swap check (controlled spillover injection) | residual_spatial_autocorrelation, relative_pairwise_ct_expression, neighbor_contamination | pending_download |
| 3 | cellxgene_census/{hlca,tabula_sapiens,dkd,immune_cell_atlas} | annotation benchmark (native CL terms, no manual crosswalk tables); spatial-anno-metrics ontology_consistency | ontology_consistency, inter_sample_consistency, relative_pairwise_ct_expression | pending_download |
| 4 | task_ist_preprocessing/allen_brain_cell_atlas_merfish_combined | annotation benchmark (non-cancer multi-section hierarchical MERFISH); imputation benchmark; spatial_organization_consistency | spatial_organization_consistency, residual_spatial_autocorrelation, relative_pairwise_ct_expression | pending_download |
| 5 | tenx_visium/post_xenium/{human_colon_cancer,mouse_brain} | imputation benchmark (cross-section WTA); same-section panel separability (real vs simulated) | relative_pairwise_ct_expression, residual_spatial_autocorrelation, moran_signal | pending_download |

## Per-Contract Verdicts

| dataset | metric | perturbation | role | verdict | v0 | v_end | rho | first_level |
|---|---|---|---|---|---|---|---|---|
| xenium_breast_golden | ext.macro_f1 | label_noise | must_detect | pass | 1.0000 | 0.4432 | 1.0000 | 0.1000 |
| xenium_breast_golden | ext.macro_f1 | merge | must_detect | pass | 1.0000 | 0.2500 | 1.0000 | 1.0000 |
| xenium_breast_golden | ext.macro_f1 | type_swap | must_detect | pass | 1.0000 | 0.3333 | 1.0000 | 1.0000 |
| xenium_breast_golden | ext.macro_f1 | nuisance_shuffle | must_be_invariant | pass | 1.0000 | 1.0000 | 0.0000 | nan |
| xenium_breast_golden | ext.macro_f1 | coord_shuffle | must_be_invariant | pass | 1.0000 | 1.0000 | 0.0000 | nan |
| xenium_breast_golden | fid.mean_auc | label_noise | must_detect | pass | 0.9605 | 0.6616 | 1.0000 | 0.1000 |
| xenium_breast_golden | fid.mean_auc | nuisance_shuffle | must_be_invariant | pass | 0.9605 | 0.9605 | 0.0000 | nan |
| xenium_breast_golden | fid.mean_auc | coord_shuffle | must_be_invariant | pass | 0.9605 | 0.9605 | 0.0000 | nan |
| xenium_breast_golden | fid.mean_auc | over_split | known_caveat | caveat_observed | 0.9605 | 0.9592 | 0.5000 | nan |
| xenium_breast_golden | pmp | label_noise | must_detect | pass | 0.6704 | 0.3684 | 1.0000 | 0.1000 |
| xenium_breast_golden | pmp | nuisance_shuffle | must_be_invariant | pass | 0.6704 | 0.6704 | 0.0000 | nan |
| xenium_breast_golden | pmp | coord_shuffle | must_be_invariant | pass | 0.6704 | 0.6704 | 0.0000 | nan |
| xenium_breast_golden | pmp | spillover | known_caveat | caveat_observed | 0.6704 | 0.5788 | 1.0000 | 0.1000 |
| xenium_breast_golden | resid_moran.explained | label_noise | must_detect | pass | 0.5007 | 0.0976 | 1.0000 | 0.1000 |
| xenium_breast_golden | resid_moran.explained | merge | must_detect | pass | 0.5007 | 0.2704 | 1.0000 | 1.0000 |
| xenium_breast_golden | resid_moran.explained | nuisance_shuffle | must_be_invariant | pass | 0.5007 | 0.5007 | 0.8660 | nan |
| xenium_breast_golden | resid_moran.explained | infiltration | known_caveat | caveat_observed | 0.5007 | 0.5364 | -0.4000 | 0.2500 |
| xenium_breast_golden | resid_moran.explained | coord_shuffle | known_caveat | caveat_observed | 0.5007 | nan | 0.0000 | nan |
| xenium_breast_golden | resid_moran.residual | label_noise | must_detect | pass | 0.0496 | 0.0896 | 1.0000 | 0.1000 |
| xenium_breast_golden | resid_moran.residual | merge | must_detect | pass | 0.0496 | 0.0724 | 1.0000 | 1.0000 |
| xenium_breast_golden | resid_moran.residual | spillover | must_detect | pass | 0.0496 | 0.1334 | 1.0000 | 0.1000 |
| xenium_breast_golden | resid_moran.residual | nuisance_shuffle | must_be_invariant | pass | 0.0496 | 0.0496 | -0.8660 | nan |
| xenium_breast_golden | resid_moran.residual | infiltration | must_not_degrade | pass | 0.0496 | 0.0268 | -1.0000 | nan |
| xenium_breast_golden | resid_moran.residual | coord_shuffle | known_caveat | caveat_observed | 0.0496 | -0.0015 | -1.0000 | nan |
| yerly2022_bcc_split | ext.macro_f1 | label_noise | must_detect | pass | 1.0000 | 0.4805 | 1.0000 | 0.1000 |
| yerly2022_bcc_split | ext.macro_f1 | merge | must_detect | pass | 1.0000 | 0.4545 | 1.0000 | 1.0000 |
| yerly2022_bcc_split | ext.macro_f1 | type_swap | must_detect | pass | 1.0000 | 0.5556 | 1.0000 | 1.0000 |
| yerly2022_bcc_split | ext.macro_f1 | nuisance_shuffle | must_be_invariant | pass | 1.0000 | 1.0000 | 0.0000 | nan |
| yerly2022_bcc_split | ext.macro_f1 | coord_shuffle | must_be_invariant | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | iv.integrated | label_noise | must_detect | pass | 0.8453 | 0.4626 | 1.0000 | 0.1000 |
| yerly2022_bcc_split | iv.integrated | over_split | must_detect | pass | 0.8453 | 0.7845 | 1.0000 | 1.0000 |
| yerly2022_bcc_split | iv.integrated | nuisance_shuffle | must_be_invariant | pass | 0.8453 | 0.8456 | -0.8660 | nan |
| yerly2022_bcc_split | iv.integrated | coord_shuffle | must_be_invariant | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | iv.integrated | infiltration | known_caveat | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | fid.mean_auc | label_noise | must_detect | pass | 0.9600 | 0.6972 | 1.0000 | 0.1000 |
| yerly2022_bcc_split | fid.mean_auc | nuisance_shuffle | must_be_invariant | pass | 0.9600 | 0.9600 | 0.0000 | nan |
| yerly2022_bcc_split | fid.mean_auc | coord_shuffle | must_be_invariant | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | fid.mean_auc | over_split | known_caveat | caveat_observed | 0.9600 | 0.9574 | 0.5000 | nan |
| yerly2022_bcc_split | pmp | label_noise | must_detect | pass | 0.8179 | 0.4219 | 1.0000 | 0.1000 |
| yerly2022_bcc_split | pmp | nuisance_shuffle | must_be_invariant | pass | 0.8179 | 0.8179 | 0.0000 | nan |
| yerly2022_bcc_split | pmp | coord_shuffle | must_be_invariant | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | pmp | spillover | known_caveat | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | nmp | label_noise | must_detect | pass | 0.9483 | 0.3750 | 1.0000 | 0.1000 |
| yerly2022_bcc_split | nmp | type_swap | must_detect | pass | 0.9483 | 0.3528 | 1.0000 | 1.0000 |
| yerly2022_bcc_split | nmp | nuisance_shuffle | must_be_invariant | pass | 0.9483 | 0.9477 | 0.8660 | nan |
| yerly2022_bcc_split | nmp | coord_shuffle | must_be_invariant | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | nmp | merge | known_caveat | caveat_observed | 0.9483 | 0.7302 | 1.0000 | 1.0000 |
| yerly2022_bcc_split | rel_pw_ct | label_noise | must_detect | pass | 0.8190 | 0.6986 | 1.0000 | 0.2500 |
| yerly2022_bcc_split | rel_pw_ct | type_swap | must_detect | pass | 0.8190 | 0.4496 | 1.0000 | 1.0000 |
| yerly2022_bcc_split | rel_pw_ct | nuisance_shuffle | must_be_invariant | pass | 0.8190 | 0.8190 | 0.0000 | nan |
| yerly2022_bcc_split | rel_pw_ct | coord_shuffle | must_be_invariant | skipped | nan | nan | nan | nan |
| yerly2022_bcc_split | rel_pw_ct | panel_growth | known_caveat | caveat_observed | 0.8190 | 0.8190 | -1.0000 | nan |
| merfish_bregma4_allen_proxy | ext.macro_f1 | label_noise | must_detect | pass | 1.0000 | 0.4087 | 1.0000 | 0.1000 |
| merfish_bregma4_allen_proxy | ext.macro_f1 | merge | must_detect | pass | 1.0000 | 0.4545 | 1.0000 | 1.0000 |
| merfish_bregma4_allen_proxy | ext.macro_f1 | type_swap | must_detect | pass | 1.0000 | 0.5556 | 1.0000 | 1.0000 |
| merfish_bregma4_allen_proxy | ext.macro_f1 | nuisance_shuffle | must_be_invariant | pass | 1.0000 | 1.0000 | 0.0000 | nan |
| merfish_bregma4_allen_proxy | ext.macro_f1 | coord_shuffle | must_be_invariant | pass | 1.0000 | 1.0000 | 0.0000 | nan |
| merfish_bregma4_allen_proxy | resid_moran.explained | label_noise | must_detect | pass | 0.4994 | 0.1194 | 1.0000 | 0.1000 |
| merfish_bregma4_allen_proxy | resid_moran.explained | merge | must_detect | pass | 0.4994 | 0.4161 | 1.0000 | 1.0000 |
| merfish_bregma4_allen_proxy | resid_moran.explained | nuisance_shuffle | must_be_invariant | pass | 0.4994 | 0.5013 | -0.8660 | nan |
| merfish_bregma4_allen_proxy | resid_moran.explained | infiltration | known_caveat | caveat_observed | 0.4994 | 0.5457 | -0.8000 | nan |
| merfish_bregma4_allen_proxy | resid_moran.explained | coord_shuffle | known_caveat | caveat_observed | 0.4994 | 0.8093 | -0.8000 | 0.2500 |
| merfish_bregma4_allen_proxy | resid_moran.residual | label_noise | must_detect | pass | 0.0467 | 0.0822 | 1.0000 | 0.1000 |
| merfish_bregma4_allen_proxy | resid_moran.residual | merge | must_detect | pass | 0.0467 | 0.0545 | 1.0000 | 1.0000 |
| merfish_bregma4_allen_proxy | resid_moran.residual | spillover | must_detect | pass | 0.0467 | 0.5435 | 1.0000 | 0.1000 |
| merfish_bregma4_allen_proxy | resid_moran.residual | nuisance_shuffle | must_be_invariant | pass | 0.0467 | 0.0467 | -0.8660 | nan |
| merfish_bregma4_allen_proxy | resid_moran.residual | infiltration | must_not_degrade | pass | 0.0467 | 0.0261 | -1.0000 | nan |
| merfish_bregma4_allen_proxy | resid_moran.residual | coord_shuffle | known_caveat | caveat_observed | 0.0467 | 0.0003 | -1.0000 | nan |
| merfish_bregma4_allen_proxy | spatial_org.corr | coord_shuffle | must_detect | pass | 0.8316 | -0.0628 | 1.0000 | 0.5000 |
| merfish_bregma4_allen_proxy | spatial_org.corr | nuisance_shuffle | must_be_invariant | pass | 0.8316 | 0.8316 | -0.8660 | nan |
| merfish_bregma4_allen_proxy | spatial_org.corr | merge | known_caveat | caveat_observed | 0.8316 | 0.8384 | -0.5000 | nan |
| merfish_bregma4_allen_proxy | spatial_org.offdiag | coord_shuffle | must_detect | pass | 0.5945 | -0.0721 | 1.0000 | 0.5000 |
| merfish_bregma4_allen_proxy | spatial_org.offdiag | label_noise | must_detect | pass | 0.5945 | 0.5350 | 1.0000 | 0.1000 |
| merfish_bregma4_allen_proxy | spatial_org.offdiag | nuisance_shuffle | must_be_invariant | pass | 0.5945 | 0.5945 | -0.8660 | nan |
| merfish_bregma4_allen_proxy | spatial_org.offdiag | merge | known_caveat | caveat_observed | 0.5945 | 0.6878 | -0.5000 | nan |
