| edge case | internal_validity | isolated_label_scores | clisi | ontology_consistency | inter_sample_consistency | marker_program_fidelity | marker_gene_overlap | purity | nmp | coexpression_similarity | ncp | neighbor_contamination | spatial_coherence | neighborhood_sanity | per_cell_confidence | label_stability | signal_qc | external | calibration_novelty |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| one_label | NaN silhouette,neighborhood_purity,orbital_medoid,ward_propmatch,avg_similarity,davies_bouldin,calinski_harabasz,integrated | NaN asw,f1 | NaN clisi | NaN rho,rwr_rho | NaN iscg,iscl,integrated | NaN mean_auc,mean_cohens_d | NaN mean_overlap | ok | NaN nmp | ok | ok | NaN neighbor_contamination | NaN mean_coherence,pas,chaos | NaN n_suspicious | ok | ok | ok | ok | ok |
| singleton_type | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| two_types | ok | ok | ok | NaN rho,rwr_rho | ok | NaN mean_auc,mean_cohens_d | NaN mean_overlap | NaN pmp | NaN nmp | NaN by_type | ok | NaN neighbor_contamination | ok | ok | ok | ok | ok | ok | ok |
| empty_marker_set | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| markers_off_panel | ok | ok | ok | ok | ok | NaN mean_auc,mean_cohens_d | NaN mean_overlap | NaN crisp_purity,mecr,pmp | ok | ok | ok | NaN neighbor_contamination | ok | ok | NaN mean_margin,mean_entropy | ok | ok | ok | ok |
| zero_count_cells | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| unused_category | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| integer_labels | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ≠ ece |
| nan_labels | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | ok | 💥 ValueError | 💥 ValueError | 💥 ValueError | ok | ok | ok | 💥 ValueError | ok |
| no_control_probes | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | NaN moran_gap,frac_above_ctrl,snr | ok | ok |
| tiny_30_cells | ok | ok | ok | NaN rho,rwr_rho | NaN iscg,iscl,integrated | ok | ok | ok | NaN nmp | NaN by_type | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| duplicate_coords | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| all_zero_marker_gene | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| dense_matrix | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| unicode_labels | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| extra_marker_type | ok | ok | ok | ok | ok | ok | ok | ≠ crisp_purity,mecr,pmp | ok | ok | ok | ok | ok | ok | ≠ mean_margin,mean_entropy | ≠ flip_rate | ok | ok | ok |
| no_spatial_coords | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | NaN neighbor_contamination | NaN mean_coherence,pas,chaos | NaN n_suspicious | ok | ok | NaN moran_gap,frac_above_ctrl | NaN spari | ok |
| lognorm_X_no_counts | ok | ok | ok | ok | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | 💥 ValueError | ok | ok | 💥 ValueError | 💥 ValueError | 💥 ValueError | ok | 💥 ValueError |

Errors:

- **nan_labels** / `internal_validity`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `isolated_label_scores`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `clisi`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `ontology_consistency`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `inter_sample_consistency`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `marker_program_fidelity`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `marker_gene_overlap`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `purity`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `nmp`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `coexpression_similarity`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `neighbor_contamination`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `spatial_coherence`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `neighborhood_sanity`: ValueError: adata.obs['label'] has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') be
- **nan_labels** / `external`: ValueError: the label array has 30 missing (NaN/None) labels of 600: drop those cells or give them a label (e.g. 'Unassigned') befor
- **lognorm_X_no_counts** / `inter_sample_consistency`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `marker_program_fidelity`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `marker_gene_overlap`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `purity`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `nmp`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `coexpression_similarity`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `ncp`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `neighbor_contamination`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `per_cell_confidence`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `label_stability`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `signal_qc`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
- **lognorm_X_no_counts** / `calibration_novelty`: ValueError: layer 'counts' is not in adata.layers (none) and X is not raw counts (negative or non-integer values): pass layer= namin
