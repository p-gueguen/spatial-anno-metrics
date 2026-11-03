| edge case | internal_validity | inter_sample_consistency | marker_program_fidelity | marker_gene_overlap | purity | nmp | ncp | neighbor_contamination | spatial_coherence | neighborhood_sanity | per_cell_confidence | label_stability | signal_qc | external |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| one_label | NaN silhouette,neighborhood_purity,orbital_medoid,ward_propmatch,avg_similarity,davies_bouldin,calinski_harabasz,integrated | NaN iscg,iscl,integrated | NaN mean_auc,mean_cohens_d | NaN mean_overlap | ok | ok | ok | NaN neighbor_contamination | ok | 💥 ValueError | ok | ok | ok | ok |
| singleton_type | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| two_types | ok | ok | NaN mean_auc,mean_cohens_d | NaN mean_overlap | NaN pmp | ok | ok | NaN neighbor_contamination | ok | ok | ok | ok | ok | ok |
| empty_marker_set | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| markers_off_panel | ok | ok | NaN mean_auc,mean_cohens_d | NaN mean_overlap | NaN pmp | ok | ok | NaN neighbor_contamination | ok | ok | ok | ok | ok | ok |
| zero_count_cells | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| unused_category | ok | ok | ok | ok | ok | ok | ok | ok | ok | ≠ n_suspicious | ok | ok | ok | ok |
| integer_labels | ok | ok | ok | ok | ok | ≠ nmp | ok | ok | ok | ok | ok | ok | ok | ok |
| nan_labels | 💥 TypeError | 💥 TypeError | ok | ok | 💥 TypeError | 💥 TypeError | ok | 💥 TypeError | 💥 TypeError | 💥 ValueError | ok | ok | ok | 💥 ValueError |
| no_control_probes | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | NaN moran_gap,frac_above_ctrl,snr | ok |
| tiny_30_cells | ok | NaN iscg,iscl,integrated | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | NaN moran_gap | ok |
| duplicate_coords | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| all_zero_marker_gene | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | NaN moran_gap | ok |
| dense_matrix | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| unicode_labels | ok | ok | ok | ok | ok | ≠ nmp | ok | ok | ok | ok | ok | ok | ok | ok |
| extra_marker_type | ok | ok | ok | ok | ≠ crisp_purity,mecr,pmp | ok | ok | ok | ok | ok | ≠ mean_margin,mean_entropy | ok | ok | ok |
| no_spatial_coords | ok | ok | ok | ok | ok | ok | ok | 💥 KeyError | 💥 KeyError | 💥 KeyError | ok | ok | 💥 KeyError | ok |
| lognorm_X_no_counts | ok | ≠ iscg,integrated | ≠ mean_auc,mean_cohens_d | ok | ≠ pmp | ≠ nmp | ok | ≠ neighbor_contamination | ok | ok | ≠ mean_margin,mean_entropy | ok | ≠ moran_gap,frac_above_ctrl,snr,detection_entropy,tx_per_area | ok |

Errors:

- **one_label** / `neighborhood_sanity`: ValueError: Expected at least `2` clusters, found `1`.
- **nan_labels** / `internal_validity`: TypeError: '<' not supported between instances of 'str' and 'float'
- **nan_labels** / `inter_sample_consistency`: TypeError: '<' not supported between instances of 'str' and 'float'
- **nan_labels** / `purity`: TypeError: '<' not supported between instances of 'str' and 'float'
- **nan_labels** / `nmp`: TypeError: '<' not supported between instances of 'str' and 'float'
- **nan_labels** / `neighbor_contamination`: TypeError: '<' not supported between instances of 'str' and 'float'
- **nan_labels** / `spatial_coherence`: TypeError: '<' not supported between instances of 'str' and 'float'
- **nan_labels** / `neighborhood_sanity`: ValueError: Found `NaN` values in `adata.obs['label']`; every cell needs a cluster.
- **nan_labels** / `external`: ValueError: Input contains NaN
- **no_spatial_coords** / `neighbor_contamination`: KeyError: 'Spatial basis `spatial` not found in `adata.obsm`.'
- **no_spatial_coords** / `spatial_coherence`: KeyError: 'Spatial basis `spatial` not found in `adata.obsm`.'
- **no_spatial_coords** / `neighborhood_sanity`: KeyError: 'Spatial basis `spatial` not found in `adata.obsm`.'
- **no_spatial_coords** / `signal_qc`: KeyError: 'Spatial basis `spatial` not found in `adata.obsm`.'
