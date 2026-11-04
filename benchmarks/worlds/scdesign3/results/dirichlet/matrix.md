| metric | label_noise<br>*label* | sibling_confusion<br>*label* | type_swap<br>*label* | patch_error<br>*label* | over_split<br>*label* | merge<br>*label* | spillover<br>*data* | depth_loss<br>*data* | ambient<br>*data* | batch<br>*data* | coord_shuffle<br>*data* | coexpression<br>*null* | rare_type<br>*null* | infiltration<br>*null* | panel_growth<br>*null* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `cal.ece` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↓ |  |  |  |  |  |
| `clisi` | ✅ 0.05 |  |  |  | 〰️ | ❌ |  |  |  | ↓ |  |  |  |  |  |
| `coex.by_type` | ✅ 0.4 |  |  |  | ❌ | ❌ |  |  |  | ↓ |  |  |  |  |  |
| `coex.global` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↓ |  |  |  |  |  |
| `conf.mean_entropy` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↑ |  |  |  |  |  |
| `conf.mean_margin` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↑ |  |  |  |  |  |
| `crisp_purity` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↑ |  |  |  |  |  |
| `ext.ari` | ✅ 0.05 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `ext.avg_bio` | ✅ 0.05 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `ext.comp_l1` | ✅ 0.2 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `ext.ecs` | ✅ 0.05 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `ext.hier_acc` | ✅ 0.1 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `ext.macro_f1` | ✅ 0.05 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `ext.spari` | ✅ 0.05 |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `fid.mean_auc` | ✅ 0.05 |  |  |  | ✅ 1 | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `fid.mean_cohens_d` | ✅ 0.05 |  |  |  | ❌ | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `isc.integrated` | ✅ 0.1 |  |  |  | ✅ 2 | 〰️ |  |  |  | ✅ 0.5 |  |  |  |  |  |
| `isc.iscg` | ✅ 0.05 |  |  |  | 〰️ | 〰️ |  |  |  | ✅ 0.25 |  |  |  |  |  |
| `isc.iscl` | ✅ 0.4 |  |  |  | ✅ 2 | ❌ |  |  |  | ❌ |  |  |  |  |  |
| `iso.asw` | ✅ 0.05 |  |  |  | ✅ 1 | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `iso.f1` | ✅ 0.05 |  |  |  | ✅ 1 | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `iv.avg_similarity` | ✅ 0.05 |  |  |  | ❌ | ✅ 1 |  |  |  | ↓ |  |  |  |  |  |
| `iv.calinski_harabasz` | ✅ 0.05 |  |  |  | ✅ 1 | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `iv.davies_bouldin` | ✅ 0.1 |  |  |  | ✅ 1 | ❌ |  |  |  | · |  |  |  |  |  |
| `iv.integrated` | ✅ 0.05 |  |  |  | ✅ 1 | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `iv.neighborhood_purity` | ✅ 0.05 |  |  |  | ✅ 1 | ⛔ |  |  |  | ↓ |  |  |  |  |  |
| `iv.orbital_medoid` | ✅ 0.1 |  |  |  | ✅ 1 | ❌ |  |  |  | ↓ |  |  |  |  |  |
| `iv.silhouette` | ✅ 0.05 |  |  |  | ✅ 1 | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `iv.ward_propmatch` | ✅ 0.1 |  |  |  | ✅ 2 | ❌ |  |  |  | 〰️ |  |  |  |  |  |
| `mecr` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↑ |  |  |  |  |  |
| `mgo.mean_overlap` | ✅ 0.4 |  |  |  | ⛔ | ❌ |  |  |  | ↓ |  |  |  |  |  |
| `ncp` | ✔ |  |  |  | ✔ | ✔ |  |  |  | · |  |  |  |  |  |
| `ncp.moldia` | ✔ |  |  |  | ✔ | ✔ |  |  |  | 〰️ |  |  |  |  |  |
| `neighbor_contamination` | ✅ 0.05 |  |  |  | 〰️ | 〰️ |  |  |  | 〰️ |  |  |  |  |  |
| `nmp` | ✅ 0.4 |  |  |  | ❌ | 〰️ |  |  |  | · |  |  |  |  |  |
| `nov.frac_novel` | ✔ |  |  |  | ✔ | ✔ |  |  |  | · |  |  |  |  |  |
| `ns.n_suspicious` | ❌ |  |  |  | ✅ 1 | ❌ |  |  |  | · |  |  |  |  |  |
| `onto.rho` | ❌ |  |  |  | ❌ | 〰️ |  |  |  | · |  |  |  |  |  |
| `onto.rwr_rho` | ❌ |  |  |  | ❌ | 〰️ |  |  |  | · |  |  |  |  |  |
| `pmp` | ✅ 0.05 |  |  |  | 〰️ | 〰️ |  |  |  | ↓ |  |  |  |  |  |
| `pmp.no_dict_frac` | ❌ |  |  |  | ✅ 1 | ✅ 1 |  |  |  | · |  |  |  |  |  |
| `sc.chaos` | ✅ 0.05 |  |  |  | ✅ 1 | ⛔ |  |  |  | · |  |  |  |  |  |
| `sc.mean_coherence` | ✅ 0.1 |  |  |  | ✅ 1 | ⛔ |  |  |  | · |  |  |  |  |  |
| `sc.pas` | ✅ 0.1 |  |  |  | ✅ 1 | ⛔ |  |  |  | · |  |  |  |  |  |
| `sig.detection_entropy` | · |  |  |  | · | · |  |  |  | ↓ |  |  |  |  |  |
| `sig.frac_above_ctrl` | ✔ |  |  |  | ✔ | ✔ |  |  |  | · |  |  |  |  |  |
| `sig.moran_gap` | ✔ |  |  |  | ✔ | ✔ |  |  |  | 〰️ |  |  |  |  |  |
| `sig.snr` | ✔ |  |  |  | ✔ | ✔ |  |  |  | 〰️ |  |  |  |  |  |
| `sig.tx_per_area` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↑ |  |  |  |  |  |
| `sig.zero_fraction` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↓ |  |  |  |  |  |
| `stab.flip_rate` | ✔ |  |  |  | ✔ | ✔ |  |  |  | ↑ |  |  |  |  |  |
