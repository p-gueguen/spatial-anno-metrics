| verdict (cells with a claim) | toy (5 seeds) | SRTsim (3 seeds) | toy (3 seeds) |
|---|---|---|---|
| detects | 134 | 112 | 134 |
| BLIND | 61 | 74 | 58 |
| FALSE-ALARM | 22 | 16 | 21 |
| WRONG-WAY | 3 | 6 | 3 |
| erratic | 35 | 41 | 38 |
| stable | 280 | 282 | 276 |

133 of 765 (scenario, metric) cells change verdict (23 to detects/stable, 43 away from it, 67 other); 3 of them also differ between the toy at 5 and at 3 seeds.

| scenario | metric | toy | SRTsim | toy v0 → vmax | SRTsim v0 → vmax | seeds? |
|---|---|---|---|---|---|---|
| ambient | `sig.moran_gap` | BLIND | detects | 0.0934 → 0.0857 | 0.0922 → 0.0528 | |
| batch | `isc.iscl` | BLIND | detects | 1 → 1 | 0.893 → 0.676 | |
| coexpression | `crisp_purity` | FALSE-ALARM | stable | 0.813 → 0.685 | 0.196 → 0.182 | |
| coexpression | `iso.asw` | FALSE-ALARM | stable | 0.701 → 0.684 | 0.585 → 0.585 | |
| coexpression | `iv.silhouette` | FALSE-ALARM | stable | 0.404 → 0.369 | 0.153 → 0.151 | |
| coexpression | `mecr` | FALSE-ALARM | stable | 0.0417 → 0.0653 | 0.284 → 0.298 | |
| coexpression | `onto.rho` | FALSE-ALARM | stable | 0.255 → -0.0534 | 0.598 → 0.556 | |
| coexpression | `onto.rwr_rho` | FALSE-ALARM | stable | 0.236 → -0.0321 | 0.538 → 0.499 | |
| coexpression | `pmp` | FALSE-ALARM | stable | 0.995 → 0.95 | 0.628 → 0.617 | |
| coord_shuffle | `ns.n_suspicious` | erratic | detects | 0 → 5 | 1 → 4.67 | |
| infiltration | `sc.mean_coherence` | FALSE-ALARM | stable | 0.92 → 0.734 | 0.452 → 0.436 | |
| infiltration | `sc.pas` | FALSE-ALARM | stable | 0.0223 → 0.159 | 0.537 → 0.539 | |
| label_noise | `ext.comp_l1` | BLIND | detects | 0 → 0.0897 | 0 → 0.254 | |
| patch_error | `clisi` | BLIND | detects | 1 → 1 | 0.961 → 0.857 | |
| patch_error | `isc.iscl` | BLIND | detects | 1 → 1 | 0.893 → 0.753 | |
| rare_type | `coex.global` | FALSE-ALARM | stable | 0.978 → 0.974 | 0.979 → 0.978 | |
| rare_type | `iso.asw` | FALSE-ALARM | stable | 0.701 → 0.641 | 0.585 → 0.583 | |
| rare_type | `iso.f1` | FALSE-ALARM | stable | 1 → 0.955 | 0.744 → 0.693 | |
| rare_type | `iv.integrated` | FALSE-ALARM | stable | 0.878 → 0.855 | 0.725 → 0.726 | |
| rare_type | `iv.silhouette` | FALSE-ALARM | stable | 0.404 → 0.31 | 0.153 → 0.153 | |
| rare_type | `iv.ward_propmatch` | FALSE-ALARM | stable | 1 → 0.942 | 0.748 → 0.714 | |
| spillover | `coex.by_type` | BLIND | detects | 0.967 → 0.967 | 0.967 → 0.966 | |
| type_swap | `ext.comp_l1` | BLIND | detects | 0 → 0.0672 | 0 → 0.481 | |
| ambient | `clisi` | flat | worse | 1 → 1 | 0.961 → 0.929 | |
| ambient | `isc.iscl` | flat | worse | 1 → 1 | 0.893 → 0.625 | |
| ambient | `iv.avg_similarity` | worse | flat | 0.804 → 0.748 | 0.603 → 0.551 | |
| ambient | `iv.calinski_harabasz` | worse | flat | 714 → 498 | 330 → 249 | |
| ambient | `iv.neighborhood_purity` | flat | worse | 1 → 0.999 | 0.795 → 0.742 | |
| ambient | `iv.orbital_medoid` | flat | worse | 1 → 1 | 0.818 → 0.747 | |
| ambient | `iv.silhouette` | worse | flat | 0.404 → 0.318 | 0.153 → 0.121 | |
| ambient | `mgo.mean_overlap` | flat | erratic | 1 → 1 | 0.653 → 0.75 | |
| ambient | `pmp` | worse | flat | 0.995 → 0.918 | 0.628 → 0.508 | |
| ambient | `sig.tx_per_area` | flat | better | 3 → 3.42 | 1.05 → 1.76 | |
| ambient | `stab.flip_rate` | flat | worse | 0 → 2.78e-05 | 0.0609 → 0.105 | |
| batch | `clisi` | flat | worse | 1 → 1 | 0.961 → 0.916 | |
| batch | `conf.mean_entropy` | flat | better | 0.0214 → 0.0555 | 0.653 → 0.457 | |
| batch | `crisp_purity` | flat | better | 0.813 → 0.883 | 0.196 → 0.305 | yes |
| batch | `ext.avg_bio` | worse | flat | 0.901 → 0.871 | 0.859 → 0.849 | |
| batch | `fid.mean_auc` | flat | worse | 0.995 → 0.986 | 0.897 → 0.842 | |
| batch | `iv.neighborhood_purity` | flat | worse | 1 → 0.994 | 0.795 → 0.715 | |
| batch | `mecr` | flat | better | 0.0417 → 0.0256 | 0.284 → 0.212 | |
| batch | `nov.frac_novel` | worse | flat | 0 → 0.998 | 0 → 0 | |
| batch | `sig.moran_gap` | better | flat | 0.0934 → 0.402 | 0.0922 → 0.116 | |
| depth_loss | `coex.by_type` | flat | worse | 0.967 → 0.965 | 0.967 → 0.949 | yes |
| depth_loss | `isc.iscl` | flat | worse | 1 → 1 | 0.893 → 0.554 | |
| depth_loss | `ncp.moldia` | flat | worse | 1 → 1 | 0.986 → 0.934 | |
| depth_loss | `neighbor_contamination` | flat | better | 0.000292 → 0.000333 | 0.023 → 0.0102 | |
| depth_loss | `sig.frac_above_ctrl` | erratic | flat | 0.866 → 0.712 | 0.665 → 0.648 | |
| depth_loss | `stab.flip_rate` | worse | flat | 0 → 0.0462 | 0.0609 → 0.0585 | |
| label_noise | `ns.n_suspicious` | BLIND | erratic | 0 → 0 | 1 → 1.67 | |
| label_noise | `onto.rho` | erratic | BLIND | 0.255 → 0.465 | 0.598 → 0.66 | |
| label_noise | `onto.rwr_rho` | erratic | BLIND | 0.236 → 0.436 | 0.538 → 0.609 | |
| merge | `isc.integrated` | WRONG-WAY | BLIND | 0.937 → 0.978 | 0.485 → 0.564 | |
| merge | `isc.iscg` | WRONG-WAY | BLIND | 0.937 → 0.978 | 0.525 → 0.578 | |
| merge | `iv.neighborhood_purity` | BLIND | WRONG-WAY | 1 → 1 | 0.795 → 0.864 | |
| merge | `mgo.mean_overlap` | BLIND | erratic | 1 → nan | 0.653 → nan | |
| merge | `neighbor_contamination` | BLIND | erratic | 0.000292 → nan | 0.023 → nan | |
| merge | `ns.n_suspicious` | BLIND | WRONG-WAY | 0 → 0 | 1 → 0 | |
| merge | `pmp` | BLIND | erratic | 0.995 → nan | 0.628 → nan | |
| merge | `sc.chaos` | BLIND | WRONG-WAY | 0.036 → 0.0356 | 0.0522 → 0.0415 | |
| over_split | `mgo.mean_overlap` | BLIND | erratic | 1 → 1 | 0.653 → 0.861 | |
| over_split | `neighbor_contamination` | BLIND | erratic | 0.000292 → 0.000249 | 0.023 → 0.00966 | |
| panel_growth | `sig.detection_entropy` | better | flat | 0.902 → 0.943 | 0.929 → 0.919 | |
| panel_growth | `sig.frac_above_ctrl` | flat | worse | 0.866 → 0.818 | 0.665 → 0.222 | |
| panel_growth | `sig.zero_fraction` | better | flat | 0.573 → 0.383 | 0.834 → 0.814 | |
| patch_error | `onto.rho` | BLIND | WRONG-WAY | 0.255 → 0.149 | 0.598 → 0.79 | |
| patch_error | `onto.rwr_rho` | BLIND | WRONG-WAY | 0.236 → 0.202 | 0.538 → 0.768 | |
| sibling_confusion | `isc.iscl` | BLIND | erratic | 1 → 1 | 0.893 → 0.996 | |
| sibling_confusion | `ns.n_suspicious` | erratic | BLIND | 0 → 1 | 1 → 1 | |
| sibling_confusion | `sc.mean_coherence` | erratic | BLIND | 0.92 → 0.883 | 0.452 → 0.484 | |
| sibling_confusion | `sc.pas` | erratic | BLIND | 0.0223 → 0.0571 | 0.537 → 0.479 | |
| spillover | `cal.ece` | worse | erratic | 0.000566 → 0.0093 | 0.0889 → 0.127 | |
| spillover | `clisi` | flat | worse | 1 → 1 | 0.961 → 0.877 | |
| spillover | `fid.mean_cohens_d` | better | flat | 6.14 → 7.58 | 2.21 → 2.12 | |
| spillover | `isc.integrated` | flat | worse | 0.937 → 0.91 | 0.485 → 0.284 | |
| spillover | `isc.iscl` | flat | worse | 1 → 1 | 0.893 → 0.604 | |
| spillover | `iso.f1` | worse | flat | 1 → 0.954 | 0.744 → 0.602 | |
| spillover | `iv.integrated` | flat | worse | 0.878 → 0.861 | 0.725 → 0.646 | |
| spillover | `iv.neighborhood_purity` | flat | worse | 1 → 0.983 | 0.795 → 0.661 | |
| spillover | `iv.orbital_medoid` | flat | worse | 1 → 0.996 | 0.818 → 0.662 | |
| spillover | `iv.silhouette` | worse | flat | 0.404 → 0.364 | 0.153 → 0.107 | |
| spillover | `iv.ward_propmatch` | worse | erratic | 1 → 0.956 | 0.748 → 0.621 | |
| spillover | `ncp` | erratic | BLIND | 1 → 0.94 | 0.997 → 0.996 | |
| spillover | `onto.rho` | better | flat | 0.255 → 0.619 | 0.598 → 0.661 | |
| spillover | `onto.rwr_rho` | better | flat | 0.236 → 0.586 | 0.538 → 0.602 | |
| spillover | `sig.moran_gap` | flat | better | 0.0934 → 0.0974 | 0.0922 → 0.154 | |
| spillover | `stab.flip_rate` | flat | worse | 0 → 0.000944 | 0.0609 → 0.0986 | |
| type_swap | `clisi` | BLIND | erratic | 1 → 1 | 0.961 → 0.961 | |
| type_swap | `iso.f1` | erratic | BLIND | 1 → 1 | 0.744 → 0.742 | |
| type_swap | `iv.ward_propmatch` | erratic | BLIND | 1 → 1 | 0.748 → 0.748 | |
| merge | `iso.asw` | detects | erratic | 0.701 → 0.627 | 0.585 → 0.563 | |
| merge | `iso.f1` | detects | BLIND | 1 → 0.852 | 0.744 → 0.738 | |
| merge | `iv.integrated` | detects | erratic | 0.878 → 0.779 | 0.725 → 0.7 | |
| merge | `iv.orbital_medoid` | detects | erratic | 1 → 0.867 | 0.818 → 0.823 | yes |
| merge | `iv.silhouette` | detects | erratic | 0.404 → 0.25 | 0.153 → 0.11 | |
| merge | `iv.ward_propmatch` | detects | BLIND | 1 → 0.849 | 0.748 → 0.714 | |
| over_split | `clisi` | detects | erratic | 1 → 0.909 | 0.961 → 0.879 | |
| over_split | `fid.mean_cohens_d` | detects | erratic | 6.14 → 5.54 | 2.21 → 2.02 | |
| over_split | `iso.f1` | detects | BLIND | 1 → 0.778 | 0.744 → 0.596 | |
| over_split | `iv.ward_propmatch` | detects | BLIND | 1 → 0.782 | 0.748 → 0.644 | |
| panel_growth | `clisi` | stable | FALSE-ALARM | 1 → 1 | 0.961 → 0.903 | |
| panel_growth | `coex.global` | stable | drifts-up | 0.978 → 0.979 | 0.979 → 0.987 | |
| panel_growth | `isc.integrated` | stable | FALSE-ALARM | 0.937 → 0.925 | 0.485 → 0.187 | |
| panel_growth | `isc.iscg` | stable | FALSE-ALARM | 0.937 → 0.925 | 0.525 → 0.266 | |
| panel_growth | `isc.iscl` | stable | FALSE-ALARM | 1 → 1 | 0.893 → 0.582 | |
| panel_growth | `iv.davies_bouldin` | stable | FALSE-ALARM | 1.06 → 0.824 | 2.23 → 3.16 | |
| panel_growth | `iv.neighborhood_purity` | stable | FALSE-ALARM | 1 → 1 | 0.795 → 0.706 | |
| panel_growth | `iv.orbital_medoid` | stable | FALSE-ALARM | 1 → 1 | 0.818 → 0.7 | |
| panel_growth | `iv.ward_propmatch` | stable | FALSE-ALARM | 1 → 1 | 0.748 → 0.532 | |
| panel_growth | `ncp.moldia` | stable | drifts-up | 1 → 1 | 0.986 → 1 | |
| panel_growth | `neighbor_contamination` | stable | drifts-up | 0.000292 → 2.27e-05 | 0.023 → 0.00439 | |
| patch_error | `coex.by_type` | detects | erratic | 0.967 → 0.964 | 0.967 → 0.963 | |
| patch_error | `iso.f1` | detects | BLIND | 1 → 0.941 | 0.744 → 0.643 | |
| patch_error | `iv.ward_propmatch` | detects | BLIND | 1 → 0.951 | 0.748 → 0.653 | |
| patch_error | `neighbor_contamination` | detects | BLIND | 0.000292 → 0.0361 | 0.023 → 0.0234 | |
| patch_error | `pmp` | detects | BLIND | 0.995 → 0.9 | 0.628 → 0.571 | |
| rare_type | `iv.neighborhood_purity` | stable | drifts-up | 1 → 0.992 | 0.795 → 0.851 | |
| rare_type | `iv.orbital_medoid` | stable | FALSE-ALARM | 1 → 0.99 | 0.818 → 0.809 | |
| sibling_confusion | `ext.hier_acc` | detects | BLIND | 1 → 0.945 | 1 → 0.961 | |
| sibling_confusion | `iso.asw` | detects | BLIND | 0.701 → 0.667 | 0.585 → 0.579 | |
| sibling_confusion | `iso.f1` | detects | BLIND | 1 → 0.907 | 0.744 → 0.683 | |
| sibling_confusion | `iv.avg_similarity` | detects | BLIND | 0.804 → 0.764 | 0.603 → 0.595 | |
| sibling_confusion | `iv.calinski_harabasz` | detects | BLIND | 714 → 567 | 330 → 318 | |
| sibling_confusion | `iv.integrated` | detects | BLIND | 0.878 → 0.827 | 0.725 → 0.721 | |
| sibling_confusion | `iv.neighborhood_purity` | detects | BLIND | 1 → 0.944 | 0.795 → 0.815 | |
| sibling_confusion | `iv.silhouette` | detects | BLIND | 0.404 → 0.318 | 0.153 → 0.143 | |
| sibling_confusion | `iv.ward_propmatch` | detects | BLIND | 1 → 0.932 | 0.748 → 0.757 | |
| sibling_confusion | `neighbor_contamination` | detects | BLIND | 0.000292 → 0.0527 | 0.023 → 0.0228 | |
| sibling_confusion | `nmp` | detects | BLIND | 1 → 0.968 | 0.978 → 0.956 | |
| sibling_confusion | `pmp` | detects | BLIND | 0.995 → 0.941 | 0.628 → 0.614 | |
| spillover | `ncp.moldia` | detects | BLIND | 1 → 0.962 | 0.986 → 0.977 | |
| spillover | `pmp` | detects | BLIND | 0.995 → 0.973 | 0.628 → 0.533 | |
| type_swap | `mgo.mean_overlap` | detects | BLIND | 1 → 0.667 | 0.653 → 0.514 | |
