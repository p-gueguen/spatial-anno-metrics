# spatial (wp/spatial) - fragment for the integrator

## CHANGELOG

### Changed (breaking)

- `spatial_coherence` PAS now follows Shang & Zhou 2022 (doi 10.1038/s41467-022-34879-1, SpatialPCA
  `fx_PAS`): abnormal = the label differs from >= 6 of the 10 nearest neighbours.
  `n_neighs=15 -> 10`, `pas_threshold=0.2 -> 0.4`, comparison `frac < pas_threshold -> frac <= pas_threshold`.
  `spatial_coherence(adata, n_neighs=15, pas_threshold=0.2)` reproduces 0.6.x except for cells with exactly 20 %
  same-label neighbours, which now count as abnormal. On the Xenium breast section (167,780 cells, marker-argmax
  labels) the per-cell fractions are identical and PAS moves 0.153952 -> 0.205150, exactly the 8,590 cells
  (0.051198) at 0.2. Synthetic sweep (3 seeds): only `sc.pas` and `sc.mean_coherence` move (129 rows each,
  36 of 38 metrics identical, edge cases identical). The new PAS reads higher and responds more: clean layout
  0.000 -> 0.024 (border cells), `label_noise` 0.4 0.344 -> 0.553, `coord_shuffle` 1.0 0.512 -> 0.976,
  `over_split` 3 0.013 -> 0.291; it stays flat across `patch_error` (0.023-0.025) and falls with `merge`
  (0.024 -> 0.011).
- `neighbor_contamination`: `n_source_target` is NaN, not 0, for a pair with no evaluable target cell (as in
  SegTraQ); the frame is float.
- `neighbor_contamination`: a nested marker-set pair (one set inside the other, e.g. `T` in `CD4 T`) is
  skipped as hierarchical, as `mecr` skips it; 0.6.x counted the subtype's private markers as spillover. On the
  Xenium section with `T = {CD3E, CD3D}` and `CD4 T = {CD3E, CD3D, CD4}`, 86 cells change from a value to NaN
  and the `CD4 T -> T` entry becomes NaN; nothing else changes.
- `spatial_coherence`, `neighborhood_sanity`, `neighbor_contamination`: a spatial graph cached by another call
  with other settings (another `n_neighs`, `moran_signal`'s graph, a user's squidpy graph) is no longer
  reused silently; the requested kNN graph is built (and not cached over the existing one). A graph put in
  `obsp['spatial_connectivities']` by hand (no `uns['spatial_neighbors']`) is still used as it is.

### Added

- `spatial_coherence`: `chaos` (Shang & Zhou 2022 CHAOS: mean distance to the nearest same-label cell on
  per-axis standardised coordinates; NaN without `obsm['spatial']`) and `graph` in the result.
- `graph=` keyword (`None` | `"delaunay"` | `"obsp"`) on the three spatial functions; `"delaunay"` is
  SegTraQ's default graph, `"obsp"` forces the graph in `obsp`.
- `neighbor_contamination`: `baseline`, `n_baseline`, `excess` and `p_value` source x target frames (the
  same share in target cells with no source neighbour, the raw value minus it, and a one-sided permutation p
  over the target cells; `n_perm=199`, `seed=0`, keyword-only), `no_dict_frac` / `labels_without_markers`
  (also on its skips), `graph`, `n_perm`.
- `neighbor_contamination` accepts SegTraQ's marker format `{type: {"positive": [...], "negative": [...]}}`
  (mixable with plain lists): a target's listed negatives replace the implicit "every other type's marker".
- `thresholds.THRESHOLDS["spatial"]["contamination_n_perm"]` (199).

### Fixed

- `neighbor_contamination` and `spatial_coherence` crashed (`AttributeError: ... 'tocsc'`) on a dense
  ndarray in `obsp['spatial_connectivities']`.
- `neighbor_contamination` with int labels and an int-keyed `marker_sets` gave all-NaN and a 0x0 matrix;
  keys are now matched as strings.
- `neighbor_contamination` no longer densifies cells x marker genes (two dense matrices, one float, one
  bool): sparse throughout. 400,000 cells x 300 markers: peak traced memory 1.88 GB -> 0.14 GB,
  10.2 s -> 3.5 s with `n_perm=0` (9.8 s with the default 199 permutations).
- Docstrings: spatial coherence / PAS / CHAOS detect random label noise only (blind to a compact mislabelled
  patch, rise when types merge, fall for correctly labelled infiltrating cells).

## README rows

| `spatial_coherence(adata, label_key, n_neighs=10, pas_threshold=0.4)` | per-cell fraction of spatial neighbors sharing the label + dataset **PAS** (Percentage of Abnormal Spots: label differs from >= 6 of 10 neighbours) and **CHAOS** (mean distance to the nearest same-label cell), both as in Shang & Zhou 2022. Detects *random* label noise only: blind to a compact mislabelled patch, better when types are merged, worse for correctly labelled infiltrating cells. Down-weight only | Shang & Zhou 2022 (SpatialPCA); SpatialScribe QC L6 |
| `neighborhood_sanity(adata, label_key)` | flags types whose strongest neighborhood-enrichment partner is a *different* type (a mislabeling smell, e.g. two types forced together) via a squidpy permutation null | SpatialScribe QC L6 |
| `neighbor_contamination(adata, marker_sets, label_key)` | per-cell share of counts on markers of the **adjacent** cell types (not its own) + directed source → target frames: raw spillover, a **baseline** from same-type cells with no such neighbour, the **excess** over it and a permutation **p-value**, which separate spillover from real co-expression. Accepts SegTraQ's `{"positive", "negative"}` marker format; `graph="delaunay"` for SegTraQ's graph; reports labels without a marker set | SegTraQ `contamination_strength` |

## Catalog

Section **3e. Spatial-context metrics (spatial only)** of `docs/cell-annotation-quality-metrics.md`, replace
the first two bullets and the neighbor-contamination bullet, and the caveat, with:

- **Spatial coherence** = mean fraction of a cell's spatial neighbours (k = 10) sharing its label
  (`spatial.py` `spatial_coherence`).
- **PAS** (Percentage of Abnormal Spots; Shang & Zhou 2022, [10.1038/s41467-022-34879-1](https://doi.org/10.1038/s41467-022-34879-1))
  = % cells whose label differs from at least 6 of their 10 nearest neighbours (same-label fraction <= 0.4).
  0.6.x used < 20 % of 15 neighbours, which reads systematically lower.
- **CHAOS** (same paper) = mean distance from each cell to its nearest same-label cell, on per-axis
  standardised coordinates. Lower = spatially smoother labels.
- **Neighbor contamination** (`spatial.py` `neighbor_contamination`, from SegTraQ's `contamination_strength`)
  — per cell, the share of counts on markers of the cell types *touching it* that its own type lacks (with a
  touching cell of that type detecting them), plus a directed source → target matrix. Reference-free;
  separates spillover from neighbors from a mislabel, which PMP / NMP cannot. The raw value does not separate
  spillover from genuine co-expression: on a synthetic two-block grid, a B marker every A cell truly expresses
  scores 0.037 and real border spillover 0.041. The `baseline` (same-type cells with no such neighbour: 0.033
  vs 0.001) and the `excess` over it (0.004 vs 0.040) do; the one-sided permutation `p_value` is 0.001 for
  spillover and 0.19 for co-expression with no border excess (it is small, 0.003, when border cells carry
  even a +0.003 real excess: with large n read the excess, not the p). Nested marker sets (`T` in `CD4 T`)
  are skipped, as in MECR.

**Caveat (load-bearing):** coherence, PAS and CHAOS detect *random* label noise only. They are blind to a
compact mislabelled patch, they *rise* when two neighbouring types are merged into one label (0.904 → 0.961 in
the synthetic `merge` scenario, mean coherence 0.918 → 0.966 and PAS 0.024 → 0.011 from 0 to 3 merges), and
they fall for correctly labelled infiltrating cells (`infiltration`: PAS 0.024 → 0.163). PAS stays at
0.023-0.025 whatever the size of a compact mislabelled patch (`patch_error`). Use them to **down-weight**,
never to hard-delete or to rank annotations by correctness.

Also in the master table (row "Spatial coherence / **PAS**"): add CHAOS; the §3e source line "STEAM (spatial
coherence / PAS)" should cite Shang & Zhou 2022 for PAS and CHAOS.

## Roadmap

- **#2** PAS definition -> **done**. Defaults 10 / 0.4 inclusive (`spatial.py` `spatial_coherence`,
  `thresholds.py` spatial); `tests/test_spatial.py::test_pas_defaults_are_shang_zhou_2022` (matches a
  brute-force 10-NN ">= 6 differ" count), `::test_pas_threshold_is_inclusive`. CHAOS added:
  `::test_chaos_coherent_vs_shuffled` (vs a by-hand implementation of the paper's formula),
  `::test_chaos_nan_without_coordinates`. Old values: `n_neighs=15, pas_threshold=0.2` reproduce 0.6.x exactly
  except the cells at exactly 0.2 (Xenium: 8,590 / 167,780).
- **#19** neighbor_contamination baseline -> **done** (baseline, excess, permutation p).
  `::test_neighbor_contamination_baseline_separates_spillover_from_coexpression` (two-block probe: raw 0.037
  vs 0.041, baseline 0.033 vs 0.001, excess 0.004 vs 0.040, spillover p = 1/1000, co-expression without a
  border excess p = 0.19), `::test_contamination_p_value_is_one_sided`. On the Xenium section (10 types,
  default kNN 6) 75 of 100 pairs have p <= 0.01: with thousands of border cells the test is very powerful
  and spatially clustered cells make it optimistic, so the excess size is the readout. The
  SPLIT secondary-signal idea and CellSPA's distance curve stay **open**.
- **#8a** neighborhood_sanity single label / unused categories -> **done** (foundation; existing tests
  `test_spatial_label_metrics_skip_a_single_label`, `test_neighborhood_sanity_leaves_the_label_column_alone`).
  Positive control added: `::test_neighborhood_sanity_flags_two_types_forced_adjacent`.
- **#8c** dictionary coverage -> **done for neighbor_contamination** (`no_dict_frac`,
  `labels_without_markers`; `::test_neighbor_contamination_reports_dictionary_coverage`).
- **#8f** randomness-detector docstrings -> **done** (module docstring and `spatial_coherence`).
- **#15** dense hotspots -> **partial**: `neighbor_contamination` is sparse (1.88 -> 0.14 GB at 400k cells x
  300 markers); the other hotspots belong to other modules.
- WF1 review defects (spatial.py): dense obsp crash, `n_source_target` 0 vs NaN, int keys, nested sets ->
  fixed; OR across sources, graph direction (asymmetric graph), matrix mean -> pinned by tests
  (`::test_neighbor_contamination_or_across_sources`, `::test_neighbor_contamination_graph_direction`,
  `::test_neighbor_contamination_matrix_is_the_mean`); 11/11 mutants of the new code caught by
  tests/test_spatial.py (OR overwritten, graph transposed, median, nested counted, strict PAS, int keys,
  baseline over all cells, unscaled CHAOS, cache reused whatever its settings, two-sided p).

## tests/test_thresholds.py DEFAULTS row (for the integrator)

    ("spatial", "neighbor_contamination", "n_perm"): ("spatial", "contamination_n_perm"),
