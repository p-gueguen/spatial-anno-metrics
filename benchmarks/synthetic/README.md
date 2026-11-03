# Synthetic sweep - which metrics measure what they claim

`sim.py` simulates a section with known truth: six types in tissue stripes (T_CD4 | T_CD8 | B | Myeloid |
Fibro | Epi), 2,400 cells, 150 counts per cell. Each type has private markers, and the two T subtypes
share two lineage markers. There are 60 background genes and 20 negative-control probes, and six samples
as horizontal bands. A matched dissociated reference comes from the same world. `sweep.py` damages that
section along 15 axes at 4-5 strengths × 5 seeds, runs every public metric (49 values from 19 adapters), and
scores each trajectory against what the metric should do. It also runs 18 edge cases and four dedicated
checks.

```bash
pip install -e ".[spatial,dev]"
python benchmarks/synthetic/sweep.py --seeds 5 --workers 32     # ~1 min on 32 cores
```

**Results (0.7.0, 5 seeds):**
- [`results/matrix.md`](results/matrix.md): the verdict per metric and scenario.
- [`results/verdicts.csv`](results/verdicts.csv): value at baseline and at maximum strength, Spearman ρ,
  the threshold, and the first strength that moves the metric.
- [`results/values.csv`](results/values.csv): the raw values, with the number of `UserWarning`s per adapter.
- [`results/edge_cases.md`](results/edge_cases.md): what each metric returns on the edge cases.
- [`results/check_*.csv`](results/): the four dedicated checks.

The run also writes `errors.csv` (every adapter error with its traceback) and `warnings.csv`; they hold
local paths and are not committed.

## Scenarios

| kind | scenario | strength | what changes |
|---|---|---|---|
| labels wrong | `label_noise` | 0-40 % | random cells get a random wrong type |
| | `sibling_confusion` | 0-75 % | T_CD8 called T_CD4 (right lineage, wrong subtype) |
| | `type_swap` | 0-100 % | Fibro ↔ Epi exchanged (at 100 %: names swapped, partition intact) |
| | `patch_error` | 0-50 % | a compact patch of Fibro at the Epi border called Epi |
| | `over_split` | 0-3 types | a type split into two random halves |
| | `merge` | 0-3 pairs | two types given one label (the third pair is the T siblings) |
| data degraded | `spillover` | 0-40 % | each cell takes that share of counts from its 6 spatial neighbours |
| | `depth_loss` | 0-94 % | binomial thinning |
| | `ambient` | 0-0.2 | Poisson background on every gene, controls included |
| | `batch` | log-sd 0-1.5 | per-sample gene-wise expression shift |
| | `coord_shuffle` | 0-100 % | coordinates permuted among cells |
| nothing wrong | `coexpression` | 0-40 % | Myeloid truly expresses Fibro markers (also in the reference) |
| | `rare_type` | 15 → 1 % | T_CD8 becomes rare |
| | `infiltration` | 0-100 % | T_CD8 cells scattered through the tissue |
| | `panel_growth` | +0-960 genes | uninformative genes added (counts on the old genes unchanged) |

**Expected behaviour:**
- When labels are wrong, every label-aware metric should get worse, and label-blind metrics should stay
  flat.
- When data is degraded, only the metrics aimed at that damage are held to it, for example purity under
  spillover and signal QC under ambient noise. Other metrics are reported without a verdict.
- When nothing is wrong, label-aware and purity metrics should stay flat, and so should the two
  annotator checks (`cal.ece`, `nov.frac_novel`).

**Adapters added in 0.7.0:** `isolated_label_scores` (`iso.asw`, `iso.f1`), `clisi`, `ontology_consistency`
(`onto.rho`, on a 10-term toy ontology that mirrors the sim hierarchy: split halves map to their type's term,
a merged pair to its lowest common ancestor), `coexpression_similarity` (`coex.global`, and `coex.by_type`
with `label_key`), `spatial_ari` (`ext.spari`), CHAOS (`sc.chaos`, from `spatial_coherence`), PMP's
dictionary coverage (`pmp.no_dict_frac`), and `calibration_novelty`: a logistic classifier trained on the
reference annotates the section, scored by `calibration_error` against the truth (`cal.ece`) and by the
share of section cells `conformal_novelty` flags (`nov.frac_novel`, reference thinned to the section's
per-cell depths). Neither of the last two reads `obs['label']`, and no main scenario holds a novel type,
so both are held to "stay flat" when only labels change or nothing is wrong; the novel-type case is in
`check_conformal.csv`.

**What counts as a move:** a change of more than 3 × the seed SD at strength 0, and also more than 5 % of
the metric's range across the whole sweep.

**Legend:**
- ✅ detects (with the first strength at which it moves)
- ❌ blind
- ⛔ moves the wrong way
- 〰️ moves, but not monotonically
- ⚠️ false alarm
- ✔ stable as it should be
- ↓ / ↑ / · moves worse / better / not at all, where no claim is made
- 💥 (edge cases) the adapter raised

## Since 0.6.0

Counted over the cells that carry a claim (expectation "down" or "flat"):

| | 0.6.0 (38 metrics) | 0.7.0, same 38 metrics | 0.7.0, 11 new metrics | 0.7.0, all 49 |
|---|---|---|---|---|
| ✅ detects | 106 | 108 | 25 | 133 |
| ❌ blind | 44 | 40 | 18 | 58 |
| ⚠️ false alarm | 13 | 13 | 8 | 21 |
| ⛔ wrong way | 3 | 3 | 0 | 3 |
| 〰️ erratic | 22 | 24 | 8 | 32 |
| 💥 edge-case raises | 13 (9 deep `TypeError` / `ValueError`) | 24 (0 deep) | 6 (0 deep) | 30 (0 deep) |

Six cells of the 38 shared metrics changed verdict; the other 564 did not. A rerun of the same commit
reproduces all 570 shared cells exactly.

- **`nmp` now catches spillover** (❌ → ✅ 0.4). The default is the Salas et al. score: 1.000 → 0.969 at
  40 % spillover; the 0.6 raw score went 0.0013 → 0.0079, under its threshold. It still stays flat under
  real co-expression (1.000 → 1.000).
- **PAS follows Shang & Zhou** (10 neighbours, abnormal at ≤ 40 % same-label): it now catches
  `over_split` (❌ → ✅ 1; 0.022 → 0.288, was 0 → 0.012), and `sibling_confusion` and `type_swap` move
  (❌ → 〰️). The clean baseline is 0.022, not 0.
- **SNR no longer reads the `eps` floor:** 5.12 on the clean section (was 30.36), 4.91 at an ambient of 0.02.
  Its `batch` trend is now 〰️ (no claim).
- **`sig.frac_above_ctrl` under `depth_loss`:** ↓ → 〰️ (no claim); 0.866 → 0.712 at 94 % (was 0.439).
- **Marker metrics disclose labels with no marker set.** At `merge` 3 every label is unknown to the
  dictionary: `pmp`, `nmp`, `fid.*`, `mgo` and `neighbor_contamination` now skip (NaN) where 0.6 returned a
  clean-looking value (`nmp` 0, which read as perfect). `pmp.no_dict_frac` reports the share of affected cells:
  0.23 / 0.42 / 0.58 for `over_split` 1-3 and 0.42 / 0.70 / 1.00 for `merge` 1-3.
- **Edge cases fail loudly instead of deep inside or silently** (see below): NaN labels raise one
  `ValueError` naming the column and the count; log-normalised `X` without a counts layer raises instead
  of returning different numbers; a single label is a skip, not a crash; an unused category no longer
  changes `neighborhood_sanity`; and `nmp` on a label vocabulary the reference does not share is a skip, not 0.

## Findings

**Each metric is blind to some kind of error, so read them as a panel, never one alone.**

1. **No single metric detects all six label errors**, external metrics aside.
   - Internal validity (silhouette, `integrated`, ...) and the isolated-label ASW and F1 catch 5/6. They miss
     the name swap, which no partition metric can see.
   - `pmp` and `nmp` catch 4/6. All marker metrics (`pmp`, `nmp`, fidelity, `neighbor_contamination`) miss
     `over_split` and `merge`: the new labels have no entry in `marker_sets`, so the score is computed on the
     remaining types. Since 0.7 this is reported: `pmp.no_dict_frac` catches both at strength 1, so `pmp`
     read with its coverage catches 6/6.
   - cLISI catches 2/6 (`label_noise` from 10 %, `over_split` at 3 types). It is a median over cells, capped
     at 1: it stays 1.000 under 75 % sibling confusion (T_CD8 is 15 % of cells), a 50 % patch error, a full
     swap and every merge.
2. **Spatial coherence, PAS and CHAOS see random errors only.**
   - A compact mislabelled patch (50 % of Fibro called Epi) leaves coherence at 0.920 vs 0.920, PAS at 0.023
     vs 0.022 and CHAOS at 0.0360 vs 0.0360.
   - Merging types *raises* coherence (0.920 → 0.967), and lowers PAS (0.022 → 0.009) and CHAOS
     (0.0360 → 0.0356).
   - Correctly labelled infiltrating cells lower coherence (0.920 → 0.734) and raise PAS (0.022 → 0.159) and
     CHAOS (0.036 → 0.044).
   - Random noise is caught early: PAS 0.022 → 0.077 at 5 % label noise, CHAOS 0.036 → 0.052.
3. **ISC regime calls are right for over-splitting only.**
   - Over-split: the split halves get `merge` in 5/5 seeds (ISCL −0.01 on average over the halves). Correct.
   - Merge: merged labels get `consistent` in 5/5, and ISCG *rises* (0.938 → 0.978). It should say `split`.
   - 20 % random label noise reads `split` for 21 of 30 type-seeds, and a batch effect (log-sd 1.0) for
     30 of 30.
   - ISCL stays ≥ 0.90 under 40 % label noise and at 1.000 under every batch strength; on this fixture it
     only reacts to over-splitting.
4. **Real co-expression is flagged as contamination** by every reference-free purity metric:
   - `crisp_purity` 0.81 → 0.70 and `mecr` 0.042 → 0.064, already at 5 % co-expression;
   - `pmp` 0.996 → 0.950;
   - `nmp` stays flat (1.000), because the reference carries the co-expression.
5. **`crisp_purity` is 0.81 on perfectly clean data.** The cause is an off-target rate of 0.1 % of
   the marker rate, summed over 20 markers from other lineages, since the rule is "≥ 1 count". That baseline
   depends on the panel and on depth, so absolute CRISP values do not compare across sections; read its
   `vs_null` instead.
6. **Contamination:** `neighbor_contamination` moves with spillover (0.0003 → 0.036 at 40 %) and stays at 0.011
   or below under co-expression. `nmp` falls 1.000 → 0.969, `ncp` moves only at 40 % (1.00 → 0.94, 〰️).
   Co-expression similarity is blind to it (`coex.global` 0.978 → 0.977, `coex.by_type` 0.967 → 0.967): in
   stripes, a cell's six neighbours are mostly its own type, so the bleed adds same-type counts and leaves
   gene-gene correlations almost unchanged.
7. **`marker_gene_overlap` stays at 1.0 at 40 % label noise**, because the top DE genes survive random
   noise. It reacts only to a type swap (1.00 → 0.67 at 75 %).
8. **Partition metrics cannot see a name swap.** ARI, spARI, ECS, silhouette, the isolated-label scores and
   ISC return to their clean values when Fibro and Epi are fully swapped (〰️); cLISI never moves. Marker metrics, macro-F1
   and `hier_acc` catch it. `composition_accuracy` is blind to errors that keep the proportions balanced.
   spARI tracks ARI and gives a compact patch error partial credit for distance (0.826 vs ARI 0.795 at 50 %).
9. **Signal QC:**
   - SNR falls 5.12 → 4.91 at an ambient of 0.02 per gene (✅ 0.02).
   - Moran's gap barely moves under ambient noise (0.093 → 0.086, threshold 0.033), and the fraction of genes
     above the control p95 does not move (0.866 → 0.866).
   - Both behave correctly under coordinate shuffling; depth loss is caught by `tx_per_area` (3.00 → 0.18)
     and `zero_fraction` (0.57 → 0.92).
10. **Silhouette-type metrics move when composition or panel changes**, even with correct labels. Do not
    compare them across sections with different composition or panels.
    - Rare T_CD8 (1 %): silhouette 0.404 → 0.310, `integrated` 0.878 → 0.855, Ward-PropMatch 1.00 → 0.94,
      isolated-label ASW 0.701 → 0.641 (a macro average weighs the rare type as much as the others) and F1
      1.00 → 0.95.
    - 960 uninformative genes: silhouette 0.404 → 0.496, Calinski-Harabasz 715 → 1,005, average similarity
      0.804 → 0.853, isolated-label ASW 0.701 → 0.747. These rise (look better) but not monotonically
      (60 genes dips), so they are scored as false alarms.
    - Co-expression: silhouette 0.404 → 0.369 and isolated-label ASW 0.701 → 0.684 as Myeloid moves towards
      Fibro.
11. **The annotator checks behave.** `cal.ece` and `nov.frac_novel` are flat in every label and null scenario
    and under coordinate shuffling (✔ 22/22). ECE grows with depth loss (0.0006 → 0.078 at 94 %), batch
    (→ 0.008) and spillover (→ 0.009).
    The novelty share stays 0 under depth loss, as it should with a depth-matched reference, and rises with
    batch (0 → 0.88 at log-sd 1.0) and ambient (0 → 0.20 at 0.2): a section whose expression has shifted
    from the reference reads as novel.

### Remaining ❌ / ⚠️ / ⛔, cell by cell

Values are means over 5 seeds, baseline → maximum strength; "threshold" is the move rule above.

**Blind by construction (the metric cannot see that error; read it beside one that can):**
- `ext.comp_l1` on `label_noise` (0 → 0.090, threshold 0.100) and `type_swap` (0 → 0.067): errors that keep
  the proportions balanced.
- `pmp.no_dict_frac` on `label_noise`, `sibling_confusion`, `type_swap`, `patch_error` (0 throughout): it
  only measures labels missing from the dictionary, and these errors use dictionary labels.
- Marker metrics on `over_split` and `merge`: `pmp` 0.996 → 0.995 then a skip, `nmp` 1.000 → 1.000 then a
  skip, `fid.mean_auc` 0.995 → 0.991 then a skip, `mgo` 1.00, `neighbor_contamination` 0.0003 → 0.0002 then
  a skip. Pair them with `no_dict_frac` (finding 1). `coex.by_type` is blind to them too (0.967 → 0.966,
  skip at `merge` 3) because it scores only labels present in the reference, and it does not report the share
  of cells it dropped.
- cLISI on `sibling_confusion`, `type_swap`, `patch_error`, `merge` (1.000 throughout): a median capped at 1.
- `iv.neighborhood_purity` on `merge` (1.000 throughout) and `iv.avg_similarity` on `over_split`
  (0.804 → 0.804): merging only makes neighbourhoods purer, and random halves have the same within-label
  similarity as the whole.
- `isc.iscl` on `batch`, `merge`, `patch_error`, `sibling_confusion` (1.000 throughout): the reciprocal
  classification still finds every label's best match.
- Spatial label metrics on `patch_error` (coherence, PAS, CHAOS flat, finding 2), CHAOS on `merge`
  (0.0360 → 0.0356) and `sibling_confusion` (0.0360 → 0.0380, threshold 0.0027), PAS on `merge`
  (0.022 → 0.009). A compact or merged error is still spatially coherent.
- `coex.global` and `coex.by_type` on `spillover` (finding 6).

**Real responses below the move rule (the metric moves the right way, too little):**
- `iv.davies_bouldin` on `sibling_confusion` (1.06 → 1.56), `patch_error` (1.06 → 1.36) and `merge`
  (1.06 → 1.81): monotone, but the threshold is 0.80 because 40 % label noise stretches its range to 3.24.
- `isc.iscg` / `isc.integrated` on `sibling_confusion` (0.938 → 0.923), `patch_error` (0.938 → 0.934) and
  ISCG on `over_split` (0.938 → 0.917), threshold 0.030-0.034; ISCG on `type_swap` returns to 0.938 (name
  swap).
- `fid.mean_auc` (0.995 → 0.986) and `fid.mean_cohens_d` (6.14 → 5.79, threshold 0.47) on
  `sibling_confusion`: the mislabelled T_CD8 cells still carry the two T lineage markers in T_CD4's set.
- `onto.rho` on `sibling_confusion` (0.25 → 0.15), `patch_error` (0.25 → 0.15) and `over_split`
  (0.25 → 0.33, the wrong way). Its seed SD is 0.08 over 15 label pairs (threshold 0.23), and its clean
  value is only 0.25: in this world every type has its own private markers, so the centroids are nearly
  equidistant and only the T pair follows the ontology. It cannot detect anything here; it needs a real
  embedding with lineage structure. Under label noise it rises (0.25 → 0.47 at 40 %, 〰️).

**Real but undesired responses:**
- ⛔ ISCG and `isc.integrated` rise under `merge` (0.938 → 0.978) and ⛔ spatial coherence rises
  (0.920 → 0.967): a merged label is more consistent and more coherent (findings 2-3).
- ⚠️ `crisp_purity`, `mecr`, `pmp` under `coexpression` (finding 4).
- ⚠️ coherence, PAS and CHAOS under `infiltration` (finding 2).
- ⚠️ silhouette, `integrated`, Ward-PropMatch, average similarity, Calinski-Harabasz and the isolated-label
  scores under `rare_type`, `panel_growth` and `coexpression` (finding 10).
- ⚠️ `coex.global` (0.978 → 0.974) and `coex.by_type` (0.967 → 0.960) under `rare_type`: whole-section
  correlations depend on the type mix (the reference is balanced), and at 1 % T_CD8 has 24 cells, just above
  `min_cells` = 20, so its per-type correlations are noisy (seed SD 0.004).
- ⚠️ `onto.rho` under `coexpression` (0.25 → −0.07 at 5 %): Myeloid's centroid moves towards Fibro, an
  immune type next to a structural one. That is what the metric reports (the layout disagrees with the
  ontology); here the biology is real.
- ❌ `ns.n_suspicious` on `label_noise`, `merge`, `patch_error` (0 suspicious types throughout); it detects
  only `over_split` (2 / 4 / 6).
- ❌ `mgo.mean_overlap` on `label_noise`, `sibling_confusion`, `patch_error` (1.00 throughout, finding 7).
- ❌ `sig.moran_gap` and `sig.frac_above_ctrl` on `ambient` (finding 9).

## Dedicated checks

- **Conformal sets (α = 0.1):**
  - Calibrated on a labelled 30 % of the section, coverage is 0.90.
  - Calibrated on the reference, it over-covers at the same depth (0.976). At a quarter of the depth it
    under-covers (0.733), and the misses are *empty* sets (26.7 % of known cells), so known cells are called
    "novel".
  - For a 2 % rare type, marginal coverage of T_CD8 is 0.796; Mondrian restores it to 0.899. Clustered
    Mondrian does not (0.796).
  - With B left out of the reference (a truly novel type), every novel cell gets an empty set (1.00), and so
    do 2.2 % of known cells at full depth (25 % at quarter depth). Conformal novelty p-values single the
    novel cells out (AUC 0.999, 97 % flagged, 0.7 % of known cells flagged, FDP 0.05 at a 0.10 target), but
    only against a depth-matched reference: the reference as it is (400 counts per cell vs 150) flags 86-100 %
    of known cells.
  - The reference-trained classifier is well calibrated on the section (ECE 0.001; 0.006 at a quarter of
    the depth).
- **`panel_resolvability`:**
  - Macro-F1 falls monotonically with depth and with markers removed: 1.00 → 0.35.
  - It names T_CD4 as T_CD8's confuser in 61 of the 62 runs that name one.
  - The tiers now use the chance-corrected F1: with *no* private markers at depth 50 (macro-F1 0.72)
    `frac_resolvable` is 0.07 (was 1.0 in 0.6.0). At depth 200 every panel is resolvable.
- **`deconvolution_metrics`:** monotone. R² falls 1.00 → 0.70 at 50 % noise and → −0.20 at pure noise; AUPR
  1.00 → 0.92 → 0.25; SSIM 1.00 → 0.82 → 0.01.
- **ISC regimes:** see finding 3.

## Edge cases ([`results/edge_cases.md`](results/edge_cases.md))

- **No deep crashes.** Every 💥 is a deliberate, named input error:
  - NaN labels: all 14 label-reading adapters raise `ValueError: adata.obs['label'] has 30 missing (NaN/None)
    labels of 600 ...` (0.6.0: 6 `TypeError`s from a sort and 2 other `ValueError`s, while fidelity and overlap
    scored the NaNs as a type).
  - Log-normalised `X` with no `counts` layer: the 12 count-based adapters raise `ValueError` asking for
    `layer=`. In 0.6.0 ISC, fidelity, PMP, NMP, `neighbor_contamination`, confidence and signal QC returned
    different numbers with no warning.
  - Missing coordinates: the four spatial adapters raise squidpy's `KeyError`.
- **Skips instead of made-up values:** a single label (`neighborhood_sanity` crashed in 0.6.0), markers off
  the panel (`crisp_purity` 1.0, `mecr` 0.0, confidence margin 0 / entropy 1 in 0.6.0), and a label vocabulary the
  reference does not share (integer codes, renamed or merged labels: `nmp` was 0) now return NaN with a reason.
- **Still changes values:** a marker set for an absent type whose markers sit at ambient level lowers
  CRISP / MECR / PMP and shifts the confidence margin and entropy.
- **Handled cleanly:** singleton type, zero-count cells, an unused categorical level (fixed since 0.6.0),
  duplicate coordinates, dense input (identical to sparse), unicode labels, and an all-zero marker gene.
- Two adapter-side NaNs are not metric problems: `cal.ece` on integer or renamed labels (the classifier's
  classes are the reference's names), and `ext.spari` without coordinates.

## Limits

- One simple world: stripes, Poisson counts, 6 types and crisp markers (a marker is 1,000× its off-target
  rate). Real panels are noisier, so detection thresholds here are optimistic.
- A metric that fails here fails a case it should handle. One that passes has only cleared the easy
  version. `onto.rho` is not exercised here at all (see above).
- Real-data checks (#7 and the P2 validation harness in `docs/ROADMAP.md`) are still needed.
