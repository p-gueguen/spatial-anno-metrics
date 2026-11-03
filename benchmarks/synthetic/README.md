# Synthetic sweep - which metrics measure what they claim

`sim.py` simulates a section with known truth: six types in tissue stripes (T_CD4 | T_CD8 | B | Myeloid |
Fibro | Epi), 2,400 cells, 150 counts per cell. Each type has private markers, and the two T subtypes
share two lineage markers. There are 60 background genes and 20 negative-control probes, and six samples
as horizontal bands. A matched dissociated reference comes from the same world. `sweep.py` damages that
section along 15 axes at 4-5 strengths × 5 seeds, runs every public metric, and scores each trajectory
against what the metric should do. It also runs 18 edge cases and four dedicated checks.

```bash
pip install -e ".[spatial,dev]"
python benchmarks/synthetic/sweep.py --seeds 5 --workers 32     # ~1 min on 32 cores
```

**Results (5 seeds):**
- [`results/matrix.md`](results/matrix.md): the verdict per metric and scenario.
- [`results/verdicts.csv`](results/verdicts.csv): value at baseline and at maximum strength, Spearman ρ,
  the threshold, and the first strength that moves the metric.
- [`results/values.csv`](results/values.csv): the raw values.
- [`results/edge_cases.md`](results/edge_cases.md): what each metric returns on the edge cases.
- [`results/check_*.csv`](results/): the four dedicated checks.

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
- When nothing is wrong, label-aware and purity metrics should stay flat.

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

## Findings

**Each metric is blind to some kind of error, so read them as a panel, never one alone.**

1. **No single metric detects all six label errors**, external metrics aside.
   - Internal validity (silhouette, `integrated`, ...) catches 5/6. It misses the name swap, which no
     partition metric can see.
   - `pmp` and `nmp` catch 4/6.
   - All marker metrics (`pmp`, `nmp`, fidelity, `neighbor_contamination`) miss `over_split` and `merge`.
     A label with no entry in `marker_sets` is skipped silently, so the score is computed on the remaining
     types and looks clean.
2. **Spatial coherence and PAS see random errors only.**
   - A compact mislabelled patch (50 % of Fibro called Epi) leaves coherence at 0.903 vs 0.904 and PAS at 0.
   - Merging types *raises* coherence (0.904 → 0.961), and PAS stays 0.
   - Correctly labelled infiltrating cells lower coherence (0.904 → 0.724) and raise PAS (0 → 0.083).
   - At the 20 % threshold, a random half-split (each half has about 50 % same-label neighbours) leaves PAS
     at 0.012.
3. **ISC regime calls are right for over-splitting only.**
   - Over-split: the split halves get `merge` in 5/5 seeds (ISCL 0.01). Correct.
   - Merge: merged labels get `consistent` in 5/5, and ISCG *rises* (0.936 → 0.978). It should say
     `split`.
   - 20 % random label noise and a batch effect (log-sd 1.0) both read as `split`.
   - ISCL stays ≥ 0.90 under 40 % label noise; on this fixture it only reacts to over-splitting.
4. **Real co-expression is flagged as contamination** by every reference-free purity metric:
   - `crisp_purity` 0.81 → 0.69 and `mecr` 0.042 → 0.065, already at 5 % co-expression;
   - `pmp` 0.995 → 0.950;
   - `nmp` stays flat (0.001), because the reference carries the co-expression.
5. **`crisp_purity` is 0.81 on perfectly clean data.** The cause is an off-target rate of 0.1 % of
   the marker rate, summed over 20 markers from other lineages, since the rule is "≥ 1 count". That baseline
   depends on the panel and on depth, so absolute CRISP values do not compare across sections.
6. **`neighbor_contamination` moves with spillover** (0 → 0.036 at 40 %) and stays at 0.011 or below
   under co-expression. `nmp` hardly moves under spillover (0.001 → 0.008), and `ncp` only at 40 %
   (1.00 → 0.94).
7. **`marker_gene_overlap` stays at 1.0 at 40 % label noise**, because the top DE genes survive random
   noise. It reacts only to a type swap.
8. **Partition metrics cannot see a name swap.** ARI, ECS, silhouette and ISC return to their clean values
   when Fibro and Epi are fully swapped. Marker metrics, macro-F1 and `hier_acc` catch it.
   `composition_accuracy` is blind to errors that keep the proportions balanced.
9. **Signal QC:**
   - SNR is tied to the control probes: with near-silent controls it reads 30 (the `eps` floor). It falls
     to 5.6 at an ambient of 0.02 per gene.
   - Moran's gap barely moves under ambient noise (0.093 → 0.086).
   - Both behave correctly under coordinate shuffling and depth loss.
10. **Silhouette-type metrics move when composition changes** (a rare type, or a bigger panel), even with
    correct labels. Do not compare them across sections with different composition.

## Dedicated checks

- **Conformal sets (α = 0.1):**
  - Calibrated on a labelled 30 % of the section, coverage is 0.90.
  - Calibrated on the reference, it over-covers at the same depth (0.976). At a quarter of the depth it
    under-covers (0.733), and the misses are *empty* sets, so known cells are called "novel".
  - For a 2 % rare type, marginal coverage is 0.796; Mondrian restores it to 0.899.
- **`panel_resolvability`:**
  - Macro-F1 falls monotonically with depth and with markers removed: 1.00 → 0.35.
  - It names T_CD4 as T_CD8's confuser in every seed.
  - The tiers are lenient: `frac_resolvable` = 1.0 at depth 50 with *no* private markers (macro-F1 0.72).
- **`deconvolution_metrics`:** monotone. R² falls 1.00 → 0.70 at 50 % noise and → −0.20 at pure noise.
- **ISC regimes:** see finding 3.

## Edge cases ([`results/edge_cases.md`](results/edge_cases.md))

- **Crashes:**
  - NaN labels crash 8 of 14 metric groups with a `TypeError` / `ValueError` deep inside. One message at
    the input would do.
  - `neighborhood_sanity` crashes on a single label.
- **Missing coordinates:** the spatial metrics raise a clear squidpy `KeyError`.
- **Silent wrong values on log-normalised `X` with no `counts` layer:** ISC, fidelity, PMP, NMP,
  `neighbor_contamination`, confidence and every signal-QC metric return different numbers with no
  warning.
- **Label vocabulary:** labels that do not match the reference's names (integer codes, renamed labels)
  give `nmp` = 0. This is disclosed only in `types_without_negatives`.
- **An unused categorical level** changes `neighborhood_sanity`'s count of suspicious types.
- **A marker set for an absent type** whose markers sit at ambient level lowers CRISP / MECR / PMP.
- **Handled cleanly:** singleton type, zero-count cells, duplicate coordinates, dense input (identical to
  sparse), unicode labels, and an all-zero marker gene.

## Limits

- One simple world: stripes, Poisson counts, 6 types and crisp markers (a marker is 1,000× its off-target
  rate). Real panels are noisier, so detection thresholds here are optimistic.
- A metric that fails here fails a case it should handle. One that passes has only cleared the easy
  version.
- Real-data checks (#7 and the P2 validation harness in `docs/ROADMAP.md`) are still needed.
