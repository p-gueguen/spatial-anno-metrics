# Cell Annotation Quality Metrics — a consolidated reference

**Status:** authoritative design reference (the metrics catalog) · **Scope:** cell-type annotation quality for scRNA-seq **and** imaging spatial transcriptomics (Xenium / MERSCOPE / CosMx) · **Audience:** anyone evaluating "is this cell-type annotation any good?" — in a package, a custom report, or a benchmark.

> **📦 Reading this inside the `spatial-anno-metrics` package.** This catalog originated in the
> **SpatialScribe** project and is vendored here as the package's reference. Two mappings to keep in mind:
> - **Code paths:** where §6 and the master table cite `src/spatialscribe/analysis/eval_metrics.py` /
>   `signal_qc.py`, the code in *this* package lives at `src/spatial_anno_metrics/eval_metrics.py` /
>   `signal_qc.py` and is imported as `spatial_anno_metrics.eval_metrics` / `.signal_qc` (top-level
>   re-exports: `import spatial_anno_metrics as sam; sam.internal_validity(...)`). As of **v0.3.0**,
>   many metrics the tables attribute to SpatialScribe modules or mark "—" now ship here:
>   - **§3c/§3f purity** (`purity.py`): `crisp_purity`, `mecr`, `pmp`, `nmp`, `ncp`.
>   - **§3e spatial label quality** (`spatial.py`): `spatial_coherence` (+ PAS), `neighborhood_sanity`.
>   - **§3d per-cell confidence** (`confidence.py`): `per_cell_confidence` (marker margin + entropy),
>     `label_stability` (subsampling flip-rate, annotator-agnostic via a callback).
>   - **§3a / §3g / external** (`eval_metrics.py`): Davies-Bouldin + Calinski-Harabasz (inside
>     `internal_validity`), `marker_gene_overlap`, `avg_bio` (AvgBIO), and `conformal_prediction_sets`.
>   - **§3b ISC** (`eval_metrics.py`, v0.5.0): `inter_sample_consistency` returns the Garnica et al. 2026 pair
>     ISCG / ISCL, their product, and a per-label `isc_regime` (consistent / merge / split / remove).
>
>   So the "In code" column citing `purity.py` / `spatial.py` / `annotate.py` (margin/entropy/stability)
>   or "—" for those metrics is superseded by the package. Metrics still attributed to `annotate.py`
>   (the confidence *funnel*), `panel_check.py`, and `methods.py` remain SpatialScribe-pipeline-only.
> - **Companion docs:** relative links to `cell-annotation-qc.md`, `annotation_qc_thresholds.yaml`, and
>   `annotation-method-selection.md` refer to sibling docs in the SpatialScribe repo, not this one.

> **Why this document exists, and how it relates to the others.** We already have [`cell-annotation-qc.md`](cell-annotation-qc.md), which is excellent — but it is a **spatial pipeline QC funnel** (six ordered gates: segmentation → counts → contamination → panel adequacy → confidence → spatial coherence) scoped Xenium-first and oriented toward *per-cell gating and abstention inside the SpatialScribe pipeline*. It does **not** consolidate two whole families of metric: (a) **reference-based / external-validation** metrics used when a ground truth or trusted reference exists (accuracy, F1, ARI, NMI, kappa, biological fidelity), and (b) **internal cluster-validity + inter-sample-consistency** metrics used when there is no ground truth (silhouette variants, neighborhood purity, medoid/Ward concordance, inter-sample consistency — the family the Carmona lab's [**scTypeEval**](https://github.com/carmonalab/scTypeEval) packages). This document is the **catalog**: every metric worth using to *evaluate* an annotation, organized by what evidence you have, what each measures, whether it needs a reference, and where it lives (or should live) in code. Use this to *choose a battery*; use `cell-annotation-qc.md` for the spatial per-cell funnel and [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml) for the tunable numbers; use [`annotation-method-selection.md`](annotation-method-selection.md) to choose *which annotator to run*.

**The one thing to internalize:** in real FGCZ projects there is **almost never a ground truth**. So the **internal / reference-free** metrics (§3) are the default regime, and the **external** metrics (§2) apply only to benchmarks with labels. Do not reach for accuracy/F1 on a delivery — you have no truth to compute them against.

---

## Table of contents

1. [Decision guide — which metrics to compute](#1-decision-guide)
2. [External / reference-based metrics (ground truth or trusted reference)](#2-external)
3. [Internal / reference-free metrics (the default regime)](#3-internal)
   - 3a. [Cluster-validity / separation (scTypeEval + classics)](#3a-cluster-validity)
   - 3b. [Inter-sample consistency — ISC (scTypeEval)](#3b-isc)
   - 3c. [Marker / contamination purity (reference-free)](#3c-purity)
   - 3d. [Per-cell confidence / uncertainty](#3d-confidence)
   - 3e. [Spatial-context metrics (spatial only)](#3e-spatial)
   - 3f. [Reference-anchored specificity (needs a matched profile, not labels)](#3f-anchored)
   - 3g. [Marker-program fidelity (needs marker sets, not truth)](#3g-marker-fidelity)
   - 3h. [Signal-quality QC (annotation-independent; SpatialQM)](#3h-signal-qc)
4. [What we standardize on — default metric batteries](#4-batteries)
5. [Master metric table (quick reference)](#5-master-table)
6. [Implementation status in SpatialScribe (`eval_metrics.py`)](#6-implementation)
7. [scTypeEval — what it adds and how to fold it in](#7-sctypeeval)
8. [Cross-cutting rules & anti-patterns](#8-rules)
9. [References](#9-references)

---

<a name="1-decision-guide"></a>
## 1. Decision guide — which metrics to compute

Pick the column that matches what you actually have. This is the whole document in one table; the rest is detail.

| You have… | Regime | Compute (in priority order) |
|---|---|---|
| **Ground-truth labels** (manual gold standard, or a benchmark dataset) | External (§2) | **macro-F1** + **balanced accuracy** + **ARI** + per-class F1; add **AvgBIO / program concordance** because accuracy ≠ biological fidelity. Never rank methods on plain accuracy alone. |
| **A trusted reference atlas** but no per-cell truth | External-ish (§2.4) + Internal | Label-transfer native confidence (SingleR `delta`/`pruned`, CellTypist prob, RCTD `spot_class`/weight margin) as a quality proxy; **NCP/NMP** (§3f) if the reference is matched scRNA/snRNA. |
| **Biological replicates** (≥5 samples with ≥10 cells per sample × label) but no truth | Internal (§3b) | **Inter-sample consistency on two axes**: ISCG (cross-sample compactness) + ISCL (cross-sample replicability) → a merge / split / remove regime per label; their product ranks annotators or reference atlases without truth (R² > 0.98 with macro-F1). Pair with cluster-validity (§3a). |
| **Just labels + an embedding / expression matrix** (the common case) | Internal (§3a, §3d) | **Silhouette + 2label-silhouette + neighborhood purity** (scTypeEval); **marker-score margin + entropy**; **subsampling label-stability**. |
| **Marker gene sets** (curated per lineage) | Internal (§3c) | **MECR + CRISP purity + PMP + on-target marker specificity** — reference-free contamination readouts. |
| **Spatial coordinates** | Internal (§3e) | **Spatial coherence / PAS**, neighborhood enrichment sanity. Down-weight, never delete (rare infiltrating cells are legitimately incoherent). |
| **A targeted panel** (Xenium/CosMx) | Internal (§3c, §3f) + panel gate | Everything above **plus panel adequacy** (can the panel even resolve the type — [`cell-annotation-qc.md` §7](cell-annotation-qc.md)). |

> Most SpatialScribe deliveries land in rows 4–7 simultaneously (labels + embedding + markers + coords, no truth). The **default spatial battery** in [§4](#4-batteries) is exactly that combination.

---

<a name="2-external"></a>
## 2. External / reference-based metrics (ground truth or trusted reference)

Use these **only** when you have labels you trust as truth — a manual gold standard, a held-out benchmark, or (with caveats) a reference atlas. They answer "how close is the predicted labeling to the reference labeling?"

### 2.1 Per-cell agreement

| Metric | Definition | Range / dir | When to use / gotcha |
|---|---|---|---|
| **Accuracy** | fraction of cells whose predicted label == truth | 0–1, higher better | Dominated by the majority class; **misleading under class imbalance** (a "call everything tumor" classifier scores high). Almost never report alone. |
| **Balanced accuracy** | mean of per-class recall | 0–1, higher better | Corrects majority-class dominance; the right "overall accuracy" number when types are imbalanced (they always are). |
| **Cohen's κ** | agreement corrected for chance | −1…1, higher better | Good single-number agreement that accounts for label prevalence; κ>0.8 strong, 0.6–0.8 moderate. |

### 2.2 Per-class quality (precision / recall / F1)

- **Precision**(t) = TP / (TP+FP) — of cells called t, how many really are t.
- **Recall**(t) = TP / (TP+FN) — of true-t cells, how many were found.
- **F1**(t) = harmonic mean of the two.
- **macro-F1** = unweighted mean of per-class F1 → **weights rare types equally** (the metric that exposes "great on tumor, useless on the 200 T cells"). **weighted-F1** weights by class size (closer to accuracy).

**Rule:** for cell typing, **macro-F1 + per-class F1 table** is the headline, because the failure that matters (dropping a rare immune type) is invisible in accuracy and weighted-F1. This is exactly the RCTD-drops-T-cells failure the SPLIT comparison surfaced.

### 2.3 Partition agreement (label sets as clusterings)

| Metric | Measures | Range | Notes |
|---|---|---|---|
| **ARI** (Adjusted Rand Index) | pairwise co-assignment agreement, chance-corrected | −1…1 (0 = random, 1 = identical) | The standard partition-agreement metric; robust, label-name-independent (doesn't need matched label vocab). |
| **NMI / AMI** (Normalized / Adjusted Mutual Information) | shared information between partitions | 0–1 | NMI is intuitive; **AMI** is the chance-corrected version — prefer AMI when cluster counts differ. |
| **ECS** (Element-Centric Similarity) | cell-level label consistency vs a reference labeling, via an affinity-matrix formulation (Gates & Ahn 2019; α=0.9) | 0–1 | **Resolution-agnostic** — designed to compare labelings with *different numbers of clusters*, where ARI/NMI are unfair. Zhu 2026's headline structure-metric; per-cell, so it also localizes disagreement. Add it whenever prediction/truth granularity differ. |
| **V-measure** (homogeneity + completeness) | homogeneity = each cluster is one type; completeness = each type in one cluster | 0–1 | Decomposes *why* a partition disagrees — over- vs under-splitting. |
| **spARI** (spatially-aware ARI) | ARI that credits spatial arrangement | — | Spatial benchmarks only; needs a ground-truth spatial-domain partition. |

Granularity trap: ARI/NMI between a fine prediction and a coarse truth (or vice-versa) is low **by construction**, not because the annotation is wrong. Match granularity, or evaluate hierarchically.

### 2.4 Reference native-confidence as a quality proxy (no hard truth, trusted reference)

When you ran a reference-based annotator, its own posterior is a *soft* quality signal (not ground truth, but informative):

- **SingleR** `delta` (score gap to next label) and `pruned.labels` (NA = abstained).
- **CellTypist** posterior probability + top1−top2 **margin**; `majority_voting`.
- **RCTD** `spot_class` (`singlet`/`doublet_*`/`reject`) + first-vs-second **weight margin**.
- **Pseudobulk-reference correlation**: Pearson/Spearman of the cell (or cluster) profile vs each reference type; low max-correlation → out-of-distribution / novel.

These also appear in the spatial funnel ([`cell-annotation-qc.md` §8.1](cell-annotation-qc.md)); listed here so the external-regime picture is complete.

### 2.5 Biological fidelity (accuracy ≠ fidelity)

The load-bearing finding of the Zhu et al. 2026 spatial benchmark (this paper — [`sources/2026.06.16.732716v1.full.pdf`](sources/2026.06.16.732716v1.full.pdf); also in [`annotation-method-selection.md`](annotation-method-selection.md)): **high classification accuracy does not guarantee preserved biology** — rankings by fidelity frequently *differ* from accuracy rankings, so accuracy alone "may substantially overestimate" real quality. Complement accuracy with:

- **AvgBIO** (Average Biological Conservation, from scIB) — a summary of structure-preservation metrics, each oriented larger-is-better and scaled to [0,1], then averaged. **Exact form used in Zhu 2026:** `AvgBIO = mean(ARI, NMI, ASW_scaled)` with `ASW_scaled = (ASW+1)/2` (ASW = Average Silhouette Width on the labeled embedding). Reduces over-reliance on any single metric.
- **Marker-gene overlap** — overlap between each method's data-derived per-type DE genes (top 10 / 20 / 30) and the reference cell-type marker genes: do the labels reproduce the *known* markers?
- **ssGSEA program concordance + AUC-ROC on enrichment** — score each cell against ground-truth marker/pathway sets with ssGSEA, then per type compute the **AUC-ROC** of that enrichment score ranking predicted-positive above predicted-negative cells. Directly asks: does a method's label assignment track the cells' actual marker-program enrichment? (Reference-free given marker sets — see [§3g](#3g-marker-fidelity).)
- **Cohen's d (score separability)** — standardized effect size of how cleanly a method's binary type call separates target vs non-target cells on the marker-enrichment score.

**Rule:** in any benchmark, report at least one fidelity metric next to accuracy/F1.

### 2.6 Hierarchy-, composition- & robustness-aware external metrics

The Zhu 2026 framework adds three axes plain accuracy misses:

- **Hierarchical annotation accuracy** (the "accuracy tree") — accuracy scored *along the cell-type lineage hierarchy*, so a call that gets the **major class right but the subtype wrong** earns partial credit instead of a flat miss. This is the principled fix for the granularity trap (§8): a coarse-vs-fine mismatch should not be penalized like a lineage error. Pair with **subtype-level accuracy** (accuracy restricted to fine subtypes) to see *where* resolution breaks down.
- **Cell-type composition / proportion accuracy** — compare predicted vs true **cell-type proportions**. Two settings: (i) **hard labels → global proportions** per section/condition (`composition_accuracy`: L1 / Pearson / JSD on the length-K proportion vectors) — catches methods that get per-cell labels roughly right but distort the population mix (the paper tracks this across a kidney injury–repair time course); (ii) **soft per-spot proportion matrices** (spots × cell-types), i.e. **deconvolution** output — scored with the **OpenProblems `task_spatial_decomposition`** metrics: **R²** (`r2_score(true, pred, multioutput='uniform_average')`, per cell type averaged; ↑) and **JSD** (mean over cell types of `scipy.jensenshannon(true[:,k], pred[:,k])`, axis=0; ↓), plus **RMSE** and per-spot JSD. Use (ii) to score a deconvolver's proportion output (RCTD/TACCO/cell2location weights) against known mixtures / pseudospots (`eval_metrics.deconvolution_metrics`).
- **Robustness to reference volume** — re-run with the reference **down-sampled** (e.g. 25/50/75%) and track accuracy / macro-F1 degradation. A method whose accuracy collapses as the reference shrinks is fragile in the common low-reference regime. Distinct from **query** subsampling ([§3d](#3d-confidence)), which needs no reference.

*(Practical fourth axis in the paper — **scalability**: wall-clock time + peak memory. Not a quality metric, but the reason GraphST OOM'd on 1.37 M cells; track it when method choice must scale.)*

---

<a name="3-internal"></a>
## 3. Internal / reference-free metrics (the default regime)

No ground truth. These test whether a labeling is **self-consistent** — with the embedding, the markers, the neighbors, and across samples. This is what you compute on a real delivery.

<a name="3a-cluster-validity"></a>
### 3a. Cluster-validity / separation (scTypeEval + classics)

Treat each annotated cell type as a cluster and ask: is it **cohesive** (cells of a type are similar) and **separated** (distinct from other types)? Computed in an embedding (PCA/latent) or in a pseudobulk/expression space.

| Metric | Definition | Range / dir | Source |
|---|---|---|---|
| **Silhouette** | per cell, (b−a)/max(a,b) where a = mean within-type distance, b = mean nearest-other-type distance; averaged per type | −1…1, higher better | classic (Rousseeuw 1987); **scTypeEval** |
| **2label-silhouette** | silhouette variant comparing a cell's own type vs **all others pooled** (one-vs-rest) — more stable with many types | −1…1, higher | **scTypeEval** |
| **Neighborhood purity** | fraction of a cell's *k* nearest neighbors (in the embedding) sharing its label | 0–1, higher | **scTypeEval** (`NeighborhoodPurity`); classic kNN purity |
| **Orbital-medoid** | fraction of a type's cells closer to their **own type medoid** than to any other type's medoid | 0–1, higher | **scTypeEval** |
| **Ward-PropMatch** | proportion of a type that falls in its **dominant Ward cluster** (agreement between the labels and an unsupervised Ward clustering) | 0–1, higher | **scTypeEval** |
| **Average-similarity** | within-type similarity relative to between-type similarity | higher | **scTypeEval** |
| **Davies–Bouldin** | mean over types of (within-scatter / between-separation) to the worst-confused neighbor | ≥0, **lower** better | classic; candidate, not in scTypeEval |
| **Calinski–Harabasz** | between/within variance ratio (F-like) | higher | classic; candidate |
| **Graph modularity / connectivity** | modularity of the label partition on the kNN graph; connectivity penalty for split types | higher modularity better | classic; candidate |

**scTypeEval's key framing** (Garnica, Andreatta & Carmona 2026, benchmark over 31 datasets - §3b): a **local** ISC (silhouette in reciprocal-classification-match space, ISCL) catches **over-partitioning**, a **global** ISC (2label-silhouette in pseudobulk-cosine space, ISCG) catches **signal degradation**; the paper integrates them as the **product** (the package README says geometric mean - same ranking, different scale). The important idea for us: *neighborhood purity / silhouette computed on the labels are a direct, reference-free "is this type real?" signal* — under-used in FGCZ spatial reports, and cheap.

Interpreting on sparse panels: silhouette on a raw sparse Xenium matrix is noisy — compute it on the **PCA/latent embedding** (or on pseudobulk profiles for a per-type read), never on raw counts, and remember the `maxRank=150` sparsity caveat when scores feed in.

<a name="3b-isc"></a>
### 3b. Inter-sample consistency — ISC (scTypeEval)

**The distinctive scTypeEval contribution**, now a preprint: Garnica, Andreatta & Carmona, *Evaluating cell type annotations in single-cell omics in the absence of ground truth*, bioRxiv 2026-06-12, DOI [10.64898/2026.06.10.731285](https://doi.org/10.64898/2026.06.10.731285) (scTypeEval is on Bioconductor; benchmark code [carmonalab/ISC_benchmark_reproducibility](https://github.com/carmonalab/ISC_benchmark_reproducibility); data Zenodo [10.5281/zenodo.18921437](https://doi.org/10.5281/zenodo.18921437)). Premise: a label is useful only if it captures a molecular pattern that is **reproducible across biological replicates**. ISC turns replicates into a reference-free quality signal, one score per label.

**How it is computed** (paper Methods; `inter_sample_consistency` follows this at pseudobulk level):
1. One observation per (sample, cell type): pseudobulk = aggregate counts, then log-normalise. A type is kept only with **≥ 5 samples of ≥ 10 cells** (`min_samples = 5`, `min_cells = 10`). 2,000 scran HVGs blocked by sample; 30 PCs when a reduction is used; TCR / Ig / Y-chromosome genes blacklisted.
2. A (sample × type) by (sample × type) **dissimilarity matrix**: pseudobulk Euclidean / cosine / Pearson; single-cell Wasserstein-2; or **reciprocal classification match (RCM)** - for every ordered sample pair, train SingleR (classic mode, pruned labels) on one sample's profiles, predict the other, swap; distance 0 only for a *mutual best hit*, else 1 (a score-based variant keeps the SingleR scores).
3. A **consistency metric** on that matrix with labels as clusters, averaged per label: silhouette, 2-label silhouette (b(i) = mean distance to *all* other labels pooled), neighborhood purity (k = 5), Ward-D2 clustering agreement, nearest-medoid agreement, average similarity.

**The two axes and the four regimes** (paper Fig. 5):

| | high ISCL (replicable) | low ISCL |
|---|---|---|
| **high ISCG** (compact) | **consistent** - keep | **over-partitioned** - compact but not distinguishable across samples (e.g. proliferating vs non-proliferating monocytes) → **merge** |
| **low ISCG** | **heterogeneous / too coarse** - matchable but no strong shared signature (e.g. a pan-B-cell label) → **split** when biology justifies it | **inconsistent** - multiplets, low-quality labels → **remove / collapse** |

- **ISCG** (global) = 2-label silhouette on pseudobulk **cosine** distance - the most sensitive detector of signal degradation.
- **ISCL** (local) = silhouette on the **RCM** distance - top-ranked metric overall, the detector of over-partitioning, and the least affected by batch / condition shifts (relative cross-sample relationships, not absolute distances).
- **Integrated** = ISCG × ISCL (paper; the package README uses the geometric mean).
- **Empirical guide thresholds** ISCG = 0.74, ISCL = 0.45 (median − 1 IQR over a 9-dataset meta-study, Fig. 5c). The reproducibility notebook derives 0.76 / 0.56 with median − 1 SD; use the preprint's numbers, and treat both as *pragmatic guides, not cut-offs* - the authors state no universal cutoff exists.
- **Refinement is testable**: merging low-ISCL labels should raise ISCL *without moving* ISCG; splitting a low-ISCG label should raise ISCG *without lowering* ISCL; unmodified labels must not move (they did not, in the HLCA case below).

**What the benchmark found** (31 datasets from 17 studies, 6 tissues, 7-15 samples each; 6 metrics × 6 dissimilarities = 36 combinations; 8 semi-synthetic tasks; task weights 0.2 signal degradation, 0.2 over-partitioning, 0.1 each for granularity, complexity, n-samples, n-cells, batch, biological perturbation):
- Overall ranking: **silhouette | RCM 0.91** > 2-label silhouette | RCM 0.88 > silhouette | pseudobulk-cosine 0.87. Silhouette-based metrics give the best sensitivity / robustness balance.
- Marginal effects: 2-label silhouette drives task 1 (degradation); silhouette drives task 2 (over-partitioning); RCM is best for task 2 and for batch / condition robustness (tasks 7-8); cosine and Euclidean pseudobulk distances are best for task 1. Hence the two-axis pair above.
- Most metrics are robust to granularity, cellular complexity, number of samples and cells per type. **Exception: neighborhood purity breaks when the number of samples approaches k** (k = 5) - do not use it with few replicates.
- Wasserstein wins some tasks but costs orders of magnitude more time and memory for no overall gain; pseudobulk methods run in < 1 s and < 100 MB per dataset.
- Pearson ≈ cosine and average-similarity ≈ 2-label silhouette (highly correlated; dropped from the figure).
- Method **rankings saturate at ~10 datasets** (median Spearman between disjoint dataset subsets 0.95, flat thereafter) - a metric benchmark on fewer datasets has an unstable leaderboard.
- Classical internal validation at single-cell resolution (cLISI aggregated per type) is **weakly and inconsistently** concordant with ISC: within-dataset separability does not measure cross-sample reproducibility.

**ISC predicts supervised label-transfer performance without ground truth** (14 classifiers incl. a random negative control; 21 intra-dataset sample-wise splits with 30 % of samples as reference / 70 % as query; 12 inter-dataset pairs): mean integrated ISC of the *predicted* labels vs mean macro-F1 gives **R² = 0.982 (intra) and 0.984 (inter)**; either component alone correlates less. ISC of the *training* annotation is associated with test F1 (a consistent reference is a prerequisite) but is not sufficient - some classifiers score low F1 on highly consistent references. Use: rank annotators, or choose a reference atlas, by ISC when no truth exists.

**Meta-study of author labels** (9 atlases, coarse and fine): on average **64 % of labels consistent, 5 % low/low, 6 % low ISCG only, 5 % low ISCL only, 3 % not computable** (< 5 samples or < 10 cells). 5 of 7 coarse and 7 of 8 fine annotation sets contained at least one inconsistent label. **HLCA** (~500 k cells, 486 individuals, 49 datasets, 6-expert consensus labels): club, multiciliated, bronchial- and nasal-goblet cells were low-ISC because nose / airway / parenchyma subpopulations were pooled under one label (ELF3, FXYD3 nose-high; SCGB1A1 airway/lung-high); splitting by tissue raised both ISCs, merging the two goblet labels (high ISCG, low ISCL) raised ISCL selectively.

**Limits the paper states**: needs ≥ 5 biological replicates and enough cells per label; large batch effects lower ISC (RCM least affected); **individual-specific populations such as malignant cells score low regardless of annotation accuracy** - the benchmark's own preprocessing excluded cancer-cell labels, so do not read low ISC on tumour cells as a bad label; interpretation depends on granularity, tissue complexity and cohort composition. The package README adds: in a high-granularity dataset (e.g. only CD8 T subsets) use **ISCL alone**.

**In this package**: `inter_sample_consistency(adata, label_key, sample_key, min_cells=10, min_samples=5, counts_layer=None)` returns `iscg`, `iscl`, `integrated`, per-type dicts, a per-type `regime` from `isc_regime()`, and `dropped_types` (labels below the minimums are reported, not folded in as 1-cell pseudobulks). RCM uses a nearest-profile Spearman classifier in place of SingleR - same reciprocal-best-hit rule, so expect the same ranking and a different absolute scale. ISCL on a targeted spatial panel (a few hundred genes) is not validated by the paper; run the over-partitioning positive control (`tests/test_eval_metrics.py`, task-2 replica) on your own data before trusting a threshold.

This is genuinely additive to `cell-annotation-qc.md` (which has no multi-sample notion) and to FGCZ practice: when a report has replicates, ISC answers "is this cluster a real cell type or a batch/over-clustering artifact?" without any atlas.

<a name="3c-purity"></a>
### 3c. Marker / contamination purity (reference-free)

Do a cell's transcripts belong to its assigned type, or are they a spillover mixture? These need only **marker gene sets** (no reference atlas). Detailed in [`cell-annotation-qc.md` §6](cell-annotation-qc.md); summarized here because they are core annotation-quality metrics, and several are already implemented in [`compare_rctd_tacco.py`](../../scripts/split/compare_rctd_tacco.py).

| Metric | Definition | Dir | Implemented |
|---|---|---|---|
| **MECR** (Mutually Exclusive Co-expression Rate) | over disjoint-lineage marker pairs, `#(both>0)/#(either>0)`, averaged | lower cleaner | `purity.py` (mecr); SPLIT compare |
| **CRISP purity** | per-cell: impure if it detects markers of ≥2 disjoint lineages | higher cleaner | `purity.py` (`crisp_purity`) |
| **PMP** (Positive Marker Purity) | per-cell fraction of counts from the **assigned** type's markers | higher cleaner | `purity.py` (`pmp`) |
| **On-target marker specificity** | on-target fraction of a lineage's marker expression: on/(on+off) mean, bounded [0,1] | higher cleaner | `compare_rctd_tacco.py` (`_specificity`) |
| **On-target signal retention** | fraction of a lineage's own-marker counts kept after a purification step | higher = less over-purified | `compare_rctd_tacco.py` |

These are the metrics the RCTD-vs-TACCO-into-SPLIT comparison uses, and they are the reference-free way to *quantify contamination* — the dominant iST annotation-quality problem.

<a name="3d-confidence"></a>
### 3d. Per-cell confidence / uncertainty

Per-cell "how sure are we about *this* label?" Detailed in [`cell-annotation-qc.md` §8 + §10](cell-annotation-qc.md); the metrics:

- **Marker-score margin** = top1 − top2 marker score. **The margin, not the absolute score, is the confidence.**
- **Entropy** of the class-score distribution (near log K → the classifier is guessing).
- **Reference posterior max-prob** (CellTypist-style) — where a reference exists.
- **Conformal label-set size** — distribution-free: set size 1 = confident, >1 = ambiguous, 0 = novel. The cleanest formalization of "abstain."
- **Annotation stability under subsampling** — drop ~20% of a cell's transcripts, re-annotate, measure label-flip rate; reference-free, evidence-backed (MERFISH: ~20% dropout flips 10–15% of labels). Implemented in `annotate.py` (`annotation_stability`).

<a name="3e-spatial"></a>
### 3e. Spatial-context metrics (spatial only)

Does the label make **spatial sense**? Blind to non-spatial metrics. Detailed in [`cell-annotation-qc.md` §9](cell-annotation-qc.md):

- **Spatial coherence** = mean fraction of a cell's *k* (~15) spatial neighbors sharing its label (`spatial.py` `spatial_coherence`).
- **PAS** (Proportion of Abnormal Spots) = % cells with <20% same-label neighbors — a dataset-level annotation-noise summary.
- **Neighborhood enrichment / co-occurrence** (squidpy) vs a permutation null — biologically impossible adjacencies flag systematic errors.

**Caveat (load-bearing):** rare infiltrating cells are legitimately spatially incoherent. Use these to **down-weight**, never to hard-delete.

<a name="3f-anchored"></a>
### 3f. Reference-anchored specificity (needs a matched profile, not per-cell labels)

Between purely internal and fully external: these need a **matched scRNA/snRNA reference profile** but not per-cell ground truth. From Salas et al. 2025:

- **NCP** (Negative Co-expression Purity) — % of gene pairs **not** co-expressed in the reference that stay non-co-expressed in situ. 0–1, higher = better specificity; commercial platforms mean **>0.8**.
- **NMP** (Negative Marker Purity) — % of reads from "negative markers" (genes a type should NOT express) found in cells assigned to that type. Lower cleaner; the most direct readout of segmentation-driven mis-assignment, and the recommended metric for **comparing segmentations**.

> Terminology trap (kept from `cell-annotation-qc.md`): **NCP** (negative *co-expression* purity, reference-anchored) ≠ the 10x **negative *control-probe* rate** (background binding). Different things; keep separate.

<a name="3g-marker-fidelity"></a>
### 3g. Marker-program fidelity (needs marker sets, not ground truth)

A powerful family the Zhu 2026 benchmark formalizes: given curated marker/program sets, you can measure whether the **labels track the biology without per-cell ground truth** — so these run on real deliveries, not just benchmarks. They are the §2.5 fidelity metrics in their reference-free form:

- **ssGSEA enrichment** (GSVA) — score each cell against each type's marker/pathway set.
- **AUC-ROC on enrichment** — per type, how well the method's predicted-positive cells are ranked above the rest by that type's enrichment score (1 = perfect ranking, 0.5 = chance). A clean, threshold-free "do the cells labeled *T cell* actually carry the highest T-cell-program score?" This **generalizes the on-target specificity** in [`compare_rctd_tacco.py`](../../scripts/split/compare_rctd_tacco.py) into a proper ranking metric.
- **Cohen's d** — standardized target-vs-non-target separation on the enrichment score (effect size, not just AUC).
- **Marker-gene overlap** — overlap of a label's data-derived top-N DE genes with its curated markers.

Use these when you have marker sets but no truth (the usual FGCZ case), as the **fidelity** complement to the **purity** metrics in [§3c](#3c-purity): purity asks "is the cell contaminated?", fidelity asks "does the label match the cell's program?".

<a name="3h-signal-qc"></a>
### 3h. Signal-quality QC (annotation-independent; SpatialQM)

Upstream of annotation quality is *signal* quality: does the section carry real, structured biology above background? These are **label-independent** - they score the data, not the labels - and are ported from the **SpatialQM** R package (Center for Spatial OMICs; NOT the identically-named Plummer/Spatial-Touchstone tool - see §9). Implemented in `analysis/signal_qc.py`, **squidpy-compatible** and fast:

- **Moran's I gene-vs-control** (`getMorans`) - the flagship. Real genes are **spatially autocorrelated** (Moran's I well above 0); negative-control probes are spatially random (I ≈ 0). The **separation** between the gene and control I distributions (and the fraction of genes above the controls' 95th percentile) is a reference-free spatial signal-quality score. Computed via squidpy's `sq.gr.spatial_autocorr` on the `sq.gr.spatial_neighbors` graph - parallel, scales to 1e6 cells.
- **Signal-to-noise** (`getMeanSignalRatio` / `getMaxRatio`) - per-cell log2(mean gene expr / mean negative-control expr) and log2(max gene / control mean); a cell whose signal doesn't clear the control floor is untypeable.
- **Sparsity** (`getSparsity`) - zero fraction of the count matrix (dataset + per cell).
- **Detection entropy** (`getEntropy`) - per-cell normalized Shannon entropy over detected genes; low → one gene dominates (rRNA / probe artifact / extreme specialization) - the continuous version of the "max single-gene fraction" flag in [`cell-annotation-qc.md` §5](cell-annotation-qc.md).
- **Transcript density per area** (`getTxPerArea`) - counts / cell area (needs `cell_area`).

These gate whether annotation should be trusted *at all* on a section (they sit alongside Layer 0-3 of the spatial funnel); run them before the label-quality metrics above. `signal_qc.run_signal_qc(adata)` returns the headline.

---

<a name="4-batteries"></a>
## 4. What we standardize on — default metric batteries

The authoritative recommendation: don't compute everything; compute the right battery for the scenario. Each battery is ordered (compute the first few always; add the rest if cheap/available).

**A. Spatial delivery, no reference (the typical Xenium/CosMx report):**
1. **Panel adequacy** (can the panel resolve the type at all — [`panel_check.py`](../../src/spatialscribe/analysis/panel_check.py))
2. **MECR + CRISP purity + PMP** (reference-free contamination)
3. **Marker-score margin + entropy** (per-cell confidence)
4. **Spatial coherence / PAS** (down-weight only)
5. **Subsampling label-stability**
6. *(if ≥5 replicates)* **Inter-sample consistency (ISCG + ISCL, regime per label)**
7. *(reporting)* section **annotatability headline**: % PASS / WARN / abstain + top reasons

**B. Spatial with a matched snRNA reference:** battery A **plus** RCTD `spot_class` + weight-margin, **NCP/NMP**, label-transfer posterior margin. (This is the highest-quality regime; reference choice matters more than method — see [`annotation-method-selection.md`](annotation-method-selection.md).)

**C. scRNA-seq, no ground truth:** **silhouette + 2label-silhouette + neighborhood purity** (scTypeEval internal) + **ISCG + ISCL across samples** (≥5 replicates) + **marker-score margin** + **AUC-ROC on marker enrichment** (§3g fidelity). Add Davies–Bouldin / Ward-PropMatch to catch over-clustering.

**D. Benchmark with ground-truth labels (the Zhu 2026 template):** **macro-F1 + per-class F1 + balanced accuracy** (classification) + **ARI/AMI + ECS** (partition, ECS when granularity differs) + **AvgBIO** (structure conservation) + **≥1 fidelity metric** (AUC-ROC on marker enrichment / marker-gene overlap) + **hierarchical accuracy** (so subtype-vs-lineage errors are scored differently) + **composition/proportion accuracy** if the question is compositional/dynamic; report **scalability** (time/mem) when method choice must scale. Never rank on accuracy alone; never on a single metric.

**Three cross-cutting rules (apply to every battery):**
1. **Margin > absolute score** — a high score not separated from the runner-up is ambiguous, not confident.
2. **Keep ≥1 method-independent metric** — if you validate an RNA annotation with a reference-derived silver label, that metric favors the reference arm by construction; MECR/CRISP (marker-only) is the independent check.
3. **Calibrate/normalize locally** — thresholds don't transfer across tissue/panel/platform (10x's own guidance); prefer distributional/relative variants; min-max scaling across a few methods inflates tiny gaps (annotate raw values).

---

<a name="5-master-table"></a>
## 5. Master metric table (quick reference)

Regime: **E** = external (needs truth) · **I** = internal (reference-free) · **A** = reference-anchored (needs a matched profile, not labels). Needs: L=labels, R=reference atlas, M=marker sets, S=biological replicates (samples), XY=spatial coords, Z=molecule/z data.

| Metric | Regime | Needs | Measures | Dir | In code |
|---|---|---|---|---|---|
| Accuracy / balanced accuracy | E | L+truth | per-cell correctness | high | `eval_metrics.py` |
| Precision / recall / **macro-F1** (+ micro=acc, weighted) | E | L+truth | per-class quality; rare-type sensitivity | high | `eval_metrics.py` |
| **Hierarchical accuracy** · subtype-level accuracy | E | L+truth+hierarchy | partial credit along the lineage tree | high | `eval_metrics.py` |
| **ARI** / AMI / NMI / V-measure | E | L+truth | partition agreement | high | `eval_metrics.py` |
| **ECS** (Element-Centric Similarity) | E | L+truth | resolution-agnostic, cell-level partition agreement | high | `eval_metrics.py` |
| **AvgBIO** = mean(ARI,NMI,ASW_scaled) · **ASW** | E | L+truth+embed | biological-structure conservation | high | — |
| **Composition / proportion accuracy** | E | L+truth | predicted vs true cell-type proportions (hard-label global) | high | `eval_metrics.py` |
| **Deconvolution props** (R²/JSD/RMSE, OpenProblems) | E | proportion matrices | soft per-spot proportions vs known mixtures | R²↑/JSD↓ | `eval_metrics.py` |
| **Robustness to reference volume** | E | L+truth+R | accuracy decay as reference shrinks | stable | — |
| Scalability (runtime · peak memory) | — | — | practical cost, not quality | lower | bench scripts |
| SingleR delta · CellTypist prob · RCTD spot_class | E/I | R | native confidence proxy | high | `methods.py` join |
| **Silhouette / 2label-silhouette** | I | L+embedding | cohesion vs separation | high | `eval_metrics.py` |
| **Neighborhood purity** | I | L+embedding | kNN label agreement | high | `eval_metrics.py` |
| Orbital-medoid · Ward-PropMatch · Avg-similarity | I | L+embedding | cluster–label concordance | high | `eval_metrics.py` |
| Davies–Bouldin · Calinski–Harabasz | I | L+embedding | separation | DB low / CH high | — |
| **ISCG** (2-label silhouette \| pseudobulk cosine) | I | L+S (≥5 samples) | cross-replicate compactness (signal degradation) | high; guide 0.74 | `eval_metrics.py` |
| **ISCL** (silhouette \| reciprocal-classification match) | I | L+S (≥5 samples) | cross-replicate replicability (over-partitioning) | high; guide 0.45 | `eval_metrics.py` |
| **MECR** | I | M | cross-lineage co-expression | low | `purity.py`, SPLIT compare |
| **CRISP purity** | I | M | per-cell lineage purity | high | `purity.py` |
| **PMP** | I | M | own-marker transcript fraction | high | `purity.py` |
| On-target specificity / retention | I | M | marker specificity after purification | high | `compare_rctd_tacco.py` |
| **AUC-ROC on marker enrichment** · Cohen's d | I | L+M | do labels track the marker-program enrichment | high | `eval_metrics.py` |
| **Marker-gene overlap** (top-N DE vs markers) | I | L+M | labels reproduce known markers | high | — |
| Marker-score **margin** · entropy | I | L+M | per-cell ambiguity | margin high / entropy low | `annotate.py` |
| Conformal set size | I | L+scores | abstention (ambiguous/novel) | =1 | — |
| Subsampling label-stability | I | counts | detection-robustness | low flip | `annotate.py` |
| Spatial coherence / **PAS** | I | L+XY | spatial label agreement | high / PAS low | `spatial.py` |
| Neighborhood enrichment sanity | I | L+XY | impossible adjacencies | — | `spatial.py` |
| **Moran's I gene-vs-control** | I | XY+control | spatial signal vs background noise | high separation | `signal_qc.py` |
| Signal-to-noise · sparsity · detection entropy · tx/area | I | counts (±control/area) | label-independent signal quality | mixed | `signal_qc.py` |
| **NCP** / **NMP** | A | R+M | reference-anchored specificity | NCP high / NMP low | — |
| Panel adequacy (coverage/confusable) | I | M+panel | can the panel resolve the type | green | `panel_check.py` |
| ovrlpy VSI / 3D-coherence | I | Z | vertical doublets/overlap | high | `qc.py` join (H5) |

---

<a name="6-implementation"></a>
## 6. Implementation status in SpatialScribe

**Contamination + confidence + spatial** (reference-free): `purity.py` (MECR, CRISP, PMP, NMP/NCP guards), `annotate.py` (marker margin, entropy path, subsampling stability, penalized-posterior confidence + abstention), `panel_check.py` (panel adequacy + identifiability AUC), `spatial.py` (spatial coherence, PAS, neighborhood sanity), and `scripts/split/compare_rctd_tacco.py` (on-target specificity, marker purity, on-target retention, MECR panel).

**Metric battery** — `analysis/eval_metrics.py` (**built, TDD**: `tests/test_eval_metrics.py`). Covers exactly the families this doc exposed:

```python
# reference-free internal validity (scTypeEval family) on an embedding + labels
internal_validity(adata, label_key="cell_type", embedding="X_pca", k=30) -> dict
    # {silhouette, silhouette_2label, neighborhood_purity, orbital_medoid,
    #  ward_propmatch, avg_similarity, integrated (geometric mean)}
inter_sample_consistency(adata, label_key, sample_key, min_cells=10, min_samples=5) -> dict  # ISCG / ISCL / integrated / regime
marker_program_fidelity(adata, label_key, marker_sets) -> dict     # §3g: per-type AUC-ROC + Cohen's d
external_scores(pred, truth) -> dict     # balanced_acc, macro/weighted/per-class F1, ARI/AMI/NMI, kappa, ECS
element_centric_similarity(a, b, alpha=0.9) -> float               # exact vectorized ECS (Gates & Ahn)
hierarchical_accuracy(pred, truth, hierarchy, partial=0.5) -> dict # partial credit along the lineage tree
composition_accuracy(pred, truth) -> dict                         # hard-label global proportions: L1 / Pearson / JSD
deconvolution_metrics(true_prop, pred_prop) -> dict               # soft per-spot proportion matrices: R2 / JSD / RMSE (OpenProblems)
annotation_quality(adata, label_key, marker_sets=None, sample_key=None) -> dict   # reference-free battery
```

`internal_validity` runs the O(N²) silhouettes on a `subsample` for scalability and falls back to a quick PCA when no embedding is present; `external_scores` is a thin `sklearn.metrics` wrapper plus the exact ECS closed form. The reference-free battery (`annotation_quality`) is **wired into `qc.run_funnel`** (guarded, non-breaking) so every labeled section's funnel headline carries `annotation_quality.internal_validity` + `marker_fidelity`; `external_scores` / `hierarchical_accuracy` / `composition_accuracy` / `deconvolution_metrics` are for benchmark scripts with ground truth. Tested against well-separated vs shuffled fixtures (clean labeling → purity≈1, high integrated; shuffled → near-chance).

**Signal-quality QC** — `analysis/signal_qc.py` (**built, TDD**: `tests/test_signal_qc.py`; §3h) ports the label-independent SpatialQM metrics, squidpy-compatible: `moran_signal` (Moran's I gene-vs-control via `sq.gr.spatial_autocorr`), `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, and `run_signal_qc` (headline). Standalone (run before annotation to gate whether a section is trustworthy at all); not auto-wired into the funnel because Moran's I over the full panel is heavier than the per-cell metrics.

> **Two scTypeEval implementations — pick the right one.** `eval_metrics.internal_validity` is an **embedding-space approximation** built for a fast, always-on `qc.run_funnel` headline: it computes the family on the PCA embedding (subsampled) and does **not** implement scTypeEval's dissimilarity spaces (the pseudobulk-level ISCG / ISCL pair, with a Spearman-classifier `recip_classif:Match`, lives in `inter_sample_consistency` since v0.5.0), and its `orbital_medoid` / `avg_similarity` / `ward_propmatch` use classic embedding-space definitions rather than the exact `IntClusVal.R` ones. A **byte-exact port** (dissimilarity spaces + `recip_classif:Match` + exact IntClusVal formulas + per-sample pseudobulk, TDD-verified) lives in [`scripts/norm_benchmark/sctypeeval_metrics.py`](../../scripts/norm_benchmark/sctypeeval_metrics.py) — use it when you need scTypeEval's recommended local(recip-match)/global(pseudobulk-cosine) combo or paper parity; use `eval_metrics.py` for the production funnel. Empirically the two **rank annotations/normalizations identically but differ in absolute scale**, so the approximation is fine for a QC headline. See [`scripts/norm_benchmark/README.md`](../../scripts/norm_benchmark/README.md) — the **normalization × annotation-quality benchmark** on the 100k breast Xenium 5K section, which anchors to this catalog and reports both implementations plus the implemented contamination/confidence/spatial families across 9 normalizations. Headline finding: normalization is **~invariant** for the contamination metrics (§3c, they read raw counts) but has a **strong** effect on internal separability (§3a/b), the abstention rate (§3d), and above all the clustering the annotation is built on.

---

<a name="7-sctypeeval"></a>
## 7. scTypeEval — what it adds and how to fold it in

[**scTypeEval**](https://github.com/carmonalab/scTypeEval) (Carmona lab; Bioconductor) is the R package behind the ISC preprint (§3b; Garnica, Andreatta & Carmona 2026). Inputs: count matrix + metadata / Seurat / SCE; feature selection by HVGs (scran), marker genes, or a custom list; dissimilarities `Pseudobulk:{Euclidean,Cosine,Pearson}`, `WasserStein`, `recip_classif:{Match,Score}`; consistency metrics `silhouette`, `2label_silhouette`, `NeighborhoodPurity`, `ward_PropMatch`, `Orbital_medoid`, `Average_similarity`; optional PCA.

**Recommended combo (benchmarked, §3b):** ISCL = `silhouette | recip_classif:Match`, ISCG = `2label_silhouette | Pseudobulk:Cosine`, integrated = product (README: geometric mean). Defaults: `Log1p`, 2,000 HVGs, 30 PCs, `min_samples = 5`, `min_cells = 10`, `knn_graph_k = 5`, `hclust_method = "ward.D2"`.

**How it is folded in:** `inter_sample_consistency` (§3b) carries a pseudobulk-level port of both axes plus the regime call; `internal_validity` carries the metric family on a cell embedding (§6 note). Call the R package as a subprocess only when byte-parity with the paper (SingleR-based RCM, Wasserstein) is required.

**Benchmark-design lessons worth reusing** (from `ISC_benchmark_reproducibility`; they apply to any metric benchmark, including ours):
- Score each task with a bounded 0-1 curve-shape statistic, not a raw metric value: **monotonicity** = fraction of consecutive perturbation steps that move the expected way; **constant-fit** = 1 − RMSD / 0.5 (0.5 is the largest RMSD a [0,1] metric can have, so the score spans the full range); **consistency drop** = ISC(complete) − ISC(split), clipped to [0,1]. Then min-max across methods and weight tasks explicitly (0.2 / 0.2 / 0.1 × 6 here).
- Perturb only the labels or the sampling, never the features: the HVG list is computed once on the unperturbed data and held fixed across perturbation levels, and the rate-1 baseline is computed once and reused by every task.
- Three replicates per perturbation level; the over-partitioning and cell-downsampling targets are the **most abundant** label, so the effect is not confounded with cell number.
- Test ranking stability before trusting a leaderboard: draw disjoint dataset subsets of size k, rank methods within each, Spearman-correlate the two rankings, then fit the IQR of that correlation against k (GAM); here it plateaus at k ≈ 10 of 31.
- Decompose the leaderboard into marginal (residual) effects of each component (metric vs dissimilarity) - that is how the two-axis result was found.
- Reading their outputs: the two-level tasks have **opposite polarity** in code - task 2 scores the raw drop (larger = more sensitive) while tasks 7-8 score 1 − drop (larger = more robust); and the granularity / complexity task keys are labelled inconsistently between the YAML descriptions, the saturation config and the figure notebook (the implementation is: label-subset removal = complexity, hierarchical merging = granularity). Both carry weight 0.1, so the overall ranking is unaffected.

---

<a name="8-rules"></a>
## 8. Cross-cutting rules & anti-patterns

- **No ground truth is the default.** Reach for §3, not §2, on real deliveries. Report external metrics only on benchmarks.
- **No single metric is an oracle.** Every real evaluation is a *battery* (§4). Accuracy alone hides rare-type failure; silhouette alone can be gamed by over-merging; MECR alone can be gamed by deleting transcripts (pair it with on-target retention — the SPLIT comparison does exactly this).
- **Accuracy ≠ biological fidelity** (Zhu 2026): always pair accuracy/F1 with a fidelity metric in benchmarks.
- **Margin > absolute score** for confidence.
- **Keep one method-independent metric** in every battery (marker-only MECR/CRISP) to avoid circularity with reference-anchored ones.
- **Spatial coherence & similar metrics penalize rare biology** — down-weight, never delete.
- **Thresholds are local** — tune per tissue/panel/platform; prefer distributional variants; annotate raw values (don't min-max a 3-method bar chart).
- **Granularity must match** for partition metrics (ARI/NMI) — fine-vs-coarse is low by construction, not by error. Fixes: use **ECS** (resolution-agnostic) for the partition score and **hierarchical accuracy** (§2.6) so a right-lineage/wrong-subtype call earns partial credit instead of a flat miss.
- **Internal metrics are consistency, not correctness** — a confidently-wrong-but-self-consistent annotation (e.g. a whole cluster mislabeled) scores well on silhouette/ISC. Internal metrics catch *incoherence*, not systematic bias; only a reference or orthogonal modality (protein) catches the latter.
- **ISC needs ≥ 5 replicates and ≥ 10 cells per (sample, label)** (Garnica 2026). Below that, report a label as *not computable*, never as 0 or as "inconsistent": a 1-cell pseudobulk is noise, not a replicate.
- **Read ISC on two axes, not one number.** Low ISCL alone → merge, low ISCG alone → split, both low → remove; the product is for ranking annotators or references, the quadrant is for fixing labels.
- **Individual-specific populations (malignant cells, patient-private clones) have low ISC by construction** - exclude them from an ISC-based quality verdict, as the benchmark's own preprocessing excluded cancer-cell labels.

---

<a name="9-references"></a>
## 9. References

**This catalog builds on the FGCZ research docs** — see them for the spatial-funnel detail and verified sources:
- [`cell-annotation-qc.md`](cell-annotation-qc.md) — the six-layer spatial QC funnel; source for Layers 0–6, contamination metrics (MECR/CRISP/PMP/NCP/NMP), confidence/abstention, spatial coherence, and all verified primary sources (Salas *Nat Methods* 2025; Plummer/Spatial-Touchstone *Nat Biotechnol* 2025; Mitchel *Nat Genet* 2025; SPLIT *Nat Methods* 2026).
- [`annotation-method-selection.md`](annotation-method-selection.md) — which annotator to run; source for the accuracy≠fidelity finding and the 20-method spatial benchmark (Zhu et al., bioRxiv 2026; AvgBIO / program concordance).
- [`annotation_qc_thresholds.yaml`](annotation_qc_thresholds.yaml) — the tunable pass/warn/fail numbers consumed by the code.

**scTypeEval (the metrics folded in here):**
- Carmona lab, `scTypeEval` — reference-free annotation-quality assessment (silhouette / 2label-silhouette / neighborhood purity / Orbital-medoid / Ward-PropMatch / average-similarity; inter-sample consistency; local/global/integrated framing). GitHub: [github.com/carmonalab/scTypeEval](https://github.com/carmonalab/scTypeEval) `[tool; benchmarked in the preprint below]`. Lab: [github.com/carmonalab](https://github.com/carmonalab).
- **Garnica J, Andreatta M, Carmona SJ.** "Evaluating cell type annotations in single-cell omics in the absence of ground truth," bioRxiv 2026-06-12, DOI [10.64898/2026.06.10.731285](https://doi.org/10.64898/2026.06.10.731285) `[preprint, not yet peer-reviewed]`. Source for ISCG / ISCL / RCM, the 36-combination × 31-dataset benchmark, the 0.74 / 0.45 guide thresholds, the four refinement regimes, and the ISC ↔ F1 R² > 0.98 result (§3b), plus the benchmark-design lessons in §7. Code: [carmonalab/ISC_benchmark_reproducibility](https://github.com/carmonalab/ISC_benchmark_reproducibility) (`targets` pipelines, `renv`-pinned); data: Zenodo [10.5281/zenodo.18921437](https://doi.org/10.5281/zenodo.18921437).

**Classic metric definitions (external + internal validity):**
- **ARI** — Hubert & Arabie, *J. Classification* 1985. **AMI/NMI** — Vinh, Epps & Bailey, *JMLR* 2010.
- **Silhouette** — Rousseeuw, *J. Comput. Appl. Math.* 1987. **Davies–Bouldin** 1979; **Calinski–Harabasz** 1974.
- **F1 / precision / recall / balanced accuracy / Cohen's κ** — standard classification metrics (`sklearn.metrics`).
- **spARI** — spatially-aware ARI (spatial-domain benchmarks).
- **Deconvolution proportion metrics** — **OpenProblems `task_spatial_decomposition`** ([github.com/openproblems-bio/task_spatial_decomposition](https://github.com/openproblems-bio/task_spatial_decomposition)): `r2` = `r2_score(uniform_average)` over cell types; `jsd` = mean per-cell-type `scipy.spatial.distance.jensenshannon` (axis=0). Ported into `eval_metrics.deconvolution_metrics` (+ RMSE, per-spot JSD).

**Benchmark-derived metrics — Zhu et al. 2026 (the paper read for this update):**
- Zhu, Hu, … Meltzer & Zhou, "Benchmarking cell type annotation in spatial transcriptomics: resolving cellular hierarchies, biological fidelity, and dynamic cell states," bioRxiv 2026-06-16, DOI [10.64898/2026.06.16.732716](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full) · local [`sources/2026.06.16.732716v1.full.pdf`](sources/2026.06.16.732716v1.full.pdf). Source for the **ECS**, **AvgBIO = mean(ARI,NMI,ASW_scaled)**, **hierarchical/subtype accuracy**, **composition/proportion accuracy**, **reference-volume robustness**, **AUC-ROC-on-enrichment / Cohen's d / marker-gene-overlap** fidelity metrics folded in here, and the accuracy≠fidelity finding.
- **ECS** — Gates & Ahn, "Element-centric clustering comparison," *Sci. Rep.* 2019 (resolution-agnostic, affinity-matrix formulation).
- **AvgBIO / ASW-based conservation** — Luecken et al., "Benchmarking atlas-level data integration" (scIB), *Nat. Methods* 2022.
- **ssGSEA** — Barbie et al., *Nature* 2009; **GSVA** — Hänzelmann et al., *BMC Bioinformatics* 2013. **Cohen's d** — Cohen 1988.

**Signal-quality QC (SpatialQM, §3h — the metrics ported into `signal_qc.py`):**
- **SpatialQM** — Center for Spatial OMICs, [github.com/Center-for-Spatial-OMICs/SpatialQM](https://github.com/Center-for-Spatial-OMICs/SpatialQM). R package; source for Moran's I gene-vs-control (`getMorans`), signal-to-noise (`getMeanSignalRatio` / `getMaxRatio`), sparsity, detection entropy, tx-per-area, and its own `getMECR`. `[tool]`
- ⚠ **Name collision:** this is **NOT** the *other* "SpatialQM" (Plummer et al., Spatial Touchstone, *Nat Biotechnol* 2025, GEO GSE277080) cited in `cell-annotation-qc.md` §17 for tissue-specific baselines. Two different tools, same acronym — keep them distinct.
- Moran's I — Moran 1950; computed here via squidpy `sq.gr.spatial_autocorr` (Palla et al., *Nat Methods* 2022).

**Spatial / contamination metric sources** (full detail + DOIs in [`cell-annotation-qc.md` §17](cell-annotation-qc.md)): Salas et al. *Nat Methods* 2025 (NCP, NMP, 3D-coherence); Hartman & Satija 2024 (MECR); Center for Spatial OMICs (CRISP); Cable et al. *Nat Biotechnol* 2022 (RCTD); Marconato et al. *Nat Methods* 2026 (SPLIT); STEAM (spatial coherence / PAS).

> **Provenance note.** The spatial/contamination/confidence metrics and their evidence are carried over (with cross-links, not duplicated) from `cell-annotation-qc.md`, which was built from an adversarial-verification research pass. The **new** material here — the external-validation family (§2), the internal cluster-validity + ISC family (§3a–3b), the marker-program fidelity family (§3g), the decision guide (§1), the batteries (§4), the scTypeEval integration (§7), and the ECS / AvgBIO-form / hierarchical / composition / reference-volume-robustness metrics (§2.5–2.6) — is synthesized from standard metric definitions, the scTypeEval repository, and a **direct reading of the Zhu et al. 2026 benchmark Methods** (`sources/2026.06.16.732716v1.full.pdf`: exact ECS/ASW/AvgBIO/AUC-ROC/Cohen's-d formulas transcribed from its Methods section). Classic-metric citation years are working knowledge and should be link-verified before any external publication.
