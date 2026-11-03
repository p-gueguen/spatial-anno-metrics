# Cell Annotation Quality Metrics — a consolidated reference

**Status:** authoritative design reference (the metrics catalog) · **Scope:** cell-type annotation quality for scRNA-seq **and** imaging spatial transcriptomics (Xenium / MERSCOPE / CosMx) · **Audience:** anyone evaluating "is this cell-type annotation any good?" — in a package, a custom report, or a benchmark.

> **📦 Reading this inside the `spatial-anno-metrics` package.** This catalog originated in the
> [**SpatialScribe**](https://github.com/p-gueguen/spatial-scribe) project and is vendored here as the
> package's reference. The §5 master table's "In code" column names the module of *this* package
> (`src/spatial_anno_metrics/`, imported as `import spatial_anno_metrics as sam`); rows marked
> "SpatialScribe only" live in the pipeline, not here. §6 describes the SpatialScribe integration
> (`qc.run_funnel`) and keeps its original module paths. Links to `cell-annotation-qc.md`,
> `annotation_qc_thresholds.yaml` and `annotation-method-selection.md` point to the public
> SpatialScribe repo.

> **Why this document exists, and how it relates to the others.** We already have [`cell-annotation-qc.md`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md), which is excellent — but it is a **spatial pipeline QC funnel** (six ordered gates: segmentation → counts → contamination → panel adequacy → confidence → spatial coherence) scoped Xenium-first and oriented toward *per-cell gating and abstention inside the SpatialScribe pipeline*. It does **not** consolidate two whole families of metric: (a) **reference-based / external-validation** metrics used when a ground truth or trusted reference exists (accuracy, F1, ARI, NMI, kappa, biological fidelity), and (b) **internal cluster-validity + inter-sample-consistency** metrics used when there is no ground truth (silhouette variants, neighborhood purity, medoid/Ward concordance, inter-sample consistency — the family the Carmona lab's [**scTypeEval**](https://github.com/carmonalab/scTypeEval) packages). This document is the **catalog**: every metric worth using to *evaluate* an annotation, organized by what evidence you have, what each measures, whether it needs a reference, and where it lives (or should live) in code. Use this to *choose a battery*; use `cell-annotation-qc.md` for the spatial per-cell funnel and [`annotation_qc_thresholds.yaml`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/annotation_qc_thresholds.yaml) for the tunable numbers; use [`annotation-method-selection.md`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/annotation-method-selection.md) to choose *which annotator to run*.

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
   - 3i. [Ontology label-relation consistency (needs a CL id map, not truth)](#3i-ontology)
4. [What we standardize on — default metric batteries](#4-batteries)
5. [Master metric table (quick reference)](#5-master-table)
6. [Implementation status (package modules)](#6-implementation)
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
| **A targeted panel** (Xenium/CosMx) | Internal (§3c, §3f) + panel gate | Everything above **plus panel adequacy** (can the panel even resolve the type — [`cell-annotation-qc.md` §7](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md)). |

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
| **spARI** (spatially aware ARI) | Rand index / ARI where a pair the two labelings disagree on scores f(d) = α·exp(−d²) (together in prediction, apart in truth) or h(d) = α·(1−exp(−d²)) (apart in prediction, together in truth), coordinates min-max scaled; chance-adjusted with the authors' expectation (Yan, Feng & Luo, Biometrics 2025, doi 10.1093/biomtc/ujaf127; α = 0.8) | ≤1 | A spatially near mislabel costs less than a far one; α = 0 is ARI exactly. O(n²) pairs, so `spatial_ari` scores a seeded 5000-cell subset by default. Asymmetric in (truth, prediction). |

Granularity trap: ARI/NMI between a fine prediction and a coarse truth (or vice-versa) is low **by construction**, not because the annotation is wrong. Match granularity, or evaluate hierarchically.

### 2.4 Reference native-confidence as a quality proxy (no hard truth, trusted reference)

When you ran a reference-based annotator, its own posterior is a *soft* quality signal (not ground truth, but informative):

- **SingleR** `delta` (score gap to next label) and `pruned.labels` (NA = abstained).
- **CellTypist** posterior probability + top1−top2 **margin**; `majority_voting`.
- **RCTD** `spot_class` (`singlet`/`doublet_*`/`reject`) + first-vs-second **weight margin**.
- **Pseudobulk-reference correlation**: Pearson/Spearman of the cell (or cluster) profile vs each reference type; low max-correlation → out-of-distribution / novel.

These also appear in the spatial funnel ([`cell-annotation-qc.md` §8.1](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md)); listed here so the external-regime picture is complete.

### 2.5 Biological fidelity (accuracy ≠ fidelity)

The load-bearing finding of the Zhu et al. 2026 spatial benchmark (this paper — [bioRxiv full text](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full); also in [`annotation-method-selection.md`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/annotation-method-selection.md)): **high classification accuracy does not guarantee preserved biology** — rankings by fidelity frequently *differ* from accuracy rankings, so accuracy alone "may substantially overestimate" real quality. Complement accuracy with:

- **AvgBIO** (Average Biological Conservation, from scIB) — a summary of structure-preservation metrics, each oriented larger-is-better and scaled to [0,1], then averaged. **Exact form used in Zhu 2026:** `AvgBIO = mean(ARI, NMI, ASW_scaled)` with `ASW_scaled = (ASW+1)/2` (ASW = Average Silhouette Width on the labeled embedding). Reduces over-reliance on any single metric.
- **Marker-gene overlap** — overlap between each method's data-derived per-type DE genes (top 10 / 20 / 30) and the reference cell-type marker genes: do the labels reproduce the *known* markers?
- **ssGSEA program concordance + AUC-ROC on enrichment** — score each cell against ground-truth marker/pathway sets with ssGSEA, then per type compute the **AUC-ROC** of that enrichment score ranking predicted-positive above predicted-negative cells. Directly asks: does a method's label assignment track the cells' actual marker-program enrichment? (Reference-free given marker sets — see [§3g](#3g-marker-fidelity).)
- **Cohen's d (score separability)** — standardized effect size of how cleanly a method's binary type call separates target vs non-target cells on the marker-enrichment score.

**Rule:** in any benchmark, report at least one fidelity metric next to accuracy/F1.

### 2.6 Hierarchy-, composition- & robustness-aware external metrics

The Zhu 2026 framework adds three axes plain accuracy misses:

- **Hierarchical annotation accuracy** (the "accuracy tree") — accuracy scored *along the cell-type lineage hierarchy*, so a call that gets the **major class right but the subtype wrong** earns partial credit instead of a flat miss. Credit is by distance on the label tree: the hierarchical F-measure of Kiritchenko et al. 2006 from the lowest-common-ancestor depth. An honest coarse call (the parent) earns 2/3 and a wrong sibling subtype 1/2, so the score prefers 'T cell' over a confident wrong 'CD8 T' for a CD4 T cell. A flat 0.5 for both cannot tell them apart. Report **accuracy on coarse (non-leaf) calls** next to it: a method that hedges to the parent is right more often at lower resolution. This is the principled fix for the granularity trap (§8): a coarse-vs-fine mismatch should not be penalized like a lineage error. Pair with **subtype-level accuracy** (accuracy restricted to fine subtypes) to see *where* resolution breaks down.
- **Cell-type composition / proportion accuracy** — compare predicted vs true **cell-type proportions**. Two settings: (i) **hard labels → global proportions** per section/condition (`composition_accuracy`: L1 / Pearson / JSD on the length-K proportion vectors) — catches methods that get per-cell labels roughly right but distort the population mix (the paper tracks this across a kidney injury–repair time course); (ii) **soft per-spot proportion matrices** (spots × cell-types), i.e. **deconvolution** output — scored with the **OpenProblems `task_spatial_decomposition`** metrics: **R²** (`r2_score(true, pred, multioutput='uniform_average')`, per cell type averaged; ↑) and **JSD** (mean over cell types of `scipy.jensenshannon(true[:,k], pred[:,k])`, axis=0; ↓), plus **RMSE** and per-spot JSD. Use (ii) to score a deconvolver's proportion output (RCTD/TACCO/cell2location weights) against known mixtures / pseudospots (`external.deconvolution_metrics`). Presence/absence **AUPR** (Spotless, Sang-aram et al. eLife 2024: predicted proportions as scores for 'type present in spot') and per-type **SSIM** after max-scaling (Li et al. Nat Methods 2022) are also returned. To compare several deconvolvers on one dataset, `deconvolution_rank_score` averages per-metric ranks (Spotless: RMSE, JSD, AUPR; Li et al.'s accuracy score: PCC, SSIM, RMSE, JSD).
- **Decomposition purity** (`decomposition.decomposition_purity`): 1 - the largest normalised deconvolution weight per cell (RCTD full mode). Higher = more of the cell's expression assigned outside its dominant type. Spatial Touchstone used this as a segmentation-error readout. It depends on the reference; the code follows SegTraQ `rctd_contamination` (MIT).
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
| **Isolated-label ASW / F1** | for labels seen in the fewest batches (every label without batches): the label's mean silhouette (rescaled to 0-1), and its best F1 against any cluster of a clustering sweep. Macro over labels, so a small label counts fully | 0–1, higher | **scib** (Luecken et al. 2022); `isolated_label_scores` (F1 sweeps Ward cuts, not Leiden) |
| **cLISI** | inverse Simpson index of labels in each cell's perplexity-30 kNN neighbourhood; (n_labels − median) / (n_labels − 1) | 0–1, higher | Harmony LISI; **scib**; `clisi` |
| **Graph modularity / connectivity** | modularity of the label partition on the kNN graph; connectivity penalty for split types | higher modularity better | classic; candidate |

**scTypeEval's key framing** (Garnica, Andreatta & Carmona 2026, benchmark over 31 datasets - §3b): a **local** ISC (silhouette in reciprocal-classification-match space, ISCL) catches **over-partitioning**, a **global** ISC (2label-silhouette in pseudobulk-cosine space, ISCG) catches **signal degradation**; the paper integrates them as the **product** (the package README says geometric mean - same ranking, different scale). The important idea for us: *neighborhood purity / silhouette computed on the labels are a direct, reference-free "is this type real?" signal* — under-used in FGCZ spatial reports, and cheap.

Interpreting on sparse panels: silhouette on a raw sparse Xenium matrix is noisy — compute it on the **PCA/latent embedding** (or on pseudobulk profiles for a per-type read), never on raw counts, and remember the `maxRank=150` sparsity caveat when scores feed in.

`internal_validity` computes neighbourhood purity, orbital-medoid and average similarity on every cell and the O(N^2) parts (silhouettes, Ward-PropMatch, Davies-Bouldin, Calinski-Harabasz) on a `subsample` (default 2,000; 5,000 inside `annotation_quality`). An undefined silhouette (one label or only singletons in the subsample) is NaN with a reason, never 0.

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

**In this package**: `inter_sample_consistency(adata, label_key, sample_key, *, min_cells=10, min_samples=5, layer=None, reduction="pca", n_pcs=30, rcm_genes=None)` returns `iscg`, `iscl`, `integrated` (each a mean over labels, as the paper), per-type dicts, the pseudobulk-weighted means (`iscg_profile_mean`, `iscl_profile_mean`), a per-type `regime` from `isc_regime()` flagged `regime_status = "experimental"`, and `dropped_types` (labels below the minimums are reported, not folded in as 1-cell pseudobulks). Parity with scTypeEval 1.0.0 on a Yerly 2022 subsample (`tests/goldens/isc_yerly2022`): ISCG equals scTypeEval's default 30-PC ISCG within 1.4e-6 (irlba vs exact SVD) and, with `reduction=None`, its `reduction = FALSE` ISCG within 1.3e-14; given scTypeEval's RCM distances, ISCL equals its ISCL within 1e-16; the RCM rule (0 for a mutual best hit, 0.5 read as 1 when either label is absent from the other sample, no match for a sample with < 2 profiles) is scTypeEval's. The classifier differs: nearest profile by Spearman over all genes (or `rcm_genes`) in place of SingleR with fine-tuning and pruning. On the authors' labels both make the same mutual hits on 1,128 / 1,128 pseudobulk pairs; on a randomly over-split label the two halves are near-identical and 30 / 1,891 pairs differ (28 with `rcm_genes="markers"`), so per-half ISCL is classifier noise (up to 0.30 apart) while the mean over labels agrees within 0.0002. Do not run ISCL on toy data with flat, low-noise background genes: Spearman ranks there let one foreign cell lift its type's genes above the background, and ISCL collapses (8 types x 64 genes: 0.92 -> 0.14 at 10 % relabelling; 1.00 with a graded per-gene baseline). Regime calls are experimental: on the synthetic sweep only over-splitting reads right (`merge`, 6 / 6); merged labels read `consistent` (Fibro+Epi ISCG 0.984, ISCL 1.000), 20 % relabelling reads `split` for 13 / 18 labels, and a per-sample batch effect `split` for 18 / 18. ISCL on a targeted spatial panel (a few hundred genes) is not validated by the paper; run the over-partitioning positive control (`tests/test_isc.py`, task-2 replica) on your own data before trusting a threshold.

This is genuinely additive to `cell-annotation-qc.md` (which has no multi-sample notion) and to FGCZ practice: when a report has replicates, ISC answers "is this cluster a real cell type or a batch/over-clustering artifact?" without any atlas.

<a name="3c-purity"></a>
### 3c. Marker / contamination purity (reference-free)

Do a cell's transcripts belong to its assigned type, or are they a spillover mixture? These need only **marker gene sets** (no reference atlas). Detailed in [`cell-annotation-qc.md` §6](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md); summarized here because they are core annotation-quality metrics, and several are already implemented in [`compare_rctd_tacco.py`](https://github.com/p-gueguen/spatial-scribe/blob/main/scripts/split/compare_rctd_tacco.py).

| Metric | Definition | Dir | Implemented |
|---|---|---|---|
| **MECR** (Mutually Exclusive Co-expression Rate) | over lineage marker pairs, `#(both>0)/#(either>0)`, averaged; each pair is judged on the markers the two sets do not share (CD3E in both CD4 T and CD8 T sets is not co-detection), nested sets (T vs CD4 T) are skipped; reported with `mecr_null` and `vs_null = mecr / mecr_null` | lower cleaner | `purity.py` (mecr); SPLIT compare |
| **CRISP purity** | per-cell: impure if it detects markers of ≥2 lineages; purity = 1 − N_impure / N_marker-positive (Rose et al. 2026). Our variant detects at >0 counts (CRISP default: >2) and judges each lineage pair on its unshared markers - not comparable to published CRISP scores. Reported with a detection-matched null and `vs_null` = observed / null impurity | higher cleaner | `purity.py` (`crisp_purity`) |
| **PMP** (Positive Marker Purity) | per-cell fraction of counts from the **assigned** type's markers | higher cleaner | `purity.py` (`pmp`) |
| **On-target marker specificity** | on-target fraction of a lineage's marker expression: on/(on+off) mean, bounded [0,1] | higher cleaner | `compare_rctd_tacco.py` (`_specificity`) |
| **On-target signal retention** | fraction of a lineage's own-marker counts kept after a purification step | higher = less over-purified | `compare_rctd_tacco.py` |

**The null, and what it cannot do.** At the ≥ 1-count rule, CRISP purity is 0.80 on a perfectly clean synthetic section (off-target counts at 0.1 % of the marker rate) and falls with depth alone, so each of `crisp_purity` and `mecr` is reported against a detection-matched null: every marker gene keeps its number of detecting cells, reassigned to random cells (20 draws). `vs_null` = observed co-detection / null co-detection reads 0.22 (CRISP) and 0.12 (MECR) on the clean section, rises with spillover (0.36 / 0.18 at 0.4) and ambient RNA (0.98 / 0.70 at 0.2), and stays stable under label changes. Genuine co-expression raises it as much as spillover (0.35 / 0.17 when Myeloid truly expresses Fibroblast markers at 40 %): no reference-free metric can tell the two apart. NMP with a matched reference can, because the reference carries the co-expression: Salas NMP stays ≥ 0.9996 in that scenario while it drops to 0.968 at spillover 0.4 (`benchmarks/synthetic`, 3 seeds).

These are the metrics the RCTD-vs-TACCO-into-SPLIT comparison uses, and they are the reference-free way to *quantify contamination* — the dominant iST annotation-quality problem.

<a name="3d-confidence"></a>
### 3d. Per-cell confidence / uncertainty

Per-cell "how sure are we about *this* label?" Detailed in [`cell-annotation-qc.md` §8 + §10](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md); the metrics:

- **Marker-score margin** = top1 − top2 marker score. **The margin, not the absolute score, is the confidence.**
- **Entropy** of the class-score distribution (near log K → the classifier is guessing).
- **Reference posterior max-prob** (CellTypist-style) — where a reference exists.
- **Conformal label-set size**: distribution-free. Set size 1 = confident, >1 = ambiguous, 0 = novel; it is the cleanest formalization of "abstain". Per-class (Mondrian) thresholds keep coverage on rare types. With many classes and few calibration cells each, clustered Mondrian (Ding et al. 2023) pools classes whose score distributions look alike. CovGap (mean |per-class coverage - target|) is the number to compare them on. `conformal_prediction_sets` in `conformal.py`.
- **Calibration error (ECE / class-wise ECE)** (Guo et al. 2017): before reading a softmax or an RCTD weight as a probability, check that cells given 0.8 are right about 80 % of the time. `calibration_error` in `conformal.py`.
- **Conformal novelty p-values** (Bates et al. 2023): an empty conformal set is a weak novelty signal, because depth loss empties sets for known types too. Synthetic check, B left out of the reference, 3 seeds:
  - At 1/4 section depth, 18 % of known cells get an empty set, against 100 % of the novel B cells.
  - The novelty p-values (PCA reconstruction error, BH 10 %, reference thinned to the section's per-cell depths) separate B at AUC 0.99. They flag 1.2 % of known cells and 85 % of B cells, at a false-discovery proportion of 0.10.
  - Without the depth matching, 87-100 % of known cells are flagged at every depth.

  `conformal_novelty` in `conformal.py`.
- **Consensus entropy with ontology-path voting** (popV, Ergen et al. 2024): when several annotators label the same cells, each label also votes for its ancestors, and the consensus is the deepest label with a majority. The entropy of the votes, projected to that depth, measures disagreement without penalising subtype splits. `consensus_entropy` in `consensus.py`.
- **Annotation stability under subsampling** — drop ~20% of a cell's transcripts, re-annotate, measure label-flip rate; reference-free, evidence-backed (MERFISH: ~20% dropout flips 10–15% of labels). Implemented in `annotate.py` (`annotation_stability`).

<a name="3e-spatial"></a>
### 3e. Spatial-context metrics (spatial only)

Does the label make **spatial sense**? Blind to non-spatial metrics. Detailed in [`cell-annotation-qc.md` §9](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md):

- **Spatial coherence** = mean fraction of a cell's spatial neighbours (k = 10) sharing its label (`spatial.py` `spatial_coherence`).
- **PAS** (Percentage of Abnormal Spots; Shang & Zhou 2022, [10.1038/s41467-022-34879-1](https://doi.org/10.1038/s41467-022-34879-1)) = % cells whose label differs from at least 6 of their 10 nearest neighbours (same-label fraction <= 0.4). 0.6.x used < 20 % of 15 neighbours, which reads systematically lower.
- **CHAOS** (same paper) = mean distance from each cell to its nearest same-label cell, on per-axis standardised coordinates. Lower = spatially smoother labels.
- **Neighborhood enrichment / co-occurrence** (squidpy) vs a permutation null — biologically impossible adjacencies flag systematic errors.
- **Neighbor contamination** (`spatial.py` `neighbor_contamination`, from SegTraQ's `contamination_strength`) — per cell, the share of counts on markers of the cell types *touching it* that its own type lacks (with a touching cell of that type detecting them), plus a directed source → target matrix. Reference-free; separates spillover from neighbors from a mislabel, which PMP / NMP cannot. The raw value does not separate spillover from genuine co-expression: on a synthetic two-block grid, a B marker every A cell truly expresses scores 0.037 and real border spillover 0.041. The `baseline` (same-type cells with no such neighbour: 0.033 vs 0.001) and the `excess` over it (0.004 vs 0.040) do; the one-sided permutation `p_value` is 0.001 for spillover and 0.19 for co-expression with no border excess (it is small, 0.003, when border cells carry even a +0.003 real excess: with large n read the excess, not the p). Nested marker sets (`T` in `CD4 T`) are skipped, as in MECR.

**Caveat (load-bearing):** coherence, PAS and CHAOS detect *random* label noise only. They are blind to a compact mislabelled patch, they *rise* when two neighbouring types are merged into one label (0.904 → 0.961 in the synthetic `merge` scenario, mean coherence 0.918 → 0.966 and PAS 0.024 → 0.011 from 0 to 3 merges), and they fall for correctly labelled infiltrating cells (`infiltration`: PAS 0.024 → 0.163). PAS stays at 0.023-0.025 whatever the size of a compact mislabelled patch (`patch_error`). Use them to **down-weight**, never to hard-delete or to rank annotations by correctness.

<a name="3f-anchored"></a>
### 3f. Reference-anchored specificity (needs a matched profile, not per-cell labels)

Between purely internal and fully external: these need a **matched scRNA/snRNA reference profile** but not per-cell ground truth. From Salas et al. 2025:

- **NCP** (Negative Co-expression Purity) — % of gene pairs **not** co-expressed in the reference that stay non-co-expressed in situ. 0–1, higher = better specificity; commercial platforms mean **>0.8**. `purity.py` `ncp` follows the paper's wording (fraction of pairs); the authors' code averages a clipped excess instead (`tests/goldens/coexpr_txsim/SOURCE.md`).
- **Co-expression similarity** (txsim) — `1 − |ρ_reference − ρ_section| / 2` per gene pair (Pearson on log-normalised counts), averaged; optionally within each cell type. In `purity.py` (`coexpression_similarity`), matched to txsim in `tests/goldens/coexpr_txsim`. On the synthetic sweep it barely moves with spillover (0.979 → 0.978 at 0.4); prefer NCP / NMP for contamination.
- **NMP** (Negative Marker Purity) — built on the reads from "negative markers" (genes a type should NOT express) found in cells assigned to that type; the most direct readout of segmentation-driven mis-assignment, and the recommended metric for **comparing segmentations**. In `purity.py`, `nmp` computes it as Salas et al. define it and txsim implements it (`negative_marker_purity_reads`, matched to ≤ 1e-18 in `tests/goldens/nmp_txsim`): per type, the type-balanced share of each gene's mean count; a (gene, type) pair is negative when the reference type holds < 0.5 % of the share; NMP = 1 − the section's excess share on those pairs, scaled so 0 = the mean over 20 random label permutations and 1 = no more than the reference (the paper's scaling; its baseline comes from one binning reference, ours from the section itself). Higher cleaner. Because the reference carries genuine co-expression, it is the purity metric to trust where cell types really share markers. `score="raw"` keeps the v0.6.0 per-cell score: the fraction of a cell's counts on genes the reference detects in ≤ 5 % of its type's cells (SegTraQ `markers_from_reference`'s cut) and ≥ 25 % of another type's (a package choice: SegTraQ uses 0.1 plus a DE step, a 10-cell minimum and drops MT-/RPL/RPS genes); cells it cannot score are NaN. The earlier bottom-decile definition picked near-silent genes on large panels: 10% fibroblast spillover into T cells scored 0.004, same as clean, on a synthetic 5,000-gene panel.

> Terminology trap (kept from `cell-annotation-qc.md`): **NCP** (negative *co-expression* purity, reference-anchored) ≠ the 10x **negative *control-probe* rate** (background binding). Different things; keep separate.

<a name="3g-marker-fidelity"></a>
### 3g. Marker-program fidelity (needs marker sets, not ground truth)

A powerful family the Zhu 2026 benchmark formalizes: given curated marker/program sets, you can measure whether the **labels track the biology without per-cell ground truth** — so these run on real deliveries, not just benchmarks. They are the §2.5 fidelity metrics in their reference-free form:

- **ssGSEA enrichment** (GSVA) — score each cell against each type's marker/pathway set.
- **AUC-ROC on enrichment** — per type, how well the method's predicted-positive cells are ranked above the rest by that type's enrichment score (1 = perfect ranking, 0.5 = chance). A clean, threshold-free "do the cells labeled *T cell* actually carry the highest T-cell-program score?" This **generalizes the on-target specificity** in [`compare_rctd_tacco.py`](https://github.com/p-gueguen/spatial-scribe/blob/main/scripts/split/compare_rctd_tacco.py) into a proper ranking metric.
- **Cohen's d** — standardized target-vs-non-target separation on the enrichment score (effect size, not just AUC).
- **Marker-gene overlap** — overlap of a label's data-derived top-N DE genes with its curated markers.

Use these when you have marker sets but no truth (the usual FGCZ case), as the **fidelity** complement to the **purity** metrics in [§3c](#3c-purity): purity asks "is the cell contaminated?", fidelity asks "does the label match the cell's program?".

<a name="3h-signal-qc"></a>
### 3h. Signal-quality QC (annotation-independent; SpatialQM)

Upstream of annotation quality is *signal* quality: does the section carry real, structured biology above background? These are **label-independent** - they score the data, not the labels - and are ported from the **SpatialQM** R package (Plummer et al., Center for Spatial Omics, St. Jude; released with the Spatial Touchstone study - see §9). Implemented in `signal_qc.py`, **squidpy-compatible** and fast:

- **Which features are controls.** They are flagged by feature type, never by name: in the 10x Xenium breast
  panel, 8 of the 28 negative-control probes are named `antisense_<GENE>`, and SpatialQM's name regexes count them
  as genes (`getMaxRatio` 4.333 vs 4.001, FDR 0.593 % vs 0.664 %). Rule: `var['control']` if present, else
  `var['feature_type'] != "Gene Expression"`.
- **Moran's I gene-vs-control** (`getMorans`). Real genes are spatially autocorrelated and controls sit at
  I ~ 0; the separation of the medians is the score. Features with undefined I (no count) are left out and
  counted. With shuffled coordinates the separation collapses to ~0 (test). Two methods:
  - `"squidpy"` (default): raw counts, binary 6-NN graph, all controls.
  - `"spatialqm"`: SpatialQM's recipe, reproduced to 5e-15. Its control median rests on the few cells with a
    probe count (150 of 2000 on the crop), so it sits near E[I] = -1/(n-1).
- **Signal-to-noise** (`getMeanSignalRatio` / `getMaxRatio`):
  - Default: per cell, log2(mean gene / mean control), defined only where the cell has a control count. 81.9 %
    of the full Xenium section's cells and 88.8 % of the crop's have none. The result therefore reports
    `n_cells_with_controls` and the pooled ratio, and never an eps floor (ROADMAP 8g).
  - `method="spatialqm"`: SpatialQM's log10 dataset scalars.
- **Sparsity** (`getSparsity`): zero fraction, matches to 3e-16.
- **Detection entropy:** per-cell normalised entropy by default; `method="spatialqm"` = `getEntropy` (bits of
  the count-value histogram of the whole matrix).
- **Transcripts per area / per cell** (`getTxPerArea`, `getTxPerCell`): mean (SpatialQM's) and median.
- **Transcript QC** (`transcript_qc`, #21): FTC and the negative-codeword FDR from a transcripts table.
  SpatialQM divides its FDR by 100; we report the fraction.

These gate whether annotation should be trusted *at all* on a section (they sit alongside Layer 0-3 of the spatial funnel); run them before the label-quality metrics above. `signal_qc.run_signal_qc(adata)` returns the headline.

<a name="3i-ontology"></a>
### 3i. Ontology label-relation consistency (needs a label -> CL id map, not ground truth)

When the labels are Cell Ontology terms, their relations are a prior on the embedding: a T cell should sit
nearer a B cell than a keratinocyte. `ontology.ontology_consistency` takes the labels' embedding centroids and
the `is_a` path length between their CL terms (the shortest walk up to a common ancestor, read from a
user-supplied OBO file with the standard library). It then reports the Spearman correlation over label pairs
(`rho`) with a one-sided Mantel permutation p. Each label also gets a leave-one-out `loo_delta`, the rise in
`rho` without it. The label with the largest rise is `most_inconsistent`, and `nearest_embedding` vs
`nearest_ontology` shows why: a "macrophage" whose nearest centroid is epithelial while its nearest ontology
term is the monocyte.

Read it as a layout score, not a verdict on any label. Raw CL path lengths are coarse: on cl-basic 2026-06-08,
macrophage-epithelial cell and macrophage-T cell are both 4 steps. Siblings that a good embedding separates
also lower `rho`. With K labels the p cannot go below 1/K!.

On the Xenium FFPE breast section, with marker-argmax labels for 10 types (12,031 confident cells, PCA 30):

- clean labels: `rho` 0.28 (p 0.08);
- planting the "macrophage" label on 1,500 luminal cells: `rho` 0.12, and `most_inconsistent` = macrophage;
- a shuffled label -> CL map: `rho` -0.16.

The order is right, but the effect is modest. ROADMAP #25 cites Genome Biol 2025,
[10.1186/s13059-025-03781-6](https://doi.org/10.1186/s13059-025-03781-6), for the idea. That paper's ontology
metric is scGraph-OntoRWR (a random walk with restart over the ontology graph). `ontology_consistency` is a
different, own statistic (`is_a` shortest path, Spearman, Mantel test), not a port of it.

---

<a name="4-batteries"></a>
## 4. What we standardize on — default metric batteries

The authoritative recommendation: don't compute everything; compute the right battery for the scenario. Each battery is ordered (compute the first few always; add the rest if cheap/available).

**A. Spatial delivery, no reference (the typical Xenium/CosMx report):**
1. **Panel adequacy** (can the panel resolve the type at all — [`panel_check.py`](https://github.com/p-gueguen/spatial-scribe/blob/main/src/spatialscribe/analysis/panel_check.py))
2. **MECR + CRISP purity + PMP** (reference-free contamination)
3. **Marker-score margin + entropy** (per-cell confidence)
4. **Spatial coherence / PAS** (down-weight only)
5. **Subsampling label-stability**
6. *(if ≥5 replicates)* **Inter-sample consistency (ISCG + ISCL, regime per label)**
7. *(reporting)* section **annotatability headline**: % PASS / WARN / abstain + top reasons

**B. Spatial with a matched snRNA reference:** battery A **plus** RCTD `spot_class` + weight-margin, **NCP/NMP**, label-transfer posterior margin. (This is the highest-quality regime; reference choice matters more than method — see [`annotation-method-selection.md`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/annotation-method-selection.md).)

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
| Accuracy / balanced accuracy | E | L+truth | per-cell correctness | high | `external.py` (`external_scores`) |
| Precision / recall / **macro-F1** (+ micro=acc, weighted) | E | L+truth | per-class quality; rare-type sensitivity | high | `external.py` (`external_scores`) |
| **Hierarchical accuracy** · subtype-level accuracy | E | L+truth+hierarchy | credit by LCA depth (parent call 2/3, wrong sibling 1/2) + coarse-call accuracy | high | `external.py` (`hierarchical_accuracy`) |
| **spARI** | E | L+truth+coords | ARI with distance-weighted credit for disagreeing pairs | ≤1, ↑ | `external.py` (`spatial_ari`) |
| **ARI** / AMI / NMI / V-measure | E | L+truth | partition agreement | high | `external.py` (`external_scores`) |
| **ECS** (Element-Centric Similarity) | E | L+truth | resolution-agnostic, cell-level partition agreement | high | `external.py` |
| **AvgBIO** = mean(ARI,NMI,ASW_scaled) · **ASW** | E | L+truth+embed | biological-structure conservation | high | `external.py` (`avg_bio`) |
| **Composition / proportion accuracy** | E | L+truth | predicted vs true cell-type proportions (hard-label global) | high | `external.py` |
| **Deconvolution props** (R²/JSD/RMSE, OpenProblems; AUPR, SSIM) | E | proportion matrices | soft per-spot proportions vs known mixtures | R²↑/JSD↓/AUPR↑/SSIM↑ | `external.py` (`deconvolution_metrics`, `deconvolution_rank_score`) |
| Decomposition purity | I | deconvolution weights | 1 − largest normalised weight per cell (segmentation / reference-gap readout) | low | `decomposition.py` (`decomposition_purity`) |
| **Robustness to reference volume** | E | L+truth+R | accuracy decay as reference shrinks | stable | — |
| Scalability (runtime · peak memory) | — | — | practical cost, not quality | lower | bench scripts |
| SingleR delta · CellTypist prob · RCTD spot_class | E/I | R | native confidence proxy | high | SpatialScribe only (`methods.py` join) |
| **Silhouette / 2label-silhouette** | I | L+embedding | cohesion vs separation | high | `internal.py` |
| **Neighborhood purity** | I | L+embedding | kNN label agreement | high | `internal.py` |
| Orbital-medoid · Ward-PropMatch · Avg-similarity | I | L+embedding | cluster–label concordance | high | `internal.py` |
| Davies–Bouldin · Calinski–Harabasz | I | L+embedding | separation | DB low / CH high | `internal.py` (`internal_validity`) |
| **ISCG** (2-label silhouette \| pseudobulk cosine) | I | L+S (≥5 samples) | cross-replicate compactness (signal degradation) | high; guide 0.74 | `isc.py` |
| **ISCL** (silhouette \| reciprocal-classification match) | I | L+S (≥5 samples) | cross-replicate replicability (over-partitioning) | high; guide 0.45 | `isc.py` |
| **MECR** | I | M | cross-lineage co-expression | low | `purity.py` |
| **CRISP purity** | I | M | per-cell lineage purity | high | `purity.py` |
| **PMP** | I | M | own-marker transcript fraction | high | `purity.py` |
| On-target specificity / retention | I | M | marker specificity after purification | high | SpatialScribe only (`compare_rctd_tacco.py`) |
| **AUC-ROC on marker enrichment** · Cohen's d | I | L+M | do labels track the marker-program enrichment | high | `fidelity.py` |
| **Marker-gene overlap** (top-N DE vs markers) | I | L+M | labels reproduce known markers | high | `fidelity.py` (`marker_gene_overlap`) |
| Marker-score **margin** · entropy | I | L+M | per-cell ambiguity | margin high / entropy low | `confidence.py` (`per_cell_confidence`) |
| Conformal set size | I | L+scores | abstention (ambiguous/novel) | =1 | `conformal.py` (`conformal_prediction_sets`) |
| Calibration error (ECE / class-wise ECE) | E | scores+truth | is confidence a probability? | low | `conformal.py` (`calibration_error`) |
| Conformal novelty p-value | I | R+features | novel / out-of-reference cells, FDR-controlled | high p = known | `conformal.py` (`conformal_novelty`) |
| Consensus entropy (ontology voting) | I | several L (+hierarchy) | agreement across annotators | low | `consensus.py` (`consensus_entropy`) |
| Subsampling label-stability | I | counts | detection-robustness | low flip | `confidence.py` (`label_stability`) |
| Spatial coherence / **PAS** / CHAOS | I | L+XY | spatial label agreement | high / PAS, CHAOS low | `spatial.py` |
| Neighborhood enrichment sanity | I | L+XY | impossible adjacencies | — | `spatial.py` |
| Neighbor contamination | I | L+XY+M | spillover from adjacent types | low | `spatial.py` |
| **Moran's I gene-vs-control** | I | XY+control | spatial signal vs background noise | high separation | `signal_qc.py` |
| Signal-to-noise · sparsity · detection entropy · tx/area · tx/cell | I | counts (±control/area) | label-independent signal quality | mixed | `signal_qc.py` |
| Transcript QC (FTC, negative-codeword FDR) | I | transcripts table | share of transcripts in cells, false-discovery rate of codewords | FTC high / FDR low | `signal_qc.py` (`transcript_qc`) |
| **NCP** / **NMP** | A | R+M | reference-anchored specificity | NCP high / NMP low | `purity.py` (`ncp`, `nmp`) |
| Panel adequacy (coverage/confusable) | I | M+panel | can the panel resolve the type; tiers on chance-corrected F1 | resolvable | `panel.py` (`panel_resolvability`) |
| Ontology label-relation consistency | I | L+emb+CL map | do label centroids follow Cell Ontology distances | high rho, low p | `ontology.py` (`ontology_consistency`) |
| ovrlpy VSI / 3D-coherence | I | Z | vertical doublets/overlap | high | SpatialScribe only (`qc.py` join) |

---

<a name="6-implementation"></a>
## 6. Implementation status

**In SpatialScribe — contamination + confidence + spatial** (reference-free): `purity.py` (MECR, CRISP, PMP, NMP/NCP guards), `annotate.py` (marker margin, entropy path, subsampling stability, penalized-posterior confidence + abstention), `panel_check.py` (panel adequacy + identifiability AUC), `spatial.py` (spatial coherence, PAS, neighborhood sanity, neighbor contamination), and `scripts/split/compare_rctd_tacco.py` (on-target specificity, marker purity, on-target retention, MECR panel).

**Metric battery** — this package, one module per family: `internal.py`, `isc.py`, `fidelity.py`, `external.py`, `panel.py`, `conformal.py`, `consensus.py`, `ontology.py`, `battery.py` (**built, TDD**: `tests/test_<module>.py`). Until 0.7.0 these lived in one `eval_metrics.py`, which is now a deprecated shim that re-exports the old names. Covers the families this doc exposed:

```python
# reference-free internal validity (scTypeEval family) on an embedding + labels
internal_validity(adata, label_key="cell_type", embedding="X_pca", n_neighbors=30) -> dict
    # {silhouette, silhouette_2label, neighborhood_purity, orbital_medoid,
    #  ward_propmatch, avg_similarity, davies_bouldin, calinski_harabasz,
    #  integrated (geometric mean of the defined parts), nan_reasons}
isolated_label_scores(adata, label_key, *, batch_key=None) -> dict  # scib isolated-label ASW / F1
clisi(adata, label_key, *, n_neighbors=90) -> dict                   # scib cLISI
inter_sample_consistency(adata, label_key, sample_key, *, min_cells=10, min_samples=5,
                         reduction="pca", n_pcs=30) -> dict  # ISCG / ISCL / integrated / regime
marker_program_fidelity(adata, label_key, marker_sets) -> dict     # §3g: per-type AUC-ROC + Cohen's d
external_scores(*, y_true, y_pred) -> dict     # balanced_acc, macro/weighted/per-class F1, ARI/AMI/NMI, kappa, ECS
element_centric_similarity(a, b, alpha=0.9) -> dict                # exact vectorized ECS (Gates & Ahn)
hierarchical_accuracy(*, y_true, y_pred, hierarchy, partial=0.5, credit="lca") -> dict  # credit by LCA depth
spatial_ari(*, y_true, y_pred, coords, alpha=0.8) -> dict           # spRI / spARI
composition_accuracy(*, y_true, y_pred) -> dict                    # hard-label global proportions: L1 / Pearson / JSD
deconvolution_metrics(*, true_prop, pred_prop) -> dict             # soft per-spot proportions: R2 / JSD / RMSE / AUPR / SSIM
deconvolution_rank_score(results, metrics=("rmse", "jsd_per_spot", "aupr")) -> dict
conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha=0.1) -> dict
calibration_error(proba, classes, y_true) -> dict                  # ECE / class-wise ECE
conformal_novelty(X_query, X_ref, *, method="pca", fdr=0.1) -> dict  # novelty p-values, BH
consensus_entropy(adata, label_keys, *, hierarchy=None) -> dict    # ontology-path voting
ontology_consistency(adata, label_key, obo_path, label_to_cl) -> dict  # §3i
annotation_quality(adata, label_key, marker_sets=None, sample_key=None) -> dict   # reference-free battery
```

`internal_validity` runs the O(N²) silhouettes on a `subsample` for scalability and falls back to a quick PCA when no embedding is present; `external_scores` is a thin `sklearn.metrics` wrapper plus the exact ECS closed form. The reference-free battery (`annotation_quality`) is **wired into `qc.run_funnel`** (guarded, non-breaking) so every labeled section's funnel headline carries `annotation_quality.internal_validity` + `marker_fidelity`; `external_scores` / `hierarchical_accuracy` / `composition_accuracy` / `deconvolution_metrics` are for benchmark scripts with ground truth. Tested against well-separated vs shuffled fixtures (clean labeling → purity≈1, high integrated; shuffled → near-chance).

**Signal-quality QC** — `signal_qc.py` (**built, TDD**: `tests/test_signal_qc.py`; §3h) ports the label-independent SpatialQM metrics, squidpy-compatible: `moran_signal` (Moran's I gene-vs-control via `sq.gr.spatial_autocorr`, or SpatialQM's `getMorans` with `method="spatialqm"`), `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, `tx_per_cell`, `transcript_qc` (FTC and FDR from a transcripts table), and `run_signal_qc` (headline). `decomposition.py` adds `decomposition_purity` (§2.6). Standalone (run before annotation to gate whether a section is trustworthy at all); not auto-wired into the funnel because Moran's I over the full panel is heavier than the per-cell metrics.

> **Two scTypeEval implementations — pick the right one.** `internal_validity` is an **embedding-space approximation** built for a fast, always-on `qc.run_funnel` headline: it computes the family on the PCA embedding (subsampled) and does **not** implement scTypeEval's dissimilarity spaces (the pseudobulk-level ISCG / ISCL pair, with a Spearman-classifier `recip_classif:Match`, lives in `inter_sample_consistency` since v0.5.0), and its `orbital_medoid` / `avg_similarity` / `ward_propmatch` use classic embedding-space definitions rather than the exact `IntClusVal.R` ones. A **byte-exact port** (dissimilarity spaces + `recip_classif:Match` + exact IntClusVal formulas + per-sample pseudobulk, TDD-verified) lives in `scripts/norm_benchmark/sctypeeval_metrics.py` (SpatialScribe, not public) — use it when you need scTypeEval's recommended local(recip-match)/global(pseudobulk-cosine) combo or paper parity; use `internal.py` for the production funnel. Empirically the two **rank annotations/normalizations identically but differ in absolute scale**, so the approximation is fine for a QC headline. See `scripts/norm_benchmark/README.md` (SpatialScribe, not public) — the **normalization × annotation-quality benchmark** on the 100k breast Xenium 5K section, which anchors to this catalog and reports both implementations plus the implemented contamination/confidence/spatial families across 9 normalizations. Headline finding: normalization is **~invariant** for the contamination metrics (§3c, they read raw counts) but has a **strong** effect on internal separability (§3a/b), the abstention rate (§3d), and above all the clustering the annotation is built on.

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
- **Internal metrics on a label-trained embedding are circular** (scANVI, a fine-tuned foundation model): the
  embedding was optimised to separate these labels, so silhouette / kNN purity / cLISI / ISC score the training,
  not the annotation. Compute them on a label-agnostic embedding (PCA of the counts). `internal_validity(...,
  embedding_label_trained=True)` flags it. (Metric Mirages, doi 10.1101/2024.04.02.587824, preprint.)
- **cLISI hides a small label, and the silhouette does not say which label is wrong**: an 8-cell
  label planted inside another type (188 cells) leaves cLISI at 1.00 (a median over cells); the silhouette drops
  (0.86 -> 0.53) without naming the label, while that label's isolated-label ASW is 0.53 (0.98 when it stands
  apart) and its F1 0.32 (1.00). Read `isolated_label_scores` when rare types matter.
- **ISC needs ≥ 5 replicates and ≥ 10 cells per (sample, label)** (Garnica 2026). Below that, report a label as *not computable*, never as 0 or as "inconsistent": a 1-cell pseudobulk is noise, not a replicate.
- **Read ISC on two axes, not one number.** Low ISCL alone → merge, low ISCG alone → split, both low → remove; the product is for ranking annotators or references, the quadrant is for fixing labels.
- **Individual-specific populations (malignant cells, patient-private clones) have low ISC by construction** - exclude them from an ISC-based quality verdict, as the benchmark's own preprocessing excluded cancer-cell labels.

---

<a name="9-references"></a>
## 9. References

**This catalog builds on the FGCZ research docs** — see them for the spatial-funnel detail and verified sources:
- [`cell-annotation-qc.md`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md) — the six-layer spatial QC funnel; source for Layers 0–6, contamination metrics (MECR/CRISP/PMP/NCP/NMP), confidence/abstention, spatial coherence, and all verified primary sources (Salas *Nat Methods* 2025; Plummer/Spatial-Touchstone *Nat Biotechnol* 2025; Mitchel *Nat Genet* 2026; SPLIT *Nat Methods* 2026).
- [`annotation-method-selection.md`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/annotation-method-selection.md) — which annotator to run; source for the accuracy≠fidelity finding and the 20-method spatial benchmark (Zhu et al., bioRxiv 2026; AvgBIO / program concordance).
- [`annotation_qc_thresholds.yaml`](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/annotation_qc_thresholds.yaml) — the tunable pass/warn/fail numbers consumed by the code.

**scTypeEval (the metrics folded in here):**
- Carmona lab, `scTypeEval` — reference-free annotation-quality assessment (silhouette / 2label-silhouette / neighborhood purity / Orbital-medoid / Ward-PropMatch / average-similarity; inter-sample consistency; local/global/integrated framing). GitHub: [github.com/carmonalab/scTypeEval](https://github.com/carmonalab/scTypeEval) `[tool; benchmarked in the preprint below]`. Lab: [github.com/carmonalab](https://github.com/carmonalab).
- **Garnica J, Andreatta M, Carmona SJ.** "Evaluating cell type annotations in single-cell omics in the absence of ground truth," bioRxiv 2026-06-12, DOI [10.64898/2026.06.10.731285](https://doi.org/10.64898/2026.06.10.731285) `[preprint, not yet peer-reviewed]`. Source for ISCG / ISCL / RCM, the 36-combination × 31-dataset benchmark, the 0.74 / 0.45 guide thresholds, the four refinement regimes, and the ISC ↔ F1 R² > 0.98 result (§3b), plus the benchmark-design lessons in §7. Code: [carmonalab/ISC_benchmark_reproducibility](https://github.com/carmonalab/ISC_benchmark_reproducibility) (`targets` pipelines, `renv`-pinned); data: Zenodo [10.5281/zenodo.18921437](https://doi.org/10.5281/zenodo.18921437).

**Classic metric definitions (external + internal validity):**
- **ARI** — Hubert & Arabie, *J. Classification* 1985. **AMI/NMI** — Vinh, Epps & Bailey, *JMLR* 2010.
- **Silhouette** — Rousseeuw, *J. Comput. Appl. Math.* 1987. **Davies–Bouldin** 1979; **Calinski–Harabasz** 1974.
- **F1 / precision / recall / balanced accuracy / Cohen's κ** — standard classification metrics (`sklearn.metrics`).
- **spARI** — Yan, Feng & Luo, "Spatially aware adjusted Rand index for evaluating spatial transcriptomics clustering", *Biometrics* 2025, doi 10.1093/biomtc/ujaf127 (preprint 10.1101/2025.03.25.645156).
- **Hierarchical precision / recall** — Kiritchenko et al. 2006 (Canadian AI 2006, LNCS 4013).
- **Deconvolution proportion metrics** — **OpenProblems `task_spatial_decomposition`** ([github.com/openproblems-bio/task_spatial_decomposition](https://github.com/openproblems-bio/task_spatial_decomposition)): `r2` = `r2_score(uniform_average)` over cell types; `jsd` = mean per-cell-type `scipy.spatial.distance.jensenshannon` (axis=0). Ported into `external.deconvolution_metrics` (+ RMSE, per-spot JSD).
- **Spotless** — Sang-aram et al., *eLife* 2024, doi 10.7554/elife.88431 (presence/absence AUPR, rank aggregation).
- **Li et al.**, *Nat Methods* 2022, doi 10.1038/s41592-022-01480-9 (per-type SSIM, accuracy score).

**Benchmark-derived metrics — Zhu et al. 2026 (the paper read for this update):**
- Zhu, Hu, … Meltzer & Zhou, "Benchmarking cell type annotation in spatial transcriptomics: resolving cellular hierarchies, biological fidelity, and dynamic cell states," bioRxiv 2026-06-16, DOI [10.64898/2026.06.16.732716](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full) · local [bioRxiv full text](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full). Source for the **ECS**, **AvgBIO = mean(ARI,NMI,ASW_scaled)**, **hierarchical/subtype accuracy**, **composition/proportion accuracy**, **reference-volume robustness**, **AUC-ROC-on-enrichment / Cohen's d / marker-gene-overlap** fidelity metrics folded in here, and the accuracy≠fidelity finding.
- **ECS** — Gates & Ahn, "Element-centric clustering comparison," *Sci. Rep.* 2019 (resolution-agnostic, affinity-matrix formulation).
- **AvgBIO / ASW-based conservation** — Luecken et al., "Benchmarking atlas-level data integration" (scIB), *Nat. Methods* 2022.
- **ssGSEA** — Barbie et al., *Nature* 2009; **GSVA** — Hänzelmann et al., *BMC Bioinformatics* 2013. **Cohen's d** — Cohen 1988.

**Signal-quality QC (SpatialQM, §3h — the metrics ported into `signal_qc.py`):**
- **SpatialQM** — Plummer et al., "Standardized metrics for assessment and reproducibility of imaging-based spatial transcriptomics datasets" (Spatial Touchstone), *Nat Biotechnol* 2025, [10.1038/s41587-025-02811-9](https://doi.org/10.1038/s41587-025-02811-9); code [github.com/Center-for-Spatial-OMICs/SpatialQM](https://github.com/Center-for-Spatial-OMICs/SpatialQM) (St. Jude Center for Spatial Omics). R package; source for Moran's I gene-vs-control (`getMorans`), signal-to-noise (`getMeanSignalRatio` / `getMaxRatio`), sparsity, detection entropy, tx-per-area, and its own `getMECR`. `[peer-reviewed]`
- Moran's I — Moran 1950; computed here via squidpy `sq.gr.spatial_autocorr` (Palla et al., *Nat Methods* 2022).
- **SegTraQ** — [github.com/LazDaria/SegTraQ](https://github.com/LazDaria/SegTraQ) (MIT); `sp.rctd_contamination` is the source of `decomposition_purity`.
- **CRISP** - Rose, J. R., Rose, E. S., Assumpção, J. A. F., Pathak, H., Peck, H. E., Sasser, L. E., Patel, C. J., Vanover, D. & Santangelo, P. J. "CRISP enables comparisons of image-based spatial transcriptomic segmentation quality across ten organs", *bioRxiv* 2026, [10.64898/2026.04.16.718947](https://doi.org/10.64898/2026.04.16.718947); code [github.com/santangelo-lab/CRISP](https://github.com/santangelo-lab/CRISP) (R `crispRSeg`, Python `crispySeg`, MIT; Santangelo lab, Emory University and Georgia Tech). Source of the CRISP purity definition, built on MECR (Hartman & Satija 2024): a cell is impure when markers of ≥2 mutually exclusive cell types pass the detection threshold, and `purity = 1 - N_impure / N_marker-positive`. `crisp_purity` keeps that rule and formula with two changes: a marker counts as detected at >0 counts (CRISP: `> threshold`, default 2, compound markers summed first, threshold recommended at about half the mean marker expression), and each lineage pair is judged on the markers the two sets do not share (nested pairs skipped; CRISP expects curated disjoint lists). On disjoint marker sets it equals CRISP run at threshold 0; its values are not comparable to published CRISP scores. `[preprint]`
- **Segmentation admixture** - Mitchel, J., Gao, T., Petukhov, V., Cole, E. & Kharchenko, P. V. "Impact and correction of segmentation errors in spatial transcriptomics", *Nat Genet* 58, 434-444 (2026), [10.1038/s41588-025-02497-4](https://doi.org/10.1038/s41588-025-02497-4); code [github.com/kharchenkolab/cellAdmix](https://github.com/kharchenkolab/cellAdmix). Re-analysing imaging-based data from multiple tissues and platforms, it shows that segmentation errors put other cell types' marker transcripts into cells and that this admixture dominates cell-state analyses (niche differential expression, neighbour influence, ligand-receptor); cellAdmix removes it by factorizing subcellular molecular neighbourhoods, "analogous to doublet filtering". It is the phenomenon MECR / CRISP / PMP read out; transcript-level correction is out of scope here. `[peer-reviewed]`
- **txsim** - [github.com/theislab/txsim](https://github.com/theislab/txsim) (commit `bfdfaa90`; no licence file, so run for the goldens and never copied); cited by OpenProblems as [10.1101/2023.02.13.528102](https://doi.org/10.1101/2023.02.13.528102). Reference implementation of Salas et al.'s NMP (`negative_marker_purity_reads`) and of `coexpression_similarity`; `nmp` (default) and `coexpression_similarity` re-implement both from their definitions and match them on the synthetic goldens. `[tool]`

**Spatial / contamination metric sources** (full detail + DOIs in [`cell-annotation-qc.md` §17](https://github.com/p-gueguen/spatial-scribe/blob/main/docs/research/cell-annotation-qc.md)): Salas et al. *Nat Methods* 2025 (NCP, NMP, 3D-coherence); Hartman & Satija 2024 (MECR); Rose et al. *bioRxiv* 2026, [10.64898/2026.04.16.718947](https://doi.org/10.64898/2026.04.16.718947) (CRISP); Mitchel et al. *Nat Genet* 2026, [10.1038/s41588-025-02497-4](https://doi.org/10.1038/s41588-025-02497-4) (segmentation admixture, cellAdmix); Cable et al. *Nat Biotechnol* 2022 (RCTD); Bilous et al. *Nat Methods* 2026, [10.1038/s41592-026-03089-8](https://doi.org/10.1038/s41592-026-03089-8) (SPLIT); Shang & Zhou, *Nat Commun* 2022, [10.1038/s41467-022-34879-1](https://doi.org/10.1038/s41467-022-34879-1) (PAS, CHAOS); STEAM (spatial coherence).

> **Provenance note.** The spatial/contamination/confidence metrics and their evidence are carried over (with cross-links, not duplicated) from `cell-annotation-qc.md`, which was built from an adversarial-verification research pass. The **new** material here — the external-validation family (§2), the internal cluster-validity + ISC family (§3a–3b), the marker-program fidelity family (§3g), the decision guide (§1), the batteries (§4), the scTypeEval integration (§7), and the ECS / AvgBIO-form / hierarchical / composition / reference-volume-robustness metrics (§2.5–2.6) — is synthesized from standard metric definitions, the scTypeEval repository, and a **direct reading of the Zhu et al. 2026 benchmark Methods** (`sources/2026.06.16.732716v1.full.pdf`: exact ECS/ASW/AvgBIO/AUC-ROC/Cohen's-d formulas transcribed from its Methods section). Classic-metric citation years are working knowledge and should be link-verified before any external publication.
