"""AQI: the single Annotation Quality Index headline, with its advisory lineage-coverage helpers."""
from __future__ import annotations

import numpy as np

from .fidelity import marker_program_fidelity
from .panel import panel_resolvability

__all__ = ["annotation_quality_index"]


# --------------------------------------------------------------------------- #
# Advisory coverage helpers + the AQI headline
# --------------------------------------------------------------------------- #
def _lineage_coverage(adata, labels, marker_sets, abstained, *,
                      min_frac: float = 0.03, min_cells: int = 50, margin: float = 0.25):
    """Advisory: which marker-set lineages does the section's EXPRESSION strongly support that the
    annotation never assigned? Catches a COLLAPSED lineage - the AQI index cannot see one (purity +
    marker-fidelity score only the labels present; adequacy scores the reference), so a whole missing
    lineage costs the score nothing. Confirmed on human colon (Stromal + T/NK dropped, yet AQI 0.73).

    Reference-free: z-score each on-panel marker gene, take the per-cell marker-argmax lineage over its
    CONFIDENT calls (top-minus-second >= ``margin``), and flag any lineage that is the argmax for a real
    fraction of cells (``>= min_cells`` or ``>= min_frac``) yet is absent from the annotation labels.
    Returns ``{expected, annotated_on_axis, missing, missing_frac}`` or ``None`` (<2 lineages on panel).
    Only meaningful when ``marker_sets`` spans the tissue's lineages (not just the present labels)."""
    import scipy.sparse as sp

    var = {str(g): i for i, g in enumerate(map(str, adata.var_names))}
    lin_idx = {L: [var[g] for g in map(str, genes) if g in var]
               for L, genes in (marker_sets or {}).items()}
    lin_idx = {L: ii for L, ii in lin_idx.items() if ii}
    if len(lin_idx) < 2:
        return None
    allidx = sorted({i for ii in lin_idx.values() for i in ii})
    Xs = adata.X
    sub = (Xs[:, allidx].toarray() if sp.issparse(Xs) else np.asarray(Xs)[:, allidx]).astype(float)
    mu, sd = sub.mean(0), sub.std(0)
    sd[sd == 0] = 1.0
    Z = (sub - mu) / sd
    pos = {g: k for k, g in enumerate(allidx)}
    Ls = list(lin_idx)
    S = np.vstack([Z[:, [pos[i] for i in lin_idx[L]]].mean(1) for L in Ls]).T   # cells x lineages
    order = np.argsort(-S, axis=1)
    top = S[np.arange(len(S)), order[:, 0]]
    second = S[np.arange(len(S)), order[:, 1]] if S.shape[1] > 1 else np.full(len(S), -np.inf)
    conf = (top - second) >= margin
    arg = np.array([Ls[i] for i in order[:, 0]])
    n = max(1, int(conf.sum()))
    counts = {L: int(((arg == L) & conf).sum()) for L in Ls}
    expected = [L for L in Ls if counts[L] >= min_cells or counts[L] / n >= min_frac]
    annotated = {str(x) for x in np.unique(labels) if str(x) not in abstained}
    missing = [L for L in expected if L not in annotated]
    missing_frac = float(np.mean(conf & np.isin(arg, missing))) if missing else 0.0
    return {"expected": sorted(expected), "annotated_on_axis": sorted(annotated & set(Ls)),
            "missing": sorted(missing), "missing_frac": round(missing_frac, 3)}


def _reference_coverage(reference, ref_key, labels, abstained, *, tau: float = 0.03):
    """Advisory (reference-anchored): which lineages ABUNDANT in the reference (>= ``tau`` of ref cells)
    have NO matching annotation label? This covers the marker-based coverage blind spot on OFF-PANEL-marker
    panels - the UM liver met (melanoma/hepatocyte markers off the immuno-oncology panel, but present in the
    reference): the marker flag can't see them, this one can. Name-normalized loose match (a label using a
    different WORD for the same lineage can false-flag), so it is a "verify these" advisory, NEVER scored.
    Returns ``{reference_abundant, missing_vs_reference, n_ref_types_abundant}`` or ``None``."""
    import re as _re

    obs = getattr(reference, "obs", None)
    if reference is None or obs is None or ref_key not in obs.columns:
        return None

    def _norm(x):
        return _re.sub(r"[^a-z0-9]", "", str(x).lower())

    vc = obs[ref_key].astype(str).value_counts(normalize=True)
    ref_abundant = [str(t) for t, f in vc.items() if float(f) >= tau]
    lab_norm = {_norm(x) for x in {str(v) for v in labels} - set(map(str, abstained)) if x}
    lab_norm.discard("")

    def _present(reft):
        r = _norm(reft)
        return bool(r) and any(r == la or (len(r) >= 4 and (r in la or la in r)) for la in lab_norm)

    missing = sorted([t for t in ref_abundant if not _present(t)])
    return {"reference_abundant": sorted(ref_abundant), "missing_vs_reference": missing,
            "n_ref_types_abundant": len(ref_abundant)}


def annotation_quality_index(adata, label_key: str, marker_sets: dict[str, list[str]] | None,
                             internal_validity_result: dict | None = None, *,
                             reference=None, ref_key: str = "cell_type", panel_genes=None,
                             median_depth: float | None = None, method_label_cols=None,
                             abstention_labels=None, provenance: dict | None = None,
                             cohort_result: dict | None = None,
                             min_type_n: int = 50, p_softmin: float = -4.0) -> dict:
    """AQI - one Annotation Quality Index in [0, 1].

    ``AQI = w_coh * min(A, soft_min_{p=-4}({C, M}))`` - the soft-min of the ACHIEVED-quality terms,
    capped by the panel-adequacy ceiling: **C** contamination/purity (macro median PMP x retention;
    panel-size invariant), **M** marker-program fidelity (macro one-vs-rest AUC), capped by **A** panel x
    depth adequacy (base-rate-normalized depth-matched macro F1 - a resolvable CEILING: you cannot type
    better than the panel resolves, but a rich panel does not by itself mean the labels are right).
    Coherence **H** is a bounded <=15% multiplier (self-consistency can never manufacture score);
    completeness is a side flag (honest abstention is not low quality). Every component is macro over
    predicted types with ``n >= min_type_n`` and any UNDEFINED term is DROPPED (never set to 0).

    Cross-method AGREEMENT **G** (>=3 voters) is NOT a cross-section term: its absolute level tracks
    REFERENCE quality, not accuracy, and does not transfer across sections (2026-07 validation). It is
    reported under ``abstention`` as a WITHIN-section per-cell signal (which cells to grey out), and
    ``regime`` (``with_agreement`` vs ``index_only``) says whether that signal is available.

    It is an INDEX (section-comparable, monotone-in-quality), **not** ``P(correct)`` - the no-transfer
    finding forbids a universal reference-free accuracy curve. Ground-truth accuracy is a caller concern
    (:func:`external_scores`), not this index.

    Returns ``{aqi, aqi_core, regime, components:{A,C,G,M,H,completeness}, abstention:{signal,n_voters,
    available}, coverage:{expected,annotated_on_axis,missing,missing_frac,reference_abundant,
    missing_vs_reference,n_ref_types_abundant}, cohort:{consistency,per_type,n_profiles}, w_coh, active_set,
    argmin, n_ge50_types, flags, provenance}``. ``coverage`` is an ADVISORY warning (never folded into the
    score): lineages the annotation dropped entirely - the index cannot see a collapsed lineage (C/M score
    only present labels; A scores the reference), confirmed on colon. It is MARKER-anchored (``missing``,
    on-panel markers) AND, when a reference is threaded, REFERENCE-anchored (``missing_vs_reference``, ref-
    abundant lineages with no matching label - covers the off-panel-marker blind spot, UM liver met).
    ``cohort`` reports inter-sample label consistency when ``cohort_result`` is passed (also not scored - its
    level does not transfer, like G). Pure/guarded; QC-safe.
    """
    from .purity import pmp as _pmp

    labels = adata.obs[label_key].astype(str).to_numpy()
    abst = {str(x) for x in (abstention_labels or [])}
    vc = adata.obs[label_key].astype(str).value_counts()
    big = [str(t) for t, n in vc.items() if int(n) >= min_type_n and str(t) not in abst]
    big_mask = np.isin(labels, big)

    def _c01(x):
        return float(min(max(float(x), 0.0), 1.0))

    # A - panel x depth adequacy: base-rate-null-normalized macro F1 at the section's own depth.
    A = None
    frac_resolvable = None
    if reference is not None and panel_genes is not None:
        try:
            pr = panel_resolvability(reference, ref_key, panel_genes,
                                     target_depth=float(median_depth) if median_depth else 50.0)
            mf, nt = pr.get("macro_f1"), pr.get("n_types")
            frac_resolvable = pr.get("frac_resolvable")
            if mf is not None and nt and int(nt) >= 2:
                null = 1.0 / int(nt)
                A = _c01((float(mf) - null) / (1.0 - null))
        except Exception:  # pragma: no cover - defensive
            A = None

    # C - contamination/purity: min(macro median PMP, retention). Deleting transcripts -> NaN PMP ->
    # retention falls, so it cannot be gamed by dropping counts. Panel-size invariant (lineage-marker denom).
    C = None
    if marker_sets:
        try:
            _pmp(adata, marker_sets, label_key)                     # writes obs['pmp']
            pv = np.asarray(adata.obs["pmp"], dtype=float)
            meds, defined = [], 0
            for t in big:
                col = pv[labels == t]
                ok = col[~np.isnan(col)]
                if ok.size:
                    meds.append(float(np.median(ok)))
                defined += int((~np.isnan(col)).sum())
            if meds and big_mask.any():
                C = _c01(min(float(np.mean(meds)), defined / int(big_mask.sum())))
        except Exception:  # pragma: no cover - defensive
            C = None

    # G - cross-method agreement (the only reference-free correctness term); needs >=3 diverse voters.
    G = None
    cols = [c for c in (method_label_cols or []) if c in adata.obs.columns]
    if len(cols) >= 3:
        try:
            votes = np.vstack([adata.obs[c].astype(str).to_numpy() for c in cols])   # voters x cells
            agree = (votes == labels[None, :]).mean(0)              # frac of voters backing the label
            per = [float(agree[labels == t].mean()) for t in big if (labels == t).any()]
            if per:
                G = _c01(float(np.mean(per)))
        except Exception:  # pragma: no cover - defensive
            G = None

    # M - marker-program fidelity: macro one-vs-rest AUC over n>=50 types, rescaled 2*(auc-0.5).
    M = None
    if marker_sets:
        try:
            mpf = marker_program_fidelity(adata, label_key, marker_sets)
            aucs = [v["auc"] for t, v in mpf.get("per_type", {}).items()
                    if t in big and v.get("auc") is not None and int(v.get("n", 0)) >= min_type_n]
            if aucs:
                M = _c01(2.0 * (float(np.mean(aucs)) - 0.5))
        except Exception:  # pragma: no cover - defensive
            M = None

    # H - coherence: bounded <=15% multiplier only (coherence != correctness).
    H = internal_validity_result.get("integrated") if isinstance(internal_validity_result, dict) else None
    w_coh = float(np.clip(0.85 + 0.15 * H, 0.85, 1.0)) if isinstance(H, (int, float)) else 1.0

    # completeness - side usability flag, NOT a term (honest abstention on a shallow panel is correct).
    typed_frac = float(np.mean(~np.isin(labels, list(abst)))) if abst else 1.0
    completeness = _c01(typed_frac / frac_resolvable) if frac_resolvable else _c01(typed_frac)

    # advisory lineage-COVERAGE: which marker-supported lineages did the annotation drop ENTIRELY? The
    # index cannot see a collapsed lineage (C/M score only present labels; A scores the reference), so
    # this is a separate WARNING, never folded into the score. Confirmed on human colon (Stromal+T/NK
    # dropped yet AQI 0.73). Only fires when marker_sets spans the tissue's lineages (else missing=[]).
    coverage = None
    if marker_sets:
        try:
            coverage = _lineage_coverage(adata, labels, marker_sets, abst)
        except Exception:  # pragma: no cover - defensive
            coverage = None
    # reference-anchored coverage (cov_ref): catches a dropped compartment even when its MARKERS are
    # off-panel (the marker flag's blind spot; UM liver met). Advisory, merged into `coverage`.
    if reference is not None:
        try:
            rcov = _reference_coverage(reference, ref_key, labels, abst)
            if rcov is not None:
                coverage = {**(coverage or {}), **rcov}
        except Exception:  # pragma: no cover - defensive
            pass

    # no-dict disclosure: fraction of big cells whose label has NO on-panel marker set (undeclared - declared
    # abstention labels are already dropped from `big`). Their PMP is NaN for lack of a marker DICTIONARY (a
    # COVERAGE gap, fixable by supplying reference-derived markers), NOT for low signal - the two are conflated
    # inside `retention`. Reported so a coverage-limited C is legible instead of an opaque under-read: MERSCOPE
    # E14 curated dicts covered ~10/105 MOCA labels -> retention ~0.06 -> AQI 0.08 on a GOOD annotation; the fix
    # is to LOAD A REFERENCE (which then grounds markers per label). NOT scored - retention already caps the
    # index (a bare "forgive no-dict from retention" change over-reads; verified 2026-07-12); this only explains WHY.
    no_dict_frac = None
    if marker_sets is not None and big_mask.any():
        _panel = set(map(str, adata.var_names))
        nod = sum(int((labels == t).sum()) for t in big
                  if not any(g in _panel for g in (marker_sets or {}).get(t, [])))
        no_dict_frac = round(nod / int(big_mask.sum()), 3)

    # combine: soft-min over the ACHIEVED-quality terms {C, M}, CAPPED by the panel x depth ceiling A,
    # then the bounded coherence multiplier. This is the 2026-07 VALIDATION correction (validate_aqi.py,
    # 3 GT sections, common 6-lineage axis). Measured section-level Spearman(component, true balanced
    # accuracy): M = +1.0 (marker fidelity orders sections perfectly), C = +0.5, but A = -1.0 and
    # G = +0.5. A ANTI-orders accuracy because it is a resolvable-CEILING (high on any deep panel
    # regardless of the labels actually assigned) - so it must CAP, never drive. G's absolute level
    # tracks REFERENCE quality, not accuracy, and does not transfer across sections - so it is a WITHIN-
    # section per-cell abstention signal (reported under `abstention`), never a cross-section term. The
    # old soft_min(A,C,G,M) mixed a ceiling and a non-transferring term into the bottleneck (rho +0.5,
    # wrong section order); min(A, soft_min(C,M)) restores rho +1.0.
    core = {"C": C, "M": M}
    active = {k: v for k, v in core.items() if v is not None}
    names = {"C": "contamination", "M": "marker fidelity", "A": "panel/depth ceiling"}
    aqi_core = argmin = aqi = None
    if active:
        clipped = {k: min(max(v, 0.01), 1.0) for k, v in active.items()}
        aqi_core = _c01(float(np.mean([c ** p_softmin for c in clipped.values()])) ** (1.0 / p_softmin))
        capped = aqi_core if A is None else min(aqi_core, A)          # adequacy is a CEILING, not a driver
        argmin = ("panel/depth ceiling" if (A is not None and A < aqi_core)
                  else names[min(clipped, key=clipped.get)])
        aqi = _c01(w_coh * capped)

    return {
        "aqi": aqi, "aqi_core": aqi_core,
        "regime": "with_agreement" if G is not None else "index_only",
        "components": {"A": A, "C": C, "G": G, "M": M, "H": H, "completeness": completeness},
        # cross-method agreement is a WITHIN-section per-cell abstention signal, not a cross-section term.
        "abstention": {"signal": G, "n_voters": len(cols), "available": G is not None},
        # advisory: lineages the annotation dropped entirely (index-invisible); marker- AND reference-anchored.
        "coverage": coverage,
        # disclosure: fraction of scored cells whose label has no on-panel marker dict (coverage gap, NOT scored).
        "no_dict_frac": no_dict_frac,
        # cohort signal (across a sample_key cohort): inter-sample label consistency; REPORTED, not scored -
        # its level does not transfer (like G), and CBCL showed a per-type spread (Keratinocyte 0.51 vs
        # Myeloid -0.11). Present only when the caller computed it (a sample_key exists).
        "cohort": ({"consistency": cohort_result.get("consistency"),
                    "per_type": cohort_result.get("per_type"),
                    "n_profiles": cohort_result.get("n_profiles")}
                   if isinstance(cohort_result, dict) and "consistency" in cohort_result else None),
        "w_coh": w_coh, "active_set": list(active.keys()) + (["A"] if A is not None else []),
        "argmin": argmin, "n_ge50_types": len(big),
        "flags": {"adequacy_unknown": A is None, "low_support": len(big) < 2,
                  "abstention_available": G is not None,
                  "missing_lineages": bool(coverage and coverage.get("missing")),
                  "missing_vs_reference": bool(coverage and coverage.get("missing_vs_reference")),
                  "coverage_limited": bool(no_dict_frac is not None and no_dict_frac > 0.5)},
        "provenance": dict(provenance or {}, n_voters=len(cols)),
    }
