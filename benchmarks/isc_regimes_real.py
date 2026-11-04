"""Known-answer benchmark of the ISC regime calls on real multi-sample scRNA-seq (roadmap 8d).

    Rscript benchmarks/isc_regimes_real_export.R <rds> <ident> <sample> <data>/<name> [excluded] [extra cols]
    PYTHONPATH=src python benchmarks/isc_regimes_real.py --data <data> [--seeds 3] [--out <dir>]
    # optional, scTypeEval on the same perturbations (Yerly 2022, seed 0):
    PYTHONPATH=src python benchmarks/isc_regimes_real.py --data <data> --sctypeeval-input <dir>
    Rscript benchmarks/isc_regimes_real_sctypeeval.R <dir>
    PYTHONPATH=src python benchmarks/isc_regimes_real.py --data <data> --sctypeeval-output <dir>

Data: the ISC authors' Zenodo record 10.5281/zenodo.18921437 (CC BY 4.0), one ``.rds`` per study, exported
(2000 HVGs, scTypeEval README recipe) by ``isc_regimes_real_export.R``. The labels are the authors'
(``annotation`` / ``cell.type``), the sample is the patient, as in the authors' benchmark config.

Every scorable label of a dataset is scored after one perturbation with a known answer:

- ``none``: the authors' labels, expected ``consistent`` (an assumption: the paper finds ~64 % of atlas labels
  consistent, so a few "wrong" calls here may be right);
- ``merge:A+B``: two distinct types pooled under one label, expected ``split`` (under-partitioned);
- ``split:T``: a random half of T in every sample relabelled (the paper's over-partitioning task), both halves
  expected ``merge``;
- ``contaminate:T@f``: in every sample a fraction f of T's cells swapped with as many random other cells,
  T expected ``remove`` (at f = 1 the label is random cells; at f = 0.5 half of it is);
- ``junk@0.03``: a new label of 3 % random cells of every sample (multiplets / low quality), ``remove``;
- ``noise@f``: the paper's label-noise task (a fraction f of each sample's cells relabelled from that
  sample's label distribution), expected ``degraded``: no per-label call is right, only ``merge`` / ``split``
  count as wrong;
- ``batch@sd``: every gene of a sample multiplied by its own log-normal factor, labels untouched, expected
  ``consistent``.

Labels not targeted keep ``consistent`` as their expected call while >= 90 % of their cells are still
their own (else ``degraded``: contaminating a dominant label spreads it over every other label). Targets are drawn at random per dataset
(``numpy`` seed 0; merge pairs are drawn independently, so a pair can repeat), so the held-out dataset is
not hand-picked. Each call is graded ``correct``, ``missed``
(``consistent`` where a problem was planted), ``inconclusive`` or ``wrong`` (a different definite call), for
the paper's quadrant (``quadrant``: fixed 0.74 / 0.45 guides) and for the calibrated ``regime``. Writes
``rows.csv`` (one row per label and run) and ``summary.md``.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import spatial_anno_metrics as sam  # noqa: E402

OUT = ROOT / "benchmarks" / "isc_regimes_real" / "results"
#: name -> (export directory, label column, role). Thresholds were set on the development sets; the held-out
#: set was exported and scored only after they were fixed.
DATASETS = {
    "yerly2022": ("yerly_full", "label", "development"),
    "ganier2024": ("ganier_full", "label", "development"),
    "ganier2024_scgate": ("ganier_full", "scGate_multi", "development"),
    "joanito2022": ("joanito_full", "label", "held-out"),
}
N_MERGE, N_SPLIT, N_CONT = 4, 4, 3
GRADES = ["correct", "missed", "inconclusive", "wrong"]


def load(data: Path, name: str) -> ad.AnnData:
    """Read one export of ``isc_regimes_real_export.R`` (CSC genes x cells of the HVGs) as cells x genes."""
    d = data / name
    g, c = (int(x) for x in (d / "shape.txt").read_text().split())
    X = sp.csc_matrix((np.fromfile(d / "counts_x.bin", dtype=np.float64),
                       np.fromfile(d / "counts_i.bin", dtype=np.int32),
                       np.fromfile(d / "counts_p.bin", dtype=np.int32)), shape=(g, c)).T.tocsr().astype(np.float32)
    a = ad.AnnData(X=X, obs=pd.read_csv(d / "obs.csv", index_col=0, dtype=str),
                   var=pd.DataFrame(index=(d / "genes.txt").read_text().split()))
    a.layers["counts"] = X
    return a


# --------------------------------------------------------------------------- #
# perturbations: each returns (labels, {label: expected call})
# --------------------------------------------------------------------------- #
def _by_sample(samples):
    return [np.flatnonzero(samples == s) for s in np.unique(samples)]


def merge(lab, smp, a, b, rng=None):
    """Pool two labels into one: under-partitioning."""
    lab = lab.copy()
    lab[np.isin(lab, [a, b])] = f"{a}+{b}"
    return lab, {f"{a}+{b}": "split"}


def split(lab, smp, t, rng):
    """Relabel a random half of ``t`` in every sample (the paper's SplitCelltype at rate 0.5)."""
    lab = lab.copy()
    for ids in _by_sample(smp):
        tid = ids[lab[ids] == t]
        lab[rng.choice(tid, len(tid) // 2, replace=False)] = f"{t}.Nspl2"
    return lab, {t: "merge", f"{t}.Nspl2": "merge"}


def contaminate(lab, smp, t, f, rng):
    """Swap a fraction ``f`` of ``t``'s cells with as many random other cells, in every sample."""
    lab0, lab = lab, lab.copy()
    for ids in _by_sample(smp):
        tid, oid = ids[lab0[ids] == t], ids[lab0[ids] != t]
        n = min(int(round(f * len(tid))), len(oid))
        x, y = rng.choice(tid, n, replace=False), rng.choice(oid, n, replace=False)
        lab[x], lab[y] = lab0[y], t
    return lab, {t: "remove"}


def junk(lab, smp, f, rng):
    """A new label 'Junk' of a fraction ``f`` of every sample's cells, drawn at random."""
    lab = lab.copy()
    for ids in _by_sample(smp):
        lab[rng.choice(ids, int(round(f * len(ids))), replace=False)] = "Junk"
    return lab, {"Junk": "remove"}


def noise(lab, smp, f, rng):
    """The paper's label-noise task: a fraction ``f`` of each sample relabelled from its label distribution."""
    lab0, lab = lab, lab.copy()
    for ids in _by_sample(smp):
        pick = rng.choice(ids, int(round(f * len(ids))), replace=False)
        lab[pick] = rng.choice(lab0[ids], len(pick), replace=True)
    return lab, "degraded"


def batch_counts(X, smp, sd, rng):
    """Multiply every gene of a sample by its own log-normal factor (log-sd ``sd``), Poisson resampled."""
    u, sidx = np.unique(smp, return_inverse=True)
    F = rng.lognormal(0, sd, (len(u), X.shape[1]))
    c = sp.coo_matrix(X)
    out = sp.csr_matrix((rng.poisson(c.data * F[sidx[c.row], c.col]).astype(np.float32), (c.row, c.col)),
                        shape=X.shape)
    out.eliminate_zeros()
    return out


def scored_labels(a, col):
    """Labels with >= 5 samples of >= 10 cells, most abundant first."""
    n = pd.crosstab(a.obs[col], a.obs["sample"])
    ok = n[(n >= 10).sum(1) >= 5]
    return list(ok.sum(1).sort_values(ascending=False).index)


def jobs(lab, smp, labels, seed):
    """All perturbations of one seed; targets drawn once per dataset (seed-independent)."""
    pick = np.random.default_rng(0)
    pairs = [tuple(pick.choice(labels, 2, replace=False)) for _ in range(N_MERGE)]
    spl = list(pick.choice(labels, min(N_SPLIT, len(labels)), replace=False))
    cont = list(pick.choice(labels, min(N_CONT, len(labels)), replace=False))
    rng = np.random.default_rng(100 + seed)
    out = [("none", lab, {})] if seed == 0 else []
    if seed == 0:
        out += [(f"merge:{x}+{y}", *merge(lab, smp, x, y)) for x, y in pairs]
    out += [(f"split:{t}", *split(lab, smp, t, rng)) for t in spl]
    out += [(f"contaminate:{t}@{f}", *contaminate(lab, smp, t, f, rng)) for t in cont for f in (0.5, 1.0)]
    out += [("junk@0.03", *junk(lab, smp, 0.03, rng))]
    out += [(f"noise@{f}", *noise(lab, smp, f, rng)) for f in (0.2, 0.4)]
    return out, rng


def score(a, lab, expected, meta, rows, orig):
    """Run inter_sample_consistency on ``lab`` and append one row per scored label.

    A label that is not a target keeps ``consistent`` as its expected call only while >= 90 % of its cells
    still carry it in ``orig``; below that (a large contaminated target spreads its cells over every other
    label) it is ``degraded``.
    """
    a.obs["_lab"] = lab
    r = sam.inter_sample_consistency(a, "_lab", "sample", layer="counts")
    for t in r.get("per_type_iscg", {}):
        exp = expected if isinstance(expected, str) else expected.get(t, "consistent")
        purity = float(np.mean(orig[lab == t] == t))
        if exp == "consistent" and purity < 0.9:
            exp = "degraded"
        ev = r["regime_evidence"][t]
        rows.append(dict(**meta, label=t, expected=exp, purity=purity, iscg=r["per_type_iscg"][t], iscl=r["per_type_iscl"][t],
                         quadrant=r["quadrant"][t], regime=r["regime"][t], iscg_calibrated=ev["iscg_calibrated"],
                         self_match=ev["self_match"], partner=ev["partner"], partner_rate=ev["partner_rate"],
                         iscg_floor=r["iscg_floor"], low_global_signal=r["low_global_signal"]))


def run(data: Path, seeds: int, only=None) -> pd.DataFrame:
    """Score every dataset under every perturbation."""
    rows = []
    for name, (d, col, role) in DATASETS.items():
        if only and name not in only:
            continue
        a = load(data, d)
        a = a[a.obs[col].notna().to_numpy()].copy()
        lab = a.obs[col].astype(str).to_numpy()
        smp = a.obs["sample"].astype(str).to_numpy()
        labels = scored_labels(a, col)
        for seed in range(seeds):
            js, rng = jobs(lab, smp, labels, seed)
            for pert, lab2, exp in js:
                score(a, lab2, exp, dict(dataset=name, role=role, pert=pert, seed=seed), rows, lab)
            for sd in (0.5, 1.0):
                b = a.copy()
                b.layers["counts"] = batch_counts(a.layers["counts"], smp, sd, rng)
                score(b, lab, {}, dict(dataset=name, role=role, pert=f"batch@{sd}", seed=seed), rows, lab)
        print(name, len(rows), flush=True)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# grading and summary
# --------------------------------------------------------------------------- #
def grade(expected: str, call: str) -> str:
    """``correct`` / ``missed`` / ``inconclusive`` / ``wrong`` of one call; ``degraded`` only forbids refinements."""
    if expected == "degraded":
        return "wrong" if call in ("merge", "split") else "correct"
    if call == expected:
        return "correct"
    if call == "consistent":
        return "missed"
    return "inconclusive" if call == "inconclusive" else "wrong"


def kind(pert: str) -> str:
    """Perturbation family, with the level for contamination / noise / batch."""
    k = pert.split(":")[0].split("@")[0]
    if k in ("contaminate", "noise", "batch", "junk"):
        return f"{k}@{pert.split('@')[1]}"
    return k


def table(df: pd.DataFrame, rule: str) -> pd.DataFrame:
    """Rows: (perturbation family, expected call); columns: share of each grade and n."""
    g = df.assign(family=df.pert.map(kind), grade=[grade(e, c) for e, c in zip(df.expected, df[rule])])
    t = g.groupby(["family", "expected"]).grade.value_counts().unstack(fill_value=0).reindex(columns=GRADES,
                                                                                               fill_value=0)
    t.insert(0, "n", t.sum(1))
    return t


def per_regime(df: pd.DataFrame, rule: str) -> pd.DataFrame:
    """Accuracy per expected call (the confusion-matrix rows), over all perturbations."""
    g = df.assign(grade=[grade(e, c) for e, c in zip(df.expected, df[rule])])
    t = g.groupby("expected").grade.value_counts().unstack(fill_value=0).reindex(columns=GRADES, fill_value=0)
    t.insert(0, "n", t.sum(1))
    return t


def _md(t: pd.DataFrame) -> str:
    t = t.reset_index()
    head = "| " + " | ".join(map(str, t.columns)) + " |"
    return "\n".join([head, "|" + "---|" * len(t.columns)] + ["| " + " | ".join(map(str, r)) + " |"
                                                            for r in t.itertuples(index=False)])


def summary(df: pd.DataFrame, ste: pd.DataFrame | None = None) -> str:
    """Markdown: per-regime accuracy for both rules, per role, plus the confusion matrices."""
    out = ["# ISC regime calls on real data with a known answer", "",
           "Generated by `benchmarks/isc_regimes_real.py`; `quadrant` = the paper's fixed 0.74 / 0.45 guides, "
           "`regime` = the calibrated call. Grades: correct / missed (`consistent` where a problem was planted) / "
           "inconclusive / wrong (another definite call). `degraded` (label noise): only merge / split are wrong.",
           ""]
    for role in ("development", "held-out"):
        sub = df[df.role == role]
        if sub.empty:
            continue
        out += [f"## {role}: {', '.join(sorted(sub.dataset.unique()))}", ""]
        for rule in ("quadrant", "regime"):
            out += [f"### `{rule}`, per expected call", "", _md(per_regime(sub, rule)), ""]
        out += ["### `regime`, per perturbation", "", _md(table(sub, "regime")), ""]
        out += ["### confusion (`regime`): expected x called", "",
                _md(pd.crosstab(sub.expected, sub.regime)), ""]
        out += ["### confusion (`quadrant`): expected x called", "",
                _md(pd.crosstab(sub.expected, sub.quadrant)), ""]
    if ste is not None:
        out += ["## scTypeEval 1.0.0 on the same perturbations (Yerly 2022, seed 0)", "",
                "scTypeEval reports ISCG / ISCL only; its regime = the paper's quadrant on its scores.", "",
                _md(per_regime(ste, "quadrant")), "", _md(pd.crosstab(ste.expected, ste.quadrant)), ""]
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- #
# scTypeEval on the same perturbations
# --------------------------------------------------------------------------- #
def write_sctypeeval_input(data: Path, out: Path) -> None:
    """Write Yerly 2022 seed-0 labels (one column per perturbation) and the batch matrices for the R script."""
    out.mkdir(parents=True, exist_ok=True)
    d, col, _ = DATASETS["yerly2022"]
    a = load(data, d)
    lab, smp = a.obs[col].astype(str).to_numpy(), a.obs["sample"].astype(str).to_numpy()
    js, rng = jobs(lab, smp, scored_labels(a, col), 0)
    tab = pd.DataFrame({"cell_id": a.obs_names, "sample": smp})
    exp = {}
    for i, (pert, lab2, e) in enumerate(js):
        tab[f"p{i}"] = lab2
        exp[f"p{i}"] = (pert, e)
    for sd in (0.5, 1.0):
        X = batch_counts(a.layers["counts"], smp, sd, rng).T.tocsc()
        X.indices.astype(np.int32).tofile(out / f"batch{sd}_i.bin")
        X.indptr.astype(np.int32).tofile(out / f"batch{sd}_p.bin")
        X.data.astype(np.float64).tofile(out / f"batch{sd}_x.bin")
    for f in ("counts_i.bin", "counts_p.bin", "counts_x.bin", "shape.txt", "genes.txt"):
        (out / f).write_bytes((data / d / f).read_bytes())
    tab.to_csv(out / "labels.csv", index=False)
    pd.DataFrame([(k, p, e if isinstance(e, str) else ";".join(f"{a}={b}" for a, b in e.items()))
                  for k, (p, e) in exp.items()]).to_csv(out / "perturbations.tsv", sep="\t", header=False, index=False)


def read_sctypeeval_output(d: Path) -> pd.DataFrame:
    """Per-label scTypeEval ISCG / ISCL of each perturbation -> rows with the paper's quadrant."""
    s = pd.read_csv(d / "sctypeeval_scores.csv").dropna(subset=["iscg", "iscl"])   # NA = a dropped label
    pert = pd.read_csv(d / "perturbations.tsv", sep="\t", header=None, names=["col", "pert", "exp"],
                       keep_default_na=False)
    rows = []
    for col, p, e in pert.itertuples(index=False):
        exp = dict(x.split("=") for x in e.split(";") if x) if e != "degraded" else "degraded"
        for r in s[s.column == col].itertuples(index=False):
            ex = exp if isinstance(exp, str) else exp.get(r.label, "consistent")
            rows.append(dict(pert=p, label=r.label, expected=ex, iscg=r.iscg, iscl=r.iscl,
                             quadrant=sam.isc_regime(r.iscg, r.iscl)))
    for sd in ("0.5", "1.0"):
        for r in s[s.column == f"batch{sd}"].itertuples(index=False):
            rows.append(dict(pert=f"batch@{sd}", label=r.label, expected="consistent", iscg=r.iscg, iscl=r.iscl,
                             quadrant=sam.isc_regime(r.iscg, r.iscl)))
    return pd.DataFrame(rows)


def main(argv=None):
    """Run the benchmark (or the scTypeEval input / output steps) and write the results."""
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, required=True, help="folder with the isc_regimes_real_export.R outputs")
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--only", nargs="*", default=None)
    p.add_argument("--sctypeeval-input", type=Path, default=None)
    p.add_argument("--sctypeeval-output", type=Path, default=None)
    args = p.parse_args(argv)
    if args.sctypeeval_input:
        write_sctypeeval_input(args.data, args.sctypeeval_input)
        return
    args.out.mkdir(parents=True, exist_ok=True)
    rows_csv = args.out / "rows.csv"
    if args.sctypeeval_output:
        df = pd.read_csv(rows_csv, keep_default_na=False, na_values=[""])
        ste = read_sctypeeval_output(args.sctypeeval_output)
        ste.to_csv(args.out / "sctypeeval_rows.csv", index=False, float_format="%.6g")
    else:
        df = run(args.data, args.seeds, args.only)
        df.to_csv(rows_csv, index=False, float_format="%.6g")
        ste_csv = args.out / "sctypeeval_rows.csv"
        ste = pd.read_csv(ste_csv) if ste_csv.exists() else None
    (args.out / "summary.md").write_text(summary(df, ste))


if __name__ == "__main__":
    main()
