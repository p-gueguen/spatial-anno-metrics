"""Write the ISC parity golden input, score it with OUR inter_sample_consistency, compare with scTypeEval.

Reads make_input.R's output (counts.mtx genes x cells, genes.txt, obs.csv), writes
tests/goldens/isc_yerly2022/{counts.npz,obs.csv}, reloads them and checks the round trip, then runs
spatial_anno_metrics.inter_sample_consistency on exactly those files for both label columns and writes
ours.json. If <golden_dir>/expected.json exists (run_sctypeeval.R), prints the SOURCE.md comparison tables.

Usage: PYTHONPATH=<repo>/src python make_golden.py <make_input_out_dir> <golden_dir> <ours.json>
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from spatial_anno_metrics import inter_sample_consistency
from spatial_anno_metrics.isc import _pseudobulk_profiles, _rcm_dissimilarity

COLUMNS = ("annotation", "annotation_oversplit")


def load_golden(gdir: Path):
    z = np.load(gdir / "counts.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return X, z["genes"].astype(str), pd.read_csv(gdir / "obs.csv", dtype=str)


def purge_label(x: str) -> str:
    """scTypeEval's label clean-up (space, _, +, - -> .; commas dropped): its labels and profile ids."""
    return re.sub(r"[ _+-]", ".", x).replace(",", "")


def score(a: ad.AnnData, col: str) -> dict:
    """Ours in the raw space (= sTE reduction = FALSE), plus ISCG in our default 30-PC space (= sTE default)."""
    res = inter_sample_consistency(a, col, "patient_bcc", layer="counts", reduction=None)
    res["per_type_iscg_pca30"] = inter_sample_consistency(a, col, "patient_bcc", layer="counts")["per_type_iscg"]
    P, lab, smp, _, _ = _pseudobulk_profiles(a, col, "patient_bcc", min_cells=10, min_samples=5, layer="counts")
    res["profile_ids"] = [f"{s}_{purge_label(t)}" for t, s in zip(lab, smp)]  # scTypeEval's row labels
    res["rcm"] = _rcm_dissimilarity(P, smp, lab)
    return res


def compare(col: str, res: dict, exp: dict) -> None:
    """Print the SOURCE.md tables: per label, means over labels, combined scores, RCM pair agreement."""
    e, n = exp["columns"][col], exp["names"]
    L, G, G0 = n["ISCL"], n["ISCG"], n["ISCG_no_pca"]
    def geo(g, lo):
        return np.sqrt(g * lo) if min(g, lo) >= 0 else np.nan  # README geometric mean, undefined below 0

    rows = []
    for t in sorted(res["per_type_iscg"]):
        p = e["per_type"][t]
        og, ol = res["per_type_iscg"][t], res["per_type_iscl"][t]
        rows.append([og, p[G0], res["per_type_iscg_pca30"][t], p[G], ol, p[L], og * ol, p[L] * p[G0], p[L] * p[G],
                     geo(og, ol), geo(p[G0], p[L]), np.nan if p["geometric_mean"] is None else p["geometric_mean"]])
    R = np.array(rows, dtype=float)
    print(f"\n## {col}\n| label | ISCG ours | ISCG sTE red=FALSE | ISCG ours 30-PC | ISCG sTE 30-PC | ISCL ours | "
          "ISCL sTE | G*L ours | G0*L sTE | G*L sTE | geo ours | geo0 sTE | geo sTE |")
    fmt = lambda v: "NA" if np.isnan(v) else f"{v:.6f}"  # noqa: E731
    for t, r in zip(sorted(res["per_type_iscg"]), R):
        print(f"| {t} | " + " | ".join(map(fmt, r)) + " |")
    print("| mean over labels | " + " | ".join(fmt(np.nanmean(c)) for c in R.T) + " |")
    print(f"ours headline (label means): iscg {res['iscg']:.6f} iscl {res['iscl']:.6f}; profile-weighted: iscg "
          f"{res['iscg_profile_mean']:.6f} iscl {res['iscl_profile_mean']:.6f}; "
          f"integrated {res['integrated']:.6f}; product of label means: ours "
          f"{R[:, 0].mean() * R[:, 4].mean():.6f}, sTE 30-PC {R[:, 3].mean() * R[:, 5].mean():.6f}, "
          f"sTE red=FALSE {R[:, 1].mean() * R[:, 5].mean():.6f}; geo NA labels: ours {int(np.isnan(R[:, 9]).sum())}, "
          f"sTE {int(np.isnan(R[:, 11]).sum())}")
    print(f"max |ISCG ours - sTE red=FALSE| {np.abs(R[:, 0] - R[:, 1]).max():.1e}; "
          f"max |ISCG ours 30-PC - sTE 30-PC| {np.abs(R[:, 2] - R[:, 3]).max():.1e}")
    D = e["dissimilarity"]["recip_classif:Match"]
    order = [D["labels"].index(i) for i in res["profile_ids"]]
    ref = np.array(D["matrix"], dtype=float)[np.ix_(order, order)]
    iu = np.triu_indices(len(order), 1)
    diff = np.flatnonzero(res["rcm"][iu] != ref[iu])
    print(f"RCM pairs differing ours vs SingleR: {len(diff)} of {len(iu[0])}; 0.5 entries in sTE: "
          f"{int((np.array(D['matrix']) == 0.5).sum())}")
    for k in diff:
        i, j = iu[0][k], iu[1][k]
        print(f"   {res['profile_ids'][i]} ~ {res['profile_ids'][j]}: ours {res['rcm'][i, j]:.0f}, sTE {ref[i, j]:.0f}")


def main(src: Path, gdir: Path, out: Path) -> None:
    X = scipy.io.mmread(src / "counts.mtx").T.tocsr()  # cells x genes
    X.sort_indices()
    genes = np.array((src / "genes.txt").read_text().split())
    obs = pd.read_csv(src / "obs.csv", dtype=str)
    assert X.shape == (len(obs), len(genes)) and len(set(genes)) == len(genes)
    assert np.array_equal(X.data, np.round(X.data)) and X.data.max() <= np.iinfo(np.uint16).max

    gdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(gdir / "counts.npz", data=X.data.astype(np.uint16),
                        indices=X.indices.astype(np.uint16), indptr=X.indptr.astype(np.int32),
                        shape=np.array(X.shape, dtype=np.int64), genes=genes.astype("U"))
    obs.to_csv(gdir / "obs.csv", index=False)

    Xg, genes_g, obs_g = load_golden(gdir)
    assert (Xg != X).nnz == 0 and np.array_equal(genes_g, genes) and obs_g.equals(obs)
    Y = Xg.copy()
    Y.data[0] += 1
    assert (Y != X).nnz == 1  # positive control: the equality check sees a single changed count
    for f in ("counts.npz", "obs.csv"):
        size = (gdir / f).stat().st_size
        print(f"{f}: {size / 1e6:.2f} MB")
        assert size < 5e6

    a = ad.AnnData(X=Xg.astype(np.float64), obs=obs_g.set_index("cell_id"))
    a.var_names = genes_g
    a.layers["counts"] = a.X
    exp = json.loads((gdir / "expected.json").read_text()) if (gdir / "expected.json").exists() else None
    ours = {}
    for col in COLUMNS:
        res = score(a, col)
        print(f"\n{col}: ISCG {res['iscg']:.4f}  ISCL {res['iscl']:.4f}  integrated {res['integrated']:.4f}  "
              f"n_profiles {res['n_profiles']}  dropped {res['dropped_types']}  regime {res['regime']}")
        if exp:
            compare(col, res, exp)
        del res["rcm"]
        ours[col] = res
    out.write_text(json.dumps(ours, indent=1, sort_keys=True))


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:4]))
