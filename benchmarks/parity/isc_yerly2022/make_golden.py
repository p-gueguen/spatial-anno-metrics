"""Write the ISC parity golden input and score it with OUR inter_sample_consistency.

Reads make_input.R's output (counts.mtx genes x cells, genes.txt, obs.csv), writes
tests/goldens/isc_yerly2022/{counts.npz,obs.csv}, reloads them and checks the round trip, then runs
spatial_anno_metrics.inter_sample_consistency on exactly those files and writes ours.json.

Usage: PYTHONPATH=<repo>/src python make_golden.py <make_input_out_dir> <golden_dir> <ours.json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from spatial_anno_metrics import inter_sample_consistency
from spatial_anno_metrics.eval_metrics import _pseudobulk_profiles, _two_label_silhouette


def load_golden(gdir: Path):
    z = np.load(gdir / "counts.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return X, z["genes"].astype(str), pd.read_csv(gdir / "obs.csv", dtype=str)


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
    res = inter_sample_consistency(a, "annotation", "patient_bcc", counts_layer="counts")

    # ISCG in scTypeEval's default space (README recipe: Pseudobulk:Cosine with reduction = TRUE, i.e. cosine
    # on the first 30 PCs of the centred, unscaled pseudobulk matrix). Exact SVD here; scTypeEval uses irlba.
    P, lab, _, _ = _pseudobulk_profiles(a, "annotation", "patient_bcc", min_cells=10, min_samples=5,
                                        counts_layer="counts")
    U, S, _ = np.linalg.svd(P - P.mean(0), full_matrices=False)
    E = U[:, :min(30, P.shape[0] - 1)] * S[:min(30, P.shape[0] - 1)]
    En = E / np.linalg.norm(E, axis=1, keepdims=True)
    g_pca = _two_label_silhouette(1.0 - En @ En.T, lab)
    res["per_type_iscg_pca30"] = {t: float(g_pca[lab == t].mean()) for t in np.unique(lab)}

    types = sorted(res["per_type_iscg"])
    rows = [(t, res["per_type_iscg"][t], res["per_type_iscg_pca30"][t], res["per_type_iscl"][t],
             res["per_type_integrated"][t], res["regime"][t]) for t in types]
    print(f"{'type':14s} {'ISCG':>7s} {'ISCG_pc30':>9s} {'ISCL':>7s} {'G*L':>7s}  regime")
    for t, g, gp, lo, p, r in rows:
        print(f"{t:14s} {g:7.4f} {gp:9.4f} {lo:7.4f} {p:7.4f}  {r}")
    print(f"mean: ISCG {res['iscg']:.4f}  ISCL {res['iscl']:.4f}  integrated {res['integrated']:.4f}  "
          f"n_profiles {res['n_profiles']}  dropped {res['dropped_types']}")
    out.write_text(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:4]))
