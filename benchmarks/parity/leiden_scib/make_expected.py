"""Parity values for the Leiden / k-means NMI and ARI: scib-metrics ``nmi_ari_cluster_labels_*`` vs ours.

Prints the scib-metrics values that ``tests/test_internal.py`` (``SCIB_LEIDEN`` / ``SCIB_KMEANS``) pins, and
how far our UMAP-connectivity port is from umap-learn's ``fuzzy_simplicial_set``. Not run in CI. Needs an
environment with scib-metrics 0.5.6 (BSD-3; it brings umap-learn, igraph and JAX), scikit-learn and anndata;
our package is imported from ``src/``, so it need not be installed there:

    <env>/bin/python benchmarks/parity/leiden_scib/make_expected.py

Both sides get the same exact scikit-learn kNN (self included), so only the clustering code is compared;
scib-metrics' own benchmark would use approximate pynndescent neighbours. The fixtures are the test suite's
``blobs`` (tests/conftest.py), rebuilt here with the same generator.
"""
from __future__ import annotations

import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
FIXTURES = {"overlap": {"sep": 1.5}, "overlap6": {"n_per": 60, "k": 6, "sep": 1.2, "seed": 3},
            "nmi_not_ari": {"n_per": 60, "k": 6, "sep": 1.0, "seed": 1}}          # best NMI and best ARI differ


def blobs(n_per=80, k=4, dim=10, sep=8.0, seed=0):
    """tests/conftest.py ``_blobs`` without AnnData: (float32 embedding, labels)."""
    rng = np.random.default_rng(seed)
    centers = rng.normal(0, sep, size=(k, dim))
    X = np.vstack([rng.normal(centers[t], 1.0, size=(n_per, dim)) for t in range(k)]).astype("float32")
    return X, np.repeat([f"T{t}" for t in range(k)], n_per)


def main() -> None:
    sys.path.insert(0, str(REPO / "src"))
    import anndata
    import pandas as pd
    from scib_metrics.metrics import nmi_ari_cluster_labels_kmeans
    from scib_metrics.metrics._nmi_ari import _compute_nmi_ari_cluster_labels
    from scib_metrics.nearest_neighbors import NeighborsResults
    from sklearn.neighbors import NearestNeighbors

    import spatial_anno_metrics as sam
    from spatial_anno_metrics.internal import _umap_connectivities

    print(f"scib-metrics {version('scib-metrics')}, umap-learn {version('umap-learn')}, igraph "
          f"{version('igraph')}, scikit-learn {version('scikit-learn')}, numpy {np.__version__}")
    for name, kw in FIXTURES.items():
        X, y = blobs(**kw)
        D, Idx = NearestNeighbors(n_neighbors=15).fit(X).kneighbors(X)
        nr = NeighborsResults(indices=Idx, distances=D)
        theirs = nr.knn_graph_connectivities.tocsr()
        ours = _umap_connectivities(D, Idx)
        diff = abs(theirs - ours)
        print(f"{name}: connectivities nnz {theirs.nnz} vs {ours.nnz}, max |diff| {diff.max():.2e}")

        sweep = [_compute_nmi_ari_cluster_labels(nr.knn_graph_connectivities, y, r, seed=42)
                 for r in [2 * x / 10 for x in range(1, 11)]]
        best = int(np.argmax([s[0] for s in sweep]))                   # as nmi_ari_cluster_labels_leiden
        km = nmi_ari_cluster_labels_kmeans(X, y)
        a = anndata.AnnData(X=X)
        a.obsm["X_pca"] = X
        a.obs["cell_type"] = pd.Categorical(y)
        mine = sam.leiden_label_agreement(a, seed=42)
        print(f"  scib leiden: best r={2 * (best + 1) / 10}, nmi={float(sweep[best][0])!r}, "
              f"ari={float(sweep[best][1])!r}")
        print(f"  ours leiden: best r={mine['resolution']}, nmi={mine['nmi']!r}, ari={mine['ari']!r}")
        print(f"  scib sweep nmi: {[float(s[0]) for s in sweep]}")
        print(f"  scib sweep ari: {[float(s[1]) for s in sweep]}")
        print(f"  ours sweep nmi: {[s['nmi'] for s in mine['sweep']]}")
        print(f"  scib kmeans: nmi={float(km['nmi'])!r}, ari={float(km['ari'])!r}; ours: "
              f"{ {k: sam.kmeans_label_agreement(a)[k] for k in ('nmi', 'ari')} }")


if __name__ == "__main__":
    main()
