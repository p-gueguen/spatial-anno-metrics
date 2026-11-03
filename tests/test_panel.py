"""Tests for panel.py - panel resolvability at the panel's own depth."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _panel_ref(n=120, n_genes=40, seed=0):
    """Reference where A/B have distinct on-panel markers but C/D differ ONLY off-panel (last gene)."""
    rng = np.random.default_rng(seed)
    off = n_genes - 1
    blocks, labels = [], []
    for t, hi in [("A", 0), ("B", 1), ("C", 2), ("D", 2)]:      # C,D share on-panel marker g2
        M = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
        M[:, hi] += rng.poisson(20.0, size=n)
        if t == "D":
            M[:, off] += rng.poisson(20.0, size=n)              # D's only distinct gene is OFF-panel
        blocks.append(M)
        labels += [t] * n
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = labels
    return a


def test_panel_resolvability_flags_off_panel_confusable_pair():
    """A panel that omits a type-pair's only distinguishing gene cannot resolve that pair."""
    from spatial_anno_metrics import panel as em

    a = _panel_ref()
    panel = {f"g{i}" for i in range(39)}                         # 39 shared genes; excludes g39 (D's marker)
    res = em.panel_resolvability(a, "cell_type", panel, target_depth=40)
    assert res["status"] == "ok" and res["n_types"] == 4
    assert res["per_type"]["A"]["tier"] == "resolvable"         # private on-panel markers -> resolvable
    assert res["per_type"]["B"]["tier"] == "resolvable"
    cd = {res["per_type"]["C"]["tier"], res["per_type"]["D"]["tier"]}
    assert "not_resolvable" in cd                                # identical on the panel -> confusable
    assert (res["per_type"]["C"]["confused_with"] == "D"
            or res["per_type"]["D"]["confused_with"] == "C")
    assert res["frac_resolvable"] < 1.0


def test_panel_resolvability_insufficient_overlap():
    """A reference sharing < 25 genes with the panel returns an honest skip, not a bogus score."""
    from spatial_anno_metrics import panel as em

    a = _panel_ref()
    res = em.panel_resolvability(a, "cell_type", {"zzz1", "zzz2"}, target_depth=40)
    assert res["status"] == "skipped" and "0 genes shared" in res["reason"] and res["n_shared_genes"] == 0


def _log_with_counts(a):
    """The same cells with log1p(CP10k) in X and the raw counts in layers['counts']."""
    b = a.copy()
    b.layers["counts"] = a.X.copy()
    X = np.asarray(a.X, dtype=float)
    b.X = np.log1p(X / X.sum(1, keepdims=True) * 1e4).astype("float32")
    return b


def test_panel_resolvability_thins_the_count_layer_never_log_X():
    """Thinning is load-bearing: log X must never be scored at full depth. The counts come from ``layer``."""
    from spatial_anno_metrics import panel as em

    a, panel = _panel_ref(), {f"g{i}" for i in range(39)}
    raw = em.panel_resolvability(a, "cell_type", panel, target_depth=3)
    full = em.panel_resolvability(a, "cell_type", panel, target_depth=1e6)
    assert raw["macro_f1"] < full["macro_f1"] - 0.1               # 3 counts per cell blur the markers
    b = _log_with_counts(a)
    via_layer = em.panel_resolvability(b, "cell_type", panel, target_depth=3)
    assert raw["matrix"] == "X" and via_layer["matrix"] == "layers['counts']"
    assert via_layer["macro_f1"] == raw["macro_f1"] and via_layer["per_type"] == raw["per_type"]
    b.layers["raw"] = b.layers.pop("counts")
    named = em.panel_resolvability(b, "cell_type", panel, target_depth=3, layer="raw")
    assert named["matrix"] == "layers['raw']" and named["macro_f1"] == raw["macro_f1"]
    with pytest.raises(ValueError, match="layer="):                # log X, no counts layer: no silent full depth
        em.panel_resolvability(b, "cell_type", panel, target_depth=3)
