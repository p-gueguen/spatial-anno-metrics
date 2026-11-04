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


def _background_only_ref(private: bool, n=150, n_genes=60, n_types=6, seed=0):
    """Six types that differ only by a mild per-type background programme (lognormal sigma 0.4 per gene),
    plus, when ``private``, one strong private marker each. At depth 50 the background alone gives
    macro-F1 ~0.73: above chance (1/6) but ~1 cell in 4 misassigned."""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(0, 1, n_genes)
    base *= 200 / base.sum()
    blocks, labels = [], []
    for k in range(n_types):
        mu = base * rng.lognormal(0, 0.4, n_genes)
        if private:
            mu[k] += 100
        blocks.append(rng.poisson(mu, size=(n, n_genes)))
        labels += [f"T{k}"] * n
    a = anndata.AnnData(X=np.vstack(blocks).astype("float32"))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = labels
    return a


def test_panel_resolvability_no_private_markers_is_not_resolvable():
    """ROADMAP 8h: with no private markers (macro-F1 ~0.7 at depth 50) no type is 'resolvable' any more; the
    same types WITH a private marker each stay resolvable (positive control)."""
    from spatial_anno_metrics import panel as em

    bg = _background_only_ref(private=False)
    r = em.panel_resolvability(bg, "cell_type", list(bg.var_names), target_depth=50)
    assert r["n_types"] == 6 and 0.6 < r["macro_f1"] < 0.85
    assert r["frac_resolvable"] == 0.0
    assert all(d["f1_null"] == pytest.approx(1 / 6) for d in r["per_type"].values())
    for d in r["per_type"].values():                             # chance-corrected F1: (f1 - null) / (1 - null)
        assert d["f1_adjusted"] == pytest.approx((d["f1"] - 1 / 6) / (1 - 1 / 6))
    pm = _background_only_ref(private=True)
    rp = em.panel_resolvability(pm, "cell_type", list(pm.var_names), target_depth=50)
    assert rp["frac_resolvable"] == 1.0


def test_panel_resolvability_old_tiers_reproducible():
    """The pre-0.7 tiers (raw F1 >= 0.5 / 0.3) come back when passed explicitly; F1s themselves are unchanged."""
    from spatial_anno_metrics import panel as em

    bg = _background_only_ref(private=False)
    new = em.panel_resolvability(bg, "cell_type", list(bg.var_names), target_depth=50)
    old = em.panel_resolvability(bg, "cell_type", list(bg.var_names), target_depth=50,
                                 f1_resolvable=0.5, f1_weak=0.3, null_adjusted=False)
    assert old["frac_resolvable"] == 1.0                         # the lenient verdict ROADMAP 8h describes
    assert old["macro_f1"] == new["macro_f1"]
    assert [d["f1"] for d in old["per_type"].values()] == [d["f1"] for d in new["per_type"].values()]


def test_panel_resolvability_two_types_at_chance_are_not_resolvable():
    """Two indistinguishable types score F1 ~0.5 = chance; the null-adjusted F1 is ~0, so neither is resolvable."""
    from spatial_anno_metrics import panel as em

    rng = np.random.default_rng(1)
    a = anndata.AnnData(X=rng.poisson(4.0, size=(200, 30)).astype("float32"))
    a.var_names = [f"g{i}" for i in range(30)]
    a.obs["cell_type"] = ["A"] * 100 + ["B"] * 100
    r = em.panel_resolvability(a, "cell_type", list(a.var_names), target_depth=50)
    assert r["status"] == "ok" and len(r["per_type"]) == 2
    assert all(d["f1_null"] == pytest.approx(0.5) for d in r["per_type"].values())
    assert all(abs(d["f1_adjusted"]) < 0.3 and d["tier"] == "not_resolvable" for d in r["per_type"].values())


def _imbalanced_pair(n_a=200, n_b=40, shift=0.4, seed=1, n_genes=30):
    """Majority type A (prevalence 5/6) and minority B whose rates differ by exp(+-shift) per gene."""
    rng = np.random.default_rng(seed)
    lam = np.full(n_genes, 4.0)
    lam_b = lam * np.exp(rng.choice([-shift, shift], n_genes))
    a = anndata.AnnData(X=np.vstack([rng.poisson(lam, (n_a, n_genes)),
                                     rng.poisson(lam_b, (n_b, n_genes))]).astype("float32"))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = ["A"] * n_a + ["B"] * n_b
    return a


def _rule(score, confused, ok, weak):
    """The documented tier rule; confused-fraction cuts 0.30 / 0.25."""
    return ("not_resolvable" if score < weak or confused >= 0.30
            else "weak" if score < ok or confused >= 0.25 else "resolvable")


def test_panel_resolvability_chance_correction_decides_imbalanced_tiers():
    """ROADMAP 8h: a 5:1 majority type scores raw F1 > 0.9 by its base rate alone. The chance correction (null =
    its prevalence 5/6, not 1/K) keeps it out of 'resolvable'; the minority type clears the F1 cut but loses
    a third of its cells to A, so the confusion cut makes it 'not_resolvable'."""
    from spatial_anno_metrics import panel as em

    a = _imbalanced_pair()
    r = em.panel_resolvability(a, "cell_type", list(a.var_names), target_depth=50)
    per = r["per_type"]
    assert per["A"]["f1_null"] == pytest.approx(200 / 240) and per["B"]["f1_null"] == pytest.approx(40 / 240)
    for d in per.values():
        assert d["f1_adjusted"] == pytest.approx((d["f1"] - d["f1_null"]) / (1 - d["f1_null"]))
        assert d["tier"] == _rule(d["f1_adjusted"], d["confused_frac"], 0.8, 0.5)
    A, B = per["A"], per["B"]
    assert A["f1"] > 0.9 and A["confused_frac"] < 0.25            # raw F1 alone would say 'resolvable'
    assert 0.5 <= A["f1_adjusted"] < 0.8 and A["tier"] == "weak"
    assert B["f1_adjusted"] >= 0.5 and B["confused_frac"] >= 0.3 and B["tier"] == "not_resolvable"


def test_panel_resolvability_old_tiers_reproducible_per_type():
    """Passing the pre-0.7 cuts back gives the pre-0.7 per-type tiers where the new ones differ."""
    from spatial_anno_metrics import panel as em

    a = _imbalanced_pair()
    new = em.panel_resolvability(a, "cell_type", list(a.var_names), target_depth=50)
    old = em.panel_resolvability(a, "cell_type", list(a.var_names), target_depth=50,
                                 f1_resolvable=0.5, f1_weak=0.3, null_adjusted=False)
    assert {k: d["tier"] for k, d in old["per_type"].items()} == {"A": "resolvable", "B": "not_resolvable"}
    for d in old["per_type"].values():
        assert d["tier"] == _rule(d["f1"], d["confused_frac"], 0.5, 0.3)
    assert new["per_type"]["A"]["tier"] != old["per_type"]["A"]["tier"]   # positive control: the rule changed
    assert [d["f1"] for d in old["per_type"].values()] == [d["f1"] for d in new["per_type"].values()]


def test_panel_resolvability_weak_confusion_cut_and_f1_weak_threshold():
    """Pins _CONF_WEAK = 0.25 independently of f1_resolvable = 0.80, and f1_weak = 0.50 independently of
    _CONF_NOT = 0.30."""
    from spatial_anno_metrics import panel as em

    n_genes = 30
    A_cells = np.zeros((200, n_genes), dtype="float32")
    A_cells[:, :5] = 50.0
    C_cells = np.zeros((200, n_genes), dtype="float32")
    C_cells[:, 10:15] = 50.0
    # 30 of 40 B cells (75 %) carry B's private markers (g5..g9); 10 of 40 (25.0 %) carry C's markers -> called C
    B_cells = np.zeros((40, n_genes), dtype="float32")
    B_cells[:30, 5:10] = 50.0
    B_cells[30:, 10:15] = 50.0
    a = anndata.AnnData(X=np.vstack([A_cells, B_cells, C_cells]))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = ["A"] * 200 + ["B"] * 40 + ["C"] * 200

    r = em.panel_resolvability(a, "cell_type", list(a.var_names), target_depth=250)
    B = r["per_type"]["B"]
    assert B["f1_adjusted"] == pytest.approx(59 / 70, abs=1e-6)    # 0.8429 > 0.80 (f1=6/7, f1_null=1/11)
    assert B["confused_frac"] == pytest.approx(0.25) and B["confused_with"] == "C"
    assert B["tier"] == "weak"
    assert r["per_type"]["A"]["tier"] == "resolvable" and r["per_type"]["C"]["tier"] == "resolvable"

    # D loses 20 % each to A, B, C (confused_frac = 0.20 < 0.25), with f1 = 4/7 and f1_adjusted = 3/7 = 0.4286 < 0.50
    D_cells = np.zeros((100, n_genes), dtype="float32")
    D_cells[:40, 15:20] = 50.0
    D_cells[40:60, :5] = 50.0
    D_cells[60:80, 5:10] = 50.0
    D_cells[80:100, 10:15] = 50.0
    B_clean = np.zeros((100, n_genes), dtype="float32")
    B_clean[:, 5:10] = 50.0
    a4 = anndata.AnnData(X=np.vstack([A_cells[:100], B_clean, C_cells[:100], D_cells]))
    a4.var_names = [f"g{i}" for i in range(n_genes)]
    a4.obs["cell_type"] = ["A"] * 100 + ["B"] * 100 + ["C"] * 100 + ["D"] * 100
    r4 = em.panel_resolvability(a4, "cell_type", list(a4.var_names), target_depth=250)
    D = r4["per_type"]["D"]
    assert D["confused_frac"] == pytest.approx(0.20)
    assert D["f1_adjusted"] == pytest.approx(3 / 7, abs=1e-6)      # 0.4286 in [0.20, 0.50), confused_frac < 0.25
    assert D["tier"] == "not_resolvable"
