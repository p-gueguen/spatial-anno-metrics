"""Sphinx configuration for the spatial-anno-metrics docs (ReadTheDocs builds it from .readthedocs.yaml)."""
import sys
from importlib.metadata import version as _version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

project = "spatial-anno-metrics"
author = "Paul Gueguen"
copyright = "2025-2026, Paul Gueguen"
release = _version("spatial-anno-metrics")
version = ".".join(release.split(".")[:2])

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",  # after napoleon
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

# Notebooks are rendered from their saved outputs; the build never runs them.
nb_execution_mode = "off"
myst_enable_extensions = ["colon_fence", "dollarmath"]
myst_heading_anchors = 3

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = True
typehints_document_rtype = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
}

exclude_patterns = ["_build", "ROADMAP.md", "**.ipynb_checkpoints"]

html_theme = "furo"
html_title = f"spatial-anno-metrics {release}"
