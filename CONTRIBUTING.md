# Contributing

Bug reports, metric requests and pull requests are welcome. For anything bigger than a fix, open an
issue first, so the metric and its source are agreed on before code is written.

## Development setup

```bash
git clone https://github.com/p-gueguen/spatial-anno-metrics
cd spatial-anno-metrics
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,spatial]"
```

`[spatial]` installs squidpy; without it the Moran's I and neighbourhood tests are skipped. Python 3.10
is the oldest supported version. On 3.10, `import squidpy` reaches `pkg_resources` (through spatialdata
and xarray-schema), which setuptools 82 and later no longer ship: `python -m venv` seeds an older
setuptools, a uv venv does not, so there add `"setuptools<82"`.

## Tests and lint

```bash
pytest -q
ruff check src tests benchmarks
```

CI runs both on Python 3.10 and 3.12, and a lowest-direct job that installs every direct dependency
at the floor declared in `pyproject.toml`. A new or raised dependency needs a floor that job passes.

Tests run offline, on synthetic data or committed fixtures under 5 MB each, with optional dependencies
behind `pytest.importorskip`.

## Adding a metric

1. **Cite the source.** A metric reproduces a published definition (paper, DOI) or a reference
   implementation (URL, version or commit, licence). MIT or BSD code may be ported, with a comment
   naming where it came from. GPL code may be run to produce expected values, never copied in.
2. **Write the test first** and watch it fail. Give it a positive control (a synthetic case where the
   metric must move: shuffled labels, injected contamination) and a negative control (a case where it
   must not). When a reference implementation exists, commit a small golden under
   `tests/goldens/<name>/`: its inputs, its outputs, and a `SOURCE.md` with the upstream name, version
   or commit, licence and URL.
3. **Put constants in `src/spatial_anno_metrics/thresholds.py`**, in your module's dict, with the
   value, its source and whether it is provisional. No bare cut-offs inside the metric.
4. **Docstrings follow numpydoc**: `ruff check --select D` on the files you touched.
5. **Add a CHANGELOG fragment**, `changes/<short-name>.md`: the CHANGELOG entry, plus the README row
   and the catalog entry (`docs/cell-annotation-quality-metrics.md`) for a new metric. Do not edit
   `CHANGELOG.md`, the README or the catalog directly; fragments are folded in at release time.
