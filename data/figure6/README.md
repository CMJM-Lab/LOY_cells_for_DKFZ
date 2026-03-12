# Figure 6 data layout

Tracked content:
- reference/: fixed reference metadata used for reproducibility checks
- derived/: plot-ready CSV/TSV/JSON outputs intended for sharing

Untracked runtime content:
- runtime/downloads/: cached GEO archives
- runtime/raw/: prepared MTX/TSV inputs
- runtime/processed/: generated h5ad cache

The runtime subtree is intentionally excluded from git and can be regenerated from notebook 01.
