# Figure 6 A549 Notes

This note tracks the current A549 / Figure 6 adaptation that was staged in `dkfz_upload` and then merged into this repository structure.

## Scope

The current merge-ready Figure 6 content lives under:

- `code/figure6/`
- `data/figure6/`
- `docs/figure6/`
- `figures/figure6/`

Only the validated main Figure 6 workflow is treated as merge-ready in this pass:

1. `code/figure6/notebooks/01_A549_Core_Pipeline.ipynb`
2. `code/figure6/notebooks/02_A549_LOY_TreeMet_Analysis.ipynb`

Supplementary Figure S6 outputs are not yet treated as stable merge targets because the `highly_meta` and `weakly_meta` thresholds still need to be refined.

## Provenance

The currently merged Figure 6 assets were validated with a raw-to-figure rerun on March 12, 2026:

1. notebook 01 regenerated `A549_integrated_final.h5ad` from locally staged raw inputs
2. the regenerated object passed the committed fingerprint check
3. notebook 02 regenerated the current main Figure 6 tables and plots from that fresh processed object

Reviewer-facing Figure 6 main-panel outputs are currently exported as both `pdf` and `png`.

## Runtime policy

Large runtime inputs and caches are intentionally not tracked in git. This includes:

- downloaded GEO archives
- prepared raw MTX/TSV staging files
- processed `*.h5ad` files

Only code, small reference files, plot-ready derived tables, and reviewer-facing figure exports are intended to be committed.

## Supporting docs

Figure-specific methods, data staging notes, and merge rationale are documented under `docs/figure6/`.
