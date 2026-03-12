# Figure 6 merge rationale

This note records what should be merged from `dkfz_upload` into `dkfz_skeleton`, what should be excluded, and why.

## Goal

`dkfz_upload` is a staging subtree used to adapt and validate the Figure 6 pipeline before the validated contents are merged into the Git-managed repository `dkfz_skeleton`.

The merge target is therefore:

- keep reviewer-facing, reproducible, small-footprint assets
- exclude runtime caches and large regenerable intermediates

## What should be merged

These categories are intended to move into `dkfz_skeleton`:

- `code/figure6/`
  notebooks, helper scripts, environment file, and figure-specific README
- `data/figure6/reference/`
  fixed reference metadata such as the committed fingerprint JSON
- `data/figure6/derived/`
  plot-ready CSV/TSV/JSON outputs used for review and figure regeneration checks
- `docs/figure6/`
  methods, data staging notes, and this merge rationale
- `figures/figure6/`
  reviewer-facing figure exports in `pdf` and `png`

## What should not be merged

These paths should remain local runtime/cache material and should not enter git history:

- `data/figure6/runtime/downloads/`
- `data/figure6/runtime/raw/`
- `data/figure6/runtime/processed/A549_integrated_final.h5ad`
- `__pycache__/`
- `.ipynb_checkpoints/`

## Rationale for excluding runtime files

The runtime layer contains either downloaded source archives, prepared local working copies, or large regenerable intermediates.

In particular, the processed `A549_integrated_final.h5ad` file should not be committed because:

- it is large and regenerable from notebook 01
- committing it would bloat repository history
- this would recreate the same large-blob problem already observed in the upstream DKFZ GitLab repository
- provenance is better preserved through code, methods, fingerprints, and documented execution order than through storing the cache itself in git

## Provenance standard for merge readiness

For Figure 6, merge readiness means:

1. notebook 01 can regenerate `data/figure6/runtime/processed/A549_integrated_final.h5ad` from `data/figure6/runtime/raw/`
2. the regenerated processed object passes the committed fingerprint check
3. notebook 02 can regenerate the current Figure 6 tables and figure files from that freshly rebuilt processed object
4. the exported tables and figures are internally consistent

This condition has been verified for the current Figure 6 outputs in `dkfz_upload`.

At the current stage, this validated status applies to:

- notebook 01 core pipeline outputs
- notebook 02 main Figure 6 outputs

The supplementary notebook 03 and the current `03_Supplementary` derived/figure outputs should still be treated as provisional because the next planned step is threshold refinement for the highly metastatic and weakly metastatic categories.

## Repository strategy

`dkfz_upload` should remain a non-repository staging subtree.

`dkfz_skeleton` should remain the only formal git repository in this workspace for the DKFZ-style publication package.

The recommended workflow is:

1. validate work inside `dkfz_upload`
2. copy only the merge-approved files into `dkfz_skeleton`
3. review the resulting diff inside `dkfz_skeleton`
4. commit there

Creating a second independent git repository inside `dkfz_upload` is not recommended because it would add parallel history management without improving the final publication package.
