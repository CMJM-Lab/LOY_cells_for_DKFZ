# Figure 6 / S6 Pipeline

Run the notebooks in order:
1. notebooks/01_A549_Core_Pipeline.ipynb
2. notebooks/02_A549_LOY_TreeMet_Analysis.ipynb
3. notebooks/03_A549_FigS6_Supplementary.ipynb
4. notebooks/04_A549_FigS6H_Threshold_Exploration.ipynb (optional exploration)

Expected repository layout:
- code lives under code/figure6/
- runtime data lives under data/figure6/runtime/
- derived tables live under data/figure6/derived/
- exported PDFs and PNGs live under figures/figure6/

Notebook 01 downloads GEO accession GSE161363 and writes A549_integrated_final.h5ad to data/figure6/runtime/processed/. Notebooks 02 and 03 read that cache and generate publication outputs. Notebook 04 reads the canonical Figure S6H clone export from notebook 03 and performs threshold sweeps without changing the production figures.

The current S6H threshold observations and candidate rationale are recorded in `docs/figure6/FigS6H_threshold_notes.md`.

The `figure6` directory is the shared analysis module for both the main Figure 6 panels and the related Figure S6 supplementary panels. Figure S6 outputs are maintained under the explicit `03_FigS6_Supplementary` submodule to avoid confusion with a separate top-level figure number. Threshold exploration lives in `04_FigS6H_Threshold_Exploration` and is intentionally separated from the production notebook.

For Figure 6, the current validated provenance is:
- notebook 01 was re-executed from runtime/raw to regenerate the processed h5ad
- the regenerated object passed the committed fingerprint check
- notebook 02 was then re-executed from that regenerated h5ad to refresh the current Figure 6 outputs
