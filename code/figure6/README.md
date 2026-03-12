# Figure 6 / S6 Pipeline

Run the notebooks in order:
1. notebooks/01_A549_Core_Pipeline.ipynb
2. notebooks/02_A549_LOY_TreeMet_Analysis.ipynb
3. notebooks/03_A549_Supplementary.ipynb

Expected repository layout:
- code lives under code/figure6/
- runtime data lives under data/figure6/runtime/
- derived tables live under data/figure6/derived/
- exported PDFs and PNGs live under figures/figure6/

Notebook 01 downloads GEO accession GSE161363 and writes A549_integrated_final.h5ad to data/figure6/runtime/processed/. Notebooks 02 and 03 read that cache and generate publication outputs.

For Figure 6, the current validated provenance is:
- notebook 01 was re-executed from runtime/raw to regenerate the processed h5ad
- the regenerated object passed the committed fingerprint check
- notebook 02 was then re-executed from that regenerated h5ad to refresh the current Figure 6 outputs
