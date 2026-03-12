# Data staging (runtime; not versioned)

The publication notebooks under code/figure6/notebooks/ download and stage GEO data under data/figure6/runtime/ so the pipeline can run end-to-end without depending on the larger exploratory project tree.

Expected runtime layout after running notebook 01:
- data/figure6/runtime/downloads/GSE161363*.tar (cached GEO archive from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161363&format=file)
- data/figure6/runtime/downloads/GSE161363_RAW/ (extracted GEO files)
- data/figure6/runtime/raw/M5k/ and data/figure6/runtime/raw/LM0/ (prepared MTX/TSV inputs for Scanpy)
- data/figure6/runtime/processed/A549_integrated_final.h5ad (derived cache consumed by notebooks 02 and 03)

Tracked companion data:
- data/figure6/reference/ stores fixed reference metadata such as the fingerprint JSON.
- data/figure6/derived/ stores plot-ready tables and small JSON/TSV outputs intended for sharing.

Notes:
- Raw downloads and large derived caches are intentionally not tracked in git.
- For a clean re-download, delete only:
  - data/figure6/runtime/downloads/
  - data/figure6/runtime/raw/
- To force regeneration of the processed AnnData cache, delete:
  - data/figure6/runtime/processed/
