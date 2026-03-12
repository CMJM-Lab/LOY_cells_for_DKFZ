# Methods (Publication-ready)

## Overview and data source
We reanalyzed single-cell RNA-seq and lineage-tracing data from an A549 mouse xenograft metastasis model (GEO accession GSE161363; Ref. 36). Analyses were performed in Python 3.10.12 using Scanpy 1.11.0 and AnnData 0.11.3. Exact package versions are provided in code/figure6/environment.yml, and all code to reproduce the results is included as reviewer-ready notebooks under code/figure6/notebooks/.

## Data staging and inputs
Raw GEO files were downloaded and cached under data/figure6/runtime/downloads/ and extracted to data/figure6/runtime/raw/ as Scanpy-readable inputs (MTX plus barcode/gene tables). Notebook 01 writes the processed AnnData cache to data/figure6/runtime/processed/, and downstream derived tables are exported to data/figure6/derived/ while publication PDFs are written to figures/figure6/.

## Provenance of the current Figure 6 outputs
The current Figure 6 outputs in this subtree were regenerated end-to-end within dkfz_upload rather than copied from the historical 06_Dissemination snapshot. In the verified run on 2026-03-12, notebook 01 was executed from the staged raw MTX/TSV inputs under data/figure6/runtime/raw/ and forced to overwrite data/figure6/runtime/processed/A549_integrated_final.h5ad. The regenerated processed object passed the committed reference fingerprint check stored under data/figure6/reference/.

Notebook 02 was then re-executed against that freshly regenerated processed AnnData object to produce the current Figure 6 derived tables and panel exports. Accordingly, the current Figure 6 outputs should be interpreted as raw-to-figure reproductions from the local runtime/raw inputs, not as figures derived from a pre-existing staged processed h5ad cache.

## Quality control and preprocessing
Cells were filtered using the following thresholds: 200 <= n_genes <= 7000 and mitochondrial fraction < 10 percent. Genes detected in fewer than 3 cells were removed. For each dataset, counts were normalized to 10000 per cell and log1p-transformed, and highly variable genes were selected using Scanpy defaults (min_mean=0.0125, max_mean=3, min_disp=0.5).

## Integration, dimensionality reduction, and clustering
Principal component analysis was computed on highly variable genes (n_comps=100). To integrate across datasets, Harmony batch correction was applied to the PCA space using the dataset label (theta=2, max_iter=100), yielding the corrected embedding X_pca_harmony. A k-nearest-neighbor graph was built using 40 PCs (n_neighbors=15, n_pcs=40), and UMAP was computed (min_dist=0.5, spread=1.0). Leiden clustering was performed using the igraph backend (flavor=igraph, n_iterations=2) at multiple resolutions (0.3, 0.6, 0.9, 1.2).

To support exact reproducibility, random seeds were fixed for neighbor graph construction and UMAP (random_state=0) and for Leiden clustering (random_state=42), matching the original analysis notebooks.

## Sex-chromosome and housekeeping signature scoring
Sex-chromosome and housekeeping gene signatures were scored per cell using Scanpy score_genes, requiring at least 3 signature genes present in the expression matrix and using 50 control genes (ctrl_size=50). We computed three scores: Y_score, X_score, and Housekeep_score.

Cells were classified as loss of Y (LOY) or retention of Y (ROY) using a fixed Y-score threshold of -0.2. This labeling is stored in adata.obs[Y_group].

## Metastatic capacity metrics and correlations
Metastatic capacity was quantified using the scTreeMetRate field provided with the dataset; for plotting and manuscript consistency we refer to this metric as scMetRate and store adata.obs[scMetRate] = adata.obs[scTreeMetRate] when needed. Associations between Y_score and metastatic capacity were evaluated using Spearman correlation on post-implantation cells (M5k).

For the main Figure 6 panels, panels D and E are all-cell compartment summaries. Panel F follows the historical lineage-enriched LOY/ROY definition traced from 05_Replot_new: lineage groups are first classified as LOY-enriched or ROY-enriched, and per-cell scMetRate values are then plotted within those lineage-defined groups.

## Clone-level LOY/ROY enrichment and metastatic phenotype (Figure S6H)
Clone identity was defined by the lineage group label (LineageGroup). For clone-level LOY/ROY enrichment, we evaluated each clone within each metastatic site and required a minimum of 20 cells per clone-site pair to be eligible. For each eligible site, we computed the LOY and ROY fractions within the clone. A clone was classified as LOY-enriched if the LOY fraction was >= 0.75 in every eligible site; it was classified as ROY-enriched if the ROY fraction was >= 0.75 in every eligible site; otherwise it was classified as Mixed and excluded from LOY-vs-ROY comparisons.

Clone metastatic phenotype was assigned based on the median clone TreeMetRate across M5k cells, using fixed cutoffs: non-metastatic (<0.001), weakly metastatic (0.001-0.008), and highly metastatic (>=0.008). Differences in metastatic phenotype composition between LOY- and ROY-enriched clones were summarized using stacked bar plots and tested using Fisher exact tests (Highly vs Others; AnyMet vs Non).

## Figure generation and export
Figures were generated with Matplotlib, Seaborn, and Scanpy plotting utilities and saved under figures/figure6/. The current validated main Figure 6 outputs are exported as both publication-ready PDFs and companion PNG files, whereas some retained pipeline or supplementary outputs may remain PDF-only unless rerun with the same dual-format save helper. To facilitate reviewer reuse in other software, plot-ready tables (CSV, TSV) and metadata JSON files were exported to data/figure6/derived/, including the exact cutoffs, lineage-classification settings, and category orders used for each panel.
