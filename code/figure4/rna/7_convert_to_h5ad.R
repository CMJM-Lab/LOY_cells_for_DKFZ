# ---- SETUP ----
library(anndataR)
library(Seurat)
library(qs2)

ANALYSIS_PATH = "analysis/rna"
ANALYSIS_PATH_PYTHON = "analysis/scanpy"


# ---- CONVERT ----
seu <- qs_read(file.path(ANALYSIS_PATH, "2_rna_doublet_filtered.qs2"))
write_h5ad(seu, file.path(ANALYSIS_PATH_PYTHON, "2_rna_doublet_filtered.h5ad"), mode = "w")


seu <- qs_read(file.path(ANALYSIS_PATH, "3_rna_preprocessed.qs2"))
write_h5ad(seu, file.path(ANALYSIS_PATH_PYTHON, "3_rna_preprocessed.h5ad"), mode = "w")
