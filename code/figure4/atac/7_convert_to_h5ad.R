# ---- SETUP ----
library(anndataR)
library(Signac)
library(Seurat)
library(qs2)

ANALYSIS_PATH = "analysis/atac"
ANALYSIS_PATH_PYTHON = "analysis/scanpy"


# ---- CONVERT ----
seu <- qs_read(file.path(ANALYSIS_PATH, "2_atac_filtered.qs2"))
write_h5ad(seu, file.path(ANALYSIS_PATH_PYTHON, "2_atac_filtered.h5ad"), mode = "w")


seu <- qs_read(file.path(ANALYSIS_PATH, "3_atac_preprocessed.qs2"))
write_h5ad(seu, file.path(ANALYSIS_PATH_PYTHON, "3_atac_preprocessed.h5ad"), mode = "w")
