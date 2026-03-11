# ---- SETUP ----
library(dplyr)
library(ggplot2)
library(Signac)
library(Seurat)
library(qs2)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/atac"
PLOTS_PATH = "plots/atac"
PLOTS_PATH_PDF = "plots/atac/pdf"


# ---- READ ATAC SEURAT OBJECT WITH QC METRICS ----
atac <- qs_read(file.path(ANALYSIS_PATH, "2_atac_imported_qc_metrics.qs2"))


# --- LOAD QC RESULTS ----
low.qc.rna <- readLines("analysis/cells_low_qc_rna.txt")
low.qc.rna.doublet <- readLines("analysis/cells_low_qc_rna_doublets.txt")
low.qc.atac <- readLines("analysis/cells_low_qc_atac.txt")
low.qc <- unique(c(low.qc.rna, low.qc.rna.doublet, low.qc.atac))
writeLines(low.qc, "analysis/cells_low_qc.txt")

all.cells <- readLines("analysis/cells_raw.txt")
high.qc <- setdiff(all.cells, low.qc)
writeLines(high.qc, "analysis/cells_high_qc.txt")


# ---- SUBSET AND SAVE SEURAT ATAC OBJECT ----
atac <- subset(atac, cells = high.qc)
qs_save(atac, file.path(ANALYSIS_PATH, "2_atac_filtered.qs2"))
