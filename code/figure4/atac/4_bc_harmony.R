# --- SETUP ----
library(Signac)
library(Seurat)
library(dplyr)
library(qs2)
library(ggplot2)
library(harmony)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/atac"
PLOTS_PATH = "plots/atac"
PLOTS_PATH_PDF = "plots/atac/pdf"


# ---- LOAD FILES ----
seu <- qs_read(file.path(ANALYSIS_PATH, "3_atac_preprocessed.qs2"))

# ---- BATCH CORRECTION VIA HARMONY ----
seu <- RunHarmony(
  object = seu,
  group.by.vars = "sample.name",
  reduction.use = "lsi",
  assay.use = "ATAC",
  project.dim = FALSE
)
DepthCor(seu, reduction = "harmony")
# Harmony components 1 and 2 correlate too strongly with sequencing depth
seu <- RunUMAP(seu, reduction = 'harmony', dims = 3:30)

bc_harmony_lsi <- DimPlot(
  seu,
  group.by = "sample.name",
  reduction = "harmony",
  dims = c(3, 4)
)
bc_harmony_umap <- DimPlot(seu, group.by = "sample.name")
ggsave("4_bc_harmony_lsi.png", bc_harmony_lsi, path = PLOTS_PATH)
ggsave("4_bc_harmony_lsi.pdf", bc_harmony_lsi, path = PLOTS_PATH_PDF)
ggsave("4_bc_harmony_umap.png", bc_harmony_umap, path = PLOTS_PATH)
ggsave("4_bc_harmony_umap.pdf", bc_harmony_umap, path = PLOTS_PATH_PDF)


# ---- SAVE INTEGRADET OBJECT ----
qs_save(seu, file.path(ANALYSIS_PATH, "4_atac_harmony.qs2"))
