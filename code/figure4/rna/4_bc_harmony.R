# ---- SETUP ----
library(ggplot2)
library(Seurat)
library(SeuratWrappers)
library(qs2)

ANALYSIS_PATH = "analysis/rna"
PLOTS_PATH = "plots/rna"


# ---- BATCH CORRECTION ----
rna <- qs_read(file.path(ANALYSIS_PATH, "3_rna_preprocessed.qs2"))
rna <- split(rna, f = rna$orig.ident, layers = c("counts", "data"))
rna <- IntegrateLayers(
  object = rna,
  method = HarmonyIntegration,
  new.reduction = "harmony"
)
rna <- JoinLayers(rna)


# ---- PCA COMPARISON ----
bc_harmony_pca <- DimPlot(rna, reduction = "pca", group.by = "orig.ident") +
  DimPlot(rna, reduction = "harmony", group.by = "orig.ident")
ggsave("4_bc_harmony_pca.png", bc_harmony_pca, path = PLOTS_PATH, width = 12)


# --- UMAP ----
rna <- FindNeighbors(rna, reduction = "harmony", dims = 1:30)
rna <- RunUMAP(rna, dims = 1:30, reduction = "harmony")
rna <- FindClusters(rna, resolution = 0.25)
bc_harmony_umap <- DimPlot(rna, reduction = "umap", group.by = "orig.ident") +
  DimPlot(rna, reduction = "umap", group.by = "condition")
ggsave("4_bc_harmony_umap.png", bc_harmony_umap, path = PLOTS_PATH, width = 12)


# ---- SAVE INTEGRATED rnaRAT OBJECT ----
qs_save(rna, file.path(ANALYSIS_PATH, "4_rna_harmony.qs2"))
