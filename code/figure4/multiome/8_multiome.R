# ---- SETUP ----
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH_RNA = "analysis/rna"
ANALYSIS_PATH_ATAC = "analysis/atac"
PLOTS_PATH = "plots/multiome"
PLOTS_PATH_PDF = "plots/multiome/pdf"


# ---- LOAD OBJECTS -----
rna <- qs_read(file.path(ANALYSIS_PATH_RNA, "3_rna_preprocessed.qs2"))
atac <- qs_read(file.path(ANALYSIS_PATH_ATAC, "3_atac_preprocessed.qs2"))

rna[["ATAC"]] <- atac[["ATAC"]]
rna[["lsi"]] <- atac[["lsi"]]
rna[["umap.ATAC"]] <- atac[["umap"]]
rna <- AddMetaData(rna, metadata = atac@meta.data)

rm(atac)


# ---- RUN WNN and UMAP ----
rna <- FindMultiModalNeighbors(
  object = rna,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 3:30),
  modality.weight.name = "RNA.weight"
)

rna <- RunUMAP(
  object = rna,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

pl_umap_multiome_by_sample <- DimPlot(
  rna,
  reduction = "wnn.umap",
  group.by = "sample.name",
  label = TRUE
) +
  theme(plot.title = element_blank()) +
  NoLegend()
pl_umap_multiome_by_condition <- DimPlot(
  rna,
  reduction = "wnn.umap",
  group.by = "condition"
) +
  scale_color_manual(values = color.condition) +
  theme(plot.title = element_blank()) +
  NoLegend()
ggsave(
  "8_pl_umap_multiome_by_sample.png",
  pl_umap_multiome_by_sample,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_by_sample.pdf",
  pl_umap_multiome_by_sample,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_by_condition.png",
  pl_umap_multiome_by_condition,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_by_condition.pdf",
  pl_umap_multiome_by_condition,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)


# ---- FEATURE PLOTS WITH NEW UMAP ----
# UTY UMAP
pl_umap_multiome_UTY <- FeaturePlot(rna, "UTY", reduction = "wnn.umap") +
  theme(plot.title = element_blank())
ggsave(
  "8_pl_umap_multiome_UTY.png",
  pl_umap_multiome_UTY,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_UTY.pdf",
  pl_umap_multiome_UTY,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# LOX UMAP
pl_umap_multiome_LOX <- FeaturePlot(rna, "LOX", reduction = "wnn.umap") +
  theme(plot.title = element_blank())
ggsave(
  "8_pl_umap_multiome_LOX.png",
  pl_umap_multiome_LOX,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_LOX.pdf",
  pl_umap_multiome_LOX,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# THY1 UMAP
pl_umap_multiome_THY1 <- FeaturePlot(rna, "THY1", reduction = "wnn.umap") +
  theme(plot.title = element_blank())
ggsave(
  "8_pl_umap_multiome_THY1.png",
  pl_umap_multiome_THY1,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_THY1.pdf",
  pl_umap_multiome_THY1,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# CDH2 UMAP
pl_umap_multiome_CDH2 <- FeaturePlot(rna, "CDH2", reduction = "wnn.umap") +
  theme(plot.title = element_blank())
ggsave(
  "8_pl_umap_multiome_CDH2.png",
  pl_umap_multiome_CDH2,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "8_pl_umap_multiome_CDH2.pdf",
  pl_umap_multiome_CDH2,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)


# ---- SAVE MULTIOME OBJECT ----
qs_save(rna, file.path("analysis/8_multiome.qs2"))
