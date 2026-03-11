# --- SETUP ----
library(Signac)
library(Seurat)
library(dplyr)
library(qs2)
library(ggplot2)
library(patchwork)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/atac"
PLOTS_PATH = "plots/atac"
PLOTS_PATH_PDF = "plots/atac/pdf"


# ---- LOAD FILES ----
seu <- qs_read(file.path(ANALYSIS_PATH, "2_atac_filtered.qs2"))


# ---- LSI ----
seu <- RunTFIDF(seu)
seu <- FindTopFeatures(seu, min.cutoff = "q25")
seu <- RunSVD(seu)

DepthCor(seu)
# LSI components 1 and 2 correlate too strongly with sequencing depth

pp_lsi <- DimPlot(
  seu,
  reduction = "lsi",
  group.by = c("sample.name", "condition"),
  dims = c(3, 4)
) +
  scale_color_manual(values = color.condition)
ggsave("3_pp_lsi.png", pp_lsi, path = PLOTS_PATH, width = 16)
ggsave("3_pp_lsi.pdf", pp_lsi, path = PLOTS_PATH_PDF, width = 16)


# ---- UMAP ----
seu <- RunUMAP(seu, reduction = "lsi", dims = 3:30)

pp_umap_by_sample <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "sample.name",
  label = TRUE
) +
  theme(plot.title = element_blank()) +
  NoLegend()
pp_umap_by_condition <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "condition"
) +
  scale_color_manual(values = color.condition) +
  theme(plot.title = element_blank()) +
  NoLegend()
ggsave(
  "3_pp_umap_by_sample.png",
  pp_umap_by_sample,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pp_umap_by_sample.pdf",
  pp_umap_by_sample,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)
ggsave(
  "3_pp_umap_by_condition.png",
  pp_umap_by_condition,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pp_umap_by_condition.pdf",
  pp_umap_by_condition,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)


# ---- GENOMIC REGIONS
Y_GENES <- c("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY")
X_GENES <- c("RPS4X", "DDX3X", "KDM6A", "KDM5C", "EIF1AX")

# THY1
pl_THY1_frequency <- CoveragePlot(seu, region = "11-119384984-119464985") &
  scale_fill_manual(values = color.sample.name)
ggsave(
  "3_pl_THY1_frequency.png",
  pl_THY1_frequency,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_THY1_frequency.pdf",
  pl_THY1_frequency,
  path = PLOTS_PATH_PDF,
  width = 12
)

# LOX
pl_LOX_frequency <- CoveragePlot(seu, region = "5-122038284-122118285") &
  scale_fill_manual(values = color.sample.name)
ggsave(
  "3_pl_LOX_frequency.png",
  pl_LOX_frequency,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_LOX_frequency.pdf",
  pl_LOX_frequency,
  path = PLOTS_PATH_PDF,
  width = 12
)

# CDH2
pl_CDH2_frequency <- CoveragePlot(seu, region = "18-27937445-28217446") &
  scale_fill_manual(values = color.sample.name)
ggsave(
  "3_pl_CDH2_frequency.png",
  pl_CDH2_frequency,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_CDH2_frequency.pdf",
  pl_CDH2_frequency,
  path = PLOTS_PATH_PDF,
  width = 12
)

# ANPEP
pl_ANPEP_frequency <- CoveragePlot(seu, region = "15-89775400-89855401") &
  scale_fill_manual(values = color.sample.name)
ggsave(
  "3_pl_ANPEP_frequency.png",
  pl_ANPEP_frequency,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_ANPEP_frequency.pdf",
  pl_ANPEP_frequency,
  path = PLOTS_PATH_PDF,
  width = 12
)

# SEX GENES
pl_SEX_GENES_frequency <- CoveragePlot(seu, c(Y_GENES, X_GENES)) &
  scale_fill_manual(values = color.sample.name)
ggsave(
  "3_pl_SEX_GENES_frequency.png",
  pl_SEX_GENES_frequency,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_SEX_GENES_frequency.pdf",
  pl_SEX_GENES_frequency,
  path = PLOTS_PATH_PDF,
  width = 12
)


# ---- SAVE PREPROCESSED OBJECT ----
qs_save(seu, file.path(ANALYSIS_PATH, "3_atac_preprocessed.qs2"))
