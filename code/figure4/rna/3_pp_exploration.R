# ---- SETUP ----
library(dplyr)
library(ggplot2)
library(Seurat)
library(qs2)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/rna"
PLOTS_PATH = "plots/rna"
PLOTS_PATH_PDF = "plots/rna/pdf"


# ---- READ FILTERED SEURAT OBJECT ----
rna <- qs_read(file.path(ANALYSIS_PATH, "2_rna_filtered.qs2"))


# ---- NORMALIZE DATA ----
rna <- NormalizeData(rna)


# ---- FEATURE SELECTION ----
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(rna), 10)
plot <- VariableFeaturePlot(rna)
pp_var_features <- LabelPoints(
  plot = plot,
  points = top10,
  repel = TRUE,
  max.overlaps = 100
) +
  WhiteBackground()
ggsave("3_pp_var_features.png", pp_var_features, path = PLOTS_PATH)
ggsave("3_pp_var_features.pdf", pp_var_features, path = PLOTS_PATH_PDF)


# ---- SCALE DATA ----
rna <- ScaleData(rna)


# ---- PCA ----
rna <- RunPCA(rna)
pp_pca <- DimPlot(
  rna,
  reduction = "pca",
  group.by = c("sample.name", "condition")
) +
  scale_color_manual(values = color.condition)
ggsave("3_pp_pca.png", pp_pca, path = PLOTS_PATH, width = 12)
ggsave("3_pp_pca.pdf", pp_pca, path = PLOTS_PATH_PDF, width = 12)
ElbowPlot(rna, ndims = 50)


# ---- UMAP ----
rna <- FindNeighbors(rna, dims = 1:30)
rna <- RunUMAP(rna, dims = 1:30)
rna <- FindClusters(rna, resolution = 0.25, algorithm = 4, random.seed = 1)

pp_umap_by_sample <- DimPlot(
  rna,
  reduction = "umap",
  group.by = "sample.name",
  label = TRUE
) +
  theme(plot.title = element_blank()) +
  NoLegend()
pp_umap_by_condition <- DimPlot(
  rna,
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


# ---- PLOTS ----
Y_GENES <- c("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY")
X_GENES <- c("RPS4X", "DDX3X", "KDM6A", "KDM5C", "EIF1AX")
HK_GENES <- c("ACTB", "TUBB", "TUBA1B", "GAPDH", "PGK1", "ENO1", "LDHA")
EPI_GENES <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "MUC1", "MUC4")
MES_GENES <- c(
  "VIM",
  "CDH2",
  "FN1",
  "ITGA5",
  "S100A4",
  "COL1A1",
  "COL1A2",
  "COL3A1"
)
EMT_GENES <- c("SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "FOXC2")

# UTY UMAP
pl_umap_UTY <- FeaturePlot(rna, "UTY") + theme(plot.title = element_blank())
ggsave(
  "3_pl_umap_UTY.png",
  pl_umap_UTY,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pl_umap_UTY.pdf",
  pl_umap_UTY,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# LOX UMAP
pl_umap_LOX <- FeaturePlot(rna, "LOX") + theme(plot.title = element_blank())
ggsave(
  "3_pl_umap_LOX.png",
  pl_umap_LOX,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pl_umap_LOX.pdf",
  pl_umap_LOX,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# THY1 UMAP
pl_umap_THY1 <- FeaturePlot(rna, "THY1") + theme(plot.title = element_blank())
ggsave(
  "3_pl_umap_THY1.png",
  pl_umap_THY1,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pl_umap_THY1.pdf",
  pl_umap_THY1,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# CDH2 UMAP
pl_umap_CDH2 <- FeaturePlot(rna, "CDH2") + theme(plot.title = element_blank())
ggsave(
  "3_pl_umap_CDH2.png",
  pl_umap_CDH2,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pl_umap_CDH2.pdf",
  pl_umap_CDH2,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# SEX GENES UMAP
pl_umap_sex_genes <- FeaturePlot(rna, c(Y_GENES, X_GENES))
ggsave(
  "3_pl_umap_sex_genes.png",
  pl_umap_sex_genes,
  path = PLOTS_PATH,
  width = 8,
  height = 8
)
ggsave(
  "3_pl_umap_sex_genes.pdf",
  pl_umap_sex_genes,
  path = PLOTS_PATH_PDF,
  width = 8,
  height = 8
)

# SEX GENES DOTPLOT
pl_dotplot_sex_genes <- DotPlot(
  rna,
  c(Y_GENES, X_GENES),
  group.by = "sample.name",
  scale = FALSE
) +
  scale_y_discrete(limits = rev(levels(rna$sample.name))) +
  ggtitle("Y AND X GENES") +
  WhiteBackground() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  "3_pl_dotplot_sex_genes.png",
  pl_dotplot_sex_genes,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_dotplot_sex_genes.pdf",
  pl_dotplot_sex_genes,
  path = PLOTS_PATH_PDF,
  width = 12
)

# EPITHELIAL GENES DOTPLOT
pl_dotplot_epi_genes <- DotPlot(
  rna,
  c(EPI_GENES),
  group.by = "sample.name",
  scale = FALSE
) +
  scale_y_discrete(limits = rev(levels(rna$sample.name))) +
  ggtitle("Epithelial marker genes") +
  WhiteBackground() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  "3_pl_dotplot_epi_genes.png",
  pl_dotplot_epi_genes,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_dotplot_epi_genes.pdf",
  pl_dotplot_epi_genes,
  path = PLOTS_PATH_PDF,
  width = 12
)

# MESENCHYMAL GENES DOTPLOT
pl_dotplot_mes_genes <- DotPlot(
  rna,
  c(MES_GENES),
  group.by = "sample.name",
  scale = FALSE
) +
  scale_y_discrete(limits = rev(levels(rna$sample.name))) +
  ggtitle("Mesenchymal marker genes") +
  WhiteBackground() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  "3_pl_dotplot_mes_genes.png",
  pl_dotplot_mes_genes,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_dotplot_mes_genes.pdf",
  pl_dotplot_mes_genes,
  path = PLOTS_PATH_PDF,
  width = 12
)

# EMT GENES DOTPLOT
pl_dotplot_emt_genes <- DotPlot(
  rna,
  c(EMT_GENES),
  group.by = "sample.name",
  scale = FALSE
) +
  scale_y_discrete(limits = rev(levels(rna$sample.name))) +
  ggtitle("EMT marker genes") +
  WhiteBackground() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  "3_pl_dotplot_emt_genes.png",
  pl_dotplot_emt_genes,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_dotplot_emt_genes.pdf",
  pl_dotplot_emt_genes,
  path = PLOTS_PATH_PDF,
  width = 12
)

# DE GENES
de.genes <- FindMarkers(
  rna,
  ident.1 = "LOY",
  ident.2 = "ROY",
  group.by = "condition"
)
downsample.size <- round(min(table(rna$sample.name)) / 100)

top.genes.downregulated <- de.genes %>%
  filter(avg_log2FC < 0) %>%
  head(50) %>%
  rownames()
pl_heatmap_genes_downregulated <- DoHeatmap(
  subset(rna, downsample = downsample.size),
  top.genes.downregulated,
  group.by = "sample.name"
)
ggsave(
  "3_pl_heatmap_genes_downregulated.png",
  pl_heatmap_genes_downregulated,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_heatmap_genes_downregulated.pdf",
  pl_heatmap_genes_downregulated,
  path = PLOTS_PATH_PDF,
  width = 12
)

top.genes.upregulated <- de.genes %>%
  filter(avg_log2FC > 0) %>%
  head(50) %>%
  rownames()
pl_heatmap_genes_upregulated <- DoHeatmap(
  subset(rna, downsample = downsample.size),
  top.genes.upregulated,
  group.by = "sample.name"
)
ggsave(
  "3_pl_heatmap_genes_upregulated.png",
  pl_heatmap_genes_upregulated,
  path = PLOTS_PATH,
  width = 12
)
ggsave(
  "3_pl_heatmap_genes_upregulated.pdf",
  pl_heatmap_genes_upregulated,
  path = PLOTS_PATH_PDF,
  width = 12
)


# ---- SAVE PREPROCESSED SEURAT OBJECT ----
qs_save(rna, file.path(ANALYSIS_PATH, "3_rna_preprocessed.qs2"))
