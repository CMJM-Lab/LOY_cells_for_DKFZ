# ---- SETUP ----
library(dplyr)
library(ggplot2)
library(Seurat)
library(qs2)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/rna"
PLOTS_PATH = "plots/rna"
PLOTS_PATH_PDF = "plots/rna/pdf"


# ---- READ MERGED SEURAT OBJECT ----
rna <- qs_read(file.path(ANALYSIS_PATH, "1_rna_imported.qs2"))


# ---- CALCULATE QC METRICS ----
rna$percent.mt <- PercentageFeatureSet(rna, pattern = "^MT-")
rna$log1p_nFeature_RNA <- log1p(rna$nFeature_RNA)
rna$log1p_nCount_RNA <- log1p(rna$nCount_RNA)


# ---- QC PLOTS ----
qc_mt <- plot_qc_metric(rna@meta.data, "percent.mt", mad.cutoff = 3)
ggsave("2_qc_mt.png", qc_mt, path = PLOTS_PATH)
ggsave("2_qc_mt.pdf", qc_mt, path = PLOTS_PATH_PDF)

qc_nfeature <- plot_qc_metric(
  rna@meta.data,
  "log1p_nFeature_RNA",
  mad.cutoff = 5
)
ggsave("2_qc_nfeature.png", qc_nfeature, path = PLOTS_PATH)
ggsave("2_qc_nfeature.pdf", qc_nfeature, path = PLOTS_PATH_PDF)

qc_ncount <- plot_qc_metric(rna@meta.data, "log1p_nCount_RNA", mad.cutoff = 5)
ggsave("2_qc_ncount.png", qc_ncount, path = PLOTS_PATH)
ggsave("2_qc_ncount.pdf", qc_ncount, path = PLOTS_PATH_PDF)

qc_nfeat_mt <- plot_qc_scatter(
  rna@meta.data,
  "log1p_nFeature_RNA",
  "percent.mt",
  5,
  3
)
ggsave("2_qc_nfeat_mt.png", qc_nfeat_mt, path = PLOTS_PATH, width = 10)
ggsave("2_qc_nfeat_mt.pdf", qc_nfeat_mt, path = PLOTS_PATH_PDF, width = 10)

qc_ncount_mt <- plot_qc_scatter(
  rna@meta.data,
  "log1p_nCount_RNA",
  "percent.mt",
  5,
  3
)
ggsave("2_qc_ncount_mt.png", qc_ncount_mt, path = PLOTS_PATH, width = 10)
ggsave("2_qc_ncount_mt.pdf", qc_ncount_mt, path = PLOTS_PATH_PDF, width = 10)


# ---- ASSESS LOW QUALITY CELLS ----
low.qc.mt <- get_low_qc_cells(rna@meta.data, "percent.mt", 3)
low.qc.nfeature <- get_low_qc_cells(rna@meta.data, "log1p_nFeature_RNA", 5)
low.qc.ncount <- get_low_qc_cells(rna@meta.data, "log1p_nCount_RNA", 5)
low.qc.rna <- unique(c(low.qc.mt, low.qc.ncount, low.qc.nfeature))
writeLines(low.qc.rna, "analysis/cells_low_qc_rna.txt")


# ---- IDENTIFY DOUBLETS ----
# Setup
library(scDblFinder)
library(BiocParallel)

# Doublet detection
sce <- as.SingleCellExperiment(rna)
sce <- scDblFinder(sce, samples = "orig.ident", BPPARAM = MulticoreParam(8))
rna$doublet_class <- sce$scDblFinder.class
rna$doublet_score <- sce$scDblFinder.score
rm(sce)

# Plot doublets
qc_doublet <- rna@meta.data %>%
  ggplot(aes(sample.name, doublet_score, fill = sample.name)) +
  geom_violin() +
  geom_jitter(aes(color = doublet_class), size = 0.1, width = 0.1) +
  scale_color_manual(values = c("singlet" = "#00AFBB", "doublet" = "#FC4E07")) +
  theme(legend.position = "right") +
  guides(fill = "none")
ggsave("2_qc_doublet.png", qc_doublet, path = PLOTS_PATH, width = 10)
ggsave("2_qc_doublet.pdf", qc_doublet, path = PLOTS_PATH_PDF, width = 10)
low.qc.rna.doublets <- rna@meta.data %>%
  filter(doublet_class == "doublet") %>%
  rownames()
writeLines(low.qc.rna.doublets, "analysis/cells_low_qc_rna_doublets.txt")


# ---- SAVE QC METRICS OBJECT ----
qs_save(rna, file.path(ANALYSIS_PATH, "2_rna_imported_qc_metrics.qs2"))
