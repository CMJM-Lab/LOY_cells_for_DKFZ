# ---- SETUP ----
library(dplyr)
library(ggplot2)
library(Seurat)
library(qs2)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/rna"
PLOTS_PATH = "plots/rna"
PLOTS_PATH_PDF = "plots/rna/pdf"


# ---- READ SEURAT RNA OBJECT ----
rna <- qs_read(file.path(ANALYSIS_PATH, "3_rna_preprocessed.qs2"))


# --- ADD EMT ENRICHMENT GENES SCORE ----
enrichr_GSEA <- readLines("res/enrichr_GSEA.txt") # 13 genes
prerank_GSEA <- readLines("res/prerank_GSEA.txt") # 58 genes
msigdb <- readLines("res/MSigDB_Hallmark_2020.txt") # 200 genes

rna <- AddModuleScore(
  rna,
  features = list(enrichr_GSEA, prerank_GSEA, msigdb),
  name = "emt_score_"
)


# ---- PLOT EMT ENRICHMENT ----
emt_score_1 <- rna@meta.data %>%
  ggplot(aes(sample.name, emt_score_1, fill = sample.name)) +
  geom_boxplot() +
  scale_fill_manual(values = color.sample.name) +
  ylab("Enrichr_GSEA (13 genes) Score") +
  xlab("Sample")
ggsave("5_emt_score_1.png", emt_score_1, path = PLOTS_PATH)
ggsave("5_emt_score_1.png", emt_score_1, path = PLOTS_PATH_PDF)

emt_score_2 <- rna@meta.data %>%
  ggplot(aes(sample.name, emt_score_2, fill = sample.name)) +
  geom_boxplot() +
  scale_fill_manual(values = color.sample.name) +
  ylab("Prerank_GSEA (58 genes) Score") +
  xlab("Sample")
ggsave("5_emt_score_2.png", emt_score_2, path = PLOTS_PATH)
ggsave("5_emt_score_2.png", emt_score_2, path = PLOTS_PATH_PDF)

emt_score_3 <- rna@meta.data %>%
  ggplot(aes(sample.name, emt_score_3, fill = sample.name)) +
  geom_boxplot() +
  scale_fill_manual(values = color.sample.name) +
  ylab("MSigDB_Hallmark (200 genes) Score") +
  xlab("Sample")
ggsave("5_emt_score_3.png", emt_score_3, path = PLOTS_PATH)
ggsave("5_emt_score_3.png", emt_score_3, path = PLOTS_PATH_PDF)
