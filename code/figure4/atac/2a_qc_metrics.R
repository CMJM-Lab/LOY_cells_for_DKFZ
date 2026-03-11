# ---- SETUP ----
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)

source("code/figure4/helper_functions.R")

ANALYSIS_PATH = "analysis/atac"
PLOTS_PATH = "plots/atac"
PLOTS_PATH_PDF = "plots/atac/pdf"


# ---- READ MERGED SEURAT OBJECT ----
atac <- qs_read(file.path(ANALYSIS_PATH, "1_atac_imported.qs2"))


# ---- CALCULATE QC METRICS ----
atac <- NucleosomeSignal(atac)
atac <- TSSEnrichment(atac)
atac$FRiP <- atac$atac_peak_region_fragments / atac$atac_fragments
atac$blacklist_ratio <- FractionCountsInRegion(atac, blacklist_hg38_unified)
atac$log1p_nCount_ATAC <- log1p(atac$nCount_ATAC)


# ---- PLOT QC METRICS -----
qc_nCount_ATAC <- plot_qc_metric(atac@meta.data, "nCount_ATAC", mad.cutoff = 5)
ggsave("2_qc_nCount_ATAC.png", qc_nCount_ATAC, path = PLOTS_PATH)
ggsave("2_qc_nCount_ATAC.pdf", qc_nCount_ATAC, path = PLOTS_PATH_PDF)

qc_TSS.enrichment <- plot_qc_metric(
  atac@meta.data,
  "TSS.enrichment",
  mad.cutoff = 5
)
ggsave("2_qc_TSS.enrichment.png", qc_TSS.enrichment, path = PLOTS_PATH)
ggsave("2_qc_TSS.enrichment.pdf", qc_TSS.enrichment, path = PLOTS_PATH_PDF)

qc_nucleosome_signal <- plot_qc_metric(
  atac@meta.data,
  "nucleosome_signal",
  mad.cutoff = 5
)
ggsave("2_qc_nucleosome_signal.png", qc_nucleosome_signal, path = PLOTS_PATH)
ggsave(
  "2_qc_nucleosome_signal.pdf",
  qc_nucleosome_signal,
  path = PLOTS_PATH_PDF
)

qc_FRiP <- plot_qc_metric(atac@meta.data, "FRiP", mad.cutoff = 5)
ggsave("2_qc_FRiP.png", qc_FRiP, path = PLOTS_PATH)
ggsave("2_qc_FRiP.pdf", qc_FRiP, path = PLOTS_PATH_PDF)


# ---- ASSESS LOW AND HIGH QUALITY CELLS ----
low.qc.ncount <- get_low_qc_cells(atac@meta.data, "nCount_ATAC", mad.cutoff = 5)
low.qc.TSS <- get_low_qc_cells(atac@meta.data, "TSS.enrichment", mad.cutoff = 5)
low.qc.nucleosome <- get_low_qc_cells(
  atac@meta.data,
  "nucleosome_signal",
  mad.cutoff = 5
)
low.qc.frip <- get_low_qc_cells(atac@meta.data, "FRiP", mad.cutoff = 5)

low.qc.atac <- unique(c(
  low.qc.ncount,
  low.qc.TSS,
  low.qc.nucleosome,
  low.qc.frip
))
writeLines(low.qc.atac, "analysis/cells_low_qc_atac.txt")


# ---- SAVE QC METRICS OBJECT ----
qs_save(atac, file.path(ANALYSIS_PATH, "2_atac_imported_qc_metrics.qs2"))
