# ---- SETUP ---
library(qs2)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(AnnotationHub)
library(GenomeInfoDb)
library(future)

plan(multisession, workers = 8)
options(future.globals.maxSize = 50000 * 1024^2)

ANALYSIS_PATH = "analysis/atac"
DATA_PATH = "data/cellranger_out"
sample.list <- list.files(DATA_PATH)


# ---- CREATE COMMON PEAK SET ----
combined.peaks <- GRanges()
for (sample.name in sample.list) {
  path = file.path(DATA_PATH, sample.name, "atac_peaks.bed")
  peaks <- read.table(path, skip = 58, col.names = c("chr", "start", "end"))
  gr.peaks <- makeGRangesFromDataFrame(peaks)
  combined.peaks <- c(combined.peaks, gr.peaks)
}
combined.peaks <- reduce(combined.peaks)


# ---- CREATE OBJECTS OBJECTS ----
seurat.list <- list()
for (sample.name in sample.list) {
  md.path = file.path(DATA_PATH, sample.name, "per_barcode_metrics.csv")
  meta.data <- read.table(md.path, header = TRUE, sep = ",", row.names = 1)
  meta.data <- meta.data[meta.data$is_cell == 1, ]
  cells = rownames(meta.data)
  frag.path <- file.path(DATA_PATH, sample.name, "atac_fragments.tsv.gz")
  frags <- CreateFragmentObject(frag.path, cells)
  counts <- FeatureMatrix(frags, combined.peaks, cells = cells)
  counts.assay <- CreateChromatinAssay(
    counts,
    fragments = frags,
    min.cells = 20
  )
  atac <- CreateSeuratObject(
    counts.assay,
    assay = "ATAC",
    meta.data = meta.data,
    project = sample.name
  )
  seurat.list[[sample.name]] <- atac
}


# ---- MERGE OBJECTS ----
atac <- merge(
  seurat.list[[1]],
  seurat.list[2:length(seurat.list)],
  add.cell.ids = names(seurat.list),
  project = "LOY vs. ROY"
)
rm(seurat.list)


# ---- MODIFY META DATA ----
atac$condition <- ifelse(
  startsWith(atac@meta.data$orig.ident, "LOY"),
  "LOY",
  "ROY"
)

lookup <- c(
  LOY_D11 = "LOY1",
  LOY_E5b = "LOY2",
  LOY_E6b = "LOY3",
  LOY_F10 = "LOY4",
  WT_A3 = "ROY1",
  WT_C1 = "ROY2",
  WT_C8 = "ROY3",
  WT_D5 = "ROY4"
)
atac@meta.data$sample.name <- lookup[atac$orig.ident]


# ORDER OF META DATA
atac$orig.ident <- factor(
  atac$orig.ident,
  levels = c(
    'LOY_D11',
    'LOY_E5b',
    'LOY_E6b',
    'LOY_F10',
    'WT_A3',
    'WT_C1',
    'WT_C8',
    'WT_D5'
  )
)
atac$condition <- factor(atac$condition, levels = c("LOY", "ROY"))
atac$sample.name <- factor(
  atac$sample.name,
  levels = c('LOY1', 'LOY2', 'LOY3', 'LOY4', 'ROY1', 'ROY2', 'ROY3', 'ROY4')
)
Idents(atac) <- "sample.name"


# --- ANNOTATION ----
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v113")
ensdb <- ah[["AH119325"]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevels(annotations) <- mapSeqlevels(seqlevels(annotations), "Ensembl")
genome(annotations) <- "hg38"

Annotation(atac) <- annotations


# ---- SAVE SEURAT ATAC FILE ----
qs_save(atac, file.path(ANALYSIS_PATH, "1_atac_imported.qs2"))
