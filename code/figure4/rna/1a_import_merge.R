# ---- SETUP ---
library(Seurat)
library(qs2)

ANALYSIS_PATH = "analysis/rna"
DATA_PATH = "data/cellranger_out"
sample.list <- list.files(DATA_PATH)


# ---- READ FILES ---
seurat.list <- list()
for (sample.name in sample.list) {
  cellranger.data <- Read10X_h5(file.path(
    DATA_PATH,
    sample.name,
    "filtered_feature_bc_matrix.h5"
  ))
  seurat.list[sample.name] <- CreateSeuratObject(
    cellranger.data$`Gene Expression`,
    project = sample.name,
    min.cells = 20
  )
}


# ---- MERGE SEURAT OBJECTS ---
rna <- merge(
  seurat.list[[1]],
  seurat.list[2:length(seurat.list)],
  add.cell.ids = names(seurat.list),
  project = "LOY vs. ROY"
)
rna <- JoinLayers(rna)


# ---- MODIFY META DATA ----
rna[["condition"]] <- ifelse(
  startsWith(rna@meta.data$orig.ident, "LOY"),
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
rna@meta.data$sample.name <- lookup[rna$orig.ident]


# ---- ORDER OF META DATA ----
rna$orig.ident <- factor(
  rna$orig.ident,
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
rna$condition <- factor(rna$condition, levels = c("LOY", "ROY"))
rna$sample.name <- factor(
  rna$sample.name,
  levels = c('LOY1', 'LOY2', 'LOY3', 'LOY4', 'ROY1', 'ROY2', 'ROY3', 'ROY4')
)
Idents(rna) <- "sample.name"


# ---- SAVE MERGED SEURAT OBJECT ---
writeLines(Cells(rna), "analysis/cells_raw.txt")
qs_save(rna, file.path(ANALYSIS_PATH, "1_rna_imported.qs2"))
