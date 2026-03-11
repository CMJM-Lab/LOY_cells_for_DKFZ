library(dplyr)
library(pheatmap)
library(RnBeads)
library(ComplexHeatmap)
library(circlize)

## Samples as they appear in YOUR 6mA data (emt_all$sample)
samples <- c("D11", "E5b", "E6b", "F10", "A3", "C1", "C8", "D5")

## Map sample name -> file prefix on disk (note: E5B/E6B on disk)
file_ids <- c(
  D11 = "D11",
  E5b = "E5B",
  E6b = "E6B",
  F10 = "F10",
  A3  = "A3",
  C1  = "C1",
  C8  = "C8",
  D5  = "D5"
)

## Folder with *_withFrac.bed files
basedir <- "data/"

## ================================
## 1. Read promoter-level 6mA files
## ================================
read_promoter_file <- function(sample_name) {
  file_id <- file_ids[[sample_name]]
  path <- file.path(basedir, paste0(file_id, "_EMT_promoter_counts.bed"))
  message("Reading: ", path)
  
  df <- read.table(path, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "")
  colnames(df) <- c("chr", "start", "end", "gene", "Nmod", "Ncan")
  df[, 'mFrac'] <- as.numeric(df[, 'Nmod'])/(as.numeric(df[, 'Nmod'])+as.numeric(df[, 'Ncan']))
  
  # "." -> NA, then numeric
  df$Nmod[df$Nmod == "."]   <- NA
  df$Ncan[df$Ncan == "."]   <- NA
  df$mFrac[df$mFrac == "."] <- NA
  
  df$Nmod  <- as.numeric(df$Nmod)
  df$Ncan  <- as.numeric(df$Ncan)
  df$mFrac <- as.numeric(df$mFrac)
  
  df$sample <- sample_name
  df
}

emt_list <- lapply(samples, read_promoter_file)
emt_all  <- bind_rows(emt_list)

## ================================
## 2. Collapse to gene-level mFrac
## ================================
emt_gene_sample <- emt_all %>%
  dplyr::group_by(gene, sample) %>%
  dplyr::summarise(
    mFrac = mean(mFrac, na.rm=TRUE),
    n_promoters = dplyr::n(),
    .groups = "drop"
  )


## ================================
## 3. Subset to EMT genes only
## ================================
emt_gene_sample_emt <- emt_gene_sample[!is.na(emt_gene_sample$mFrac), ]
genes <- sort(unique(emt_gene_sample_emt$gene))


## ================================
## 4. Build gene x sample matrix (mFrac)
## ================================
mat_6mA <- matrix(NA_real_,
                  nrow = length(genes),
                  ncol = length(samples),
                  dimnames = list(genes, samples))

for (s in samples) {
  df_s <- emt_gene_sample_emt %>% dplyr::filter(sample == s)
  mat_6mA[df_s$gene, s] <- df_s$mFrac
}

## ================================
## 5. Z-score per gene (row-wise)
## ================================
mat_6mA_z <- t(scale(t(mat_6mA)))


## ================================
## 6. Sample annotation (ROY / LOY)
## ================================
rnbSet <- load.rnb.set("data/rnbSet_preprocessed")
pheno_info <- pheno(rnbSet)
group <- pheno_info$Group
names(group) <- pheno_info$Sample_Name  # D11, E5B, E6B, F10, A3, C1, C8, D5

# Use RnBeads sample names, but map to our 'samples' with E5b/E6b
pheno_names_for_group <- c("D11", "E5B", "E6B", "F10", "A3", "C1", "C8", "D5")

annotation_col <- data.frame(
  Group = group[pheno_names_for_group]
)
rownames(annotation_col) <- samples  # D11, E5b, E6b, F10, A3, C1, C8, D5

annotation_colors <- list(
  Group = c(
    ROY = "#ED7D31",   # orange
    LOY = "#0000C0"    # blue
  )
)

print(annotation_col)


###############################################################
## Export 6mA Z-score matrix + sample groups for Python heatmap
###############################################################

#py_outdir <- "/omics/groups/OE0219/internal/Gizem/Projects/Nanopore/A549_DiMeLo_All_clones/06_Results"
#dir.create(py_outdir, showWarnings = FALSE, recursive = TRUE)

# 2) Export sample-group mapping
sample_groups <- data.frame(
  Sample = samples,                           # D11, E5b, E6b, F10, A3, C1, C8, D5
  Group  = group[pheno_names_for_group]       # LOY / ROY
)


## ================================
## 7. Heatmap (on screen)
## ================================
heat_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(
  mat_6mA_z,
  color = heat_cols,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 7,
  show_rownames = TRUE,
  main = "EMT Promoter 6mA Methylation (Z-score)"
)


## ================================
## 8. Save PNG + PDF
## ================================
outdir <- "."

col_fun <- colorRamp2(
  breaks = c(min(mat_6mA_z), -1, 0, 1, max(mat_6mA_z)),
  colors = c("#224486", '#1414ff', "#FFFFFF", '#ff2b2b', "#862244")
)

Heatmap(mat_6mA_z,
    cluster_rows=TRUE,
    cluster_columns=FALSE,
    col=col_fun)