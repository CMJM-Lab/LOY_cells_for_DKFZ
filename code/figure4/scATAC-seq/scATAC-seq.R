.libPaths("path/to/Rlibraries/4.4.3")

library(ArchR)
library(dplyr)
library(Rcpp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(epiCHAOS)
library(ggplot2)
library(corrr)
library(S7)
library(readr)
library(biomaRt)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(tibble)
library(tidyr)
library(pheatmap)

################################################################ Load the data ################################################################
#--- for downstream analysis, the peaks-by-cells matrix is here:
mat <- readRDS("path/to/ArchRproject_rerun/peaks_matrix_rerun.Rds")


peakmat <- mat@assays@data$PeakMatrix
rownames(peakmat) <- mat@rowRanges %>% paste0()

dim(peakmat)
peakmat[1:10,1:10] # To see a summary

#--- and sample annotation is here:
mat@colData

proj <- loadArchRProject(path = "path/to/ArchRproject_rerun")

################################################################ Define the regions ################################################################

EMT_genes_GSEA_A549 <- read_csv("/path/to/20251031_EMT_genes_GSEA_A549.csv", col_names = FALSE)

# prerank_GSEA genes
prerank_GSEA <- as.character(EMT_genes_GSEA_A549[1, ])
prerank_GSEA <- prerank_GSEA[!is.na(prerank_GSEA)]
prerank_GSEA <- prerank_GSEA[-1]

# Download all the gene info from ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_biomart <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "chromosome_name", "start_position",
                 "end_position", "strand"),
  mart = mart
)

# Select only the genes of interest
genes_prerank_GSEA_biomart <- genes_biomart[genes_biomart$external_gene_name %in% prerank_GSEA, ]

# Convert to a GRanges object
genes_prerank_GSEA_GRanges <- GRanges(
  seqnames = paste0("chr", genes_prerank_GSEA_biomart$chromosome_name),
  ranges = IRanges(start = genes_prerank_GSEA_biomart$start_position,
                   end = genes_prerank_GSEA_biomart$end_position),
  strand = ifelse(genes_prerank_GSEA_biomart$strand == 1, "*", "*"),
  gene_id = genes_prerank_GSEA_biomart$ensembl_gene_id
)

mcols(genes_prerank_GSEA_GRanges)$external_gene_name <- genes_prerank_GSEA_biomart$external_gene_name


# Define the genome
genome(genes_prerank_GSEA_GRanges) <- "hg38"
genes_prerank_GSEA_GRanges


#########################################################################################################################################################################
################################################################ Examine the heterogeneity for promoters ################################################################
#########################################################################################################################################################################

promoters_prerank_GSEA <- promoters(
  genes_prerank_GSEA_GRanges,
  upstream = 2500,
  downstream = 2500
)

promoters_prerank_GSEA_visualization <- as.data.frame(promoters_prerank_GSEA)
write.csv(promoters_prerank_GSEA_visualization, "path/to/ArchRproject_rerun/promoters_prerank_GSEA_visualization.csv", row.names = FALSE)


proj <- addFeatureMatrix(
  proj,
  features = promoters_prerank_GSEA,
  matrixName = "CustomRegionsMatrixPromoters_GSEA"
)

################################################################ Compute the heterogeneity for promoter regions (GSEA_enrich) ################################################################

markersCustom <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "CustomRegionsMatrixPromoters_GSEA",
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

################################################################ Extract the gene names ################################################################

markersCustom <- rowData(markersCustom) # Convert rowData to a data.frame
markersCustom_GR <- GRanges(
    seqnames = as.character(markersCustom$seqnames),
    ranges = IRanges(start = markersCustom$start,
                     end   = markersCustom$end)
)

hits <- findOverlaps(promoters_prerank_GSEA, markersCustom_GR)

gene_names <- rep(NA, length(markersCustom_GR))
gene_names[subjectHits(hits)] <- mcols(promoters_prerank_GSEA)$external_gene_name[queryHits(hits)]
mcols(markersCustom_GR)$external_gene_name <- gene_names
markersCustom_GR

save(markersCustom_GR, file = "path/to/markersCustom_GR.RData")

#########################################################################################################################################################################
###################################################### Examine the heterogeneity for promoters using counts matrix ######################################################
#########################################################################################################################################################################

load("path/to/markersCustom_GR.RData")

#--- for downstream analysis, the peaks-by-cells matrix is here:
mat <- readRDS("path/to/ArchRproject_rerun/peaks_matrix_rerun.Rds")


peakmat <- mat@assays@data$PeakMatrix
rownames(peakmat) <- mat@rowRanges %>% paste0()

dim(peakmat)
peakmat[1:10,1:10] # To see a summary

#--- and sample annotation is here:
mat@colData

#########################################################################################################################################################################
###################################################### Apply Nima´s filtering ###########################################################################################
#########################################################################################################################################################################

# Load Nima's filtered cells
Nimas_filtering <- read.table("path/to/nima_multiome_plasticity/atac_filtered_cell_ids.txt", stringsAsFactors = FALSE)

# Replace _ by #
convert_id <- function(x) {
  sub("^([^_]+)_([^_]+)_", "\\1_\\2#", x)
}
Nimas_filtering[,1] <- convert_id(Nimas_filtering[,1])

# Select the cells in the peak matrix
filtered_mat <- mat[, rownames(mat@colData) %in% Nimas_filtering$V1]
length(rownames(filtered_mat@colData))

# Extract the GeneScoreMatrix without filtering
geneScoreMat <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix"
)
saveRDS(geneScoreMat, file = "path/to/geneScoreMat.rds")

geneScoreMat <- readRDS("path/to/geneScoreMat.rds")

# Filter the cells in the GeneScoreMatrix
filtered_geneScoreMat <- geneScoreMat[, 
  colnames(geneScoreMat) %in% Nimas_filtering$V1
]

saveRDS(filtered_geneScoreMat, file = "path/to/geneScoreMat_filteredNima.rds")

geneScoreMat <- readRDS("path/to/geneScoreMat_filteredNima.rds")

# Get the names of your genes
genes_of_interest <- markersCustom_GR$external_gene_name # 58 genes of interest

# Find the indices of the genes of interest in the matrix
rowGenes <- rowData(geneScoreMat)$name  
idx <- which(rowGenes %in% genes_of_interest) # 56 genes present
genes_not_present <- setdiff(genes_of_interest, rowGenes) # "CCN1", "CCN2" are absent


# Selct the genes of interest from the matrix
genes_matrix <- assay(geneScoreMat)[idx, ]  # rows = genes, columns = cells

# Get the cell metadata to group by "Sample"
meta <- getCellColData(proj, select = "Sample")
keep_cells <- rownames(meta) %in% Nimas_filtering$V1
meta <- meta[keep_cells, , drop = FALSE]

# Calculate average signal per group (bulk)
bulk_matrix <- t(apply(genes_matrix, 1, function(x) tapply(x, meta$Sample, mean)))
bulk_matrix

# Inspect the final matrix
head(bulk_matrix)
rownames(bulk_matrix) <- rowGenes[idx] # Asign the gene names

bulk_matrix_z <- t(apply(bulk_matrix, 1, function(x) (x - mean(x)) / sd(x))) # Z-score normalization

write.csv(
  bulk_matrix_z,
  file = "path/to/bulk_matrix_z.csv",
  quote = FALSE
)


#######################################################################################################################################################
###################################################### Plot ###########################################################################################
#######################################################################################################################################################

# Read the CSV
bulk_matrix_z <- read.csv(
  file = "path/to/bulk_matrix_z.csv",
  header = TRUE,      # the first row contains column names
  row.names = 1,      # the first column contains row names (genes)
  check.names = FALSE # keep column names exactly as they are
)

# Verify
dim(bulk_matrix_z)
head(bulk_matrix_z)

# Load pheatmap package
library(pheatmap)

# Plot the heatmap
pheatmap(
  bulk_matrix_z,
  scale = "none",                     # already Z-score normalized
  cluster_rows = FALSE,               # do not cluster genes
  cluster_cols = FALSE,               # do not cluster samples
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100), # blue = low, red = high
  fontsize_row = 6,
  fontsize_col = 8
)
