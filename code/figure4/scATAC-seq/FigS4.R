library(ArchR)

# Supplementary Figure S4 (G–J)

######################## Do the filtering ########################

Nimas_filtering <- read.table("path/to/atac_filtered_cell_ids.txt", stringsAsFactors = FALSE) # Load Nima's filtered cells

convert_id <- function(x) {sub("^([^_]+)_([^_]+)_", "\\1_\\2#", x)} # Replace _ by #

Nimas_filtering[,1] <- convert_id(Nimas_filtering[,1])

cells_proj <- getCellNames(proj)
keep_cells <- cells_proj %in% Nimas_filtering$V1

# Create a filtered project for the plots
proj_filtered <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells_proj[keep_cells],
  outputDirectory = "/path/to/ArchR_rerun/ArchRproject_filtered",
  dropCells = TRUE,
  force = TRUE
)

######################## Figures for the article

# CDH2
proj <- loadArchRProject(path = "/path/to/ArchR_rerun/ArchRproject_filtered")

pal_samples <- c(
  "WT_A3"   = "#ED7D31",
  "WT_C1"   = "#ED7D31",
  "WT_D5"   = "#ED7D31",
  "WT_C8"   = "#ED7D31",
  "LOY_D11" = "#0000C0",
  "LOY_F10" = "#0000C0",
  "LOY_E6b" = "#0000C0",
  "LOY_E5b" = "#0000C0"
)

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Sample", 
    geneSymbol = c("CDH2"), 
    upstream = 240000,
    downstream = 40000,
    pal = pal_samples
)


grid::grid.newpage()
grid::grid.draw(p$CDH2)


g <- "CDH2"

message("Saving PDF for: ", g)

pdf(file = paste0("/path/to/ArchR_rerun/Tracksplot_rerun/", g, "_Fullgene_browserTrack_cols.pdf"), width = 8, height = 6)
grid::grid.newpage()
grid::grid.draw(p[[g]])
dev.off()



# Rest of the genes
genes <- c(
  "THY1", "LOX", "ANPEP"
)

pal_samples <- c(
  "WT_A3"   = "#ED7D31",
  "WT_C1"   = "#ED7D31",
  "WT_D5"   = "#ED7D31",
  "WT_C8"   = "#ED7D31",
  "LOY_D11" = "#0000C0",
  "LOY_F10" = "#0000C0",
  "LOY_E6b" = "#0000C0",
  "LOY_E5b" = "#0000C0"
)

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Sample", 
    geneSymbol = genes, 
    upstream = 40000,
    downstream = 40000,
    pal = pal_samples
)

for(g in genes){
  message("Plotting: ", g)
  grid::grid.newpage()
  grid::grid.draw(p[[g]])
}

for(g in genes){
  message("Saving PDF for: ", g)

  pdf(file = paste0("/path/to/ArchR_rerun/Tracksplot_rerun/", g, "_browserTrack_cols.pdf"), width = 8, height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[g]])
  dev.off()
}