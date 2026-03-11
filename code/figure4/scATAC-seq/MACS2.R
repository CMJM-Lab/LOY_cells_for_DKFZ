#install.packages("data.table", lib = "path/to/Rlibraries/4.4.3", type = "source")

# Load the needed modules 

#module load GCCcore/14.1.0
#module load R/4.4.3-GCCcore-14.1.0

# Peak calling Script
.libPaths("/path/to/Rlibraries/4.4.3")

library(ArchR)
library(dplyr)
library(Rcpp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#--- load ArchR project
proj <- readRDS("/path/to/ArchR_rerun/ArchRproject_rerun/Save-ArchR-Project.rds")

#--- add group coverages for peak calling
proj <- addGroupCoverages(proj)

#--- calling peaks using macs2
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  pathToMacs2 = 'path/to/MINIFORGE_local/envs/EpiCHAOS/bin/macs2' # venv created using version python=2.7 (check MACS2 requirements)
)

#--- add a peaks matrix
proj <- addPeakMatrix(proj)

mat <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = T,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, file = "/path/to/ArchR_rerun/ArchRproject_rerun/peaks_matrix_rerun.Rds")

saveArchRProject(proj)