library(RnBeads)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(viridis)
library(ggplot2)
library(Seurat)
color_samples <- c('LOY'='#0000C0',
                'ROY'='#ED7D31')
plot_theme <- theme(text=element_text(color='black', size=10),
                    axis.ticks = element_line(color='black'),
                    axis.text = element_text(color='black', size=10),
                    axis.line.x.top = element_blank(),
                    axis.line.y.right = element_blank(),
                    panel.background=element_blank(),
                    axis.line.x.bottom = element_line(),# draw bottom
                    axis.line.y.left   = element_line())
rnbSet <- load.rnb.set("data/rnbSet_preprocessed")
diff_meth_cpgs <- read.csv('data/diffMethTable_region_cmp1_cpgislands.csv')
edb <- EnsDb.Hsapiens.v86
gencode_genes <- genes(edb, columns = c("gene_name", "gene_biotype"))
seqlevelsStyle(gencode_genes) <- "UCSC"
gencode_genes <- gencode_genes[!(grepl('AS1', gencode_genes$gene_name)), ]
gencode_genes <- gencode_genes[!(grepl('RNA', gencode_genes$gene_name)), ]
gencode_genes <- gencode_genes[!(grepl('RNU', gencode_genes$gene_name)), ]
gencode_genes <- gencode_genes[!(grepl('LINC', gencode_genes$gene_name)), ]
gencode_genes <- gencode_genes[!(grepl('RP[[:digit:]]', gencode_genes$gene_name)), ]
gencode_genes <- gencode_genes[!(grepl('MIR[[:digit:]]', gencode_genes$gene_name)), ]

top_cpgs <- diff_meth_cpgs[order(diff_meth_cpgs$combinedRank), ]
top_cpgs <- top_cpgs[!(top_cpgs$Chromosome%in%c('chrX', 'chrY')), ]
top_cpgs_gr <- makeGRangesFromDataFrame(top_cpgs)
tss <- resize(gencode_genes, width=1, fix='start')
closest_gene <- sapply(1:nrow(top_cpgs), function(x){
  distances <- distance(top_cpgs_gr[x], tss)
  if(min(distances, na.rm=TRUE)<10000){
    tss$gene_name[which.min(distances)]
  }else{
    NA
  }
})
top_cpgs <- top_cpgs[!is.na(closest_gene), ]
top_cpgs$closest_gene <- closest_gene[!is.na(closest_gene)]

EMT_gsea_list <- c(
  "ANPEP","WNT5A","THY1","GREM1","SPOCK1","SPARC","CDH2","ITGB3","GLIPR1","EDIL3",
  "SCG2","FBN1","SNAI2","COL5A2","BDNF","FZD8","COL5A1","FBN2","FN1","PLOD2",
  "CCN1","SERPINE2","SERPINE1","TNC","IGFBP4","CXCL1","LOX","COL8A2","TAGLN",
  "EFEMP2","CCN2","GJA1","MEST","SLIT3","PMP22","PCOLCE2","CXCL8","SPP1",
  "COL7A1","MATN3","CD59","DPYSL3","COL16A1","FSTL1","FERMT2","NNMT","LAMA3",
  "CADM1","SERPINH1","TPM1","DST","SFRP1","COL6A3","DAB2","SAT1","GADD45A","TPM4","CALD1"
)

#selected_top <- top_cpgs[top_cpgs$closest_gene%in%EMT_gsea_list, ]
selected_top <- top_cpgs
selected_top <- selected_top[which(abs(selected_top$mean.mean.diff)>0.2), ]
highest_diff <- sapply(unique(selected_top$closest_gene), function(gene){
  sel_data <- selected_top[(selected_top$closest_gene==gene), ]
  sel_data$id[which.max(abs(sel_data$mean.mean.diff))]
})
selected_top <- selected_top[as.character(selected_top$id)%in%as.character(highest_diff), ]
cpg_data <- meth(rnbSet, type='cpgislands', row.names=TRUE)
to_plot <- cpg_data[as.character(selected_top$id), ]
desired_order <- c("D11", "E5B", "E6B", "F10", "A3", "C1", "C8", "D5")

write.table(selected_top[, c('Chromosome', 'Start', 'End', 'closest_gene')], 'data/EMT_CGI.bed', sep='\t',quote=FALSE, row.names=FALSE, col.names=FALSE)