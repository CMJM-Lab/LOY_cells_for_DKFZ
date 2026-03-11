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
dmrs <- readRDS('data/dmr_gene_overlap.RDS')
dmrs <- as.data.frame(dmrs)
colnames(dmrs)[1:3] <- c('Chromosome', 'Start', 'End')
dmrs <- dmrs[!duplicated(dmrs$symbol),]
rownames(dmrs) <- dmrs$symbol
rnb.set.annotation(type='DMRs', regions=dmrs[, 1:3], assembly='hg38')
rnbSet <- summarize.regions(rnbSet, region.type='DMRs')

to_plot <- meth(rnbSet, 'DMRs', row.names=TRUE)
desired_order <- c("D11", "E5B", "E6B", "F10", "A3", "C1", "C8", "D5")
to_plot <- to_plot[!(grepl('AS1', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('RNA', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('RNU', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('LINC', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('RP[[:digit:]]', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('MIR[[:digit:]]', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('AC[[:digit:]]', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('RMST', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('Metazoa', row.names(to_plot))), ]
to_plot <- to_plot[!(grepl('snoMe', row.names(to_plot))), ]

to_plot <- to_plot[, desired_order]
top_100 <- abs(rowMeans(to_plot[,1:4])- rowMeans(to_plot[, 5:8]))
Heatmap(na.omit(to_plot[order(top_100, decreasing=TRUE)[1:100], ]),
    cluster_rows=TRUE,
    cluster_columns=FALSE,
    col=rev(inferno(100)))

to_plot_boxplot <- data.frame(LOY=rowMeans(to_plot[, 1:4]), ROY=rowMeans(to_plot[, 5:8]))
to_plot_boxplot <- reshape2::melt(to_plot_boxplot)
ggplot(to_plot_boxplot, aes(x=variable, y=value, fill=variable))+geom_boxplot()+scale_fill_manual(values=color_samples)+theme(axis.text.x = element_text(angle=90, hjust=1), text = element_text(size=16, color = 'black'))+NoLegend()+plot_theme+
  xlab("")+ylab('DMR methylation')