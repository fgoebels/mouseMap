#!/usr/bin/Rscript
library(gplots)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

my_palette <- colorRampPalette(c("blue",  "red"))(n = 299)
dataTable = read.table(args[1], row.names=1, header=T, check.names=F)
tmp = t(t(dataTable))

pdf(args[3])
heatmap.2(tmp, dendrogram="none", col=my_palette, margins =c(2,2), trace="none", density.info="none", Colv = "NA", Rowv = "NA", main = "", xlab=args[2], ylab="Experiments")
dev.off()
