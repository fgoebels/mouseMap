#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

dataTable <- read.table(args[1], check.names=F, header=T, row.names=1, sep = "\t")


colors = rainbow(dim(dataTable[1])) #c("yellow", "cyan", "blue", "green", "red")

#dataTable <- t(t(dataTable)/rowSums(t(dataTable)))

pdf(args[2])
par(mar=c(5.1,5.1,5.1,2.1))
barplot(t(dataTable), beside=T, cex.lab = 2, cex.axis=2, cex.main=3, xlab = args[3], ylab=args[4], las =2)
#legend("top", colnames(dataTable), fill = colors, ncol=2)
dev.off()
