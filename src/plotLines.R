#!/usr/bin/Rscript
library(lattice)
library(akima)

args <- commandArgs(trailingOnly = TRUE)


dataTable = read.table(args[1])


pdf(args[4])
plot(dataTable[,1], dataTable[,4], "l", ylim = c(0,as.integer(args[3])), xlab="Number of fractions", ylab=args[2])
lines(dataTable[,1], dataTable[,2], lty=2)
lines(dataTable[,1], dataTable[,3], lty=3)
legend("topleft", c("max", "mean", "min"), lty=c(2,1,3))
dev.off()
