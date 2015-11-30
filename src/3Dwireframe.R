#!/usr/bin/Rscript
library(lattice)
library(akima)

args <- commandArgs(trailingOnly = TRUE)

protData <- read.table(args[1], check.names=F, header=T, sep = "\t")

mycolors.trans = rgb(c(255,0,0,255,255), c(0,255,0,165,255), c(0,0,255,0,0), alpha = 70, maxColorValue = 255)

mycolors = rgb(c(255,0,0,255,255), c(0,255,0,165,255), c(0,0,255,0,0), maxColorValue = 255)

numProts = length(sapply(unique(protData[,4]), paste))

pdf(args[2])
wireframe(Counts~IEF*IEX,data=protData,group=ProtName,
          col.groups=mycolors.trans,
          scales = list(arrows=FALSE, col="black",font=10),
          xlab = list("IEF",rot=30),
          ylab = list("IEX-HPLC",rot=-30),
          zlab = list("Counts",rot=90),
          zlim = c(0,1),
          #auto.key=TRUE,
          key=list(text=list(sapply(unique(protData[,4]), paste),col=mycolors),
                   lines=list(lty=rep.int(1,numProts),col=mycolors[1:numProts])),
          par.settings = list(axis.line = list(col = "transparent")),
)
dev.off()
