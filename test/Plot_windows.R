#!/usr/local/bin/Rscript
library(HH)
library(lattice)
library(akima)

args <- commandArgs(trailingOnly = TRUE)

dataTable <- read.table(args[1], check.names=F, header=T, sep = "\t")
dataTable <- dataTable[which(dataTable[,2]>9),]
scores = colnames(dataTable)[-1:-3]

for(i in scores){
pdf(paste(c(args[2], ".", i, ".pdf"), collapse=""))
wireframe(i ~ colnames(dataTable)[3] + "Windowsize",data=dataTable
          scales = list(arrows=FALSE, col="black",font=10),
          xlab = list("Window size",rot=30),
          ylab = list(colnames(dataTable)[3],rot=-30),
          zlab = list(i,rot=90),
#         zlim = c(0,1),
#          auto.key=TRUE,
#          key=list(text=list(sapply(unique(protData[,4]), paste),col=mycolors),
#                   lines=list(lty=rep.int(1,numProts),col=mycolors[1:numProts])),
#          par.settings = list(axis.line = list(col = "transparent")),
)

#regr2.plot(dataTable[,2], xlab="Size",
#           dataTable[,3],  ylab=colnames(dataTable)[3],
#           dataTable[,i], zlab=i,
#           resid.plot=FALSE,
#           theta=140, phi=35, r=sqrt(15), ## used only in R
#           box=is.R(),
#           plot.back.planes=FALSE,
#           plot.base.points=FALSE,
#           main="Least-squares with two X-variables")
dev.off()
}

library(rgl)
plot3d(dataTable[,2], xlab="Size", dataTable[,3],  ylab=colnames(dataTable)[3], dataTable[,i], zlab=i)
x <- dataTable[,2]
y <- dataTable[,3]
z <- dataTable[,i]
fit <- lm(z ~ x + y)
coefs <- coef(fit)
a <- coefs["x"]
b <- coefs["y"]
c <- -1
d <- coefs["(Intercept)"]
planes3d(a, b, c, d, alpha=0.5)
#rgl.postscript(paste(c(args[2], ".", i, ".pdf"), collapse=""), "pdf")
