coEluteData = read.table("test/Ming_HPLC_120_CELE.dat", header=T, sep = "\t")[,3:4]

coEluteList = list()
coEluteList$positive <- coEluteData[which(coEluteData[,2]=="positive"),1]
coEluteList$negative <- coEluteData[which(coEluteData[,2]=="negative"),1]
wilcox.test(coEluteList$positive, coEluteList$negative)
ggplot(coEluteData, aes(x=score, fill=Class)) + geom_density(alpha=.3) + xlab("Co-elution score")

summary(coEluteList$positive)
summary(coEluteList$negative)

pdf("test/Ming_co_elutehist.pdf")
ggplot(coEluteData, aes(x=score, fill=Class)) + geom_density(alpha=.3) + xlab("Co-elution score")
dev.off()
