#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(GOSemSim, args[1])
library(R.utils, args[1])
library(hash, args[1])

tmp = scan(args[2], what = "character", sep = "\n")

mygomap = hash()
for(i in 1:length(tmp)){
	thisSplit = strsplit(tmp[i], "\t")[[1]]
	mygomap[[thisSplit[1]]] = thisSplit[0:-1]
}


go2fun = read.table(args[3], header=FALSE, sep = "\t")
bpGO = as.character(go2fun[go2fun[,3]=="biological_process",1])
mfGO = as.character(go2fun[go2fun[,3]=="molecular_function",1])
ccGO = as.character(go2fun[go2fun[,3]=="cellular_component",1])

getSim<-function(a,b, ontology = "MF", method = "Wang"){
	score = matrix(nrow= length(a), ncol = length(b))
	rownames(score)  = a
	colnames(score)  = b
	for( i in 1:length(a)){
		for( j in 1:length(b)){
			tmp = goSim(a[i],b[j], ont= ontology, measure = method)
			if(typeof(tmp) == "double"){
				score[i,j] = tmp
			}
		}
	}
	score[is.na(score)] <- 0
	
	return( sum(c(rowMax(score),rowMax(t(score))))/(length(a)+length(b)))
#	return(mean(score))
}

getGO<- function(gene){
	return(mygomap[[gene]])
}

getSims<-function(interactions, ontology, valGOs, method = "Wang"){
        out = 1:dim(interactions)[1]
        for(i in 1:dim(interactions)[1]){
                idA = as.character(interactions[i,1])
                idB = as.character(interactions[i,2])
		goA = intersect(valGOs, getGO(idA))
		goB = intersect(valGOs, getGO(idB))
		if(length(goA) > 0 && length(goB) > 0 ){
			out[i] = getSim(goA, goB, ontology, method)
		}else{
			out[i] = NA
		}
	}
        return(out)
}


interactions = read.table(args[4], header=FALSE)

#setOntology("MF")
print("Processing MF")
mf = getSims(interactions, "MF", mfGO)
#setOntology("BP")
print("Processing BP")
bp = getSims(interactions, "BP", bpGO)
#setOntology("CC")
print("Processing CC")
cc = getSims(interactions, "CC", ccGO)
meanSim = 1:dim(interactions)[1]
for(i in 1:dim(interactions)[1]){
	meanSim[i] = mean(c(mf[i], cc[i], bp[i]), na.rm=TRUE)
}
out = cbind(interactions,mf,bp,cc, meanSim)
#rownames(out) = paste(interactions[,1],interactions[,2], sep = "\t")
out[is.na(out)] <- "?"

write.table(out, args[5], quote=FALSE, sep="\t",col.names=FALSE, row.names=FALSE)
