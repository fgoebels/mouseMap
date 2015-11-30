library(parallel)
library(MASS)

coelutionMatN <- function(elutionMatrix, n){
	out = mclapply(1:n, function(x){print(x); calculateCoElutionMatrix(addNoiseToMat(elutionMatrix))})
	return(Reduce('+', out)/n)
}

calculateCoElutionMatrix <- function(elutionMatrix){
	numberOfGenes = length(rownames(elutionMatrix))
	out = matrix(NA, numberOfGenes, numberOfGenes)
	rownames(out) = rownames(elutionMatrix)
	colnames(out) = rownames(elutionMatrix)
	for(i in 1:(numberOfGenes-1)){
		matLine = c(rep(0,i-1), 1, getCorelation(elutionMatrix[i,], elutionMatrix[-1:-i,]))
		out[i,] = matLine
		
	}
	out[numberOfGenes,] = c(rep(0,numberOfGenes-1), 1)
	return(out)
}

getCorelation <- function(gene, geneMat){
	tmpMat = geneMat
	if(length(rownames(tmpMat)) != 0){
		tmpMat = t(tmpMat)
	}
        return(cor(gene, tmpMat, method="pearson"))
	
}

getCorelationWithNoise <- function(gene, geneMat, M){
	runs = 1000
	size = length(rownames(geneMat))
	if(size ==0){
		size = size + 1
	}
	out = rep(0, size)
	for(i in 1:runs){
		tmpGene = addNoise(gene, M)
		tmpGene = tmpGene/sum(tmpGene)
		if(size != 1){
			tmpMat = apply(geneMat, 1, addNoise, M)
			tmpMat = apply(tmpMat, 2, normalize)
		}
		else{
			tmpMat = addNoise(geneMat, M)
			tmpMat = tmpMat/sum(tmpMat)
		}
			out = out + cor(tmpGene, tmpMat, method="pearson")
	}
	return(out/runs)
}

addNoiseToMat <- function(matrix){
	modMat = apply(matrix, 2, addNoise, 1/length(colnames(matrix)))
        modMat = apply(modMat, 2, normalize)
	return(modMat)
}

addNoise <- function(gene, M){
		return(gene + rpois(length(gene), fitdistr(gene, "poisson")$estimate) + M)
}

normalize <- function(x){
	return(x/sum(x))
}
