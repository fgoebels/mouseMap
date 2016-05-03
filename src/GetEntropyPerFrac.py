#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionFiles, refF, windowSize, outF) = sys.argv[1:]
	elutionFilesFH = open(elutionFiles)
	windowSize = int(windowSize)
	outData = ['']*3
	for line in elutionFilesFH:
		line = line.rstrip()
		print "processing %s" % (line)
		reference, elutionData, scoreCalc = calcS.loadData(refF, line)
		j = 0
		name = line.split("Ce_")[1].split(".")[0]
		for resultScore in getFracEvals(elutionData.elutionMat):
			data_lines = entropyVSprecision(elutionData, reference, resultScore, windowSize)
			for i in range(len(data_lines)):
				outData[j] += "\n%s\t%i\t%s" % (name, windowSize, data_lines[i])
			j += 1
	elutionFilesFH.close()

	printTable("%s_Entropy_%i.dat" % (outF, windowSize), "Entropy", outData[0])
	printTable("%s_Prot-prob_%i.dat" % (outF, windowSize), "Prot-prob", outData[1])
	printTable("%s_Num-prots_%i.dat" % (outF, windowSize), "Num-prots", outData[2])

def printTable(outF, Scorename, lines):
	outFH = open(outF, "w")
	print >> outFH, "Source\tWindowsize\t%s\tPrecision\tRecall\tFmeasure\tAuc_pr\tAuc_roc%s" % (Scorename, lines)
	outFH.close()
	

#########################################################
######	Code for calculating entropy based heatmap ######
#	
#	outFH = open(outF + ".dat", "w")
#	print >> outFH, "Experiment_name\t%s" % ("\t".join(map(str,range(1, maxSize+1))))
#	for dataset in outData:
#		scores = outData[dataset]
#		numFracs = len(scores)
#		outline = "%s\t%s"  % (dataset, "\t".join(map(str, scores)))
#		if maxSize-numFracs > 0:
#			outline = "%s\t%s" % (outline, "\t".join(["NA"]*(maxSize-numFracs)))
#		print >> outFH, outline
#	outFH.close()

def entropyVSprecision(elutionData, reference, entopies, windowSize=20):
	out = []
	sortedFracs = np.array([i[0] for i in sorted(enumerate(entopies), key=lambda x:x[1], reverse=False)])
	entopies = np.array(entopies)
	numFracs = len(sortedFracs)
	if windowSize > numFracs: return []
	for i in range(numFracs-windowSize+1):
		thisWindow = range(i, i+windowSize)
		selectedFracs = sortedFracs[thisWindow]
		averageEntropy = np.mean(entopies[selectedFracs])
		tmpElution = copy.copy(elutionData)
		tmpElution.getSubset(selectedFracs)
		scoreCalc = calcS.CalculateCoElutionScores(tmpElution)
#		scoreCalc.calculateAllScores([calcS.Euclidiean(), calcS.Jaccard(), calcS.Pearson(), calcS.Wcc()], reference)
		scoreCalc.calculateAllScores([calcS.Jaccard()], reference)
		data, targets = scoreCalc.toSklearnData()
		clf = calcS.CLF_Wrapper(data, targets, True)
		scores =  "\t".join(map( str, clf.getValScores()))
		out.append("%f\t%s" % (averageEntropy, scores))
	return out

def getFracEvals(elutionMat, cutoff=2):
	entropies = []
	p_vals = []
	num_prots = []
	numFracs = elutionMat.shape[1]
	for frac in range(numFracs):
		fracCounts = elutionMat[:,frac]
		fracs_with_prot = set([i for i,v in enumerate(fracCounts) if v > cutoff])
		num_prots.append(len(fracs_with_prot))
		p_has_prot = len(fracs_with_prot)/len(fracCounts)
		p_vals.append(p_has_prot)
		p_not_has_prot = 1 - p_has_prot 
		if p_has_prot == 0 or p_has_prot == 1:
			entropies.append(0)
			continue
		entropy = -( p_has_prot*math.log(p_has_prot,2) + p_not_has_prot*math.log(p_not_has_prot,2))
#		fracProbs = fracCounts/np.nansum(fracCounts)
#		entropy = -np.nansum(fracProbs*np.log(fracProbs))/np.log(len(fracProbs))
		entropies.append(entropy)
	return entropies, p_vals, num_prots

def arrayToString(array):
	return "\t".join(map(str, array))

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
