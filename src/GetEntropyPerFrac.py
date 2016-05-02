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
	outData = {}
	maxSize = 0
	windowNamesFH = open(outF + "window-names.txt", "w")
	for line in elutionFilesFH:
		line = line.rstrip()
		print "processing %s" % (line)
		reference, elutionData, scoreCalc = calcS.loadData(refF, line)
		scores = getEntropy(elutionData.elutionMat)
		name = line.split("Ce_")[2].split(".")[0]
		outData[name] = scores
		maxSize = max(len(scores), maxSize)
		line_data = entropyVSprecision(elutionData, reference, scores, windowSize)
		fileName = "%s_%s_%s.dat" % (outF, windowSize, name)
		print >> windowNamesFH, "%s\t%s" % (name, fileName)
		tmp_outFH = open(fileName, "w")
		print >> tmp_outFH, line_data
		tmp_outFH.close()
	elutionFilesFH.close()
	windowNamesFH.close()
	
	outFH = open(outF + ".dat", "w")
	print >> outFH, "Experiment_name\t%s" % ("\t".join(map(str,range(1, maxSize+1))))
	for dataset in outData:
		scores = outData[dataset]
		numFracs = len(scores)
		outline = "%s\t%s"  % (dataset, "\t".join(map(str, scores)))
		if maxSize-numFracs > 0:
			outline = "%s\t%s" % (outline, "\t".join(["NA"]*(maxSize-numFracs)))
		print >> outFH, outline
	outFH.close()

def entropyVSprecision(elutionData, reference, entopies, windowSize=20):
	out = []
	sortedFracs = np.array([i[0] for i in sorted(enumerate(entopies), key=lambda x:x[1], reverse=False)])
	entopies = np.array(entopies)
	numFracs = len(sortedFracs)
	for i in range(numFracs-windowSize+1):
		thisWindow = range(i, i+windowSize)
		selectedFracs = sortedFracs[thisWindow]
		averageEntropy = np.mean(entopies[selectedFracs])
		tmpElution = copy.copy(elutionData)
		tmpElution.getSubset(selectedFracs)
		scoreCalc = calcS.CalculateCoElutionScores(tmpElution)
		scoreCalc.calculateAllScores([calcS.Euclidiean()], reference)
		data, targets = scoreCalc.toSklearnData()
		clf = calcS.CLF_Wrapper(data, targets)
		scores =  clf.getValScores()
		precision = scores[0]
		out.append("%f\t%f" % (averageEntropy, precision))
	return "\n".join(out)

def getEntropy(elutionMat, cutoff=2):
	out = []
	numFracs = elutionMat.shape[1]
	for frac in range(numFracs):
		fracCounts = elutionMat[:,frac]
		fracs_with_prot = set([i for i,v in enumerate(fracCounts) if v > cutoff])
		p_has_prot = len(fracs_with_prot)/len(fracCounts)
		p_not_has_prot = 1 - p_has_prot 
		if p_has_prot == 0 or p_has_prot == 1:
			out.append(0)
			continue
		entropy = -( p_has_prot*math.log(p_has_prot,2) + p_not_has_prot*math.log(p_not_has_prot,2))
#		fracProbs = fracCounts/np.nansum(fracCounts)
#		entropy = -np.nansum(fracProbs*np.log(fracProbs))/np.log(len(fracProbs))
		out.append(entropy)
	return out

def arrayToString(array):
	return "\t".join(map(str, array))

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
