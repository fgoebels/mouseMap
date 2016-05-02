#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionFiles, refF, direction, outF) = sys.argv[1:]
	elutionFilesFH = open(elutionFiles)
	outData = {}
	maxSize = 0
	for line in elutionFilesFH:
		line = line.rstrip()
		reference, elutionData, scoreCalc = calcS.loadData(refF, line)
		scores = removeFracs(elutionData, reference, scoreCalc, direction)
		name = line.split("Ce_")[2].split(".")[0]
		outData[name] = scores
		maxSize = max(len(scores), maxSize)
	elutionFilesFH.close()
	
	outFH = open(outF, "w")
	print >> outFH, "Experiment_name\tFraction_%s" % ("\tFraction_".join(map(str,range(1, maxSize+1))))
	for dataset in outData:
		scores = outData[dataset]
		numFracs = len(scores)
		outline = "%s\t%s"  % (dataset, "\t".join(map(str, scores)))
		if maxSize-numFracs > 0:
			outline = "%s\t%s" % (outline, "\t".join(["NA"]*(maxSize-numFracs)))
		print >> outFH, outline
	outFH.close()

def getFracsToRemove(elutionData, direction):
	numFracs = elutionData.elutionMat.shape[1]
	if direction == "left":
			return range(0,numFracs)
	elif direction == "right":
			return list(reversed(range(0,numFracs)))
	else:
		print "Direction needs to be right or left"
		sys.exit()

def removeFracs(elutionData, reference, scoreCalc, direction):
	out = []
	fracsToRemove = getFracsToRemove(elutionData, direction)
	tmpFracs = copy.copy(fracsToRemove)
	for frac in fracsToRemove:
		if len(tmpFracs) == 0: continue
		print tmpFracs
		tmpElution = copy.copy(elutionData)
		tmpElution.getSubset(tmpFracs)
		scoreCalc = calcS.CalculateCoElutionScores(tmpElution)
		scoreCalc.calculateAllScores([calcS.Euclidiean()], reference)
		data, targets = scoreCalc.toSklearnData()
		clf = calcS.CLF_Wrapper(data, targets)
		scores =  clf.getValScores()
		out.append(scores[1])
		if frac in tmpFracs:
			tmpFracs.remove(frac)
	return out

def plotHeatmap( mat, name, outDir):
	outFH = open("%s.dat" % outDir, "w")
        print >> outFH, "\t" + "\t".join(map(str, range(50)))
	for i in range(49):
		print >> outFH, "%i\t%s" % (i, "\t".join(map(str, mat[i,:])))
	outFH.close()

	cmd = "src/plotHeatmap.R %s.dat \"%s\" %s.pdf" % (outDir, name, outDir)
	os.system(cmd)

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
