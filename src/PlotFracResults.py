#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils 
import sys
import os
import math

def main():
	(resultFileF, outDir) = sys.argv[1:]

	evalScores = [(0, "Number of interactions with >.5"), (2, "Precision"), (3, "Recall"), (4, "F-measure"), (5, "auc_pr"), (6, "aurc_roc")]
	resultFileFH = open(resultFileF)
	results = {}
	for line in resultFileFH:
		line = line.rstrip()
		lineSplit = line.split("\t")
		numFracs = int(lineSplit[0])
		if numFracs not in results: results[numFracs] = []
		tmp = []
		tmp.append(int(lineSplit[1]))
		tmp.append(lineSplit[2])
		for i in range(3, len(lineSplit)):
			tmp.append(float(lineSplit[i]))
		results[numFracs].append(tmp)


	for numF in results:
		results[numF] = np.array(results[numF])

	eval(outDir + "lineplot" , results, evalScores)

	getFracCounts(outDir + ".heatmap.upper", results[100], evalScores)
	getFracCounts(outDir + ".heatmap.lower", results[100], evalScores, False)

def getFracCounts(outDir, mat, scores, upper = True):
	counts = {}
	for (rowNum, name) in scores:
		quantileLower, quantileUpper = np.percentile(map(float, mat[:,rowNum]), [25,75], axis = 0)
		for row in mat:
			fracs = map(getFracFromInt, map(int, row[1].split(",")))
			score = float(row[rowNum])
			for frac in fracs:
				if frac not in counts: counts[frac] = 0
				if upper:
					if score >= quantileUpper:  counts[frac] += 1
				else:
					if score <= quantileLower: counts[frac] += 1

		heatMap = np.matrix([[0]*48]*5)

		for frac in counts:
			counts[frac] = counts[frac]
		
				
		for (ief, iex) in sorted(counts, key = counts.get, reverse=True):
#			print "%i\t%i\t%.4f" % (ief, iex, counts[(ief, iex)])
			heatMap[ief-1, iex-1] = counts[(ief, iex)]

		heatMap = normalize_fracs(heatMap, norm_rows=True, norm_cols=True)
		outFH = open("%s.tmp.dat" % outDir, "w")
		print >> outFH, "\t" + "\t".join(map(str, range(1,49)))
		for i in range(1, 6):
#			line = "%i" % (i)
#			for j in range(1, 49):
#				line +=  "\t%i" % (counts[(i, j)]) 
#			print line
			print >> outFH, "%i\t%s" % (i, "\t".join(map(str, heatMap[i-1,:])))
		outFH.close()
		cmd = "src/plotHeatmap.R %s.tmp.dat \"%s\" %s.%i.pdf" % (outDir, name, outDir, rowNum)
		os.system(cmd)
	
def arr_norm(arr, axis=0):
    """
    axis=0: normalize each column; 1: normalize each row
    """
    mat = np.asmatrix(arr)
    return np.asarray(np.nan_to_num(mat / np.sum(mat, axis)))

def normalize_fracs(arr, norm_rows=True, norm_cols=True):
    if norm_cols:
        # Normalize columns first--seems correct for overall elution profile if
        # you're calculating correlation-type scores
        arr = arr_norm(arr, 0)
    if norm_rows:
        arr = arr_norm(arr, 1)
    return arr


def getFracFromInt(numFrac):
	return (int(math.floor(numFrac/48))+1, numFrac % 48+1)

def getIntFromFrac(numFrac):
	return int(math.floor(numFrac[0]-1)*48) + numFrac[1] - 1


def eval(outDir, results, evalScores):
	for (rowNum, name) in evalScores:
		resultEval, maxVal = getVals(results, rowNum,  float)
		tmpFH = open(outDir + ".tmp.dat", "w")
		for frac in resultEval:
			print >> tmpFH, "%i\t%s" % (frac, "\t".join(map(str, resultEval[frac])))
		tmpFH.close()
		if maxVal < 1: maxVal = 1
		cmd = "src/plotLines.R %s.tmp.dat \"%s\" %i %s.%i.pdf" % (outDir, name, maxVal, outDir, rowNum)
		os.system(cmd)

def getVals(res, row, typeCast):
	out = {}
	maxVal = 0
	for fracNum in res:
		out[fracNum] = []
		tmp = map( typeCast, res[fracNum][:,row])
		if max(tmp) > maxVal: maxVal = max(tmp)
		out[fracNum].append(max(tmp))
		out[fracNum].append(min(tmp))
		out[fracNum].append(np.mean(tmp))
		out[fracNum].append(np.std(tmp))
	return out, maxVal


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
