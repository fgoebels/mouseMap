#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os

def main():
	(elutionF, outD) = sys.argv[1:]
	print elutionF
	elutionData = calcS.ElutionData(elutionF)
	
	mat =  elutionData.elutionMat

	peptideCounts = []
	uniqProts = []
	for i in range(mat.shape[1]):
		col = mat[:,i]
		uniqProts.append(np.count_nonzero(col))
		peptideCounts.append(np.sum(col))

	peptideCounts = np.array(peptideCounts)
	uniqProts = np.array(uniqProts)
	peptideCounts = peptideCounts.reshape(5,48)
	uniqProts = uniqProts.reshape(5,48)

	plotHeatmap(peptideCounts, "Peptide counts", "%s.pcounts" % outD)

	plotHeatmap(uniqProts, "Unique Prots", "%s.puniq" % outD)

def plotHeatmap( mat, name, outDir):
	outFH = open("%s.dat" % outDir, "w")
        print >> outFH, "\t" + "\t".join(map(str, range(1,49)))
	for i in range(1, 6):
		print >> outFH, "%i\t%s" % (i, "\t".join(map(str, mat[i-1,:])))
	outFH.close()

	cmd = "src/plotHeatmap.R %s.dat \"%s\" %s.pdf" % (outDir, name, outDir)
	os.system(cmd)

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
