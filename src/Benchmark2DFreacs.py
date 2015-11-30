#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionF, refF, outD) = sys.argv[1:]
	reference, elutionData, scoreCalc = calcS.loadData(refF, elutionF)
	iexFractions = range(1, 49)
	
	out = np.array([[-1.00]*48]*48)
	for removeLeft in range(1, 49):
		tmpFracs = copy.copy(iexFractions)
		for i in range(1, removeLeft):
			if i in tmpFracs: tmpFracs.remove(i)
		for removeRight in reversed(range(removeLeft+1, 50)):
			if removeRight in tmpFracs: tmpFracs.remove(removeRight)
			print tmpFracs
			fractions = getIEXFracs(tmpFracs)
			tmpElution = copy.copy(elutionData)
			tmpElution.getSubset(fractions)
			scoreCalc = calcS.CalculateCoElutionScores(tmpElution)
			scoreCalc.calculateAllScores([calcS.Euclidiean()], reference)
			data, targets = scoreCalc.toSklearnData()
			clf = calcS.RandomForest(data, targets)
			scores =  clf.getValScores()
			out[removeLeft-1][49-removeRight] = scores[1]
			print "%i\t%i\t%.2f" % (removeLeft-1, 49-removeRight, scores[1])

	outFH = open(outD + ".iex.dat", "w")
	print >> outFH, "\t" + "\t".join(map(str, range(48)))
	for i in range(48):
		print >> outFH, "%i\t%s" % (i, "\t".join(map('{0:.2f}'.format, out[i])))
	outFH.close()


def getIEXFracs(IEXfracs):
	out = []
	for IEXfrac in IEXfracs:
		twoDfracs = [(1, IEXfrac),(2, IEXfrac),(3, IEXfrac),(4, IEXfrac),(5, IEXfrac)]
		for frac in twoDfracs:
			out.append(getIntFromFrac(frac))
	return out

def getIntFromFrac(numFrac):
        return int(math.floor(numFrac[0]-1)*48) + numFrac[1] - 1

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
