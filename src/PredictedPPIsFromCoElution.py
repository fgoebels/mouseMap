#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(scoreF, refF, elutionF, geneNameF, outF) = sys.argv[1:]
	
	geneNameFH = open(geneNameF)
	geneName= {}
	species = {}
	for line in geneNameFH:
		line = line.rstrip()
		ida, idb, spec = line.split("\t")
		if ida not in geneName: geneName[ida] = set([])
		geneName[ida].add(idb)
		species[ida] = spec
		species[idb] = spec 
	geneNameFH.close()

	toLearn, toPred = calcS.loadScoreData(scoreF, refF)
	
	rfc =  calcS.trainML(toLearn)

	ref, eluD, calc = calcS.loadData(refF, elutionF)
	
	calc.calculate2DScores(ref)
	outFH = open(outF + ".arff", "w")
	outFH.write(calc.toArffData())
	outFH.close()
	print "Calculated scores"
	
	
	data, targets = toPred.toSklearnData()
	dataL, targetsL = toLearn.toSklearnData()
	preds = rfc.predict(data)
	prots = []
	for protA, protB, label in toPred.scores:
		prots.append((protA, protB))

	outFH = open(outF, "w")
	for i in range(len(preds)):
		protA, protB = prots[i]
		if protA in geneName: geneA = ",".join(geneName[protA])
		if protB in geneName: geneB = ",".join(geneName[protB])
		spec = species[protA]
		if preds[i][1]>0.5:
			print >> outFH, "%s\t%s\t%s\t%s\t%s\t%f" % (protA, protB, geneA, geneB, spec, preds[i][1])
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
