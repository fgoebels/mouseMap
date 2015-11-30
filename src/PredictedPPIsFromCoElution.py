#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionF, refF, geneNameF, outF) = sys.argv[1:]
	
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
	reference, elutionData, scoreCalc = calcS.loadData(refF, elutionF)
	reference, elutionData, toPred = calcS.loadData(refF, elutionF)
#	toPred = copy.copy(scoreCalc)
	rfc =  calcS.trainML(reference, scoreCalc)

	toPred.calculateAllPairs()
	data, targets = toPred.toSklearnData()
	print "Calculated scores"
	preds = rfc.predict(data)
	prots = []
	for protA, protB, label in toPred.scores:
		prots.append((protA, protB))

	outFH = open(outF, "w")
	for i in range(len(preds)):
		protA, protB = prots[i]
		if protA not in species or protB not in species: continue
		if species[protA] != species[protB]: continue
		geneA = ""
		geneB = ""
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
