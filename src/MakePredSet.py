#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionF, geneNameF, outPrefix) = sys.argv[1:]
	
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

	elutionData, scoreCalc = calcS.loadEData(elutionF)
	preds = scoreCalc.getAllPairs()
	out = {}
	for protA, protB, _ in preds:
		if protA not in species or protB not in species: continue
		if species[protA] != species[protB]: continue
		if species[protA] not in out: out[species[protA]] = set()
		out[species[protA]].add("\t".join(sorted([protA, protB])))	

	for species in out:
		outFH = open("%s.%s.topred.txt" % (outPrefix, species) , "w")
		print >> outFH, "ProtA\tProtB"
		print >> outFH, "\n".join(out[species])
		outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
