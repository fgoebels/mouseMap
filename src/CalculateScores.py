#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy, inspect

def main():
	(elutionF, ontoF, gene_ass_F, ppiF, outF) = sys.argv[1:]
	
	elutionData, toPred = calcS.loadEData(elutionF)
	
	toPredPPIs = []
	ppiFH = open(ppiF)
	for line in ppiFH:
		line = line.rstrip()
		protA, protB = line.split("\t")
		toPredPPIs.append((protA, protB, "?"))

	toPred.calculate2DScores(toPredPPIs)
	out = toPred.toTable(False)

	outFH = open(outF, "w")
	print >> outFH, out
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
