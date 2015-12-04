#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy, inspect

def main():
	(elutionF, ontoF, gene_ass_F, ppiF, scores, outF) = sys.argv[1:]
	print  os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
	scores = scores.split("-")
	
	elutionData, toPred = calcS.loadEData(elutionF)
	thisScores = []
	twoDscores = toPred.twoDscores
	for score in scores:
		print score
		if score == "MatrixNorms":
			thisScores.append(calcS.MatrixNorms())
			continue
		if score == "GOSim":
			thisScores.append(calcS.GOSim(ontoF, gene_ass_F))
			continue
		for twoDscore in twoDscores:
			if score == twoDscore.name: thisScores.append(twoDscore)

	print thisScores
	toPredPPIs = []
	ppiFH = open(ppiF)
	for line in ppiFH:
		line = line.rstrip()
		protA, protB = line.split("\t")
		toPredPPIs.append((protA, protB, "?"))

	toPred.calculateAllScores(thisScores, toPredPPIs)
	out = toPred.toTable(False)

	outFH = open(outF, "w")
	print >> outFH, out
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
