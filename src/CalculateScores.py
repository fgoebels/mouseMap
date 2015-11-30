#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionF, scores, outF) = sys.argv[1:]
	
	scores = scores.split("-")
	
	elutionData, toPred = calcS.loadEData(elutionF)
	thisScores = []
	twoDscores = toPred.twoDscores
	for score in scores:
		for twoDscore in twoDscores:
			if score == twoDscore.name: thisScores.append(twoDscore)
			if score == "MatrixNorms": thisScores.append(calcS.MatrixNorms())

	toPred.calculateAllPairs(thisScores)
	out = toPred.toTable(False)

	outFH = open(outF, "w")
	print >> outFH, out
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
