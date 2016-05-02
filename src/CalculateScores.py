#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy, inspect

def main():
	(elutionF, outF) = sys.argv[1:]
	
	elutionData, scoreCalc = calcS.loadEData(elutionF)

	scoreCalc.calculateAllPairs([calcS.MutualInformation(2)])

	outFH = open(outF, "w")
	outFH.write(scoreCalc.toTable(False))
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
