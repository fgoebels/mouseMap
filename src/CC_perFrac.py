#!/usr/bin/env python
from __future__ import division
import sys
import os
import re
import subprocess
import gzip
import math
import numpy as np
import utils 

def main():
	(protCCF, fracF, outDir) = sys.argv[1:]

	allCCs = set([])	
	prot2CC = {}
	protCCFH = open(protCCF)
	for line in protCCFH:
		line = line.rstrip()
		(prot, cc) = line.split("\t")
		if prot not in prot2CC: prot2CC[prot] = set([])
		prot2CC[prot].add(cc)
		allCCs.add(cc)
	protCCFH.close()

	counts = {}
	for i in range(1,9):
		counts[i] = {}
		for cc in allCCs:
			counts[i][cc] = 0

	outFH = open(outDir +".with_cc.txt", "w")
	fracFH = open(fracF)
	print >>outFH, fracFH.readline().rstrip()
	for line in fracFH:
		line = line.rstrip()
		linesplit = line.split("\t")
		protID = linesplit[0].split("|")[1]
		ccs = ""
		if protID in prot2CC:
			fracCounts = linesplit[3:]
			ccs = ",".join(prot2CC[protID])
			for i in range(len(fracCounts)):
				for cc in prot2CC[protID]:
	#				counts[i+1][cc] += int(fracCounts[i])
					if int(fracCounts[i]) > 2: counts[i+1][cc] += 1
		print >> outFH, line + "\t" + ccs

	fracFH.close()
	outFH.close()	

	allCCs = sorted(allCCs)

	out = "FracID\t" + "\t".join(allCCs) + "\n"
	for i in range(1,9):
		out += str(i)
		for cc in allCCs:
			out += "\t%s" % (str(counts[i][cc]))
		out += "\n"
	dataFH = open(outDir + ".dat", "w")
	print >> dataFH, out.rstrip()
	dataFH.close()

if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
