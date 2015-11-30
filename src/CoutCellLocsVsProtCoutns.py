#!/usr/bin/env python
from __future__ import division
import sys
import os
import numpy as np

def main():
	(protCountsF, tissueCountsF, outF) = sys.argv[1:]

	protCounts = readCounts(protCountsF)

	tissueCounts = readCounts(tissueCountsF, False)
	
	outFH = open(outF, "w")
	print >> outFH, "\tNumberofTissue\tNumberofFractions\tTotalProtCounts"

	
	for geneName in tissueCounts:
		if geneName not in protCounts: continue
#		print >> outFH, "%s\t%i\t%i\t%i" % (geneName, np.count_nonzero(tissueCounts[geneName]), np.count_nonzero(protCounts[geneName]), sum(protCounts[geneName]))
		print >> outFH, "%s\t%s\t%i\t%i" % (geneName, "\t".join(tissueCounts[geneName]), np.count_nonzero(protCounts[geneName]), sum(protCounts[geneName]))

	outFH.close()

def readCounts(countF, mapToFloat = True):
	out = {}
	countFH = open(countF)
	countFH.readline()
	for line in countFH:
		line = line.rstrip()
		if line.startswith("#"): continue
		lineSplit = line.split("\t")
		geneName = lineSplit[0]
		if geneName not in out: out[geneName] = []
		for count in lineSplit[1:]:
			if mapToFloat: count = float(count)
			out[geneName].append(count)	
	countFH.close()
	return out

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
