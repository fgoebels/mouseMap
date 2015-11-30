#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils 
import sys

def main():
	(scoresF, rownamesF, orthmapF, outF) = sys.argv[1:]


	orthmapFH = open(orthmapF)
	orthmapFH.readline()
	orthmapping = {}
	for line in orthmapFH:
		line = line.rstrip()
		(orthID, sgdID) = line.split("\t")
		if orthID not in orthmapping: orthmapping[orthID] = set([])
		orthmapping[orthID].add(sgdID)
	orthmapFH.close()

	names = {}
	i = 0
	rownamesFH = open(rownamesF)
	rownamesFH.readline()
	for line in rownamesFH:
		line = line.rstrip()
		geneName = line
		if geneName in orthmapping and len(orthmapping[geneName])==1:
			geneName = list(orthmapping[geneName])[0]
		names[i] = geneName
		i += 1
	rownamesFH.close()
	

	coelutionData = []
	scoresFH = open(scoresF)
	for line in scoresFH:
		line = line.rstrip()
		coelutionData.append(line.split("\t"))
	scoresFH.close()

	
	outFH = open(outF, "w")
	print >> outFH, "IDA\tIDB\t%s" % (scoresF.split(".")[-1])
	out = {}	
	for i in range(len(names)):
		protA = names[i]
#		if (protA, ) not in prots: continue
		for j in range(i+1, len(names)):
			protB = names[j]
#			if (protB, ) not in prots: continue
			score = coelutionData[i][j]
			if score != coelutionData[j][i]:
				print "matrix not symetical for %i and %i" % (i,j)
				sys.exit()
			edge = tuple(sorted([protA, protB]))
			score = float(score)
			print >> outFH, "%s\t%s" % ("\t".join(edge), score)
#			if edge not in out:
#				out[edge] = score
#			elif out[edge]< score:
#				out[edge] = score
#	
#	outFH = open(outF, "w")
#	print >> outFH, "IDA\tIDB\tscore"
#	for edge in out:
#		score = out[edge]
#		print >> outFH, "%s\t%f" % ("\t".join(edge), score)
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
