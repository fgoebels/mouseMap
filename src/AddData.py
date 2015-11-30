#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils 
import sys

def main():
	(netF, toAddF, loc, force, outF) = sys.argv[1:]


	
	toAdd = {}
	toAddFH = open(toAddF)
	featureNames = toAddFH.readline()
	featureNames = featureNames.rstrip()
	featureNames = featureNames.split("\t")[2:]
	numFeatures = len(featureNames)
	for line in toAddFH:
		line = line.rstrip()
		lineSplit = line.split("\t")
		edge = tuple(sorted(lineSplit[:2]))
		features = lineSplit[2:]
		toAdd[edge] = features
	toAddFH.close()
	

	netFH = open(netF)
	outFH = open(outF, "w")
	curHead = netFH.readline()
	curHead = curHead.rstrip()
	curHead = curHead.split("\t")
	if loc == "b":
		print >> outFH, "%s\t%s\t%s" % ("\t".join(curHead[:2]), "\t".join(featureNames), "\t".join(curHead[2:]))
	else:
		print >> outFH, "%s\t%s" % ("\t".join(curHead), "\t".join(featureNames))
	for line in netFH:
		line = line.rstrip()
		lineSplit = line.split("\t")
		edge = tuple(sorted(lineSplit[:2]))
		features = "\t".join(lineSplit[2:])
		addedFeatures = "\t".join(["0"]*numFeatures)
		if edge not in toAdd and force == "F": continue
#			if features == "positive":
#				print "no match for %s" % (",".join(edge))
		if edge in toAdd: addedFeatures = "\t".join(toAdd[edge])
		edge = "\t".join(edge)
		if loc == "b":
			print >> outFH, "%s\t%s\t%s" % (edge, addedFeatures, features)
		else:
			print >> outFH, "%s\t%s\t%s" % (edge, features, addedFeatures)
	outFH.close()
	netFH.close()	

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
