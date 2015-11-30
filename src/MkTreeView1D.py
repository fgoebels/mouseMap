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
	(complexF, eluteF, outF) = sys.argv[1:]

	eluteData = {}
	eluteFH = open(eluteF)
	head = "ID\tNAME\tGWEIGHT\tBlank\t" + "\t".join(map(lambda x: "Fraction "+str(x) , range(1,121)))
	for line in eluteFH:
		line = line.rstrip()
		lineSplit= line.split("\t")
		if len(lineSplit[0].split("|")) != 3:
			print "No valid ID: %s" % lineSplit[0]
			continue
		protID = lineSplit[0].split("|")[1]	
		scores = lineSplit[1:241]
		eluteData[protID] = scores
	eluteFH.close()

	clusters = {}
	complexFH = open(complexF)
	complexFH.readline()
	for line in complexFH:
		line = line.rstrip()
		(clusterID, protID) = line.split("\t")
		clusterID = int(clusterID)
		if clusterID not in clusters: clusters[clusterID] = []
		clusters[clusterID].append(protID)
	complexFH.close()

	outFH = open(outF, "w")
	print >> outFH, head
	blankLine = "Blank\tNA\t1\t" + "\t".join(["NA"]*120)
	for cluster in sorted(clusters.keys()):
		for prots in clusters[cluster]:
			matched = ""
			for prot in prots.split(","):
				if matched != "": continue
				if prot in eluteData: matched = prot
			counts = "\t".join(["0.00"]*120)
			if matched != "":
				counts = "\t".join(eluteData[matched])
			print >> outFH, "%i\t%s\t1\tNA\t%s" % (cluster, prot, counts)
		print >> outFH, blankLine
	outFH.close()
	
if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
