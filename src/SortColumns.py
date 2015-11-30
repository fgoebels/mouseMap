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
	(inF,  outF) = sys.argv[1:]
	
	inFH = open(inF)
	header = inFH.readline()
	header = header.split("\t")
	portIDname = header[0]
	colNames = header[1:]
	i = 0
	sortedColnames = {}
	for colName in colNames:
		newIndex = re.search("(\d+).mzXML_dta", colName).group(1)
		sortedColnames[i] = int(newIndex)-1
		i += 1
	outFH = open(outF, "w")
	
	print >> outFH, "%s\t%s" % ( portIDname, "\t".join(sortArray(colNames, sortedColnames)))
	for line in inFH:
		line = line.rstrip()
		linesplit = line.split("\t")
		name = linesplit[0]
		vals = linesplit[1:]
		vals = "\t".join(sortArray(vals, sortedColnames))
		print >> outFH, "%s\t%s" % (name, vals)

def sortArray(array, mapping):
	out = [0]*len(array)
	for i in range(len(array)):
		out[mapping[i]] = array[i]
	return out

if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
