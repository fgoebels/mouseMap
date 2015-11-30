#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils 
import sys

def main():
	(dataF, outF) = sys.argv[1:]
	
	dataFH = open(dataF)
	head = dataFH.readline()
	head = head.rstrip()
	head = head.split("\t")[2:]
	outFH = open(outF, "w")
	print >> outFH, "@RELATION COFrac"
	for colName in head[:-1]:
		print >> outFH, "@Attribute %s NUMERIC" % (colName)
	print >> outFH, "@ATTRIBUTE class {positive, negative}\n@DATA"
	for line in dataFH:
		line = line.rstrip()
		lineSplit = line.split("\t")[2:]
		outline = ",".join(lineSplit)
		print >> outFH, outline
	dataFH.close()
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
