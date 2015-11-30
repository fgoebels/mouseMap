#!/usr/bin/env python
import numpy as np
import utils 
import sys

def main():
	(firstF, secondF, outF) = sys.argv[1:]

	out = {}
	firstFH = open(firstF)
	header = firstFH.readline()
	out = addData(firstFH, out)
	firstFH.close()
	secondFH = open(secondF)
	secondFH.readline()
	out = addData(secondFH, out)
	secondFH.close()
	outFH = open(outF, "w")
	print >> outFH, header
	for (idA, idB, label) in out:
		(scoreA, scoreB) = out[(idA, idB, label)]
		print >> outFH, "%s\t%s\t%f\t%f\t%s" % (idA, idB, scoreA, scoreB, label)
	outFH.close()

def addData(dataFH, out):
	for line in dataFH:
		line = line.rstrip()
		(idA, idB, scoreA, scoreB, label) = line.split("\t")
		edge = (idA, idB, label)
		if edge not in out:
			out[edge]= [float(scoreA), float(scoreB)]
		else:
			tmpscoreA = (float(scoreA) + out[edge][0])/2
			tmpscoreB = (float(scoreB) + out[edge][1])/2
			out[edge] = [tmpscoreA, tmpscoreB]
	return out	

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
