#!/usr/bin/env python
import numpy as np
import utils 
import sys
import xmlrpclib
import random

def main():
	(GOfilesF, outF) = sys.argv[1:]

	server = xmlrpclib.ServerProxy('http://funsimmat.bioinf.mpi-inf.mpg.de/xmlrpc.php')
	
	allIDs = set([])
	GOfilesFH = open(GOfilesF)
	for line in GOfilesFH:
	        line = line.rstrip()
	        allIDs.add(line)
	GOfilesFH.close()

	allIDs = list(allIDs)
	
	allPairs = {}
	for i in range(len(allIDs)):
	        for j in range(i+1, len(allIDs)):
			pair = tuple(sorted([allIDs[i], allIDs[j]]))
        	        allPairs[pair] = "NA"



	toq = set(allPairs.keys())
	print len(toq)
	outFH = open(outF, "w")
	i = 0
	for idA, idB in toq:
		result = server.Semantic.getSemSims(",".join([idA, idB]))
		if len(result) >=4:
			print >> outFH, "\t".join(result[3])
		i += 1
		if i % 100 == 0:
			print len(toq)-i
	outFH.close()

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
