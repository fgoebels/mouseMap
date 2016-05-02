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
	(oboMapF, oboQueryF, oboF, outF) = sys.argv[1:]
	
	obo = readOBO(oboF)

	oboMapFH = open(oboMapF)
	oboMapFH.readline()
	obomapping = {}
	for line in oboMapFH:
		line = line.rstrip()
#		(doid, gene) = line.split("\t")
		(gene, doid) = line.split("\t")
		if gene not in obomapping: obomapping[gene] = set([])
		obomapping[gene].add(doid)
	
	oboQueryFH = open(oboQueryF)
	outFH = open(outF, "w")
	for line in oboQueryFH:
		line = line.rstrip()
		(queryID, name) = line.split("\t")
		counts = 0
		childrens = getAll(queryID, obo)
		for site in obomapping:
			if len(obomapping[site] &  childrens)>0:
				print >> outFH, "%s\t%s" % (site, name)

def getAll(doidID, onto):
	out = set([doidID])
	tosearch = [doidID]
	while(len(tosearch)>0):
		curID = tosearch.pop(0)
		for child in onto[curID]:
			if child in out: continue
			out.add(child)
			tosearch.append(child)
	return out

def readOBO(ontoF, descending=True):
	out = {}
	ontoFH = open(ontoF)
	ontoFH.readline()
	goID = ""
	for line in ontoFH:
		line = line.rstrip()
		if line.startswith("id: "):
			goID = line[4:]
			continue
		if line == "[TERM]":
			goID = ""
		if line.startswith("is_a: "):
			isA = line[6:].split(" ! ")[0]
			if isA not in out: out[isA] = set()
			if goID not in out: out[goID] = set()
			if descending:
				out[isA].add(goID)
			else:
				out[goID].add(isA)
	ontoFH.close()
	return out

def readOntology(ontoF, descending=True):
	out = {}
	ontoFH = open(ontoF)
	ontoFH.readline()
	for line in ontoFH:
		line = line.rstrip()
		(ontoID, name, termType, isObsolete, isA, omimID) = line.split("\t")
		if isA not in out: out[isA] = set()
		if ontoID not in out: out[ontoID] = set()
		if descending:
			out[isA].add(ontoID)
		else:
			out[ontoID].add(isA)
	return out

if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
