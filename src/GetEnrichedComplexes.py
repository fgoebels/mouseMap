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
from scipy.stats import hypergeom
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def main():
	(wormF, complexF, outF) = sys.argv[1:]


#	complexes = utils.readData(complexF , np.array([2,1]), np.array([0]))
	complexes = readComplexes(complexF)

	wormFH = open(wormF)
	coFracComp = {}
	compID = 0
	data = {}
	for line in wormFH:
		mappedComplexes = set([])
		line = line.rstrip("\n")
		complexProts = line.split("\t")[0].split(",")
		coFracComp[str(compID)] = complexProts
		compID +=1
		data[str(compID)] = line
	wormFH.close()

	print "Raw corum: %i" % (len(complexes))
	print "Raw cofrac: %i" % (len(coFracComp))
	print len(getOverlapp(coFracComp, complexes))
	complexes = filterComplexes(complexes, 4, 20)
	coFracComp = filterComplexes(coFracComp, 4, 20)
	print "filtered corum: %i" % (len(complexes))
	print "Filtered cofrac: %i" % (len(coFracComp))
	print len(getOverlapp(coFracComp, complexes))
	complexes = mergeComplexes(complexes)
	coFracComp = mergeComplexes(coFracComp)
	print "Merged corum: %i" % (len(complexes))
	print "Merged cofrac: %i" % (len(coFracComp))
	matchedComp = getOverlapp(coFracComp, complexes)
	print len(matchedComp)
	outFH = open(outF, "w")
	for comp in data:
		matched = ""
		if comp in matchedComp: matched = ";".join(matchedComp[comp])
		print >> outFH, "%s\t%s" % (data[comp], matched)
	outFH.close()
	
def getOverlapp(compA, compB):
	mappedComplexes = {}
	for compa in compA:
		for compb in compB:
			if overlap(compA[compa], compB[compb])>=0.5:
				if compa not in mappedComplexes: mappedComplexes[compa] = set()
				mappedComplexes[compa].add(compb)
	return mappedComplexes
	

def con(protList, sep = ","):
	tmp = set()
	for prot in protList:
		tmp.add(prot[0])
	return sep.join(tmp)

def getAllCompProts(complexList):
	out = set()
	for comp in complexList:
		out.update(complexList[comp])
	return out

def filterComplexes(complexList, sizeL, sizeU):
	todel = set()
	for comp in complexList:
		if len(complexList[comp]) < sizeL or len(complexList[comp]) > sizeU: todel.add(comp)
	for delComp in todel:
		del complexList[delComp]
	return complexList

def readComplexes(compF):
	out = {}
	compFH = open(compF)
	i = 0
	for line in compFH:
		line = line.rstrip()
		linesplit = line.split("\t")
		out["CORUM:%s" % linesplit[0]] = linesplit[1:]
		i += 1
	compFH.close()
	return out

def mergeComplexes(complexList, cutoff= 0.8):
	merged = set()
	allComplexes = complexList.keys()
	out = {}
	for i in range(len(allComplexes)):
		toMerge = set()
		compI = allComplexes[i]
		if compI in merged: continue
		for j in range(i+1, len(allComplexes)):
			compJ = allComplexes[j]
			if compJ in merged: continue
			if overlap(complexList[compI],complexList[compJ]) >= cutoff:
				toMerge.add(compJ)
				merged.add(compJ)
		if len(toMerge)>0:
			toMerge.add(compI)
			(newName, newProts) = set(), set()
			for name in toMerge:
				newName.add(name)
				newProts.update(complexList[name])
			out[";".join(map(str, newName))] = newProts
			merged.add(compI)
		else:
			out[compI] = complexList[compI]
			merged.add(compI)

	return out


def overlap(a,b):
	tmpa = set(a)
	tmpb = set(b)
	#overlap = len(tmpa & tmpb)/len(tmpa | tmpb)
	overlap = len(tmpa & tmpb)/min(len(tmpa), len(tmpb))
	return overlap
	
if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
