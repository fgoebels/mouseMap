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
	(complexF, eluteF, outDir) = sys.argv[1:]

	eluteData = {}
	eluteFH = open(eluteF)
	eluteFH.readline()
	head = "ID\tNAME\tGWEIGHT\tBlank"
	for line in eluteFH:
		line = line.rstrip()
		lineSplit= line.split("\t")
		if len(lineSplit[0].split("|")) != 3:
			print "No valid ID: %s" % lineSplit[0]
			continue
		protID = lineSplit[0].split("|")[1]	
		scores = np.array(lineSplit[1:241]).reshape(5,48)
		eluteData[protID] = matToLines(scores, "\t\t\t")
	eluteFH.close()

	clusters = {}
	complexFH = open(complexF)
	complexFH.readline()
	for line in complexFH:
		line = line.rstrip()
		(clusterID, protIDs) = line.split("\t")
		clusterID = int(clusterID)
		if clusterID not in clusters: clusters[clusterID] = set()
		for protID in protIDs.split(","):
			if protID in eluteData: clusters[clusterID].add(protID)
	complexFH.close()


	maxSize = 0
	for clust in clusters:
		if len(clusters[clust])>maxSize: maxSize = len(clusters[clust])

	tmp = "\tF"+"\tF".join(map(str, range(48))) + "\t\t\t\t"
	head +=  tmp*maxSize
	zeroMat = matToLines([["0.00"]*48]*5, "\t\t\t")
	NAMat = matToLines([["NA"]*48]*5, "\t\t\t")

	outFH = open(outDir + ".heatmap.txt", "w")
	print >> outFH, head
	blankLine = "Blank\tNA\t1\tNA" + "\tNA"*((48*maxSize)+maxSize)
	for cluster in sorted(clusters.keys()):
		matchedProts = []
		allprots = []
		for prot in sorted(list(clusters[cluster])):
			counts = eluteData[prot]
			totalCounts = getCounts(counts)
			if totalCounts < 20: continue
			matchedProts.append(prot)
			allprots.append(counts)
		matchedProtCounts = len(matchedProts)
		if matchedProtCounts < 3: continue	
		if matchedProtCounts < 6:
			outClustFH = open(outDir + ".%s.dat" % (cluster), "w")
			print >> outClustFH, printRmat(allprots, matchedProts)
			outClustFH.close()
			print "Plotting cluster:%s" % (cluster)
			wireframCMD = "src/3Dwireframe.R %s.%s.dat %s.%s.pdf" % (outDir, cluster, outDir, cluster)
			os.system(wireframCMD)
			os.system("rm %s.%s.dat" % (outDir, cluster))
		for i in range(maxSize-len(allprots)):
			allprots.append(NAMat)
		print >> outFH, matsToString(cluster, allprots)
		print >> outFH, blankLine
	outFH.close()

def printRmat(mats, names):
	out = "IEF\tIEX\tCounts\tProtName\n"
	for i in range(len(names)):
		mat = mats[i]
		counts = getCounts(mat)
		name = names[i]
		for IEFIndex in range(len(mat)):
			row = mat[IEFIndex].rstrip()	
			row = row.split("\t")
			for IEXIndex in range(len(row)):
				out += "%i\t%i\t%s\t%s\n" % (IEFIndex+1, IEXIndex+1, float(row[IEXIndex])/counts, names[i])
	return out

def getCounts(mat):
	counts = 0
	for IEFIndex in range(len(mat)):
		row = mat[IEFIndex].rstrip()
		counts += sum(map(int, row.split("\t")))
	return counts

def matsToString(cluster, mats):
	out = []
	for rowIndex in range(len(mats[0])):
		line = "Complex_%i\tPlaceholder_%i_%i\t1\tNA" % (cluster, cluster, rowIndex)
		for prot in mats:
			line += "\t" + prot[rowIndex]
		out.append(line)
	return "\n".join(out)

def matToLines(mat, suffix=""):
	out = []
	for row in mat:
		out.append("%s\t%s" % ("\t".join(row), suffix))
	return out
	
if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
