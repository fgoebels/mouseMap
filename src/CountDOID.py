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
from scipy.stats import fisher_exact
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def main():
	(protAnnoF, oboQueryF, oboF, omim_doid_F, gene2DOIDF, outF) = sys.argv[1:]
	
	obo = readOntology(oboF, False)

	omim_doid = utils.readData(omim_doid_F, np.array([2]), np.array([0]))

	qnames = utils.readData(oboQueryF, np.array([0,1]), np.array([0]))	

	protAnnoFH = open(protAnnoF)
	protAnnoFH.readline()
	prot2DOID = {}

	omimperOBO = {}
	COUNTperQTerm = {}
	for line in protAnnoFH:
		line = line.rstrip()
		prot = line.split("\t")[0]
		for m in re.finditer('OMIM:\d{6}', line):
			omimID = m.group(0)
			if (omimID, ) not in omim_doid: continue
			for doidID in omim_doid[(omimID, )]:
				parentDOID = getAll(doidID[0], obo)
				parentDOID.add(doidID[0])
				for (queryDOID, name) in qnames:
					if queryDOID in parentDOID:
						if queryDOID not in COUNTperQTerm: COUNTperQTerm[queryDOID] = set([])
						if queryDOID not in omimperOBO: omimperOBO[queryDOID] = set([])
						COUNTperQTerm[queryDOID].add(prot)
						omimperOBO[queryDOID].add((prot, omimID))

	protAnnoFH.close()


	allProts = len(utils.readData(protAnnoF, np.array([0]), np.array([0])))
        doid2Gene= {}
        gene2DOIDFH = open(gene2DOIDF)
        for line in gene2DOIDFH:
                line = line.rstrip()
                (gene, thisDOID) = line.split("\t")
                if not thisDOID.startswith("DOID"): continue
                allParents = getAll(thisDOID, obo)
                allParents.add(thisDOID)
                for doid in allParents:
                        if doid not in doid2Gene: doid2Gene[doid] = set([])
                        doid2Gene[doid].add(gene)
        gene2DOIDFH.close()

        allGenes = 20313 #len(allGenes)

	pvals = []
	names = []
	for (doid, name) in qnames:
		if doid not in doid2Gene: continue
		if doid not in COUNTperQTerm: continue
		mat = [[allProts, allGenes],[len(COUNTperQTerm[doid]), len(doid2Gene[doid])]]
		pval = fisher_exact(mat)[1]
		pvals.append(pval)
		names.append((doid, name))
	
	stats = importr('stats')
	selNames = {}
	pvals = stats.p_adjust(FloatVector(pvals), method = 'BH')
	for i in range(len(pvals)):
#		if pvals[i]>0.05: continue
		selNames[names[i]] = pvals[i]
	
	catCounts = {}
	for selname in selNames:
		pval = selNames[selname]
		doid, name = selname
		if len(COUNTperQTerm[doid]) == 0: continue
		counts = len(COUNTperQTerm[doid])
		omimCounts = len(omimperOBO[doid])
		if counts not in catCounts: catCounts[counts] = set([])
		catCounts[counts].add("%s\t%s\t%i\t%i\t%f" % (name, doid, counts, omimCounts, pval))

	outFH = open(outF + ".dat", "w")
	print >> outFH, "Name\tDOID\tCounts\nAll sites\tNA\t%i" % (allProts)
	for counts in sorted(catCounts.keys(), reverse=True):
		print >> outFH, "\n".join(catCounts[counts])
		print "\n".join(catCounts[counts])
	outFH.close()


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
