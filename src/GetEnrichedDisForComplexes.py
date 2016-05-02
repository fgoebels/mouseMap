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

# src/GetEnrichedDisForComplexes.py test/Final_worm_complexes.go.txt data/WORM/worm_gene_uniprot.map.txt data/WORM/Human_elegans_uniprot_names.txt data/WORM/gene.omim.combined data/WORM/uniprotGeneNameMap.txt data/doid.obo.tab data/doid.omim.map data/Omim_id_name.txt data/gene2DOID.txt test/Final_worm_complexes.go-dis.txt	

def main():
	(clusterF, worm_gene2UniprotF, othologMapF, disF, human_gene2UniprotF, oboF, omim_doid_F, omim2NameF, gene2DOIDF, outF) = sys.argv[1:]
	
	
	obo = readOntology(oboF, False)

	omim_doid = utils.readData(omim_doid_F, np.array([2]), np.array([0]))
	omim_name = utils.readData(omim2NameF, np.array([0]), np.array([1]))
	orthmap = utils.readData(othologMapF, np.array([1]), np.array([0]))
	worm_gene2uniprot = utils.readData(worm_gene2UniprotF, np.array([0]), np.array([1]))
	doid2Name = utils.readData(oboF, np.array([0]), np.array([1]))	

	mapData = utils.readData(human_gene2UniprotF, np.array([1]), np.array([0]))
	disData = utils.readData(disF, np.array([0]), np.array([1]), primKeyMap= mapData)
	stats = importr('stats')

	doid2Gene={}
        gene2DOIDFH = open(gene2DOIDF)
        for line in gene2DOIDFH:
                line = line.rstrip()
                (gene, thisDOID) = line.split("\t")
                if not thisDOID.startswith("DOID"): continue
                allParents = set([])
#                allParents = getAll(thisDOID, obo)
                allParents.add(thisDOID)
                for doid in allParents:
                        if doid not in doid2Gene: doid2Gene[doid] = set([])
                        doid2Gene[doid].add(gene)
        gene2DOIDFH.close()


        allGenes = 20313 
	allProts= 1713 # TOUPDATE when new complexes are used

	outFH = open(outF, "w")
	clusterFH = open(clusterF)
	for line in clusterFH:
		line = line.rstrip("\n")
		cluster_genes = line.split("\t")[0].split(",")
		cluster_prots = set([])
		cluster_omim = set([])
		doid_cluster_counts = {}
		for gene in cluster_genes:
			cluster_prots.update(getAnnotation((gene,), orthmap, [worm_gene2uniprot]))
			cluster_omim.update(getAnnotation((gene,), disData, [worm_gene2uniprot, orthmap]))
			doidIDs = mapIDs(cluster_omim, omim_doid)
			if len(doidIDs)>0:
				for doidID in doidIDs:
					parentDOID = set([])
#					parentDOID = getAll(doidID[0], obo)
					parentDOID.add(doidID[0])
					for allIDs in parentDOID:
						if allIDs not in doid_cluster_counts: doid_cluster_counts[allIDs] = set([])
						doid_cluster_counts[allIDs].add(gene)
		pvals = []
		ids = []
		for doidID in doid_cluster_counts:
			if doidID=="DOID:0000004" or doidID=="---" :continue
			if len(doid_cluster_counts[doidID]) == 1: continue
			mat = [[allProts, allGenes],[len(doid_cluster_counts[doidID]), len(doid2Gene[doidID])]]
			pval = fisher_exact(mat)[1]
			pvals.append(pval)
			ids.append(doidID)
		pvals = stats.p_adjust(FloatVector(pvals), method = 'BH')
		enrichedDOIDs = set([])
			
		for i  in range(len(pvals)):
			if pvals[i] <= 0.05:
				enrichedDOIDs.add("%s,%s,%i,%.4f" % (ids[i], doid2Name[(ids[i],)][0], len(doid_cluster_counts[ids[i]]),pvals[i]))
		
		tmp = set()
		for mim in cluster_omim:
			if mim in omim_name:
				tmp.add((mim[0], omim_name[mim][0][0]))
		doidIDs = mapIDs(cluster_omim, omim_doid)
		cluster_omim = tmp
		print >> outFH, "%s\t%s\t%s\t%s" % (line, annoToString(cluster_prots, sepB=","), annoToString(cluster_omim), ";".join(enrichedDOIDs))

	clusterFH.close()
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
		isA = isA.split(" ")
                if descending:
			for doidID in isA:
	                        out[doidID].add(ontoID)
                else:
			for doidID in isA:
                        	out[ontoID].add(doidID)
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

def annoToString(annotation, sepA = ",", sepB=";"):
        tmp = []
        for anno in annotation:
                tmp.append(sepA.join(anno))
        return sepB.join(tmp)

def getAnnotation(prot, anno, mappings = ""):
        out = set([])
        mappedIDs = set([prot])
        if mappings != "":
                for maping in mappings:
                        mappedIDs = mapIDs(mappedIDs, maping)
        for protID in mappedIDs:
                if protID in anno:
                        out.update(anno[protID])
        return out

def mapIDs(ids, mapping):
        out = set([])
        for currentID in ids:
                if currentID in mapping:
                        out.update(mapping[currentID])
        return out


if __name__ == "__main__":
	try:
	        main()
	except KeyboardInterrupt:
	        pass
