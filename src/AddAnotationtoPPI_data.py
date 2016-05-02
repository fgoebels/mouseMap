#!/usr/bin/env python
from __future__ import division
import numpy as np
import re
import utils 
import sys

def main():
#	src/AddAnotationtoPPI_data.py data/WORM/final_PPI_ratio5.ids.txt data/WORM/worm_gene_uniprot.map.txt data/WORM/final_PPI_ratio5.gprofiler_out.map.txt data/WORM/Human_elegans_uniprot_names.txt data/WORM/uniprotGeneNameMap.txt data/WORM/human.var.txt data/WORM/gene.omim.combined data/WORM/allComplexes.tab test/test.out
	(protIDsF, worm_gene2UniprotF, goMapF, othologMapF, human_gene2UniprotF, snpF, disF, complexF, outF) = sys.argv[1:]

	protIDs = utils.readData(protIDsF, np.array([0]), np.array([0]))
	worm_gene2uniprot = utils.readData(worm_gene2UniprotF, np.array([0]), np.array([1]))
	goData = utils.readData(goMapF, np.array([0]), np.array([1,2]))

	orthmap = utils.readData(othologMapF, np.array([1]), np.array([0]))
	mapData = utils.readData(human_gene2UniprotF, np.array([1]), np.array([0]))
	snpData = utils.readData(snpF, np.array([0]), np.array([1,2,3]), mapData)
	disData = utils.readData(disF, np.array([0]), np.array([1]), primKeyMap= mapData)
	complexData = utils.readData(complexF, np.array([0]), np.array([2]))

	header =  "GeneName\tUniprotIDs\tHumanOrthologIDs\tEnriched_GO_terms\tSNP\tDisease\tComplexe"
	cats = header.split("\t")
	counts = {"wAnno" : set([])}
	for cat in cats:
		counts[cat] = set([])

	outFH = open(outF, "w")
	print >> outFH, header
	for protID in protIDs:
		line = protID[0]
		line += "\t" + annoToString(getAnnotation(protID, worm_gene2uniprot))
		line += "\t" + annoToString(getAnnotation(protID, orthmap, [worm_gene2uniprot]))
		line += "\t" + annoToString(getAnnotation(protID, goData))
		line += "\t" + annoToString(getAnnotation(protID, snpData, [worm_gene2uniprot, orthmap]))
		line += "\t" + annoToString(getAnnotation(protID, disData, [worm_gene2uniprot, orthmap]))
		line += "\t" + annoToString(getAnnotation(protID, complexData, [worm_gene2uniprot, orthmap]))
		lineSplit = line.split("\t")
		for i in range(len(lineSplit)):
			col = lineSplit[i]
			if col != "-":
				counts[cats[i]].add(protID[0])
		
		for i in range(3,len(lineSplit)):
			col = lineSplit[i]
			if col != "-":
				counts["wAnno"].add(protID[0])
		
		print >> outFH, line

	for cat in counts:
		print "%s\t%i" % (cat, len(counts[cat]))

def annoToString(annotation):
	tmp = []
	for anno in annotation:
		tmp.append(",".join(anno))
	return ";".join(tmp) 

def getAnnotation(prot, anno, mappings = ""):
	out = set([])
	mappedIDs = set([prot])
	if mappings != "":
		for maping in mappings:
			mappedIDs = mapIDs(mappedIDs, maping)
	for protID in mappedIDs:
		if protID in anno:
			out.update(anno[protID])
	if len(out) == 0: out = "-"
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
