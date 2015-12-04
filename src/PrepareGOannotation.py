#!/usr/bin/env python
import numpy as np
import utils
import CalculateCoElutionScores as calcS
import sys
import os, math, copy

def main():
	(elutionF, geneAnnoF, outF) = sys.argv[1:]
	
	prots = set([])
	elutionFH = open(elutionF)
	elutionFH.readline()
	for line in elutionFH:
		line = line.rstrip()
		prots.add(line.split("\t")[0])
	elutionFH.close()

	tofilter = ["GO:0003674", "GO:0005575", "GO:0008150"]
	out = set([])
	geneAnnoFH = open(geneAnnoF)
	for line in geneAnnoFH:
		line = line.rstrip()
		if line != "" and not line.startswith('!'):
	                line = line.split('\t')
	                term = line[4].strip()
			if term in tofilter: continue
	                gene = line[1].strip()
	                code = line[6].strip()
			syns  = line[10].strip().split("|")
			syns.append(gene)
			for name in syns:
				if name in prots:
					out.add("%s\t%s\t%s" % (name, term, code))
	geneAnnoFH.close()

	
	outFH = open(outF, "w")
	print >> outFH, "\n".join(out)
	outFH.close()
	
if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
