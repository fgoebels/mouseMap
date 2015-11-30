#!/usr/bin/env python
import sys

def main():
	(ppiF, annoF, outF) = sys.argv[1:]
	
	anno = {}
	annoFH = open(annoF)
	for line in annoFH:
		line = line.rstrip()
		(protID, geneName, species) = line.split("\t")
		anno[protID] = (geneName, species)
	
	outFH = open(outF, "w")	
	ppiFH = open(ppiF)
	for line in ppiFH:
		line = line.rstrip()
		(idA, idB, score) = line.split("\t")
		if idA not in anno or idB not in anno:
			print "No infor for %s %s" % (idA, idB)
			continue
		nameA, specA = anno[idA]
		nameB, specB = anno[idB]
		if specA != specB:
			print "Inter species"
			continue
		print >> outFH, "%s\t%s\t%s\t%s\t%s\t%s" % (idA, idB, nameA, nameB, score, specA)

	outFH.close()
	ppiFH.close()


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pas
