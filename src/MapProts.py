#!/usr/bin/env python
import numpy as np
import utils 
import sys

def main():
	(rawDataF, mappingF, outF) = sys.argv[1:]

	mapData = utils.readData(mappingF, np.array([0]), np.array([1]))	

	mapping = {}

	mappingFH = open(mappingF)
	mappingFH.readline()
	for line in mappingFH:
		line = line.rstrip()
		if len(line.split("\t")) != 2: continue
		(wormID, geneName) = line.split("\t")
		(geneName, wormID) = line.split("\t")
		
		if geneName not in mapping: mapping[geneName] = set([])
		if wormID not in mapping: mapping[wormID] = set([])
		mapping[geneName].add(wormID)
		mapping[wormID].add(geneName)

	outFH = open(outF, "w")
	dataFH = open(rawDataF)
	outFH.write(dataFH.readline())
	for line in dataFH:
		line = line.rstrip()
		lineSplit = line.split("\t")
		idA = lineSplit[0]
		mapA = mapID(idA, mapping)
#		mapB = mapID(idB, mapping)
		if mapA == "" : continue
		print >> outFH, "%s\t%s" % (mapA, "\t".join(lineSplit[1:]))
	dataFH.close()
	outFH.close()

def mapID(toMap, mapping):
		mappedID = ""
		if toMap not in mapping:
			print "No map for %s" % (toMap)
			return mappedID
		if len(mapping[toMap]) > 1 :
			print "No uniq map for %s" % (toMap)
			return mappedID
		tmpID = "".join(mapping[toMap])
#		if len(mapping[tmpID]) > 1 :
#			print "No uniq map for %s" % (tmpID)
#			return mappedID
		mappedID = tmpID
		return mappedID


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
