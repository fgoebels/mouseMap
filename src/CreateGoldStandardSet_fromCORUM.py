#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils 
import sys

def main():
	(proteinComplexes, goAnnotation, outF) = sys.argv[1:]

	complexes2prot = utils.readData(proteinComplexes, np.array([1]), np.array([0]))
	prot2GO = utils.readData(goAnnotation, np.array([1]), np.array([0]))

	complexes2prot = filterComplexes(complexes2prot, 2,40)
	print len(complexes2prot)
	complexes2prot = mergeComplexes(complexes2prot)
	print len(complexes2prot)


	prot2complexes = getProtsToComplex(complexes2prot)
	out = getGoldStandard(prot2complexes, prot2GO)
	outFH = open(outF, "w")
	outFH.write( "IDA\tIDB\tLabel\n" + "\n".join(out))
	outFH.close()
	

def filterComplexes(complexList, sizeL, sizeU):
        todel = set()
        for comp in complexList:
                if len(complexList[comp]) < sizeL or len(complexList[comp]) > sizeU:
			todel.add(comp)
        for delComp in todel:
                del complexList[delComp]
        return complexList

def getGoldStandard(prot2complexes, goAnnotation, factor=10):
	pos = set([])
	neg = set([])
	prots = prot2complexes.keys()
	for i in range(len(prots)):
		protA = prots[i]
		for j in range(i+1, len(prots)):
			protB = prots[j]
			if len(prot2complexes[protA] & prot2complexes[protB]) > 0:
				pos.add("%s\t%s\tpositive" % (protA[0], protB[0]))
				continue
			elif protB not in goAnnotation or protA not in goAnnotation:
				continue
			elif len(set(goAnnotation[protA]) & set(goAnnotation[protB])) == 0:
				neg.add("%s\t%s\tnegative" % (protA[0], protB[0]))
				continue
	neg = np.random.choice(list(neg), len(pos)*factor)
	pos.update(neg)
	return pos
		
def mergeComplexes(complexList):
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
                        if overlap(complexList[compI],complexList[compJ]) >= 0.5:
                                toMerge.add(compJ)
                                merged.add(compJ)
                if len(toMerge)>0:
                        toMerge.add(compI)
                        (newName,  newProts) = set(), set()
                        for name in toMerge:
                                newName.add(name[0])
                                newProts.update(complexList[name])
                        out[(",".join(newName),)] = newProts
                        merged.add(compI)
                else:
                        out[compI] = complexList[compI]
                        merged.add(compI)

        return out


def overlap(a,b):
        tmpa = set(a)
        tmpb = set(b)
        overlap = len(tmpa & tmpb)/len(tmpa | tmpb)
        return overlap


def getProtsToComplex(complexes2prot):
	out = {}
	for cluster in complexes2prot:
		for prot in complexes2prot[cluster]:
			if prot not in out: out[prot] = set([])
			out[prot].add(cluster)
	return out

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
