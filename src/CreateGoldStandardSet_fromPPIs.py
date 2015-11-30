#!/usr/bin/env python
from __future__ import division
import numpy as np
import utils 
import sys

def main():
	(ppiF, goAnnotation, outF) = sys.argv[1:]

	ppis = utils.readData(ppiF, np.array([0,1]), np.array([0]))
	prot2GO = utils.readData(goAnnotation, np.array([1]), np.array([0]))

	out = getGoldStandard(ppis, prot2GO)
	outFH = open(outF, "w")
	outFH.write( "IDA\tIDB\tLabel\n" + "\n".join(out))
	outFH.close()

def getGoldStandard(positivePPIs, goAnnotation):
	prots = getAllProts(positivePPIs)
	pos = set([])
	neg = set([])
	for i in range(len(prots)):
		protA = prots[i]
		for j in range(i+1, len(prots)):
			protB = prots[j]
			if (protA, protB) in positivePPIs or (protB, protA) in positivePPIs: 
				pos.add("%s\t%s\tpositive" % (protA, protB))
				continue
			elif (protB,) not in goAnnotation or (protA,) not in goAnnotation:
				continue
			elif len(set(goAnnotation[(protA,)]) & set(goAnnotation[(protB,)])) == 0:
				neg.add("%s\t%s\tnegative" % (protA, protB))
				continue
	neg = np.random.choice(list(neg), len(pos)*6)
	pos.update(neg)
	return pos

def getAllProts(protList):
	out = set()
	for protA, protB in protList:
		out.add(protA)
		out.add(protB)
	return list(out)
		
if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
