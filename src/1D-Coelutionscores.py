#!/usr/bin/env python
from __future__ import division
import numpy as np
import scipy.stats
import utils 
import sys
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from scipy.spatial import distance
import itertools
import operator
import copy
from collections import defaultdict
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support
from sklearn.cross_validation import train_test_split, cross_val_predict
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import inspect
import os

class ElutionData():
	def __init__(self, elutionProfileF):
		self.elutionMat, self.prot2Index  = self.loadElutionData(elutionProfileF)
		self.normedElutionMat = normalize_fracs(self.elutionMat)
		self.elutionMat = np.array(self.elutionMat)

	def loadElutionData(self, elutionProfileF):
		elutionProfileFH = open(elutionProfileF)
		elutionProfileFH.readline()
		elutionProfile = {}
		i = 0
		elutionMat = []
		prot2Index = {}
		for line in elutionProfileFH:
			line = line.rstrip()
			line = line.split("\t")
			protID = line[0]
			counts = np.array(map(float, line[1:]))
			elutionMat.append(counts)
			prot2Index[protID] = i
			i += 1
		elutionProfileFH.close()
		return elutionMat, prot2Index

	def normalizeCoEulitionMat(self):
		self.elutionMat = normalize_fracs(self.elutionMat)

	def getElution(self, prot, normed=False):
		if prot not in self.prot2Index:
			return float('NaN')
		tmpout = ""
		if normed: 
			tmpout = self.normedElutionMat[self.prot2Index[prot]]
		else:
			tmpout = self.elutionMat[self.prot2Index[prot]]
		return tmpout
	

	def hasProt(self, prot):
		return prot in self.prot2Index

	def printMat(self, mat):
		out = "ProtiID\tFraction_" + "\tFracion_".join(map(str, range(1,240)))
		for prot in self.prot2Index:
			index = self.prot2Index[prot]
			out += "\n%s\t%s" %  (prot, "\t".join(map(str,self.normedElutionMat[index])))
		return out

def arr_norm(arr, axis=0):
    """
    axis=0: normalize each column; 1: normalize each row
    """
    mat = np.asmatrix(arr)
    return np.asarray(np.nan_to_num(mat / np.sum(mat, axis)))

def normalize_fracs(arr, norm_rows=True, norm_cols=True):
    if norm_cols:
        # Normalize columns first--seems correct for overall elution profile if
        # you're calculating correlation-type scores
        arr = arr_norm(arr, 0)
    if norm_rows:
        arr = arr_norm(arr, 1)
    return arr

class Apex(object):
	def __init__(self):
		self.name="apex"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))


	def setAll(self, elution):
		self.apex_array = np.argmax(elution, axis=1)
		self.shape = (len(self.apex_array),len(self.apex_array))

	def calculateScore(self,a,b):
		self.setAll(np.array([a,b]))
		return self.apex_scores_toarray_fast()[0][1]

	def __getitem__(self, index):
		 return int(self.apex_array[index[0]] == self.apex_array[index[1]])

	def apex_scores_toarray_fast(smat):
		dmaxes = defaultdict(set)
		for row, mx in enumerate(smat.apex_array):
			dmaxes[mx].add(row)
		arr = np.zeros(smat.shape)
		for mx,rows in dmaxes.items():
			for r1,r2 in itertools.permutations(rows,2):
				arr[r1,r2] = 1
		return arr

class Wcc:
	def __init__(self):
		self.name="wcc"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a,b):
		rpackages.importr('wccsom')
		r_wcc = robjects.r['wcc']
		return r_wcc(robjects.FloatVector(a), robjects.FloatVector(b), 20)[0]

class Poisson:
	def __init__(self):
		self.name="poisson"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))


	def calculateScore(self, a,b):
		return traver_corr(np.asmatrix([a,b]), verbose=False)[0][1]

def traver_corr(mat, repeat=10, norm='columns', verbose=True):
# As described in supplementary information in paper.
# Randomly draw from poisson(C=A+1/M) for each cell
# where A = the observed count and M is the total
# normalize each column to sum to 1 
# then correlate, and average together for repeat tries.
	def poisson_corr(mat, iteration_display, norm):
		if verbose: print iteration_display
		M = mat.shape[1]
		C = mat + 1/M
		poisson_mat = np.matrix(np.zeros(C.shape))
		for i in range(C.shape[0]):
			for j in range(M):
				poisson_mat[i,j] = np.random.poisson(C[i,j])
		if norm=='columns':
			poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
		elif norm=='rows': # seems to make no performance difference 1/25
			poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
		corr = np.nan_to_num(np.corrcoef(poisson_mat))
		return corr
	avg_result = (reduce(operator.add, (poisson_corr(mat, i, norm=norm) for i in range(repeat))) / repeat)
	return avg_result



class Jaccard():
        def __init__(self):
                self.name="Jaccard"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a,b):
		j11 = 0
		j01 = 0
		for i in range(len(a)):
			if a[i] == 0 and b[i]==0: continue
			if a[i] > 0 and b[i] > 0 :
				j11 += 1
			else:
				j01 += 1
		if j11+j01 > 0:
			return j11/(j11+j01)
		else:
			return 0

class Pearson:
        def __init__(self):
                self.name = "Pearson"

        def getScores(self, a, b, elutionData):
                return (elutionData.getElution(a), elutionData.getElution(b))


        def calculateScore(self, a,b):
                return scipy.stats.pearsonr(a, b)[0]

class Euclidiean:
        def __init__(self):
                self.name = "Euclidiean"

        def getScores(self, a, b, elutionData):
                return (elutionData.getElution(a, normed=True), elutionData.getElution(b, normed=True))


        def calculateScore(self, a,b):
                return 1-distance.euclidean(a,b)

class CalculateCoElutionScores():
	def __init__(self, elutionData=""):
		self.elutionData = elutionData
		self.scores = {}
		self.header = ["ProtA","ProtB"]
		self.oneDscores = [Poisson(), Pearson(), Wcc(), Apex(), Jaccard(), Euclidiean()]
		

	def getAllPairs(self):
		allprots = self.elutionData.prot2Index.keys()
		allPPIs = set([])
		for i in range(len(allprots)):
			for j in range(i+1, len(allprots)):
				protA = allprots[i]
				protB = allprots[j]
				if protA == protB: continue
				if protA > protB:
					allPPIs.add((protA, protB, "?"))
				else:
					allPPIs.add((protB, protA, "?"))
					
		return list(allPPIs)
		

	def calculateAllPairs(self, scores = ""):
		if scores == "": scores = self.twoDscores
		print scores
		allprots = self.getAllPairs()
		self.calculateAllScores(scores, allprots)		

	def calculate1DScores(self, PPIs):
		self.calculateAllScores(self.oneDscores, PPIs)

	def calculateAllScores(self, scoreTypes, PPIs):
		for scoreType in scoreTypes:
			self.header = np.append(self.header, scoreType.name)
		for protA, protB, label in PPIs:
			if not self.elutionData.hasProt(protA): continue
			if not self.elutionData.hasProt(protB): continue
			if (protA, protB, label) not in self.scores: self.scores[(protA, protB, label)] = []
			for scoreType in scoreTypes:
#				print "Calculating %s scores" % (scoreType.name)
				profileA, profileB = scoreType.getScores(protA, protB, self.elutionData)
				self.scores[(protA, protB, label)] = np.append(self.scores[(protA, protB, label)], scoreType.calculateScore(profileA, profileB))

	def toTable(self, labels=True):
		out = "\t".join(self.header)
		if labels: out += "\tClass"
		for protA, protB, label in self.scores:
			out += "\n%s\t%s\t%s" % (protA, protB, "\t".join(map(str,self.scores[(protA, protB, label)])))
			if labels: out += "\t" + label
		return out
	
	def toSklearnData(self):
		data = []
		targets = []
		for idA, idB, label in self.scores:
			if label == "positive": targets.append(1)
			if label == "negative": targets.append(0)
			if label == "?": targets.append("?")
			data.append(self.scores[(idA, idB, label)])
		return np.array(data), np.array(targets)

class RandomForest:
	def __init__(self, data, targets):
		self.data = data
		self.targets = targets #label_binarize(targets, classes=[0, 1])
		self.rfc = RandomForestClassifier(n_estimators=100)
		self.rfc.fit(data, targets)
		
	def kFoldCV(self, folds=10):
		return cross_val_predict(self.rfc, self.data, self.targets, cv=10)

	def getValScores(self, folds=10):
		preds = self.kFoldCV(folds)
		precision, recall, fmeasure = precision_recall_fscore_support(self.targets, preds, pos_label= 1, average='binary')[:3]
		auc_pr = average_precision_score(self.targets, preds)
		auc_roc = roc_auc_score(self.targets, preds) 
		return [precision, recall, fmeasure, auc_pr, auc_roc]

	def getPRcurve(self, folds=10):
		preds = self.kFoldCV(folds)
		return precision_recall_curve(self.targets, preds)
	
	def predict(self, toPred):
		return self.rfc.predict_proba(toPred)


def getPRcurve(elutionData, reference, score):
        scoreCalc = CalculateCoElutionScores(elutionData)
        if type(score)!=list:
                scoreCalc.calculateAllScores([score], reference)
        else:
                scoreCalc.calculateAllScores(score, reference)
        data, targets = scoreCalc.toSklearnData()
        clf = RandomForest(data, targets)
        precision, recall, _ =  clf.getPRcurve()
        return precision, recall


def plotPRcurve(prcurves, outF):
        plt.clf()
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        cols = ['b', 'r', 'c', 'm', 'y', 'k', '#0000ff', '#005D00', '#00FF00', '#8B5D00', '#FF6600', '#FF0066', '#5c5c8a']
        for (name, precision, recall) in prcurves:
                plt.plot(recall, precision, label=name, color = cols.pop())
        plt.legend(loc="upper right", ncol = 5, fontsize=8)
        plt.savefig(outF)



def loadScoreData(scoreF, refF):
	dataRef, headerRef = readTable(refF)
	dataScores, headerScores = readTable(scoreF)
	tolearn = CalculateCoElutionScores()
	topred = CalculateCoElutionScores()
	tolearn.header = headerScores
	topred.header = headerScores
	for edge in dataScores:
		label = "?"
		scores = dataScores[edge]
		protA, protB = edge
		if edge in dataRef:
			label = dataRef[edge][0]
			tolearn.scores[(protA, protB, label)] = scores
		else:
			topred.scores[(protA, protB, label)] = scores
	return tolearn, topred
	
def readTable(tabelF):
	data = {}
	header = {}
	tabelFH = open(tabelF)
	header = tabelFH.readline()
	header = header.rstrip()
	header = header.split("\t")
	for line in tabelFH:
		line = line.rstrip()
		line = np.array(line.split("\t"))
		key = tuple(sorted(line[:2]))
		data[key] = line[2:]
	tabelFH.close()
	return data, header

def loadEData(elutionProfileF):
	elutionData = ElutionData(elutionProfileF)
	scoreCalc = CalculateCoElutionScores(elutionData)
	return elutionData, scoreCalc

def loadData(goldstandardF, elutionProfileF):
	reference = readGoldStandard(goldstandardF)
	elutionData, scoreCalc = loadEData(elutionProfileF)
	return reference, elutionData, scoreCalc

def trainML(scoreCalc):
	data, targets = scoreCalc.toSklearnData()
	clf = RandomForest(data, targets)
	return clf

def readGoldStandard(refF):
	reference = set([])
	goldstandardFH = open(refF)
	for line in goldstandardFH:
		line = line.rstrip()
		(tmpidA, tmpidB, label) = line.split("\t")
		idA, idB = sorted([tmpidA, tmpidB])
		reference.add((idA, idB, label))
	return reference

def main():
	(elutionProfilesF, refF, outF) = sys.argv[1:]
	reference, elutionData, scoreCalc = loadData(refF, elutionProfilesF)
	scores = Euclidiean()
        precision, recall = getPRcurve(elutionData, reference, scores)
	plotPRcurve([tuple(["Euclidean", precision, recall])], outF)

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pas
