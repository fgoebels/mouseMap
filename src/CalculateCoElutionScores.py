#!/usr/bin/env python
from __future__ import division
import numpy as np
import scipy.stats
import utils 
import sys
import math
import MDAnalysis.analysis.psa
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from scipy.spatial import distance
import itertools
import operator
import copy
from sklearn import cross_validation
from sklearn.cross_validation import StratifiedKFold
from collections import defaultdict
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support
from sklearn.cross_validation import train_test_split, cross_val_predict
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import inspect
import os
from sklearn import svm
from sklearn import datasets
from sklearn import metrics
import random

subfldr = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"TCSS")))
sys.path.append(subfldr)

from main import load_semantic_similarity, calculate_semantic_similarity
#from ontology import GOGraph

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
		

	def getElution(self, prot, reshape = False, normed=False):
		if prot not in self.prot2Index:
			return float('NaN')
		tmpout = ""
		if normed: 
			tmpout = self.normedElutionMat[self.prot2Index[prot]]
		else:
			tmpout = self.elutionMat[self.prot2Index[prot]]
		if reshape: tmpout = tmpout.reshape(5,48)
		return tmpout
	

	def hasProt(self, prot):
		return prot in self.prot2Index

	def getProtIndex(self, prot):
		return self.prot2Index[prot]

	def getRandomSubSet(self, size):
		fractions = np.random.choice(range(np.shape(self.elutionMat)[1]), size, replace=False)
		self.elutionMat =  self.elutionMat[:, fractions]
		return fractions

	def getSubset(self, fracs):
		self.elutionMat =  self.elutionMat[:, fracs]
		self.normedElutionMat =  self.normedElutionMat[:, fracs]

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

class GOSim(object):

	objs = "" 
	
	def __init__(self, onto_F, gene_anno_F):
		self.name = ["Sim_CC", "Sim_BP", "Sim_MF"]
#		if GOSim.objs == "":
#			GOSim.objs = load_semantic_similarity(onto_F, gene_anno_F, "C:2.4,P:3.5,F:3.3", "")

	def getScores(self, a, b, elutionData):
		return (a,b)

	def calculateScore(self, a, b):
		out = []
		domain_def = {'C':'Cellular Component', 'P':'Biological Process', 'F':'Molecular Function'}
		for domain in domain_def:
			score = GOSim.objs[domain]._semantic_similarity(a, b)[0]
			if score is None: score = 0
			out.append(score)
		return out

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
	def __init__(self, repeat=100):
		self.name="poisson-%i" % (repeat)
		self.repeat=repeat
		self.noiseMats = []

	def getScores(self, a, b, elutionData):
		if len(self.noiseMats)> 0 and elutionData.elutionMat.shape != self.noiseMats[0].shape:
			self.noiseMats = []
		
		if len(self.noiseMats)<self.repeat:
			for i in range(self.repeat):
				self.noiseMats.append(self.makenoisyMat(elutionData.elutionMat))
		return (elutionData.getProtIndex(a), elutionData.getProtIndex(b))
		

	def calculateScore(self, a,b):
		out = []
		for mat in self.noiseMats:
			profile_a = mat[a].getA1()
			profile_b = mat[b].getA1()
			mat_correlation = scipy.stats.pearsonr(profile_a, profile_b)[0]
			out.append(mat_correlation)
		return sum(out)/len(out)
#		return traver_corr(np.asmatrix([a,b]), verbose=False)[0][1]

	def makenoisyMat(self, mat):
		M = mat.shape[1]
		C = mat + 1/M
		poisson_mat = np.matrix(np.zeros(C.shape))
		for i in range(C.shape[0]):
			for j in range(M):
				poisson_mat[i,j] = np.random.poisson(C[i,j])
		poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
		return poisson_mat

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


class MutualInformation():
	def __init__(self, minCounts = 2):
		self.name="MI"
		self.minCounts = minCounts
	
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a, b):
		numFracs = len(a)
		(a_upper, a_lower) = self.getFracs(a, self.minCounts)
		(b_upper, b_lower) = self.getFracs(b, self.minCounts)
		entropy_a = self.bin_entropy(len(a_upper)/numFracs)
		entropy_b = self.bin_entropy(len(b_upper)/numFracs)
		joint_probs = np.array(map(lambda x: len(x)/numFracs, [a_upper&b_upper, a_upper&b_lower, a_lower&b_upper, a_lower&b_lower]))
		joint_entropy_a_b = self.entropy(joint_probs, 2)
		mutual_information =  entropy_a  + entropy_b - joint_entropy_a_b
		return mutual_information

	def bin_entropy(self, p):
		return self.entropy(np.array([p,1-p]))

	def entropy(self, probs, base=0):
		if base ==0: base = len(probs)
		tmp_probs = probs
		tmp_probs[tmp_probs==0] = 1
		return -sum(probs*map(lambda x: math.log(x,base), tmp_probs))

	def getFracs(self, a, cutoff):
		upper = set([i for i,v in enumerate(a) if v > cutoff])
		lower = set([i for i,v in enumerate(a) if v <= cutoff])
		return (upper, lower)

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

class Hausdorff:		
	def __init__(self):
		self.name="Hausdorff"
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, reshape=True), elutionData.getElution(b, reshape=True))

	def calculateScore(self, a,b):
		a3D = self.toPointMat(a)
		b3D = self.toPointMat(b)
		return MDAnalysis.analysis.psa.hausdorff(a3D, b3D, 240)

	def toPointMat(self, a):
		out = np.ndarray(shape=(5*48,3))
		n = 0
		for i in range(5):
			for j in range(48):
				out[n][0] = i
				out[n][1] = i
				out[n][2] = a[i,j]
				n += 1
		return out

class CorrEigenVals:
	def __init__(self):
		self.name = "CorrEigenVals"	

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, reshape=True), elutionData.getElution(b, reshape=True))


	def calculateScore(self, a, b):
		vA = np.linalg.svd(a)[1]
		vB = np.linalg.svd(b)[1]
		return scipy.stats.pearsonr(vA, vB)[0]
class Herdin:
	def __init__(self):
                self.name = "Herdin"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, reshape=True), elutionData.getElution(b, reshape=True))


	def calculateScore(self, a,b):
		return np.trace(np.dot(a,np.transpose(b)))/(np.linalg.norm(a)*np.linalg.norm(b))
class Pearson:
	def __init__(self):
		self.name = "Pearson"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))


	def calculateScore(self, a,b):
		score = scipy.stats.pearsonr(a, b)[0]
		if math.isnan(score): return 0.0
		return scipy.stats.pearsonr(a, b)[0]
class Euclidiean:
	def __init__(self):
		self.name = "Euclidiean"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, normed=True), elutionData.getElution(b, normed=True))


	def calculateScore(self, a,b):
		return 1-distance.euclidean(a,b)

class MatrixNorms:
	def __init__(self, name = ["RowMax","RowMin","2-Norm","Smallest_singular_value","Frobenius_norm","ColMax","ColMin"]):
		self.name = name

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, reshape=True, normed=True), elutionData.getElution(b, reshape=True, normed=True))

	def calculateScore(self, a,b):
		out = []
		elutionAMinusB = np.subtract(a, b)	
	        out.append( np.linalg.norm(elutionAMinusB, 1))
	        out.append( np.linalg.norm(elutionAMinusB, -1))
	        out.append( np.linalg.norm(elutionAMinusB, 2))
	        out.append( np.linalg.norm(elutionAMinusB, -2))
	        out.append( 1 - np.linalg.norm(elutionAMinusB))
	        out.append(np.linalg.norm(np.matrix.transpose(elutionAMinusB), 1))
	        out.append(np.linalg.norm(np.matrix.transpose(elutionAMinusB), -1))
		return out

class CalculateCoElutionScores():
	def __init__(self, elutionData=""):
		self.elutionData = elutionData
		self.scores = {}
		self.header = ["ProtA","ProtB"]
		self.oneDscores = [Poisson(), Pearson(), Wcc(), Apex(), Jaccard(), Euclidiean()]
#		self.twoDscores = [Euclidiean()] 
		self.twoDscores = [Pearson(), Wcc(), Apex(), Jaccard(), Euclidiean(), Herdin(), MatrixNorms(), GOSim("src/TCSS/gene_ontology.obo.txt", "Yeast/data/gene_association.tab")]
		
	def mergeScoreCalc(self, toMerge):
		numFeature_in_merge = len(toMerge.scores[toMerge.scores.keys()[0]])
		numFeature_in_self = len(self.scores[self.scores.keys()[0]])
		for edge in toMerge.scores:
			if edge in self.scores:
				self.scores[edge] = np.append(self.scores[edge], toMerge.scores[edge])
			else:
				self.scores[edge] = np.append(np.array([0]*numFeature_in_self), toMerge.scores[edge])
		for edge in self.scores:
			if edge not in toMerge.scores:
				self.scores[edge] = np.append(self.scores[edge], np.array([0]*numFeature_in_merge))

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

	def calculate2DScores(self, PPIs):
		self.calculateAllScores(self.twoDscores, PPIs)

	def calculateAllScores(self, scoreTypes, PPIs):
		for scoreType in scoreTypes:
			self.header = np.append(self.header, scoreType.name)
		for protA, protB, label in PPIs:
			if (protA, protB, label) not in self.scores: self.scores[(protA, protB, label)] = []
			if not self.elutionData.hasProt(protA) or not self.elutionData.hasProt(protB):
				self.scores[(protA, protB, label)] = np.append(self.scores[(protA, protB, label)], 0)
				continue
			for scoreType in scoreTypes:
				profileA, profileB = scoreType.getScores(protA, protB, self.elutionData)
				self.scores[(protA, protB, label)] = np.append(self.scores[(protA, protB, label)], scoreType.calculateScore(profileA, profileB))

	def toTable(self, labels=True):
		out = "\t".join(self.header)
		if labels: out += "\tClass"
		for protA, protB, label in self.scores:
			out += "\n%s\t%s\t%s" % (protA, protB, "\t".join(map(str,self.scores[(protA, protB, label)])))
			if labels: out += "\t" + label
		return out
	
	def toArffData(self):
		out = ["@RELATION COFrac"]
		for colname in self.header[2:]:
			out.append("@Attribute %s NUMERIC" % (colname))
		out.append("@ATTRIBUTE class {positive, negative}")
		out.append("@Data")
		for idA, idB, label in self.scores:
			tmp = ",".join(map(str,self.scores[(idA, idB, label)]))
			tmp += ",%s" % (label)
			out.append(tmp)
		return "\n".join(out)
		

	def toSklearnData(self):
		data = []
		targets = []
		for idA, idB, label in self.scores:
			if label == "positive": targets.append(1)
			if label == "negative": targets.append(0)
			if label == "?": targets.append("?")
			data.append(self.scores[(idA, idB, label)])
		return np.nan_to_num(np.array(data)), np.array(targets)

class CLF_Wrapper:
	def __init__(self, data, targets, forest=False):
		self.data = data
		self.targets = targets #label_binarize(targets, classes=[0, 1])
		if forest:
			self.clf = RandomForestClassifier(n_estimators=100)
		else:	
			self.clf = svm.SVC(kernel="linear", probability=True)
		self.clf.fit(data, targets)
		
	def kFoldCV(self, folds=10):
		folds = StratifiedKFold(self.targets, folds)
		return cross_val_predict(self.clf, self.data, self.targets, cv=folds)

	def getValScores(self, folds=10):
#		return cross_validation.cross_val_score(self.clf, self.data, self.targets, cv=10, scoring='f1')
		preds = self.kFoldCV(folds)
		precision = metrics.precision_score(self.targets, preds, average=None)[1]
		recall = metrics.recall_score(self.targets, preds, average=None)[1]
		fmeasure = metrics.f1_score(self.targets, preds, average=None)[1]
		auc_pr = average_precision_score(self.targets, preds)
		auc_roc = roc_auc_score(self.targets, preds) 
		return [precision, recall, fmeasure, auc_pr, auc_roc]

	def getPRcurve(self, folds=10):
		all_probas = []
		all_targets = []
		for train, test in StratifiedKFold(self.targets, folds):
			probas = self.clf.fit(self.data[train], self.targets[train]).predict_proba(self.data[test])
			all_probas.extend(probas[:,1]) # make sure that 1 is positive class in binarizied class vector
			all_targets.extend(self.targets[test])
		return precision_recall_curve(all_targets, all_probas)
#		return roc_curve(all_targets, all_probas)
	
	def predict(self, toPred):
		return self.clf.predict_proba(toPred)

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
	elutionData, scoreCalc = loadEData(elutionProfileF)
	reference = readGoldStandard(goldstandardF, [elutionData])
	return reference, elutionData, scoreCalc

def readGoldStandard(refF, elutionData, ratio=5):
	positive = set([])
	negative = set([])
	protsWithElutionProfile = set([])
	for data in elutionData:
		protsWithElutionProfile = protsWithElutionProfile | set(data.prot2Index.keys())
	reference = set([])
	goldstandardFH = open(refF)
	for line in goldstandardFH:
		line = line.rstrip()
		(tmpidA, tmpidB, label) = line.split("\t")
		if tmpidA not in protsWithElutionProfile or tmpidB not in protsWithElutionProfile: continue
		idA, idB = sorted([tmpidA, tmpidB])
		if label == "positive": positive.add((idA, idB, label))
		if label == "negative": negative.add((idA, idB, label))

	if len(positive)*ratio>len(negative):
		print "Warning: not enough negative data points in reference to create desired ratio"
		return positive | negative
	
	reference = positive | set(random.sample(negative, len(positive)*ratio))
	return reference


def benchMarkFracSize(elutionProfileF, goldstandardF, size, runs, outF):
	reference, elutionData, _ = loadData(goldstandardF, elutionProfileF)
	outFH = open(outF, "w")
	size = int(size)
	runs = int(runs)
	for i in range(runs):
		tmpElutData = copy.copy(elutionData)
		tmpElutData.normalizeCoEulitionMat()
		fractions = tmpElutData.getRandomSubSet(size)
		scoreCalc = CalculateCoElutionScores(tmpElutData)
		numPPIs =  scoreCalc.calculateAllPairs()
		scoreCalc.calculateBenchmarkScore(reference)
		data, targets = scoreCalc.toSklearnData()
		clf = CLF_Wrapper(data, targets)
		scores =  clf.getValScores()
		print >> outFH, "%i\t%i\t%s\t%s" % (size, numPPIs, ",".join(map(str, fractions)), "\t".join(map(str, scores)))
	outFH.close()

def getPRcurve(elutionData, reference, score):
	scoreCalc = CalculateCoElutionScores(elutionData)
	if type(score)!=list:
		scoreCalc.calculateAllScores([score], reference)
	else:
		scoreCalc.calculateAllScores(score, reference)
	data, targets = scoreCalc.toSklearnData()
	clf = CLF_Wrapper(data, targets)
	precision, recall, _ =  clf.getPRcurve()
	return score.name, precision, recall

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

def featurePRcurves(elutionProfileF, goldstandardF, outF):
	prCurves = []
	reference, elutionData, scoreCalc = loadData(goldstandardF, elutionProfileF)
	scores = [Poisson(), Pearson(), Wcc(), Apex(), Jaccard(), Euclidiean(), Hausdorff(), CorrEigenVals(), Herdin(), MatrixNorms(name = "Matrix norms")]
	for score in scores:
		precision, recall = getPRcurve(elutionData, reference, score)
		prCurves.append(tuple([score.name, precision, recall]))
	precision, recall = getPRcurve(elutionData, reference, scores)
	prCurves.append(tuple(["combined", precision, recall]))
	plotPRcurve(prCurves, outF)

def benchamrkElutionProfiels(profilesF, outF):
	profilesFH = open(profilesF)
	profilesFH.readline()
	prCurves = []
	for line in profilesFH:
		line = line.rstrip()
		print line
		(refF, elutionF, elutionType, name) = line.split("\t")
		elutionF = elutionF.split(",")
		reference, elutionData, scoreCalc = loadData(refF, elutionF)
		print elutionData.elutionMat.shape
#		scores = scoreCalc.oneDscores
#		if elutionType == "2D": scores = scoreCalc.twoDscores
		scores = Euclidiean()
		precision, recall = getPRcurve(elutionData, reference, scores)
		prCurves.append(tuple([name, precision, recall]))
	plotPRcurve(prCurves, outF)

def benchmarkOneProfile(elutionProfileF, expType, goldstandardF, outF):
#	featurePRcurves(elutionProfileF, goldstandardF, outF)
	reference, elutionData, scoreCalc = loadData(goldstandardF, elutionProfileF)
	if expType == "2D":
		scoreCalc.calculate2DScores(reference)
	else:
		scoreCalc.calculate1DScores(reference)
	print "Done calc features"
	data, targets = scoreCalc.toSklearnData()
	clf = CLF_Wrapper(data, targets)
	scores =  clf.getValScores()
	print scores

def main():
	refF, elutionFiles = sys.argv[1:]
	elutionFH = open(elutionFiles)
	scoreCals = []
	elutionDatas = []
	
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = ElutionData(elutionFile)
		elutionDatas.append(elutionData)

	reference = readGoldStandard(refF, elutionDatas)
	out = []
	global_all_scoreCalc = []
#	for score in [Euclidiean(), Pearson(), Wcc(), Apex(), Jaccard(), MutualInformation(2)]:
	for score in [MutualInformation(2), Euclidiean(), Pearson(), Wcc(), Apex(), Jaccard(), Poisson(10)]:
		for elutionD in elutionDatas:
			scoreCalc = CalculateCoElutionScores(elutionD)
			scoreCalc.calculateAllScores([score], reference)
			scoreCals.append(scoreCalc)
			global_all_scoreCalc.append(scoreCalc)
		all_scoreCalc = scoreCals[0]
		for i in range(1,len(scoreCals)):
			all_scoreCalc.mergeScoreCalc(scoreCals[i])
		data, targets = all_scoreCalc.toSklearnData()
		clf = CLF_Wrapper(data, targets)
		print clf.getValScores()
		precision, recall, _ =  clf.getPRcurve()
		out.append((score.name, precision, recall))
		scoreCals = []
	all_scoreCalc = global_all_scoreCalc[0]
	for i in range(1, len(global_all_scoreCalc)):
		all_scoreCalc.mergeScoreCalc(global_all_scoreCalc[i])

	data, targets = all_scoreCalc.toSklearnData()	
	clf = CLF_Wrapper(data, targets)
	print clf.getValScores()
	precision, recall, _ =  clf.getPRcurve()
	out.append(("combined", precision, recall))

	plotPRcurve(out, "test/Worm_1D_MI.pdf")
	

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pas
