#!/usr/bin/env	python

"""

 Need to map all values to somewhere in 0-1 range for PSL. A good number need to be in the .66+ range here, 
 it's a lot like paradigm
"""
import networkx as nx
import math

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--geneMap",dest="map",action="store",type="string",default=None)
parser.add_option("--inferences",dest="inferences",action="store",type="string",default=None)
parser.add_option("--truth",dest="truth",action="store",type="string",default=None)
(opts, args) = parser.parse_args()


PSEUDO_COUNT = 0.01

def parseMap(file, network_nodes=None):
	
	map = {}
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	for line in fh:
		index, name = line.rstrip().split("\t")
		map[index] = name

	fh.close()
	return map

def parseInf(file):

	inf = {}
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	for line in fh:

		if not line.startswith("SL("):
			continue

		try:
			edge, val = line.rstrip().split("\t")
			idxA, idxB = edge.split('(')[1].rstrip(')').split(", ")
			inf[(idxA, idxB)] = float(val)
		except:
			continue

	fh.close()
	return inf

def parseTruth(file):

	edges = {}
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	for line in fh:
		geneA, geneB, val = line.rstrip().split("\t")
		edges[(geneA, geneB)] = float(val)


	return edges

def getRates(inferences, truth, threshold):

	TP = 0
	FP = 0
	for (A, B) in inferences:

		if inferences[(A,B)] < threshold:
			continue

		if (A,B) in truth and truth[(A,B)] == 1:
			TP += 1
		elif (B,A) in truth and truth[(B,A)] == 1:
			TP += 1
		else:
			FP += 1
		
	return (TP, FP)	

truth = parseTruth(opts.truth)
inferences = parseInf(opts.inferences)
map = parseMap(opts.map)

total_TP = 0.0
total_FP = 0.0
for pair in truth:
	if truth[pair] == 1:
		total_TP += 1
	else:
		total_FP += 1

# generate ROC curve 
for t in [t/float(100) for t in range(100,-1,-1)]:
	TP, FP = getRates(inferences, truth, t)
	print "\t".join([str(t), str(TP/total_TP), str(FP/total_FP)])

#for (A, B) in inferences:
#	inf_val = inferences[(A,B)]
#	t_val = None
#	if (A,B) in truth:
#		t_val = truth[(A,B)]
#	elif (B,A) in truth:
#		t_val = truth[(B,A)]
##	print "\t".join([map[A], map[B], inf_val, t_val])
