#!/usr/bin/env	python

"""
Just make a test/train split for code-testing purposes
"""
import networkx as nx
import math
import random

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--trainFolder",dest="train",action="store",type="string",default="TEST/train/")
parser.add_option("--testFolder",dest="test",action="store",type="string",default="TEST/test/")
parser.add_option("--slgraph",dest="slgraph",action="store",type="string",default=None)
(opts, args) = parser.parse_args()

def parseLST(file):
	
	vals = set()
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	lineno = 1
	for line in fh:
		vals.add( line.rstrip() )

		lineno += 1

	fh.close()
	return vals

def parseVals(file, network_nodes=None):
	
	vals = set()
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	lineno = 1
	for line in fh:
		parts = line.rstrip().split("\t")
		vals.add( tuple(parts) )

		lineno += 1

	fh.close()
	return vals

def parseNet(file):

	genes = set()
	edges = {}
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	for line in fh:
		geneA, geneB = line.rstrip().split("\t")

		# no duplicates
		if (geneB, geneA) in edges or (geneA, geneB) in edges:
			continue

		genes.add( geneA )
		genes.add( geneB )
		edges[ (geneA, geneB) ] = "1.0"


	map = {}
	i = 1
	for gene in genes:
		map[gene] = str(i)
		i += 1
	
	return (map, edges)

def printGO(fh, vals, map, consider=None, symmetric=False):

	for (geneA, geneB, val) in vals:

		if consider and (geneA not in consider or geneB not in consider):
			continue

		if geneA not in map or geneB not in map:
			continue
		
		if geneA == geneB:
			continue

		fh.write("\t".join( [map[geneA], map[geneB], val] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], val] )+"\n")

	fh.close()	

def printNet(fh, edges, map, consider=None, symmetric=True):

	for (geneA, geneB) in edges:
		
		if consider and (geneA not in consider or geneB not in consider):
			continue

		if geneA not in map or geneB not in map:
			continue

		if geneA == geneB:
			continue

		val = edges[(geneA, geneB)]

		fh.write("\t".join( [map[geneA], map[geneB], val] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], val] )+"\n")

	fh.close()	

def printEL(fh, edges, map, consider=None, symmetric=True):

	for (geneA, geneB) in edges:
		
		if consider and (geneA not in consider or geneB not in consider):
			continue

		if geneA not in map or geneB not in map:
			continue

		if geneA == geneB:
			continue

		fh.write("\t".join( [map[geneA], map[geneB]] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA]] )+"\n")

	fh.close()	

def naiveDataSplit(edges, train_fraction):

	sample_size = int(len(edges) * train_fraction)
	edge_l = list(edges)
	train_indexes = random.sample(range(0, len(edge_l)), sample_size)
	training_set = {}
	test_set = {}
	for i in range(0, len(edge_l)):
		if i in train_indexes:
			training_set[edge_l[i]] = "1.0"
		else:
			test_set[edge_l[i]] = "1.0"


	return (training_set, test_set)

def splitThirds(edges):
	"""
	Split into training set, learning set, and the held-out test set
	"""
	sample_size = int(len(edges) * 0.33)
	edge_l = list(edges)
	# indexes for edge_l of training sets
	train_indexes = random.sample(range(0, len(edge_l)), sample_size)
	# all other indexes
	other = [i for i in range(0, len(edge_l)) if i not in train_indexes]
	learn_indexes = random.sample(other, sample_size)

	# training set is just the slObserved used to train the initial model
	training_set = {}
	# the learning set is used for weight learning: includes labels for all
	# known sl observed and sl observations
	learning_set = {}
	# the test set is completely held out from PSL and used for evaluation only
	test_set = {}
	for i in range(0, len(edge_l)):
		if i in train_indexes:
			training_set[edge_l[i]] = "1.0"
			learning_set[edge_l[i]] = "1.0"
		elif i in learn_indexes:
			learning_set[edge_l[i]] = "1.0"
		else:
			test_set[edge_l[i]] = "1.0"

	return (training_set, learning_set, test_set)

def addPriors(edges, categories):

	PRIOR = "0.1"
	plus_priors = edges
	for c in categories:
		for (a, b, v) in c:
			if (a,b) not in plus_priors and (b,a) not in plus_priors:
				plus_priors[(a,b)] = PRIOR

	return plus_priors

# subset the data to look at just these edges
consider_nodes = None


# get mappings, name to id
name2id, edges = parseNet(opts.slgraph)

goBP = parseVals("goBP.tab")
goCC = parseVals("goCC.tab")
goMF = parseVals("goMF.tab")

# parse PPI
ppiN2ID, ppiEdges = parseNet("ppi.sgdid.tab")

# parse PPI Kernel
ppiKernel = parseVals("ppi.kernelEdges.tab")

# parse SGD Negative Kernel
gnegKernel = parseVals("gNeg.kernelEdges.tab")

# split to train/test on the target 
#slTrain, slTest = naiveDataSplit(edges, 0.9)
slTrain, slLearn, slTest = splitThirds(edges)
# slTEST needs to include all observed edges, plus the held-out set
# as well as prior values for any edge that has enough evidence to 
# have a grounded rule for it (otherwise we'll get a runtime exception)
# add in any that are just in one of these categories
slLearn = addPriors(slLearn, [goBP, goCC, goMF, ppiKernel, gnegKernel])

out = opts.train
fh = open(out+'gene.txt', 'w')
for name in name2id:
	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id, consider_nodes)

fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id, consider_nodes)
fh = open(out+'goMF.txt', 'w')
printGO(fh, goCC, name2id, consider_nodes)

fh = open(out+'ppiEdges.txt', 'w')
printNet(fh, ppiEdges , name2id, consider_nodes)

fh = open(out+'slObserved.txt', 'w')
printNet(fh, slTrain , name2id, consider_nodes)

fh = open(out+'sl.txt', 'w')
printNet(fh, slLearn , name2id, consider_nodes)

# print just the nodes to consider for SL learning
fh = open(out+'consider.txt', 'w')
printEL(fh, slLearn , name2id, consider_nodes)

fh = open(out+'ppiKernel.txt', 'w')
printGO(fh, ppiKernel, name2id, consider_nodes)

fh = open(out+'negKernel.txt', 'w')
printGO(fh, gnegKernel , name2id, consider_nodes)

out = opts.test
fh = open(out+'gene.txt', 'w')
for name in name2id:
	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id, consider_nodes)

fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id, consider_nodes)
fh = open(out+'goMF.txt', 'w')
printGO(fh, goCC, name2id, consider_nodes)

fh = open(out+'ppiEdges.txt', 'w')
printNet(fh, ppiEdges , name2id, consider_nodes)

# add all the data for testing here, to infer new edges
fh = open(out+'slObserved.txt', 'w')
printNet(fh, slLearn, name2id, consider_nodes)

fh = open(out+'ppiKernel.txt', 'w')
printGO(fh, ppiKernel, name2id, consider_nodes)

fh = open(out+'negKernel.txt', 'w')
printGO(fh, gnegKernel , name2id, consider_nodes)

fh = open(out+'heldOutSL.tab', 'w')
printNet(fh, slTest, name2id, consider_nodes)

# infer just these values
fh = open(out+'consider.txt', 'w')
printNet(fh, slTest, name2id, consider_nodes)

