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
(opts, args) = parser.parse_args()


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
	edges = set()
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
		edges.add( (geneA, geneB) )


	map = {}
	i = 1
	for gene in genes:
		map[gene] = str(i)
		i += 1
	
	return (map, edges)

def printGO(fh, vals, map, symmetric=False):

	for (geneA, geneB, val) in vals:

		if geneA not in map or geneB not in map:
			continue
		
		if geneA == geneB:
			continue

		fh.write("\t".join( [map[geneA], map[geneB], val] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], val] )+"\n")

	fh.close()	

def printNetWeight(fh, edges, map, symmetric=True):

	for (geneA, geneB) in edges:

		if geneA not in map or geneB not in map:
			continue
		
		if geneA == geneB:
			continue


		fh.write("\t".join( [map[geneA], map[geneB], "1.0"] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], "1.0"] )+"\n")

	fh.close()	

def printNet(fh, edges, map, symmetric=True):

	for (geneA, geneB) in edges:

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
	training_set = []
	test_set = []
	for i in range(0, len(edge_l)):
		if i in train_indexes:
			training_set.append(edge_l[i])
		else:
			test_set.append(edge_l[i])


	return (training_set, test_set)

# get mappings, name to id
name2id, edges = parseNet("negGraph.tab")

goBP = parseVals("goBP.tab")
goCC = parseVals("goCC.tab")
goMF = parseVals("goMF.tab")

# parse PPI
ppiN2ID, ppiEdges = parseNet("ppi.sgdid.tab")

out = opts.train
fh = open(out+'gene.txt', 'w')
for name in name2id:
	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id)
fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id)
fh = open(out+'goMF.txt', 'w')
printGO(fh, goCC, name2id)

fh = open(out+'ppiEdges.txt', 'w')
printNet(fh, ppiEdges , name2id)

# split to train/test
slTrain, slTest = naiveDataSplit(edges, 0.9)

fh = open(out+'slObserved.txt', 'w')
printNetWeight(fh, slTrain , name2id)

fh = open(out+'sl.txt', 'w')
printNetWeight(fh, slTest , name2id)

