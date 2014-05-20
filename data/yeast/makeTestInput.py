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
parser.add_option("--subnet",dest="subnet",action="store_true",default=False)
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

def printGO(fh, vals, map, edge_universe, consider=None, symmetric=False):

	done = set()
	for (geneA, geneB, val) in vals:
		
		if (geneA, geneB) in done:
			continue		

		if (geneA, geneB) not in edge_universe and (geneB, geneA) not in edge_universe:
			continue

		if geneA not in map or geneB not in map:
			continue
		
		if geneA == geneB:
			continue
		
		done.add( (geneA, geneB) )
		done.add( (geneB, geneA) )

		fh.write("\t".join( [map[geneA], map[geneB], val] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], val] )+"\n")

	fh.close()	

def printNet(fh, edges, map, edge_universe, consider=None, symmetric=True):

	done = set()
	for (geneA, geneB) in edges:

		if (geneA, geneB) in done:
			continue		

		if (geneA, geneB) not in edge_universe and (geneB, geneA) not in edge_universe:
			continue

		if geneA not in map or geneB not in map:
			continue

		if geneA == geneB:
			continue

		val = edges[(geneA, geneB)]

		done.add( (geneA, geneB) )
		done.add( (geneB, geneA) )

		fh.write("\t".join( [map[geneA], map[geneB], val] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], val] )+"\n")

	fh.close()	

def printEL(fh, edges, map, edge_universe, consider=None, symmetric=True):

	done = set()
	for (geneA, geneB) in edges:
	
		if (geneA, geneB) in done:
			continue		

		if (geneA, geneB) not in edge_universe and (geneB, geneA) not in edge_universe:
			continue

		if geneA not in map or geneB not in map:
			continue

		if geneA == geneB:
			continue
		
		done.add( (geneA, geneB) )
		done.add( (geneB, geneA) )

		fh.write("\t".join( [map[geneA], map[geneB]] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA]] )+"\n")

	fh.close()	

def randomNegEdge(g):

	first_node = random.choice(g.nodes())
	possible_nodes = set(g.nodes())
	nbrs = g.neighbors(first_node) + [first_node]
	possible_nodes.difference_update(nbrs) 
	second_node = random.choice(list(possible_nodes))   

	return (first_node, second_node)

def splitThirds(edges, g):
	"""
	Split into training set, learning set, and the held-out test set
	"""
	sample_size = int(len(edges) * 0.33)
	edge_l = list(edges)
	# indexes for edge_l of training sets
	# take the first third
	#train_indexes = random.sample(range(0, len(edge_l)), sample_size)
	train_indexes = range(0, sample_size)
	# all other indexes
	learn_indexes = range(0, sample_size*2) 

	# training set is just the slObserved used to train the initial model
	training_set = {}
	# the learning set is used for weight learning: includes labels for all
	# known sl observed and sl observations
	learning_set = {}
	# the test set is completely held out from PSL and used for evaluation only
	test_set = {}
	for i in range(0, len(edge_l)):
		if i % 10000 == 0:
			print "..."

		if i in train_indexes:
			training_set[edge_l[i]] = "1.0"
			learning_set[edge_l[i]] = "1.0"
		elif i in learn_indexes:
			learning_set[edge_l[i]] = "1.0"
		else:
			test_set[edge_l[i]] = "1.0"

	# now select a random set of edges between nodes in 'edges'
	# that is of the same size as 'edges', split into thirds
	# assigning priors to each
	print "selecting random negative training cases..."
	for i in range(0, len(edge_l)):
		neg_edge = randomNegEdge(g)
		if i % 10000 == 0:
			print "..."

		if i in train_indexes:
			training_set[neg_edge] = "0.01"
			learning_set[neg_edge] = "0.01"
		elif i in learn_indexes:
			learning_set[neg_edge] = "0.01"
		else:
			test_set[neg_edge] = "0.01"

	print "done!"
	return (training_set, learning_set, test_set)

def getSubnet(g, min_edges=1000):

	print "Growing a subnetwork..."
	seed_node = random.choice(g.nodes())
	first_neighbors = g.neighbors(seed_node)
	first_neighbors.append(seed_node)

	# add all interconnections between this set
	h = g.subgraph(first_neighbors)
	# grow the subnetwork until we have at least min_edges
	while (len(h.edges()) < min_edges):
		seed = random.choice(h.nodes())
		nbrs = g.neighbors(seed)
		nbrs.append(seed)
		for n in h.nodes():
			nbrs.append(n)
		h = g.subgraph(nbrs)

	print "Selected a random subnetwork of "+str(len(h.edges()))+" edges..."
	return h

	
# subset the data to look at just these edges
consider_nodes = None

# get mappings, name to id
name2id, edges = parseNet(opts.slgraph)

print "Building the graph..."

g = nx.Graph()
g.add_edges_from(edges)

h = None
# test only a subnetwork if this is a test run 
if opts.subnet:
	h = getSubnet(g)
else:
	h = g

edge_universe = h.edges()
node_universe = h.nodes()

# if using a subset of all nodes, select a random
# node and grow it by first-neighbors, then add
# all inter-connections to get a maximally connected
# subgraph centered at that node.

print "Parsing data.."

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
print "Doing data split..."
slTrain, slLearn, slTest = splitThirds(edge_universe, h)

edge_universe = {}
for edge in slLearn:
	edge_universe[edge] = 1
for edge in slTest:
	edge_universe[edge] = 1
edge_universe = set(edge_universe.keys())

print len(edge_universe)

print "writing..."

out = opts.train
fh = open(out+'gene.txt', 'w')
for name in name2id:

	#if name2id[name] not in node_universe:
	#	continue
		
	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id, edge_universe)

fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id, edge_universe)
fh = open(out+'goMF.txt', 'w')
printGO(fh, goCC, name2id, edge_universe)

fh = open(out+'ppiEdges.txt', 'w')
printNet(fh, ppiEdges , name2id, edge_universe)

fh = open(out+'slObserved.txt', 'w')
printNet(fh, slTrain , name2id, edge_universe)

fh = open(out+'sl.txt', 'w')
printNet(fh, slLearn, name2id, edge_universe)

# print just the nodes to consider for SL learning
fh = open(out+'consider.txt', 'w')
printEL(fh, slLearn , name2id, edge_universe)

fh = open(out+'ppiKernel.txt', 'w')
printGO(fh, ppiKernel, name2id, edge_universe)

fh = open(out+'negKernel.txt', 'w')
printGO(fh, gnegKernel , name2id, edge_universe)

out = opts.test
fh = open(out+'gene.txt', 'w')
for name in name2id:

	if name2id[name] not in node_universe:
		continue

	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id, edge_universe)

fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id, edge_universe)
fh = open(out+'goMF.txt', 'w')
printGO(fh, goCC, name2id, edge_universe)

fh = open(out+'ppiEdges.txt', 'w')
printNet(fh, ppiEdges , name2id, edge_universe)

# add all the data for testing here, to infer new edges
fh = open(out+'slObserved.txt', 'w')
printNet(fh, slLearn, name2id, edge_universe)

fh = open(out+'ppiKernel.txt', 'w')
printGO(fh, ppiKernel, name2id, edge_universe)

fh = open(out+'negKernel.txt', 'w')
printGO(fh, gnegKernel , name2id, edge_universe)

fh = open(out+'heldOutSL.tab', 'w')
printNet(fh, slTest, name2id, edge_universe)

# infer just these values
fh = open(out+'consider.txt', 'w')
printEL(fh, slTest, name2id, edge_universe)



