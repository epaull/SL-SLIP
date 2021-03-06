#!/usr/bin/env	python

"""
Just make a test/train split for code-testing purposes
"""
import networkx as nx
import math
import random
import itertools
import os

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--output",dest="output",action="store",type="string",default="TEST")
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

def parseGRAPH(file):

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
		edges[ (geneA, geneB) ] = "SL"


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

		# take the cube root for GO mapping
		val = math.pow(float(val), 0.333)

		fh.write("\t".join( [map[geneA], map[geneB], str(val)] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], str(val)] )+"\n")

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

		fh.write("\t".join( [map[geneA], map[geneB], str(val)] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], str(val)] )+"\n")

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

def randomNegEdgeSafe(g, exclude):

	A = None
	B = None
	edge = None
	while True:
		try:
			edge = randomNegEdge(g)
			if edge in exclude or (edge[1],edge[0]) in exclude:
				continue
			else:
				break
		except:
			continue


def randomNegEdge(g):

	first_node = random.choice(g.nodes())
	possible_nodes = set(g.nodes())
	nbrs = g.neighbors(first_node) + [first_node]
	possible_nodes.difference_update(nbrs) 
	second_node = random.choice(list(possible_nodes))   

	return (first_node, second_node)

def splitData(folds, trainIdx, learnIdx, testIdx):

	# training set is just the slObserved used to train the initial model
	training_set = {}
	# the learning set is used for weight learning: includes labels for all
	# known sl observed and sl observations
	learning_set = {}
	# the test set is completely held out from PSL and used for evaluation only
	test_set = {}

	for fold in folds:
		for edge in folds[fold]:
			
			if str(fold) in trainIdx:
				training_set[edge] = folds[fold][edge]
				learning_set[edge] = folds[fold][edge]
			elif str(fold) in learnIdx:
				learning_set[edge] = folds[fold][edge]
			else:
				test_set[edge] = folds[fold][edge]
		print "done!"

	return (training_set, learning_set, test_set)

def split(edges, g):
	"""
	Split into training set, learning set, and the held-out test set
	"""
	sample_size = int(len(edges) * 0.2)
	edge_l = list(edges)
	# indexes for edge_l of training sets
	# take the first third
	#train_indexes = random.sample(range(0, len(edge_l)), sample_size)
	fold1 = range(0, sample_size)
	fold2 = range(sample_size, sample_size*2) 
	fold3 = range(sample_size*2, sample_size*3) 
	fold4 = range(sample_size*3, sample_size*4) 
	fold5 = range(sample_size*4, len(edges))

	folds = {}
	folds[1] = {}
	folds[2] = {}
	folds[3] = {}
	folds[4] = {}
	folds[5] = {}

	for i in range(0, len(edge_l)):
		if i in fold1:
			folds[1][edge_l[i]] = 1.0
		elif i in fold2:
			folds[2][edge_l[i]] = 1.0
		elif i in fold3:
			folds[3][edge_l[i]] = 1.0
		elif i in fold4:
			folds[4][edge_l[i]] = 1.0
		elif i in fold5:
			folds[5][edge_l[i]] = 1.0

	# assigning priors to each
	print "selecting random negative training cases..."
	# store a list so we don't duplicate negative cases
	neg_training_cases = set()
	for i in range(0, len(edge_l)):
		neg_edge = randomNegEdgeSafe(g, neg_training_cases)
		neg_training_cases.add( neg_edge )
		if i in fold1:
			folds[1][neg_edge] = 0.01
		elif i in fold2:
			folds[2][neg_edge] = 0.01
		elif i in fold3:
			folds[3][neg_edge] = 0.01
		elif i in fold4:
			folds[4][neg_edge] = 0.01
		elif i in fold5:
			folds[5][neg_edge] = 0.01

	print "done!"
	return folds

def getSubnet(g, min_edges=500):

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

def writeGraph(folds, file):


	fh = open(file, 'w')	
	for fold in folds:
		for (A,B) in folds[fold]:
			type = None
			val = float(folds[fold][(A,B)])
			if val == 1:
				type = "SL"
			else:	
				type = "notSL"

			fh.write("\t".join([A, B, str(fold), type])+"\n")

	fh.close()

# subset the data to look at just these edges
consider_nodes = None


# get mappings, name to id
name2id, edges = parseNet(opts.slgraph)

g = nx.Graph()
g.add_edges_from(edges)

goBP = parseVals("goBP.all.tab")
print "Parsing data.."
goCC = parseVals("goCC.all.tab")
print "Parsing data.."
goMF = parseVals("goMF.all.tab")
print "Parsing data.."

# parse PPI
ppiN2ID, ppiEdges = parseNet("ppi.sgdid.tab")

# parse PPI Kernel
ppiKernel = parseVals("ppi.kernelEdges.tab")

# parse SGD Negative Kernel
#gnegKernel = parseVals("gNeg.kernelEdges.tab")

# for each fold, print out data

print "writing data..."

i = 1
for min in range(100, 1600, 100):

	h = getSubnet(g, min)

	edge_universe = h.edges()
	node_universe = h.nodes()
	folds = split(edge_universe, h)

	node_universe = set()
	edge_universe = set()
	for fold in folds:
		for edge in folds[fold]:
			edge_universe.add( edge )
			node_universe.add( edge[0] )
			node_universe.add( edge[1] )

	for combo in [ ('1','2','3'), ('5','4','3'), ('4','3','1') ]:
		output = opts.output+"/"+str(i)
		trainDir = output+"/train/"
		testDir = output+"/test/"
		if not os.path.exists(output):
			os.mkdir(output)
		if not os.path.exists(trainDir):
			os.mkdir(trainDir)
		if not os.path.exists(testDir):
			os.mkdir(testDir)

		writeGraph(folds, output+"/subgraph_folds.txt")
	
		trainIdx = combo
		learnIdx = None
		testIdx = None
		for j in range(1, 6):
			if str(j) in trainIdx:
				continue
			learnIdx = str(j)
			break
		for j in range(1, 6):
			if str(j) in learnIdx:
				continue
			testIdx = str(j)
			break
	
		# flatten data into train/learn/test splits
	
		slTrain, slLearn, slTest = splitData(folds, trainIdx, learnIdx, testIdx)
		out = trainDir
		fh = open(out+'gene.txt', 'w')
		for name in name2id:

			if name not in node_universe:
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
	
		fh = open(out+'slObserved.txt', 'w')
		printNet(fh, slTrain , name2id, edge_universe)
	
		fh = open(out+'sl.txt', 'w')
		printNet(fh, slLearn, name2id, edge_universe)
	
		# print just the nodes to consider for SL learning
		fh = open(out+'consider.txt', 'w')
		printEL(fh, slLearn , name2id, edge_universe)
	
		fh = open(out+'ppiKernel.txt', 'w')
		printGO(fh, ppiKernel, name2id, edge_universe)
	
		#fh = open(out+'negKernel.txt', 'w')
		#printGO(fh, gnegKernel , name2id, edge_universe)
	
		out = testDir
		fh = open(out+'gene.txt', 'w')
		for name in name2id:
			if name not in node_universe:
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
		
		#fh = open(out+'negKernel.txt', 'w')
		#printGO(fh, gnegKernel , name2id, edge_universe)
		
		fh = open(out+'heldOutSL.tab', 'w')
		printNet(fh, slTest, name2id, edge_universe)
		
		# infer just these values
		fh = open(out+'consider.txt', 'w')
		printEL(fh, slTest, name2id, edge_universe)
	
		# need this for populating values	
		fh = open(out+'sl.txt', 'w')
		printNet(fh, slTest, name2id, edge_universe)

		i += 1
	
