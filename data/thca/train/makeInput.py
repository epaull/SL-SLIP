#!/usr/bin/env	python

"""

 Need to map all values to somewhere in 0-1 range for PSL. A good number need to be in the .66+ range here, 
 it's a lot like paradigm
"""
import networkx as nx
import math

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--goBP",dest="goBP",action="store",type="string",default=None)
parser.add_option("--goCC",dest="goCC",action="store",type="string",default=None)
parser.add_option("--allpairs",dest="allpairs",action="store",type="string",default=None)
parser.add_option("--network",dest="network",action="store",type="string",default=None)
(opts, args) = parser.parse_args()


PSEUDO_COUNT = 0.01

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


def parseCorr(file):

	prot2prot = set()
	expr2expr = set()
	prot2expr = set()
	expr2prot = set()
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	first = True
	for line in fh:
		if first:
			first = False
			continue
		geneA, geneB, brafR, brafP, rasR, rasP, diff, sign, signaling_type = line.rstrip().split("\t")

		if signaling_type == "Prot -> Prot":
			prot2prot.add( (geneA, geneB, brafR, rasR) )
		elif signaling_type == "Prot -> Expr":
			prot2expr.add( (geneA, geneB, brafR, rasR) )
		elif signaling_type == "Expr -> Expr":
			expr2expr.add( (geneA, geneB, brafR, rasR) )
		elif signaling_type == "Expr -> Prot":
			expr2prot.add( (geneA, geneB, brafR, rasR) )

	return (prot2prot, expr2expr, expr2prot, prot2expr)

def parseNet(file):

	genes = set()
	edges = set()
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	for line in fh:
		geneA, link, geneB = line.rstrip().split("\t")
		genes.add( geneA )
		genes.add( geneB )
		edges.add( (geneA, geneB) )


	map = {}
	i = 1
	for gene in genes:
		map[gene] = str(i)
		i += 1
	
	return (map, edges)

def printLengths(fh, vals, map):

	# linearly scale from 0-1
	max = 0
	for (geneA, geneB) in vals:
		val = vals[(geneA, geneB)]
		if abs(float(val)) > max:
			max = abs(float(val)) + PSEUDO_COUNT

	for (geneA, geneB) in vals:
		val = vals[(geneA, geneB)]

		# add pseudo-counts
		val = float(val)
		val += PSEUDO_COUNT

		fh.write("\t".join( [map[geneA], map[geneB], str(abs(float(val))/max)] )+"\n")
		fh.write("\t".join( [map[geneB], map[geneA], str(abs(float(val))/max)] )+"\n")

	fh.close()	
def printGO(fh, vals, map, symmetric=True):

	# linearly scale from 0-1
	#max = 0
	#for (geneA, geneB, val) in vals:
	#	if abs(float(val)) > max:
	#		max = abs(float(val)) + PSEUDO_COUNT

	for (geneA, geneB, val) in vals:

		# add pseudo-counts
		val = float(val)

		# ERASE
		max = 1
		if val > 0.5:
			val = 1.0
		if val > 0.3:
			val = 0.9
		if val > 0.2:
			val = 0.8
		if val > 0.1:
			val = 0.75
		if val > 0:
			val = 0.6

		val += PSEUDO_COUNT

		fh.write("\t".join( [map[geneA], map[geneB], str(abs(float(val))/max)] )+"\n")
		# print the symmetric case
		if symmetric:
			fh.write("\t".join( [map[geneB], map[geneA], str(abs(float(val))/max)] )+"\n")

	fh.close()	

def printCorr(fh, vals, map, idx):
	for (geneA, geneB, val1, val2) in vals:

		if val1 > 0.5:
			val1 = 0.95
		elif val1 > 0.25:
			val1 = 0.9
		elif val1 > 0.1:
			val1 = 0.75

		if val2 > 0.5:
			val2 = 0.95
		elif val2 > 0.25:
			val2 = 0.9
		elif val2 > 0.1:
			val2 = 0.75


		if idx == 1:
			fh.write("\t".join( [map[geneA], map[geneB], str(abs(float(val1)))] )+"\n")
			fh.write("\t".join( [map[geneB], map[geneA], str(abs(float(val1)))] )+"\n")
		else:
			fh.write("\t".join( [map[geneA], map[geneB], str(abs(float(val2)))] )+"\n")
			fh.write("\t".join( [map[geneB], map[geneA], str(abs(float(val2)))] )+"\n")
	fh.close()	

	
# get mappings, name to id
name2id, edges = parseNet(opts.network)

# make a digraph
graph = nx.Graph()
graph.add_edges_from( edges )
length = nx.all_pairs_shortest_path_length(graph)
spls = {}
for geneA in name2id:
	for geneB in name2id:
		if geneA == geneB:
			continue
		# don't duplicate
		if (geneB, geneA) in spls:
			continue
		if geneA not in length or geneB not in length[geneA]:
			spls[(geneA, geneB)] = "0"
		else:
			l = length[geneA][geneB]	
			# convert, 0-1 range, biased towards the top
			l = 1/math.sqrt(float(l))
			spls[(geneA, geneB)] = l


goBP = parseVals(opts.goBP)
goCC = parseVals(opts.goCC)
prot2prot, expr2expr, expr2prot, prot2expr = parseCorr(opts.allpairs)

out = 'braf/'
fh = open(out+'gene.txt', 'w')
for name in name2id:
	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'influences.txt', 'w')
printLengths(fh, spls, name2id)

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id)
fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id)

fh = open(out+'prot2protCOR.txt', 'w')
printCorr(fh, prot2prot, name2id, 1)
fh = open(out+'prot2exprCOR.txt', 'w')
printCorr(fh, prot2expr, name2id, 1)
fh = open(out+'expr2exprCOR.txt', 'w')
printCorr(fh, expr2expr, name2id, 1)
fh = open(out+'expr2protCOR.txt', 'w')
printCorr(fh, expr2prot, name2id, 1)
		
out = 'ras/'
fh = open(out+'gene.txt', 'w')
for name in name2id:
	fh.write(name2id[name]+"\t"+name+"\n")
fh.close()

fh = open(out+'influences.txt', 'w')
printLengths(fh, spls, name2id)

fh = open(out+'goBP.txt', 'w')
printGO(fh, goBP, name2id)
fh = open(out+'goCC.txt', 'w')
printGO(fh, goCC, name2id)

fh = open(out+'prot2protCOR.txt', 'w')
printCorr(fh, prot2prot, name2id, 2)
fh = open(out+'prot2exprCOR.txt', 'w')
printCorr(fh, prot2expr, name2id, 2)
fh = open(out+'expr2exprCOR.txt', 'w')
printCorr(fh, expr2expr, name2id, 2)
fh = open(out+'expr2protCOR.txt', 'w')
printCorr(fh, expr2prot, name2id, 2)
		
