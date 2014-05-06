#!/usr/bin/env	python

import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--map",dest="genes",action="store",type="string",default=None)
(opts, args) = parser.parse_args()

def parseMap(file, network_nodes=None):
	
	map = {}
	fh = None
	try:
		fh = open(file, 'r')
	except:
		raise Exception("Error: can't open file: "+file)

	for line in fh:
		parts = line.rstrip().split("\t")
		# maps ids to gene names
		map[parts[0]] = parts[1]

	return map

map = parseMap(opts.genes)

for line in sys.stdin:
	connection, value = line.rstrip().split("\t")		
	geneA, geneB = connection.lstrip("INFLUENCES(").rstrip(")").split(", ")
	#print "\t".join([map[geneA], map[geneB], value])
	print "\t".join([geneA, geneB, value])
				
