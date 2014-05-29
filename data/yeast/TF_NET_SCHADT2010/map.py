#!/usr/bin/env  python2.7

import sys
import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument('map_name', type=str, default=None, help="<SGDID>   <Name>  <Yeast Gene>")
args = parser.parse_args()

def parseMap(file):

	name2id = {}
	yeast2id = {}
	for line in open(file, 'r'):
		sgdid, yeast_name = line.rstrip().split("\t")
		yeast2id[yeast_name] = sgdid
	#	   try:
	#			   sgdid, yeast_name, name = line.rstrip().split("\t")	 
	#			   name2id[name] = sgdid
	#			   yeast2id[yeast_name] = sgdid
	#	   except:
	#			   continue
#
	return yeast2id

map_yeastID = parseMap(args.map_name)
for line in sys.stdin:
	geneA = None
	geneB = None
	try:
		geneA, geneB, score = line.rstrip().split("\t")
	except:
		continue
	if geneA in map_yeastID and geneB in map_yeastID:
		print map_yeastID[geneA]+"\t"+map_yeastID[geneB]+"\t"+str(score)
