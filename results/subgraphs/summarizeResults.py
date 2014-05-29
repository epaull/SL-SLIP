#!/usr/bin/env	python

import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--inf",dest="inf",action="store",type="string",default=None)
parser.add_option("--subnet",dest="subnet",action="store_true",default=False)
(opts, args) = parser.parse_args()

baseline_aucs = {}
aucs = {}
train_set_size = {}

def getAUC(file):

	fh = open(file)
	auc = None
	for line in fh:
		auc = line.rstrip().split(" ")[1]

	fh.close()			
	return auc

def parseINF(file):

	fh = open(file)
	auc = None

	#{0.0} ( ( ( CONSIDER(A, B) & SL(A, X) ) & SL(X, B) ) & #NOTEQUAL(A, B) ) >> ~( SL(A, B) ) {squared}
	#{1.0220033144408238} ( ( ( CONSIDER(A, B) & SL(A, X) ) & PPIEDGES(X, B) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{4.694577841312776} ( ( ( CONSIDER(A, B) & SL(A, X) ) & PPIKERNEL(X, B) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{0.0} ( ( CONSIDER(A, B) & GOBP(A, B) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{0.0} ( ( CONSIDER(A, B) & GOCC(A, B) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{0.0} ( ( CONSIDER(A, B) & GOMF(A, B) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{0.0} ( ( ( CONSIDER(A, B) & GOBP(A, B) ) & ~( GOMF(A, B) ) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{1.0} ( ( ( CONSIDER(A, B) & GOCC(A, B) ) & ~( GOMF(A, B) ) ) & #NOTEQUAL(A, B) ) >> SL(A, B) {squared}
	#{7.615980462131928} ( CONSIDER(A, B) & PPIKERNEL(A, B) ) >> ~( SL(A, B) ) {squared}
	#{5.238885973414747} ( CONSIDER(A, B) & PPIEDGES(A, B) ) >> ~( SL(A, B) ) {squared}

	# parse rules: starts with '{', but doesn't include 'consider' in that entry. 
	for line in fh:
		line = line.rstrip()

	fh.close()			

def countNET(file):

	fh = open(file)
	lineno = 0
	for line in fh:
		lineno += 1

	fh.close()			

	return lineno


for i in range(1, 51):
	auc_f = "aucs/"+str(i)+".auc.txt"	
	base_f = "aucs/"+str(i)+".baseline-auc.txt"	
	graph = "graphs/"+str(i)+".subgraph_folds.txt"	

	train_set_size[i] = countNET(graph)
	baseline_aucs[i] = getAUC(base_f)	
	aucs[i] = getAUC(auc_f)	


for i in range(1,51):
	try:
		print "\t".join([str(i), aucs[i], baseline_aucs[i], str(train_set_size[i])])
	except:
		continue
