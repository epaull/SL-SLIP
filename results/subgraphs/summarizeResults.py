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
inferences = {}

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
	state_p = False
	p_line = None

	weights = {}
	for line in fh:
		line = line.rstrip()

		val = None
		if state_p:
			if p_line == 10:
				# friend of friend rule
				val = line.split(" ")[0].lstrip("{").rstrip("}")
				weights[p_line] = val
			elif p_line == 11:
				# 2 hop PPI/SL
				val = line.split(" ")[0].lstrip("{").rstrip("}")
				weights[p_line] = val

			elif p_line == 12:
				# 2 hop PPI Kernel/SL
				val = line.split(" ")[0].lstrip("{").rstrip("}")
				weights[p_line] = val

			elif p_line == 18:
				# PPI Kernel predictive
				val = line.split(" ")[0].lstrip("{").rstrip("}")
				weights[p_line] = val

			elif p_line == 19:
				# PPI Edges predictive
				val = line.split(" ")[0].lstrip("{").rstrip("}")
				weights[p_line] = val
			elif p_line == 20:
				fh.close()
				return weights

			p_line += 1


		elif line == "\t\tLEARNING WEIGHTS DONE":
			state_p = True
			p_line = 1


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
	inf = "inferences/"+str(i)+".inferences.txt"	

	train_set_size[i] = countNET(graph)
	baseline_aucs[i] = getAUC(base_f)	
	aucs[i] = getAUC(auc_f)	
	inferences[i] = parseINF(inf)

print "\t".join(["Index", "AUC", "Baseline AUC", "Number Edges", "Friend of Friend", "2-Hop PPI/SL", "2-Hop PPI Kernel/SL", "PPI Kernel Predictive", "PPI Edges Predictive"])
for i in range(1,51):
	try:
		printstr = "\t".join([str(i), aucs[i], baseline_aucs[i], str(train_set_size[i])])
		for j in [10,11,12,18,19]:
			printstr += "\t"+inferences[i][j]
		print printstr
	except:
		continue

