

data:
	~/bin/python2.7 ./makeTestInput.py --slgraph fold_labels.tab --output FULL_5CV

auc:
	~/bin/python2.7 ./getAUC.py --geneMap TEST/test/GENE.txt --inferences TEST/test/out.txt --truth TEST/test/heldOutSL.tab
