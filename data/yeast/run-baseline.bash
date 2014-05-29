#!/bin/bash

for fold in 13; do

	sh chooseFold.SUBGRAPHS.sh $fold
	pushd ../../develop/sl-slip
	make k-baseline > ../../data/yeast/SUBGRAPHS/test/baseline-inferences.txt 
	popd
	cp SUBGRAPHS/test/baseline-inferences.txt SUBGRAPHS/$fold/
	~/bin/python2.7 ./getAUC.py --geneMap SUBGRAPHS/test/GENE.txt --inferences SUBGRAPHS/test/baseline-inferences.txt --truth SUBGRAPHS/test/heldOutSL.tab > SUBGRAPHS/$fold/baseline-auc.txt

done

