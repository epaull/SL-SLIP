#!/bin/bash

for fold in {50..51}; do

	sh chooseFold.sh $fold
	pushd ../../develop/sl-slip
	make k-run > ../../data/yeast/FULL_5CV/test/inferences.txt 
	popd
	cp FULL_5CV/test/inferences.txt FULL_5CV/$fold/
	~/bin/python2.7 ./getAUC.py --geneMap FULL_5CV/test/GENE.txt --inferences FULL_5CV/test/inferences.txt --truth FULL_5CV/test/heldOutSL.tab	

done

