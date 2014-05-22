
fold=$1

mkdir -p FULL_5CV/train
mkdir -p FULL_5CV/test

cp FULL_5CV/$fold/train/goBP.txt FULL_5CV/train/GOBP.txt
cp FULL_5CV/$fold/train/goCC.txt FULL_5CV/train/GOCC.txt
cp FULL_5CV/$fold/train/goMF.txt FULL_5CV/train/GOMF.txt
cp FULL_5CV/$fold/train/consider.txt FULL_5CV/train/CONSIDER.txt
cp FULL_5CV/$fold/train/gene.txt FULL_5CV/train/GENE.txt
cp FULL_5CV/$fold/train/ppiEdges.txt FULL_5CV/train/PPIEDGES.txt
cp FULL_5CV/$fold/train/sl.txt FULL_5CV/train/SL.txt
cp FULL_5CV/$fold/train/slObserved.txt FULL_5CV/train/SLOBSERVED.txt
cp FULL_5CV/$fold/train/negKernel.txt FULL_5CV/train/NEGKERNEL.txt
cp FULL_5CV/$fold/train/ppiKernel.txt FULL_5CV/train/PPIKERNEL.txt


cp FULL_5CV/$fold/test/goBP.txt FULL_5CV/test/GOBP.txt
cp FULL_5CV/$fold/test/goCC.txt FULL_5CV/test/GOCC.txt
cp FULL_5CV/$fold/test/goMF.txt FULL_5CV/test/GOMF.txt
cp FULL_5CV/$fold/test/consider.txt FULL_5CV/test/CONSIDER.txt
cp FULL_5CV/$fold/test/gene.txt FULL_5CV/test/GENE.txt
cp FULL_5CV/$fold/test/ppiEdges.txt FULL_5CV/test/PPIEDGES.txt
cp FULL_5CV/$fold/test/sl.txt FULL_5CV/test/SL.txt
cp FULL_5CV/$fold/test/slObserved.txt FULL_5CV/test/SLOBSERVED.txt
cp FULL_5CV/$fold/test/negKernel.txt FULL_5CV/test/NEGKERNEL.txt
cp FULL_5CV/$fold/test/ppiKernel.txt FULL_5CV/test/PPIKERNEL.txt
cp FULL_5CV/$fold/test/heldOutSL.tab FULL_5CV/test/heldOutSL.tab
