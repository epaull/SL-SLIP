
fold=$1

mkdir -p SUBGRAPHS/train
mkdir -p SUBGRAPHS/test

cp SUBGRAPHS/$fold/train/goBP.txt SUBGRAPHS/train/GOBP.txt
cp SUBGRAPHS/$fold/train/goCC.txt SUBGRAPHS/train/GOCC.txt
cp SUBGRAPHS/$fold/train/goMF.txt SUBGRAPHS/train/GOMF.txt
cp SUBGRAPHS/$fold/train/consider.txt SUBGRAPHS/train/CONSIDER.txt
cp SUBGRAPHS/$fold/train/gene.txt SUBGRAPHS/train/GENE.txt
cp SUBGRAPHS/$fold/train/ppiEdges.txt SUBGRAPHS/train/PPIEDGES.txt
cp SUBGRAPHS/$fold/train/sl.txt SUBGRAPHS/train/SL.txt
cp SUBGRAPHS/$fold/train/slObserved.txt SUBGRAPHS/train/SLOBSERVED.txt
cp SUBGRAPHS/$fold/train/ppiKernel.txt SUBGRAPHS/train/PPIKERNEL.txt


cp SUBGRAPHS/$fold/test/goBP.txt SUBGRAPHS/test/GOBP.txt
cp SUBGRAPHS/$fold/test/goCC.txt SUBGRAPHS/test/GOCC.txt
cp SUBGRAPHS/$fold/test/goMF.txt SUBGRAPHS/test/GOMF.txt
cp SUBGRAPHS/$fold/test/consider.txt SUBGRAPHS/test/CONSIDER.txt
cp SUBGRAPHS/$fold/test/gene.txt SUBGRAPHS/test/GENE.txt
cp SUBGRAPHS/$fold/test/ppiEdges.txt SUBGRAPHS/test/PPIEDGES.txt
cp SUBGRAPHS/$fold/test/sl.txt SUBGRAPHS/test/SL.txt
cp SUBGRAPHS/$fold/test/slObserved.txt SUBGRAPHS/test/SLOBSERVED.txt
cp SUBGRAPHS/$fold/test/ppiKernel.txt SUBGRAPHS/test/PPIKERNEL.txt
cp SUBGRAPHS/$fold/test/heldOutSL.tab SUBGRAPHS/test/heldOutSL.tab
