
fold=$1

mkdir -p TEST/train
mkdir -p TEST/test

cp TEST/$fold/train/goBP.txt TEST/train/GOBP.txt
cp TEST/$fold/train/goCC.txt TEST/train/GOCC.txt
cp TEST/$fold/train/goMF.txt TEST/train/GOMF.txt
cp TEST/$fold/train/consider.txt TEST/train/CONSIDER.txt
cp TEST/$fold/train/gene.txt TEST/train/GENE.txt
cp TEST/$fold/train/ppiEdges.txt TEST/train/PPIEDGES.txt
cp TEST/$fold/train/sl.txt TEST/train/SL.txt
cp TEST/$fold/train/slObserved.txt TEST/train/SLOBSERVED.txt
cp TEST/$fold/train/negKernel.txt TEST/train/NEGKERNEL.txt
cp TEST/$fold/train/ppiKernel.txt TEST/train/PPIKERNEL.txt


cp TEST/$fold/test/goBP.txt TEST/test/GOBP.txt
cp TEST/$fold/test/goCC.txt TEST/test/GOCC.txt
cp TEST/$fold/test/goMF.txt TEST/test/GOMF.txt
cp TEST/$fold/test/consider.txt TEST/test/CONSIDER.txt
cp TEST/$fold/test/gene.txt TEST/test/GENE.txt
cp TEST/$fold/test/ppiEdges.txt TEST/test/PPIEDGES.txt
cp TEST/$fold/test/sl.txt TEST/test/SL.txt
cp TEST/$fold/test/slObserved.txt TEST/test/SLOBSERVED.txt
cp TEST/$fold/test/negKernel.txt TEST/test/NEGKERNEL.txt
cp TEST/$fold/test/ppiKernel.txt TEST/test/PPIKERNEL.txt
cp TEST/$fold/test/heldOutSL.tab TEST/test/heldOutSL.tab
