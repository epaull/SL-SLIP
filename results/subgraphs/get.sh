

for file in `find ../../data/yeast/SUBGRAPHS/ -name baseline-auc.txt`; do

	num=`echo $file | sed -e 's/.*SUBGRAPHS\///'	-e 's/\/baseline-auc.txt//'`
	echo $num
	cp $file $num.baseline-auc.txt
done
