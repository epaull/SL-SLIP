

for file in `find ../../data/yeast/SUBGRAPHS/ -name inferences.txt`; do

	num=`echo $file | sed -e 's/.*SUBGRAPHS\///'	-e 's/\/inferences.txt//'`
	echo $num
	cp $file inferences/$num.inferences.txt
done
