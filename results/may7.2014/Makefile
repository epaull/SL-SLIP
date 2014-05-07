

corr.in:
	../mapBack.py --map gene.txt < out.txt > gene_inferences.txt
	join.pl -1 1,2 -2 1,2 gene_inferences.txt influences.txt | cut -f 3,4 > $@
	../mapBack.py --map gene.txt --names < out.txt > gene_inferences.txt

