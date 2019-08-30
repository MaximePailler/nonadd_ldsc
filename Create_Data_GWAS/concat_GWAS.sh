#!/bin/bash

cat imputed.allWhites.23104-0.0.chr1.csv > all_chr_output.csv

for i in {2..22}
do
	tail -n+2 imputed.allWhites.23104-0.0.chr${i}.csv >> all_chr_output.csv
done
