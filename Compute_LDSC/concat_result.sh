#!/bin/bash


cat ./chr_1/chr1_data_output_pos.csv > all_chr_result.csv

for i in {2..22}
do
	tail -n+2 ./chr_${i}/chr${i}_data_output_pos.csv >> all_chr_result.csv
done
