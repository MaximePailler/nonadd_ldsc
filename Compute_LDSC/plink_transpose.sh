

for chr in {1..22}
do

                /usr/local/bin/plink1.9 \
		--bfile /home/common/g1000p3_EUR/genetic_map/with_cms_chr"$chr" \
		--recode 01 \
		--transpose \
		--maf 0.05 \
		--output-missing-genotype N \
		--make-bed \
                --out /home/common/projects/nonadd_ldsc/pailler/data/chr_"$chr"/chr"$chr"_data \
		--threads 5

done

