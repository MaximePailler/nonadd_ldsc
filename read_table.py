import operator
from pandas_plink import read_plink1_bin
import pandas as pd
import numpy as np

class Read_table :
	'''Functions allowing the script to read the database and return the data in pd.DataFrame format'''

	def read_tped(path, min_maf, type_pos):
		'''
		Read the file that are in .tped format
		
		Parameter :
		---------------
		path : string
			path to .tped file without the extension
		min_maf : float
			minimal value choosen by the user
		type_pos : integer
			equal to 1 if we use the Morgans, equal to 2 if we use the base-pair coordinates

		Return :
		---------------
		gen_A : pandas.DataFrame
			DataFrame with all the genotypes in the additiv model
		pos : pandas.Series
			Series with the position of the SNP in Morgans
		chr : pandas.Series
			Series with the number of the chromosome for each SNP
		sample : pandas.Series
			Series of the Id of the individuals
		MAF : pandas.Series
			Series of the MAF for each SNP
		
		'''
		dict_gen = {}
		dict_chrom = {}
		dict_pos = {}
		dict_maf = {}

		with open(path + ".tped") as f:
			i = 1
			for row in f:
				print(i)
				i +=1
				row = row.split(" ")
				(chrom, snp) = (row[0], row[1])
				if type_pos == 1:
					pos = float(row[2])
				elif type_pos == 2:
					pos = float(row[3])
				row[-1] = row[-1].replace('\n', '')
				gen = row[4:]
				gen = [np.nan if i == 'N' else int(i) for i in gen]
				gen = list(map(operator.add, gen[::2], gen[1::2]))
				maf = min(np.nanmean(gen)/2, 1- np.nanmean(gen)/2)
				if maf > min_maf :
					dict_gen.update({snp: gen})
					dict_chrom.update({snp: chrom})
					dict_pos.update({snp: pos})
					dict_maf.update({snp: maf})

		with open(path + ".tfam") as f:
			sample = []
			for row in f:
				row = row.split(" ")
				sample.append(row[1])
	
		gen_A = pd.DataFrame(dict_gen, index = sample).T
		pos = pd.Series(dict_pos)
		chrom = pd.Series(dict_chrom)
		MAF = pd.Series(dict_maf)
		return(gen_A, pos, chrom, MAF, sample)



	def read_bed(path, type_pos):
		'''
		Read the .bed, .bim and .fam files
		
		Parameter :
		---------------
		path : string
			path to .bed, .bim. fam files without the extension. They need to be all in the same folder
		type_pos : integer
			equal to 1 if we use the Morgans, equal to 2 if we use the base-pair coordinates

		Return :
		---------------
		gen_A : pandas.DataFrame
			DataFrame with all the genotypes in the additiv model
		pos : pandas.Series
			Series with the position of the SNP in Morgans
		chr : pandas.Series
			Series with the number of the chromosome for each SNP
		sample : pandas.Series
			Series of the Id of the individuals
				
		'''
	
		data = read_plink1_bin(path + ".bed", path + ".bim", path + ".fam", verbose=False)
		if type_pos == 1:
			pos = pd.Series(np.transpose(data.variant.cm.values), index = data.variant.snp.values)
		elif type_pos == 2:
			pos = pd.Series(np.transpose(data.variant.pos.values), index = data.variant.snp.values)
		chrom = pd.Series(np.transpose(data.variant.chrom.values), index = data.variant.snp.values)
		sample = data.sample.iid
		gen_A = pd.DataFrame(data = np.transpose(data.values), index = data.variant.snp.values, columns = sample)
		return(gen_A, pos, chrom, sample)
