import pandas as pd
import numpy as np
import math
from sklearn.linear_model import LinearRegression
	
class Computation:
	'''Functions that will do all the computations during the script'''

	def calc_MAF(gen_A):
		'''
		Compute the MAF for all the variants in the database given in input

		Parameters :
		---------------
		gen_A : pandas.DataFrame
			Genotype in the additive model given the reading of .bed, .bim and .fam files
		
		Returns
		---------------
		MAF : pandas.Series
			Series of the MAF for all the SNP in the database

		'''
		MAF = pd.Series(gen_A.mean(axis=1)/2)
		MAF = MAF.apply(lambda x: min(x, 1-x))
		return MAF

	

	def calc_D_resid(gen_A, gen_D):
		'''
		Compute the values of the D_residuals for each variant in the additive model and the dominant model

		Parameters :
		---------------
		gen_A : pandas.DataFrame
			Genotype in the additive model 
		gen_D : pandas.DataFrame
			Genotype in the dominant model
		
		Returns :
		---------------
		resid : pandas.DataFrame
			DataFrame of all the D_residuals

		'''
		resid = pd.DataFrame(np.nan, index = gen_A.index.values, columns = range(gen_A.shape[1]))
		i=1
		for snp in gen_A.index.values:
			print(str(i) + " / "  + str(len(gen_A.index.values)))
			i+=1
			x = gen_A.loc[snp].values
			y = gen_D.loc[snp].values
			indexNA = [elem for elem in range(len(x)) if np.isnan(x[elem])]
			x_tmp = np.delete(x,indexNA).reshape(-1,1)
			y_tmp = np.delete(y,indexNA).reshape(-1,1)

			model = LinearRegression().fit(x_tmp,y_tmp)
			residual = y_tmp - model.predict(x_tmp)

			for elem in indexNA:
				residual = np.insert(residual, elem, np.nan)

			residual = pd.Series(residual.flatten())
			resid.loc[snp] = residual.values
		return resid



	def calc_LD_Score(gen_A, resid, pos, wind, MAF, chrom):
		'''
		Compute the LD_score with different methods (without adjustment, and 2 differents adjustments)
		
		Parameters :
		--------------
		gen_A : pandas.DataFrame
			Genotype in the additive model 
		resid : pandas.DataFrame
			Values of the D_residuals
		pos : pandas.Series
			position in Morgans for each variants
		wind : float
			Size of the window (in Morgans) that we will study around each SNP
		MAF : pandas.Series
			Value of th MAF for each variants
		chrom : pandas.Series
			Pair of chromosome corresponding to the SNP

		Returns :
		-------------
		result : pandas.DataFrame
			DataFrame gathering all the results of the study (SNP's Id, values for the LD_score, MAF, etc)
		'''
		#Initialization of the data
		LD_Score = pd.DataFrame(np.nan, index = gen_A.index.values, columns = ['l_aa_orig', 'l_ad_orig', 'l_aa', 'l_ad', 'l_aa_new', 'l_ad_new', 'l_ad_new_adjust'])
		N = len(gen_A.columns)
		dict_corr_originals = {}
		dict_corr = {}
		dict_corr_new = {}
		dict_corr_new_adjust = {}
		Nsnp = {}
		Nsnp_used = {}
		i=1
		threshold = 1/(N-1)
		for snp1 in gen_A.index.values:
			#For each variant, we keep only the variants inside the window and we initialize the LD_score to 0
			print(str(i) + " / "  + str(len(gen_A.index.values)))
			i+=1

			a = gen_A.loc[snp1]
			pos_snp = pos[snp1]

			A = gen_A[abs(pos-pos_snp) <= wind]
			D = resid[abs(pos-pos_snp) <= wind]

			(l_aa_originals, l_ad_originals) = (0, 0)
			(l_aa, l_ad) = (0, 0)
			(l_aa_new, l_ad_new) = (0, 0)
			l_ad_new_adjust = 0

			list_snp = A.index.values.tolist()
			Nsnp[snp1] = len(list_snp)
			Nsnp_used_value = len(list_snp)

			for snp2 in list_snp[:list_snp.index(snp1)]:
				#We are in the case where the correlation was already calculated. We just look for the value in the dictionaries and then delete it
				corr_A_originals = dict_corr_originals[(snp2, snp1)][0]
				corr_D_originals = dict_corr_originals[(snp2, snp1)][1]
				corr_A = dict_corr[(snp2, snp1)][0]
				corr_D = dict_corr[(snp2, snp1)][1]
				corr_A_new = dict_corr_new[(snp2, snp1)][0]
				corr_D_new = dict_corr_new[(snp2, snp1)][1]

				corr_D_new_adjust = dict_corr_new_adjust[(snp2, snp1)]
				if corr_D_new_adjust == 0:
					Nsnp_used_value -= 1

				#We delete elements of the dictionnaries we used here, we won't use them anymore
				del dict_corr_originals[(snp2, snp1)]
				del dict_corr[(snp2, snp1)]
				del dict_corr_new[(snp2, snp1)]
				del dict_corr_new_adjust[(snp2, snp1)]

				
				l_aa_originals += corr_A_originals
				l_ad_originals += corr_D_originals
				l_aa += corr_A
				l_ad += corr_D
				l_aa_new += corr_A_new
				l_ad_new += corr_D_new
				l_ad_new_adjust += corr_D_new_adjust





			for snp2 in list_snp[list_snp.index(snp1):]:
				#We are in the case where the correlation was not computed yet
				#So we calculate it and then make 2 types of adjustment
				(corr_A_originals, corr_D_originals) = (a.corr(A.loc[snp2], method = 'pearson')**2, a.corr(D.loc[snp2], method = 'pearson')**2)
				(corr_A, corr_D) = (corr_A_originals - (1 - corr_A_originals)/(N - 2), corr_D_originals - (1 - corr_D_originals)/(N - 2))
				(corr_A_new, corr_D_new) = (1 - ((N-3)*(1-corr_A_originals)/(N-2))*(1+2*(1-corr_A_originals)/(N-3.3)),1 - ((N-3)*(1-corr_D_originals)/(N-2))*(1+2*(1-corr_D_originals)/(N-3.3)))

				#If the value is missing, it's considered to be null
				if math.isnan(corr_A_originals): corr_A_originals = 0
				if math.isnan(corr_D_originals): corr_D_originals = 0
				if math.isnan(corr_A): corr_A = 0
				if math.isnan(corr_D): corr_D = 0
				if math.isnan(corr_A_new): corr_A_new = 0
				if math.isnan(corr_D_new): corr_D_new = 0

				if abs(corr_D_originals) <= threshold:
					corr_D_new_adjust = 0
					Nsnp_used_value -= 1
				else:
					corr_D_new_adjust = corr_D_new

				#We input the values calculated into dictionaries. We won't have to compute them anymore
				dict_corr_originals[(snp1, snp2)] = [corr_A_originals, corr_D_originals]
				dict_corr[(snp1, snp2)] = [corr_A, corr_D]
				dict_corr_new[(snp1, snp2)] = [corr_A_new, corr_D_new]
				dict_corr_new_adjust[(snp1, snp2)] = corr_D_new_adjust

				l_aa_originals += corr_A_originals
				l_ad_originals += corr_D_originals
				l_aa += corr_A
				l_ad += corr_D
				l_aa_new += corr_A_new
				l_ad_new += corr_D_new
				l_ad_new_adjust += corr_D_new_adjust

			LD_Score.loc[snp1] = [l_aa_originals, l_ad_originals, l_aa, l_ad, l_aa_new, l_ad_new, l_ad_new_adjust]
			Nsnp_used[snp1] = Nsnp_used_value
		
		#We regroup all the results into one big DataFrame
		Nsnp_serie = pd.Series(Nsnp)
		Nsnp_used_serie = pd.Series(Nsnp_used)
		result = pd.concat([chrom, LD_Score, Nsnp_serie, Nsnp_used_serie, MAF], axis = 1)
		return result
