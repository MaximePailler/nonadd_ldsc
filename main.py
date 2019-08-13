# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 14:17:19 2019

@author: utilisateur
"""

from read_table import Read_table
from computation import Computation
from filtration import Filtration

#-------------------Initialisation of the data---------------------------

while True:
	try:
		type_file = int(input("Which type of file do you want to read?:\n [1].tped file \n [2].bed, .bim, .fam files\n"))
	except ValueError:
		print("Your answer is not valid. Please, input 1 or 2 to continue")
		continue	

	if not(type_file in [1,2]):
		print("Your answer is not valid. Please, input 1 or 2 to continue")
		continue
	else:
		break

path = str(input("Write here the path to your file (without the extension):"))

while True:
	try:
		min_maf = float(input("Enter the minimum value of MAF you want :"))
	except ValueError:
		print("Your answer is not valid. Please, enter a numeric value")
		continue	

	if min_maf <= 0:
		print("You can't use a negative or null value. Please, enter another value")
		continue
	else:
		break

while True:
	try:
		wind = float(input("Enter the size of the window (in Morgans) :"))
	except ValueError:
		print("Your answer is not valid. Please, enter a numeric value")
		continue	

	if wind <= 0:
		print("You can't use a negative or null value. Please, enter another value")
		continue
	else:
		break

#-----------------------------------------------------------------------------------------------------------------------
#                                                     Using .tped file
#-----------------------------------------------------------------------------------------------------------------------
if type_file == 1:
	(gen_A, pos, chrom, MAF, sample) = Read_table.read_tped(path, min_maf)
	print("data, MAF and filterMAF ok")

#----------------------------------------------------------------------------------------------------------------------
#                                                 Using .bed .bim .fam files
#----------------------------------------------------------------------------------------------------------------------
elif type_file == 2:

	(gen_A, pos, chrom, sample) = Read_table.read_bed(path)
	print('data ok')

	MAF = Computation.calc_MAF(gen_A)
	print('MAFs ok')
	
#--------------------Filtration by MAF------------------------------

	gen_A = Filtration.table_filter(gen_A, MAF, min_maf)
	pos = Filtration.table_filter(pos, MAF, min_maf)
	chrom = Filtration.table_filter(chrom, MAF, min_maf)
	MAF = Filtration.table_filter(MAF, MAF, min_maf)
	print('Filter MAF ok')



#----------------------------------------------------------------------------------------------------------------------
#                                                For any type of data 
#----------------------------------------------------------------------------------------------------------------------


#--------------------Genotype with dominant model------------------------

gen_D = gen_A.replace(2, 1)
print('gen_D ok')

#---------------------Computation of D_residuals-----------------------

resid = Computation.calc_D_resid(gen_A, gen_D, sample)
print('D_resid ok')

#-------------------Computation of the variance--------------------

var_D_resid = resid.var(axis = 1, skipna = True)
print('Variance ok')

#--------------------Filtration by variance--------------------

gen_A = Filtration.table_filter(gen_A, var_D_resid, 10**(-3))
gen_D = Filtration.table_filter(gen_D, var_D_resid, 10**(-3))
pos = Filtration.table_filter(pos, var_D_resid, 10**(-3))
chrom = Filtration.table_filter(chrom, var_D_resid, 10**(-3))
MAF = Filtration.table_filter(MAF, var_D_resid, 10**(-3))
resid = Filtration.table_filter(resid, var_D_resid, 10**(-3))
var_D_resid = Filtration.table_filter(var_D_resid, var_D_resid, 10**(-3))
print('Filter variance ok')

#-------------------Computation of the LD Score------------------

result = Computation.calc_LD_Score(gen_A, resid, pos, wind, MAF)
print('Results obtained')

#Output the results to a .csv file
result.to_csv(path_or_buf = path + "_output.csv", header = ["L_aa_orig", "L_ad_orig", "L_aa", "L_ad", "L_aa_new", "L_ad_new", "L_ad_new_adjust", "Nsnp","Nsnp_used", "MAF"], index = True, sep = " ", na_rep = "Nan")



