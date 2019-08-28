import pandas as pd

df_gwas_rs = pd.read_csv("all_chr_rs.csv")
df_gwas_pos = pd.read_csv("all_chr_pos.csv")
df_result = pd.read_csv("../../all_chr_result.csv")

df1 = df_gwas_rs.merge(df_result, how = 'inner', left_on = "SNP", right_on = "SNP")
 
df2 = df_gwas_pos.merge(df_result, how = 'inner', left_on = "SNP", right_on = "Id_pos")

df = pd.concat([df1, df2], ignore_index = True)

df = df.reindex(columns = ["Id_pos", "SNP", "Chr", "pos", "L_aa_orig", "L_ad_orig", "L_aa", "L_ad", "L_aa_new", "L_ad_new", "L_ad_new_adjust", "iscores", "NBETA-23104-0.0", "NSE-23104-0.0", "PV-23104-0.0", "MAF", "Nsnp", "Nsnp_used"])

df.to_csv("all_data.csv")

