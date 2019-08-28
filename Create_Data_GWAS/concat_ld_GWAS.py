import pandas as pd

df_gwas = pd.read_csv("./all_chr_output.csv")
df_ld = pd.read_csv("../../all_chr_result.csv")

df_ld = df_ld.reindex(columns = ['Id_pos', 'SNP', 'L_aa_orig', 'L_ad_orig', 'L_aa', 'L_ad', 'L_aa_new', 'L_ad_new', 'L_ad_new_adjust', 'Nsnp', 'Nsnp_used', 'MAF'])
 
df_concat = df_gwas.merge(df_ld, how = 'left', left_on = "SNP", right_on ="Id_pos")
df_concat = df_concat.merge(df_ld, how = 'left', left_on = "SNP", right_on ="SNP")

df_concat.to_csv("./all_data.csv")
