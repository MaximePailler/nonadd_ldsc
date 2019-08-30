import pandas as pd

#Read the .csv file which contains all the GWAS results
df = pd.read_csv(filepath_or_buffer = "all_chr_output.csv", sep = " ")

#The file is splited into two dataframe file according to the way the SNP Id is written
boolean = [False]*df.shape[0]
i = 0
for elem in df["SNP"].values:
	if ":" in elem:
		boolean[i] = True
	i += 1

boolean = pd.Series(boolean)
df_SNPpos = df[boolean]
df_SNPpos["SNP"] = df_SNPpos["SNP"].str.slice(stop = -4)

df_rs = df[~boolean]

#Exportation of the two dataframe into two different .csv files
df_SNPpos.to_csv("./all_chr_pos.csv")
df_rs.to_csv("./all_chr_rs.csv")
