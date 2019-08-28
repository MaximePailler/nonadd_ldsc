import pandas as pd

df = pd.read_csv(filepath_or_buffer = "all_chr_output.csv", sep = " ")
print(df.shape)

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

print(df_SNPpos.shape)
print(df_rs.shape)

df_SNPpos.to_csv("./all_chr_pos.csv")
df_rs.to_csv("./all_chr_rs.csv")
