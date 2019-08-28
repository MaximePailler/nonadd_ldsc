import pandas as pd

for i in range(1, 23):
	path_data = "./chr_" + str(i) + "/chr" + str(i) + "_data_output.csv"
	path_pos = "./chr_" + str(i) + "/chr" + str(i) + "_data.map"
	df_data = pd.read_csv(path_data, sep = " ")
	df_pos = pd.read_csv(path_pos, sep = "\t", header = None, names = ["chr", "SNP", "Morgan", 'pos'])
	df_pos = df_pos.iloc[:, [1,3]]
	result = df_data.merge(df_pos, on = 'SNP')
	result["Id_pos"] = result['Chr'].map(str) + ":" + result['pos'].map(str)
	print(result)
	result.to_csv(path_or_buf = "./chr_" + str(i) + "/chr" + str(i) + "_data_output_pos.csv", index = False)
	
