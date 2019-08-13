class Filtration:
	"""The function which allow the filtration of the data according to their MAF and variance of D_residuals"""


	def table_filter(table_to_filter, table_filtering, min_value):
		"""Delete the data from table_to_filter where the value in tale_filtering is inferior to the value given in input
		
		Parameters:
		----------------
		table_to_filter : DataFrame or Series
			Table that will be filtered
		table_filtering : Series
			Table in which is the value that will be used as criteria
		min_value : float
			Minimal value wanted for the criteria in table_filtering

		Returns:
		---------------
		filtered_table : DataFrame or Series
			Table after being filtered

		"""
		filtered_table = table_to_filter[table_filtering > min_value]
		return filtered_table
