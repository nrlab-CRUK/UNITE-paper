

# Define motif_filter_nrow function
def motif_filter_nrow(df, colname_to_use, output_col):
	# only keep cols: id, frag_len, colname_to_use
	df = df[['id', 'frag_len', colname_to_use]]
	# change dtype of id frag_len and colname_to_use to str, use .loc
	df[colname_to_use] = df[colname_to_use].astype('str')
	prefix = 'n_' + colname_to_use
	df_dummy = pd.get_dummies(df, columns=[colname_to_use], dummy_na=False, dtype='int8', prefix=prefix)
	result = (df_dummy
		.groupby(['id', 'frag_len'], observed = False)
		.agg({output_col: 'sum'})
		.reset_index()
	)
	# sort the result by id and frag_len
	return result

def each_output_layer(output_col, input_df, ref, **kwargs):
	
	# Filter reference DataFrame based on output_col
	colname_to_use = ref.loc[ref['output_col'] == output_col, 'source_col'].values[0]
	method_to_use = ref.loc[ref['output_col'] == output_col, 'method'].values[0]
	# Calculate metrics based on method
	# sum or mean
	if method_to_use == 'sum' or method_to_use == 'mean':
		#result = input_df.groupby(['id', 'frag_len']).agg({colname_to_use: 'sum'}).reset_index()
		result = (input_df.groupby(['id', 'frag_len'], observed=False)
			.agg({colname_to_use: method_to_use})
			# rename the column name to output_col
			.rename(columns={colname_to_use: output_col})
			.reset_index()
		)
	elif method_to_use == 'motif_filter_nrow':
		result = motif_filter_nrow(input_df, colname_to_use=colname_to_use, output_col=output_col)
	else:
		raise ValueError("Undefined calculation method, please check function `each_output_layer` and `ref` DataFrame.")
	return result

