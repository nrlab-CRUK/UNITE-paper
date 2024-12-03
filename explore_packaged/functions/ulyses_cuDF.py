
import pandas as pd
import os
import sys
import time
import re
import pyreadr
from functools import reduce
import pickle


################################################################################
# parameters
################################################################################
overlap_file = '/scratchc/nrlab/wang04/ulyses_second_batch/try_ulyses_gpu/SLX-18445.SXTLI192.HiSeq4000.bam.0.1x.mrkdup.bam.ulyses_bind_olap_chunks.rds.csv'
isize_from = 20
isize_to = 500
ref_csv_file = "/home/nrlab/wang04/ulyses/explore_packaged/resources/overlap_ref.csv"
layers = ["n_isize",
          "n_motif_s1_C",
          "n_motif_s1_A",
          "n_motif_s1_G",
          "n_motif_s1_T",
          "n_neg_nm", 
          "n_pos_nm", 
          "bin_mean_GC", 
          "bin_mean_mappability"]

# make output file name, use the overlap_file name with .csv replaced by .pkl
output_file = re.sub(r'\.csv$', '.pkl', overlap_file)

################################################################################
# functions
################################################################################


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



################################################################################
# read in files 
################################################################################
# say "start to read in rds files"
print("start to read in csv files")
# read ref file
ref = pd.read_csv(ref_csv_file)
# read in csv file using pyarrow
df = pd.read_csv(overlap_file, engine='pyarrow')

# say "finished reading in rds files"
print("finished reading in csv files")

################################################################################
# main
################################################################################
# report the time start to run the script
start_time = time.time()
#print start time
print("Start time: ", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)))
# get unique values of bin_size column
bin_size_levels = df['bin_size'].unique().tolist()

# make a dict called final_result, each value is the output of each_output_layer function, key is the bin_size
final_result = {}
for bin_size in bin_size_levels:
    print("handling bin_size: ", bin_size)
	# filter df based on bin_size
    df1 = df.query('bin_size == @bin_size')
    # get unique values of id column
    id_levels = df1['id'].unique().tolist()
    # reset the category levels of id column to id_levels, ordered=True
    df1['id'] = df1['id'].astype('category')
    df1['id'] = df1['id'].cat.set_categories(id_levels, ordered=True)
    # set the frag_len column as category type, levels to range(isize_min, isize_max), ordered=True
    df1['frag_len'] = df1['frag_len'].astype('category')
    df1['frag_len'] = df1['frag_len'].cat.set_categories(range(isize_from, isize_to+1), ordered=True)
    # make a dict called result, each value is the output of each_output_layer function, key is the output_col
    result = {}
    for output_col in layers:
        result[output_col] = each_output_layer(output_col, df1, ref)
    # merge all the values in result dict to one dataframe, merge on id and frag_len, how='outer'
    # Merge all DataFrames in the dictionary using an outer join
    merged_df = reduce(lambda x, y: pd.merge(x, y, on=['id', 'frag_len'], how='outer'), result.values())
    # add merged_df to final_result dict, key is bin_size
    final_result[bin_size] = merged_df



# record the time end to run the script
end_time = time.time()
#print end time
print("End time: ", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time)))
print("Time used to run the script in seconds: ", end_time - start_time)


# save the final_result dict to pickle file
print("start to save the final_result dict to pickle file")
with open(output_file, 'wb') as f:
    pickle.dump(final_result, f)
print("finished saving the final_result dict to pickle file.")