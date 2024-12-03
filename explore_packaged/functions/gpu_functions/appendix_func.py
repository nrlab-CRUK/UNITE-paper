result = (df_dummy
	.to_pandas()
	.groupby(['id', 'frag_len'], observed = False)
	.agg({output_col: 'sum'})
	.reset_index()
)



result = motif_filter_nrow(input_df, colname_to_use=colname_to_use, output_col=output_col)





# viz etc..





# convert df2 to matrix, with id as row index, frag_len as column index, frag_len_seen as values
mat = df2.pivot(index='id', columns='frag_len', values='frag_len_seen')

# convert df2 to np array
mat_array = mat.to_numpy()

# min-max normalization

def min_max_normalization(arr):
    min_val = np.min(arr)
    max_val = np.max(arr)
    normalized_arr = (arr - min_val) / (max_val - min_val)
    return normalized_arr

# Example usage
normalized_array = min_max_normalization(mat_array)

# reverse the order of rows in normalized_array
normalized_array = normalized_array[::-1]

# plot the normalized array as heatmap, save the plot to pdf file
import seaborn as sns
import matplotlib.pyplot as plt

# use imshow to plot the normalized array, no interpolation, use gray colormap
plt.imshow(mat, cmap='gray', interpolation='none')

# show x and y ticks

plt.savefig('/scratchc/nrlab/wang04/ulyses_second_batch/try_ulyses_gpu/SLX-18445.SXTLI192.HiSeq4000.bam.0.1x.mrkdup.bam.ulyses_bind_olap_chunks_imshow_gray.pdf')

# clear the plot
plt.clf()




# save the dataframe to csv file under path: /scratchc/nrlab/wang04/ulyses_second_batch/try_ulyses_gpu
df2.to_csv('/scratchc/nrlab/wang04/ulyses_second_batch/try_ulyses_gpu/SLX-18445.SXTLI192.HiSeq4000.bam.0.1x.mrkdup.bam.ulyses_bind_olap_chunks.csv')

# print the unique id
print(df1['id'].unique())

