
import numpy as np
import pandas as pd
from skimage.transform import resize
import pickle
import os



def fragmentim(group, resize_dim = True):

	n_channel = group.channel.nunique()
	n_id = group.id.nunique()
	n_isize = group.frag_len.nunique()
	dim = (n_channel, n_id, n_isize)

	result = group['pixel'].to_numpy().reshape(dim).transpose(1, 2, 0)

	if resize_dim is True:
		result = resize(result, (224, 224, n_channel))

	return result


def df_to_tensor(csv_file, tile_size_label = 'tile_size_label', output_dir = None, resize_dim = True):

	if output_dir is None:
		output_dir = os.path.dirname(csv_file)
	
	df = pd.read_csv(csv_file)[[tile_size_label, "channel", "id", "frag_len", "pixel"]]
	df_grouped = df.groupby([tile_size_label])
	tensor = df_grouped.apply(fragmentim, resize_dim = resize_dim)

	# save the tensor to a file
	base_name = os.path.basename(csv_file)
	output_file = os.path.join(output_dir, base_name + ".tensor") 
	
	pickle.dump(tensor, open(output_file, 'wb'))

	# message the output file
	return print(f"Saved file: {output_file}")


def load_pickle(fin):
	with open(fin, 'rb') as f:
		obj = pickle.load(f)
	return obj