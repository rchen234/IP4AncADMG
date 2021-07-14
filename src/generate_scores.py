from test_latent_scores import generate_scores_bidirect,generate_scores_bidirect_m3hc
import numpy as np



datasets = ['example.txt']


for dataset in datasets:

	with open('../Instances/data/'+dataset, 'rb') as f:
	    data = np.loadtxt(f, skiprows=0)

	# visible
	observed_data = data


	c_size = 2
	single_parent_size = 3
	other_c_parent_size = 1
	file_name = '../Instances/data/score_' + dataset[:-4]
	print(file_name)


	generate_scores_bidirect(observed_data,
	                         single_c_parent_size = single_parent_size,
	                         other_c_parent_size = other_c_parent_size,
	                         c_size = c_size,
	                         file_name = file_name)