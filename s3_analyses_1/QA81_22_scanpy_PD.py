import itertools
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import bbknn
from collections import Counter


###############################
#      loading in the data
###############################
data_dir = 'path/to/data'
meta_qa22 = pd.read_csv(data_dir+'/QA22_meta_after_qc.csv',index_col=0)
meta_qa81 = pd.read_csv(data_dir+'/QA81_meta_after_qc.csv', index_col=0)

qa81_data = sc.read_csv(data_dir+'/QA81_raw_counts_after_qc.csv')
qa22_data = sc.read_csv(data_dir+'/QA22_raw_counts_after_qc.csv')

qa81_data.obs = meta_qa81
qa22_data.obs = meta_qa22

 

###############################
#  combining the 2 scanpy object
###############################

data_81_22 = qa81_data.concatenate(qa22_data)


###############################
#  filtering genes not expressed in at least 3 cells
###############################

sc.pp.filter_genes(data_81_22, min_cells=3)

###############################
#  find total number of counts
###############################

data_81_22.obs['n_counts'] = data_81_22.X.sum(axis=1)

###############################
# normalise
###############################

sc.pp.normalize_total(data_81_22, target_sum=1e4)


###############################
# natural log
###############################

sc.pp.log1p(data_81_22)


###############################
# combat with covariate     
###############################

sc.pp.combat(data_81_22) 


###############################
# saving check point
###############################
sc.write('QA81_22_scanpy_PD_object_1',data_81_22)


###############################
# highly variable genes   
###############################

sc.pp.highly_variable_genes(data_81_22,min_mean=0.05, max_mean=13, min_disp=0.1, max_disp=3)

sc.pl.highly_variable_genes(data_81_22)


###############################
# subsetting scanpy object with highly variable genes   
###############################

adata = data_81_22[:,data_81_22.var['highly_variable']]


###############################
# PCA  
###############################

sc.tl.pca(adata, svd_solver='arpack') 

sc.pl.pca(adata, color='Details')


###############################
# finding nearest neighbors
###############################

batch_key = 'batch'
data_neighbor = adata.copy()
pca = data_neighbor.obsm['X_pca']
batch_list = data_neighbor.obs[batch_key].values

bbknn_out = bbknn.bbknn_pca_matrix(pca=pca, batch_list=batch_list,approx=False, metric='euclidean')

data_neighbor.uns['neighbors'] = {}
#we'll have a zero distance for our cell of origin, and nonzero for every other neighbour computed
data_neighbor.uns['neighbors']['params'] = {'n_neighbors': len(bbknn_out[0][0,:].data)+1, 'method': 'umap'}
data_neighbor.uns['neighbors']['distances'] = bbknn_out[0]
data_neighbor.uns['neighbors']['connectivities'] = bbknn_out[1]


###############################
# diffusion map
###############################

sc.tl.diffmap(data_neighbor, n_comps=3)


###############################
# setting root cell for pseudotime calculation
###############################

data_neighbor.uns['iroot'] = 365


###############################
# calculating pseudotime
###############################

sc.tl.dpt(data_neighbor, n_dcs=3)



###############################
# saving check point
###############################
sc.write('QA81_22_scanpy_PD_object_2', data_neighbor)



