import itertools
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import bbknn
from collections import Counter

data_dir = 'path/to/data'
meta_qa22 = pd.read_csv(data_dir+'/QA22_meta_after_qc.csv',index_col=0)
meta_qa81 = pd.read_csv(data_dir+'/QA81_meta_after_qc.csv', index_col=0)

qa81_data = sc.read_csv(data_dir+'/QA81_raw_counts_after_qc.csv')
qa22_data = sc.read_csv(data_dir+'/QA22_raw_counts_after_qc.csv')

qa81_data.obs = meta_qa81
qa22_data.obs = meta_qa22


qa22_data = qa22_data[qa22_data.obs['Details'] != 'LT_72h_PD',:]
qa81_data = qa81_data[qa81_data.obs['Details'] != 'LT_24h_PD',:]

#################################################

sc.pp.filter_genes(qa22_data, min_cells=3)
sc.pp.filter_genes(qa81_data, min_cells=3)

sc.pp.normalize_total(qa22_data, target_sum=1e4)
sc.pp.normalize_total(qa81_data, target_sum=1e4)


################################################

df_counts_22 = pd.DataFrame(qa22_data.X, columns= qa22_data.var_names, index=qa22_data.obs_names)
df_counts_81 = pd.DataFrame(qa81_data.X, columns= qa81_data.var_names, index=qa81_data.obs_names)

t_22 = df_counts_22.T
t_81 = df_counts_81.T

t_22.to_csv('counts_qa22__norm10k_noLOG.csv')
t_81.to_csv('counts_qa81__norm10k_noLOG.csv')

