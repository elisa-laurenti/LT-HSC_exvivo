
import scanpy as sc
import numpy as np
import pandas as pd


adata = sc.read('slx_20048_cells_225_genes_65988_adata_qc_meta_fixed.h5ad')


meta = adata.obs

counts = adata.X

t_counts = counts.T


meta.to_csv('meta_20048.csv')

np.savetxt('counts_20048.txt', t_counts,delimiter=',')

adata.var.to_csv('adata_var_20048.csv')

