import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import os
import smqpp
import re 
from collections import Counter
from my_functions import my_function_smqpp


adata = sc.read('slx_20048_cells_225_genes_3388_logNorm_scaled_HVG_meta_fixed')



#################
#  PCA
#################
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

##################
#  neighbors
##################
sc.pp.neighbors(adata)


##################
# UMAP
##################

sc.tl.umap(adata)


##################
# saving
##################

sc.write('slx_20048_cells_225_genes_3388_logNorm_scaled_HVG_meta_fixed_pca_umap',adata)


