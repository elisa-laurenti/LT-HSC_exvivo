import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import os
import smqpp
import re 
from collections import Counter
from my_functions import my_function_smqpp

adata = sc.read('slx_20048_cells_225_genes_65988_adata_qc_meta_fixed.h5ad')

adata

sc.pp.filter_genes(adata, min_cells=1)

my_function_smqpp.my_normalise_data(adata)

adata.raw = adata

############
#  check point 1
#############
sc.write('slx_20048_cells_225_genes_40449_logNorm_adata_qc_meta_fixed',adata)


##############
#  HVG
##############

adata_hvg  = my_function_smqpp.my_tech_var(adata, useERCC=True, copy=True)

smqpp.plot_tech_var(adata_hvg)

adata_hvg = adata_hvg[:,adata_hvg.uns['varGenes']['genes']['highVar']].copy()


###############
#  regressing out effects
###############
sc.pp.regress_out(adata_hvg, ['n_counts', 'percent_mito'])



####  scale for pca
sc.pp.scale(adata_hvg)



############################
#  
############################
sc.write('slx_20048_cells_225_genes_3388_logNorm_scaled_HVG_meta_fixed',adata_hvg)

