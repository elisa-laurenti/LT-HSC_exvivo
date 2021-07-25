import numpy as np
import scanpy as sc
import pandas as pd
import os


################################################
#
################################################

adata_22 = sc.read('QA22_after_qc.h5ad')

meta_22 = adata_22.obs

counts_22 = adata_22.X

meta_22.to_csv('QA22_meta_after_qc.csv')

counts_22.to_csv('QA22_raw_counts_after_qc.csv')




################################################
#
################################################


adata_81 = sc.read('QA81_after_qc.h5ad')

meta_81 = adata_81.obs

counts_81 = adata_81.X

meta_81.to_csv('QA81_meta_after_qc.csv')

counts_81.to_csv('QA81_raw_counts_after_qc.csv')


