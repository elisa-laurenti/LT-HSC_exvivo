#%matplotlib nbagg
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import os
import smqpp
import re
import anndata

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list(name='gene_cmap', colors=['lightgrey', 'thistle', 'red', 'darkred']) 

sc.settings.set_figure_params(dpi=80, color_map='viridis', vector_friendly=False,  dpi_save=300)



meta = pd.read_csv('meta_MPB.txt',delimiter='\t')
meta.index = meta['ID']
print(meta.shape)
display(meta.head(1))


meta['con_comb'] = meta['Day'].astype(str)+'_'+meta['Cell_Type'].astype(str)+'_'+meta['Condition'].astype(str)
meta['con_comb'].value_counts().sort_index()


gtpfile = 'Homo_sapians_GRCh38_ERCC92_GFP.gtf'
ftable_loc = 'features.tsv'
smqpp.generate_feature_table(gtpfile, ftable_loc)

Indir = 'MPB_data/Count_files/'
folder_list = np.array(os.listdir(Indir))
folder_list = np.sort(folder_list)
print(folder_list)

# define QC thresholds
cutoff1 = {'nMapped (log10)': 0,
           'nNuclear (log10)': np.log10(10**5.2),
           'fGenes:nTotal': 0.3,
           'nHCGenes': 0,
           'mito:nGenes': 0.2,
           'nERCC:nMapped': 0.4
          }
cutoff2 = {'nMapped (log10)': 0,
           'nNuclear (log10)': np.log10(10**5.2),
           'fGenes:nTotal': 0.3,
           'nHCGenes': 1000,
           'mito:nGenes': 0.2,
           'nERCC:nMapped': 0.4
          }
cutoff3 = {'nMapped (log10)': 0,
           'nNuclear (log10)': np.log10(15000),
           'fGenes:nTotal': 0.1,
           'nHCGenes': 1000,
           'mito:nGenes': 0.2,
           'nERCC:nMapped': 0.2
          }
cutoff4 = {'nMapped (log10)': 0,
           'nNuclear (log10)': np.log10(100000),
           'fGenes:nTotal': 0.2,
           'nHCGenes': 1000,
           'mito:nGenes': 0.2,
           'nERCC:nMapped': 0.2
          }
cutoff = [cutoff1, cutoff2, cutoff3, cutoff4]

adata = []
for i in range(len(folder_list)):
    print(folder_list[i])
    adata_sub = smqpp.read_in_files(Indir+folder_list[i], ftable_loc, method='HTSeqcount')
    adata_sub.obs['eGFP'] = adata_sub[:,'eGFP'].X.flatten()
    adata_sub = adata_sub[:, adata_sub.var_names != 'eGFP'].copy()
    adata_sub.var_names_make_unique()
    adata_sub = smqpp.smartseq_qc(adata_sub,cutoff=cutoff[i], title=folder_list[i], MTpattern='MT-')
    adata.append(adata_sub)
    del adata_sub
    
# combine data
adata = anndata.AnnData.concatenate(adata[0], adata[1], adata[2], adata[3])
adata.obs_names = [x[:-2] for x in adata.obs_names]
# add in meta
adata.obs = pd.concat([adata.obs, meta], axis=1, join='inner')
adata


# Here are the outlier cells either empty or behave weiredly, 15 in total
# 1st list: the cells that are supposed to be empty
outlier1 = adata.obs_names[adata.obs['Condition'] == 'Empty'].values
# 2nd list: the cells that effect HVG selection
outlier2 = np.array(["SLX.12981.i705_i503", "SLX.17052.i706_i506","SLX.17052.i721_i518", "SLX.17052.i729_i518", "SLX.17053.i727_i515"])
# 3nd list: the Day0 cells are clustered into Day3 cells
outlier3 = np.genfromtxt('exCells_outliers.txt', delimiter='\t', dtype=str)
outliers = np.concatenate([outlier1, outlier2, outlier3])
print(len(outliers))


adata = adata[~np.in1d(adata.obs_names, outliers),:].copy()
adata

adata.write('MPB_raw.h5ad')


