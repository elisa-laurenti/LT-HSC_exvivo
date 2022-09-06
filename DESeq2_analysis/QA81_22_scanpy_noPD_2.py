import itertools
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import bbknn
from collections import Counter


adata = sc.read('QA81_22_scanpy_noPD_object_2.h5ad')

###################################
# calculating UMAP
###################################

sc.tl.umap(adata, n_components=3)


sc.write('QA81_22_scanpy_noPD_object_3.h5ad', adata)


###################################
# calculating force directed graph
###################################


adata = sc.read('QA81_22_scanpy_noPD_object_3.h5ad')

adata.uns['iroot']

sc.tl.draw_graph(adata, root=297)

sc.pl.draw_graph(adata)


adata.write('QA81_22_scanpy_noPD_object_4.h5ad')

####################################
# extracting 1 and 2 blob on FDG
###################################

new_index = adata.obs['CRI_identifier'].str.replace('-','.')+ '.' +adata.obs['Index'].str.replace('-','_')

adata.obs_names = new_index

plot_fa = pd.DataFrame(adata.obsm['X_draw_graph_fa'])

plot_fa.columns = ['fa1','fa2']

plot_fa['Details'] =adata.obs['Details'].values
plot_fa.index= adata.obs_names


fig, ax = plt.subplots(figsize=(10,8))
sns.scatterplot(x="fa1", y="fa2", hue="Details", 
                     data=plot_fa, ax=ax ,palette=["blue", "cornflowerblue", "green", "red"])

ax.legend(loc='best', frameon=False,fontsize=9)
ax.axhline(y=298)
ax.axvline(x=-227)

plt.show()


# === 0h_A
hour0 = plot_fa['Details'] == 'LT_0h'

greater_298 = plot_fa['fa2'] > 298

first_blob_df = plot_fa[hour0 & greater_298]

# === 0h_B
lower_298 = plot_fa['fa2'] < 298 

greater_227 =  plot_fa['fa1'] > -227

second_blob_df = plot_fa[hour0 & greater_227 & lower_298]

first_blob_df.to_csv('first_blob.csv')

second_blob_df.to_csv('second_blob.csv')

