import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import Counter
from sklearn.cluster import KMeans

meta_plot = pd.read_csv('meta_noPD_cycle_regress_pseudotime_UMAP.csv',index_col=0)

only_0h = meta_plot[meta_plot.Details == 'LT_0h']

umap = only_0h.loc[:,['UMAP1','UMAP2']]
umap.head(1)

kmeans = KMeans(n_clusters=2, random_state=0).fit(umap.values)

umap['kmean'] = kmeans.labels_
umap['kmean'] = umap['kmean'].astype(str)

fig, axs = plt.subplots(1,1,figsize=(6,6))
sns.scatterplot(x = "UMAP1", y = "UMAP2", data = umap, hue = "kmean",ax=axs) 
plt.show()

new_kmean=[]
for ind in meta_plot.index:
    if ind in umap.index:
        temp_df = umap.loc[ind,:]
        new_kmean.append(temp_df['kmean'])
    else:
        new_kmean.append('other')
meta_plot['kmean'] = new_kmean

meta_plot.to_csv('meta_noPD_cycle_regress_pseudotime_UMAP_KMeans_0h_2clusters.csv')

########################################################

new_timepoints = []

for index, row in meta_plot.iterrows():
    if row['kmean'] == 'other':
        new_timepoints.append(row['Details'])
    elif row['kmean'] == '1':
        new_timepoints.append(row['kmean'])
    elif row['kmean'] == '0':
        new_timepoints.append(row['kmean'])
        
meta_plot['timepoints'] = new_timepoints

meta_plot['timepoints'] = meta_plot['timepoints'].map({'0': 'LT_0h_B',
                             '1': 'LT_0h_A',
                             'LT_6h': 'LT_6h',
                             'LT_24h_UNTR': 'LT_24h_UNTR',
                             'LT_72h_UNTR': 'LT_72h_UNTR'})

meta_plot.to_csv('meta_noPD_cycle_regress_pseudotime_UMAP_KMeans_0h_2clusters_for_deseq.csv')

######################################################

meta['kmean'] = meta['kmean'].map({'0':'LT_0h_B','1':'LT_0h_A','other':'other'})

color_dict ={'LT_0h':'darkblue',
             'LT_6h':'cornflowerblue',
             'LT_24h_UNTR':'seagreen',
             'LT_72h_UNTR':'darkred',
             'other':'lightslategrey',
            'LT_0h_B':'C1',
            'LT_0h_A':'darkmagenta'}

fig, axs = plt.subplots(2, 2,figsize=(11,11))
axs = axs.flatten()
fig.suptitle('noPD Cycle Regress UMAP', fontsize=14)

sns.scatterplot(x = "UMAP1", y = "UMAP2", data = meta, hue = "kmean", palette=color_dict,
                ax=axs[0]) 

sns.scatterplot(x = "UMAP1", y = "UMAP2", data = meta, hue = "Details", palette=color_dict,
                ax=axs[1]) 

sns.scatterplot(x = "UMAP1", y = "UMAP2", data = meta, hue = "orig.ident", 
                ax=axs[2]) 

axs[0].set_title('Kmean 2 clusters')
axs[1].set_title('timepoints')

plt.tight_layout()
fig.subplots_adjust(top=0.94)
plt.show()