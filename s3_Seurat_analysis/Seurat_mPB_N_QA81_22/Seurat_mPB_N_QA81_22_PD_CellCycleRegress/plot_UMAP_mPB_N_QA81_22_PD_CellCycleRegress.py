import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from collections import Counter
# import smqpp
import seaborn as sns
import plotly.express as px
import umap

meta = pd.read_csv('meta_PD_QA81_22_N_mPB_CycleRegress__UMAP.csv',index_col=0)

combined_time=[]
for index, row in meta.iterrows():
    batch = row['orig.ident']
    
    if batch == 'mPB':
        new_time=row['Day'] +'_'+row['Cell_Type']+'_'+row['Condition']
        combined_time.append(new_time)
    else:
        combined_time.append(row['Details'])
        
meta['both_time'] = combined_time

meta.to_csv('meta_PD_QA81_22_N_mPB_CycleRegress__UMAP_forPlotting.csv')

##############################################

meta = pd.read_csv('meta_PD_QA81_22_N_mPB_CycleRegress__UMAP_forPlotting.csv',index_col=0)

color_dict = {'LT_0h':'blue', 
             'LT_6h':'cornflowerblue',
             'LT_24h_UNTR': 'mediumseagreen', 
             'LT_72h_UNTR': 'darkred',
             'LT_24h_PD':'blueviolet',
             'LT_72h_PD':'hotpink',
             'Day0_LT-HSC_NT':'darkorange',
             'Day3_LT-HSC_GFP+PD-':'lawngreen',
             'Day3_LT-HSC_GFP+PD+':'black'}

fig, ax = plt.subplots(1, 1, figsize=(10,6))

fig.suptitle('Seurat alignment umap mPB_PD and CB_PD', fontsize=18)

sns.scatterplot(data=meta, x="UMAP1", y="UMAP2", hue="both_time",ax=ax, palette=color_dict)

ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('UMAP1',fontsize=14)
ax.set_ylabel('UMAP2',fontsize=14)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
fig.subplots_adjust(top=0.93, left =0.2)
plt.show()

