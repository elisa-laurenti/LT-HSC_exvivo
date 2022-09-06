import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from collections import Counter
import seaborn as sns
import umap

meta = pd.read_csv('meta_PD_QA81_22_N_mPB__UMAP.csv',index_col=0)

combined_time=[]
for index, row in meta.iterrows():
    batch = row['orig.ident']
    
    if batch == 'mPB':
        new_time=row['Day'] +'_'+row['Cell_Type']+'_'+row['Condition']
        combined_time.append(new_time)
    else:
        combined_time.append(row['Details'])


meta['both_time'] = combined_time
 
meta.to_csv('meta_PD_QA81_22_N_mPB__UMAP_forPlotting.csv')

###############################################################

meta = pd.read_csv('meta_PD_QA81_22_N_mPB__UMAP_forPlotting.csv',index_col=0)

sub_meta = meta[meta['both_time']=='LT_0h']


fig, ax = plt.subplots(1, 1, figsize=(10,6))

fig.suptitle('Seurat alignment umap mPB_PD and CB_PD', fontsize=18)
sns.scatterplot(data=sub_meta, x="UMAP1", y="UMAP2", hue="both_time",ax=ax)

ax.set_xlabel('UMAP1',fontsize=14)
ax.set_ylabel('UMAP2',fontsize=14)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
fig.subplots_adjust(top=0.93, left =0.2)
plt.show()


top_4 = sub_meta[sub_meta['UMAP2']>3.1]

new_color =[]

for i in sub_meta['UMAP2'].values:
    if 3.1354 < i < 3.17278:
#         print(i)
        new_color.append(f'root_{i}')
    else:
        new_color.append('not_root')
sub_meta['new_color'] = new_color

fig, ax = plt.subplots(1, 1, figsize=(10,6))

fig.suptitle('Seurat alignment umap mPB_PD and CB_PD', fontsize=18)

sns.scatterplot(data=sub_meta, x="UMAP1", y="UMAP2", hue="new_color",ax=ax)

ax.set_xlabel('UMAP1',fontsize=14)
ax.set_ylabel('UMAP2',fontsize=14)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
fig.subplots_adjust(top=0.93, left =0.2)
plt.show()

sub_meta[sub_meta['new_color']!='not_root']

