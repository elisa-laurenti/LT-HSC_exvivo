import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

meta = pd.read_csv('meta_noPD_cycle_regress_pseudotime_UMAP.csv',index_col=0)

norm_2D = plt.Normalize(meta['pseudotime_2D'].min(), meta['pseudotime_2D'].max())
sm_2D = plt.cm.ScalarMappable(cmap="viridis", norm=norm_2D)
sm_2D.set_array([])


fig, axs = plt.subplots(1, 1,figsize=(12,12),sharey=False)

fig.suptitle('noPD Cycle Regress UMAP Monocle pseudotime', fontsize=14)


D2 = axs.scatter(meta["UMAP1"], meta["UMAP2"],
                     c=meta["pseudotime_2D"], cmap="viridis")


fig.colorbar(D2, ax=axs,pad=0.0001)


axs.set_title('monocle pseudotime 2D')

axs.set_xlabel('UMAP1')
axs.set_ylabel('UMAP2')


plt.tight_layout()
fig.subplots_adjust(top=0.92)
plt.show()