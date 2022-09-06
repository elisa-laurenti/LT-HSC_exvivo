import scanpy as sc
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from scipy.stats import shapiro

meta22 = pd.read_csv('meta_qa22_only__scEntropy_scores.csv',index_col=0)

meta81 = pd.read_csv('meta_qa81_only__scEntropy_scores.csv',index_col=0)

meta81['old_Entropy_score'] = meta81['Entropy_score'].values

meta22['old_Entropy_score'] = meta22['Entropy_score'].values

meta22.Entropy_score = meta22.Entropy_score - (np.mean(meta22.Entropy_score) - np.mean(meta81.Entropy_score))

frames = [meta22, meta81]
result = pd.concat(frames)

####################################################################

fig = plt.figure(figsize=(6,6))

fig.suptitle('batch adjusted scEntropy norm10k no log qa22_81', fontsize=14)
ax = fig.add_subplot(1,1,1)

plt.ylim(4, 5.6)

sns.boxplot(x="Details", y="Entropy_score", data=result, ax=ax)

ax.set_xlabel('incubation time')

plt.tight_layout()
fig.subplots_adjust(top=0.93)
plt.show()

