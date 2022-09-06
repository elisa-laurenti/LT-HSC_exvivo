import itertools
import pickle

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from collections import Counter

import my_scEntropy


print('loading in data')

meta = pd.read_csv('QA22_meta_after_qc.csv',index_col=0)
counts = pd.read_csv('counts_qa22__norm10k_noLOG.csv',index_col=0)

print(Counter(counts.columns == meta.index))

print(meta.shape)

print(counts.shape)


print('calculating entropy')

scEntropy_arr = my_scEntropy.scEntropy(counts, option='RCSA')


print(scEntropy_arr.shape)


print('saving results to meta')

meta['Entropy_score'] = scEntropy_arr

print(meta.shape)

meta.to_csv('meta_qa22_only__scEntropy_scores.csv')


print('all done')
