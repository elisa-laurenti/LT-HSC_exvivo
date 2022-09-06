import pandas as pd

raw_counts_22 = pd.read_csv('QA22_raw_counts.csv',index_col=0)

raw_counts_81 = pd.read_csv('QA81_raw_counts.csv',index_col=0)

ERCC_22 = raw_counts_22.iloc[:,['ERCC-' in x for x in raw_counts_22.columns]]

ERCC_81 = raw_counts_81.iloc[:,['ERCC-' in x for x in raw_counts_81.columns]]

ERCC_22.to_csv('ercc_qa22_92.csv')

ERCC_81.to_csv('ercc_qa81_92.csv')
