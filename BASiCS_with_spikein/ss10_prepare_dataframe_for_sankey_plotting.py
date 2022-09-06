import pandas as pd
import numpy as np
from collections import Counter, defaultdict

timepoint_unique = pd.read_csv('correct_unqiue_BASiCS_DV_genes_at_each_timepoint.csv',index_col=0)

cont_down = pd.read_csv('Continuous_down.csv')
cont_up = pd.read_csv("Continuous_up.csv")
trans_up = pd.read_csv("Transient_up.csv")
trans_down = pd.read_csv("Transient_down.csv")
up_post_6 = pd.read_csv("up_post_6h.csv")


all_genes_trends = list(cont_down.geneName.values)+list(cont_up.geneName.values)+list(trans_up.geneName.values)+list(trans_down.geneName.values)+list(up_post_6.geneName.values)

dict_trend_genes = {'cont_down':cont_down, 'cont_up':cont_up, 
                    'trans_up':trans_up, 'trans_down':trans_down,
                   'up_host_6':up_post_6}

dict_for_df = defaultdict(list)

for timepoint in timepoint_unique.columns:
    temp_genes = list(timepoint_unique[timepoint].values)
    just_genes = list(filter(('--').__ne__, temp_genes))

    
    total_in_trends = [i for i in just_genes if i in all_genes_trends]
    print(len(total_in_trends))
    
    for key, value in dict_trend_genes.items():
        dict_for_df['Source'].append(timepoint)
        dict_for_df['Target'].append(key)
        
        in_trend = [i for i in just_genes if i in value.geneName.values]
        
        percentage = round((len(in_trend) / len(total_in_trends))*100, 1)
        
        dict_for_df['Value'].append(percentage)
    
plotting_df = pd.DataFrame(dict_for_df)

plotting_df.to_csv('dataframe_for_sankey_plotting_BASiCS_DV_N_geneTrends.csv')

