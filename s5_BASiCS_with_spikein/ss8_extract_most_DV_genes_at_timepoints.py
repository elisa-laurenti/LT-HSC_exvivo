import pandas as pd
import numpy as np
import glob
from collections import Counter
import pickle

file_list = glob.glob('BASiCs_DE_results_separate_full_csv_ERCC_Regression__noFilter/*residual_over-dispersion*.csv')

dict_for_df={}

for file in file_list:
    file_name = file.split('/')[1]
    temp_name1 = file_name.replace('__residual_over-dispersion_ERCC_Regression__noFilter.csv','')
#     temp_name2 = temp_name1[3:]
    print(file_name)
    print(temp_name1)
#     print(temp_name2)
    temp_file = pd.read_csv(file,index_col=0)
    
    print(temp_file.shape)
#     print()
    dict_for_df[temp_name1] = temp_file


with open('all_timepoints_residual_over-dispersion_dataframes.pkl', 'wb') as file:
      
    # A new file will be created
    pickle.dump(dict_for_df, file)

###################################################

dict_unique_genes={}

for timepoint in ['0h','6h','24h','72h']:
    temp_list = []
    for key in dict_for_df.keys():
        if timepoint in key:
            temp_list.append(key)
    print(timepoint)
    print(temp_list)
    
    temp_genes= []
    for comparison in temp_list:
        temp_df = dict_for_df[comparison]
        filtered_df = temp_df[temp_df['ResultDiffResDisp'] == f'{timepoint}+']
#         print(filtered_df['ResultDiffResDisp'].value_counts())
        temp_genes = temp_genes + list(filtered_df['GeneName'].values)
    print(len(temp_genes))
    print(len(set(temp_genes)))
    
    dict_unique_genes[timepoint] = list(set(temp_genes))
    
for key, value in dict_unique_genes.items():
    if len(value) == 1092:
        pass
    else:
        number = 1092 - len(value)
        new_list = ['--']*number
        
        dict_unique_genes[key] = value + new_list
df_unique = pd.DataFrame(dict_unique_genes)

df_unique.to_csv('unqiue_BASiCS_DV_genes_at_each_timepoint.csv')

#################################################

new_dict = dict()

new_dict['GeneName'] = dict_for_df['de_0h__24h']['GeneName'].values
new_dict['ResDisp_0h'] = dict_for_df['de_0h__24h']['ResDisp1'].values

new_dict['ResDisp_6h'] = dict_for_df['de_6h__24h']['ResDisp1'].values

new_dict['ResDisp_24h'] = dict_for_df['de_24h__72h']['ResDisp1'].values

new_dict['ResDisp_72h'] = dict_for_df['de_0h__72h']['ResDisp2'].values

res_df = pd.DataFrame(new_dict)

timepoints = ['0h','6h','24h','72h']
index_dict = {'0h':0,'6h':1,'24h':2,'72h':3}

unqiue_df=dict()

for timepoint in timepoints:
    temp_series = list(df[timepoint].values)
    
    new_gene_list=[]
    
    for gene in temp_series:
        
        if gene == '--':
            break
            
        temp_gene_row = res_df[res_df['GeneName'] == gene]
        
        temp_list = [temp_gene_row['ResDisp_0h'].values[0],
                     temp_gene_row['ResDisp_6h'].values[0],
                     temp_gene_row['ResDisp_24h'].values[0],
                    temp_gene_row['ResDisp_72h'].values[0]]
        
        max_value = max(temp_list)
        
        index = temp_list.index(max_value)
        
        if index == index_dict[timepoint]:
            new_gene_list.append(gene)
    
    print(timepoint)
    print(len(new_gene_list))
    unqiue_df[timepoint] = new_gene_list


for key, value in unqiue_df.items():
    if len(value) == 864:
        pass
    else:
        number = 864 - len(value)
        new_list = ['--']*number
        
        unqiue_df[key] = value + new_list
        
df_unique = pd.DataFrame(unqiue_df)


df_unique.to_csv('correct_unqiue_BASiCS_DV_genes_at_each_timepoint.csv')


