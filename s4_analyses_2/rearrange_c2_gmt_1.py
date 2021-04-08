import pandas as pd
import numpy as np
from collections import Counter

##################
##  download from http://www.gsea-msigdb.org/gsea/downloads.jsp
##################

with open('c2.all.v7.1.symbols.gmt') as f:
    lines = [line.rstrip() for line in f]


name_website_dict ={}
name_gene_list_dict ={}

for line in lines:
    items = line.split('\t')
    
    name = items[0]
    website = items[1]
    

    items.pop(0)
    items.pop(0)
    
    name_website_dict[name] = website
    
    name_gene_list_dict[name] = items


#########################
#  saving the pathways
#########################

save_dir='gene_set_gmt/c2_all_v7_1_symbols_csv'

for key, value in name_gene_list_dict.items():
    temp_dict = {key:value}
    temp_df = pd.DataFrame(temp_dict)
    
    file_name = f'{save_dir}/{key}.csv'
    temp_df.to_csv(file_name)

