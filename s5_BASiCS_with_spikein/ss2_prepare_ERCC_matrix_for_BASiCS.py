import pandas as pd
import numpy as np
from collections import Counter

thermo_fiser = pd.read_csv('cms_095046_1.txt',index_col=0,sep='\t')
thermo_fiser.index = thermo_fiser['ERCC ID']

ercc_22 = pd.read_csv('ercc_qa22_92.csv',index_col=0)
ercc_81 = pd.read_csv('ercc_qa81_92.csv',index_col=0)

dictionary = dict(zip(thermo_fiser['ERCC ID'], thermo_fiser['concentration in Mix 1 (attomoles/ul)']))

def calculate_molecule(conc):
    volume = 0.1
    dilution_factor = 300000
    
    molecule = conc * (10e-18) * (6.022 * 10e23)* volume * dilution_factor
    
    return molecule

dict_ercc={'SpikeID':[],'SpikeInput':[]}

for ercc_index in ercc_22.index:
    concentration = dictionary[ercc_index]
    molecule = calculate_molecule(concentration)
    dict_ercc['SpikeID'].append(ercc_index)
    dict_ercc['SpikeInput'].append(molecule)
    
df_ercc = pd.DataFrame(dict_ercc)

df_ercc.to_csv('ERCC_92_num_molecules_for_BASiCS.csv')



