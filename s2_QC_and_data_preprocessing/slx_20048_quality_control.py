import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import os
import smqpp
import re 
from collections import Counter


smqpp.generate_feature_table('Homo_sapiens.GRCh38.81_and_ERCC92.gtf',
                                           'geneTable_ensembl_geneName.tsv')

Indir = 'star_featureCounts_output'
ftable_loc = 'geneTable_ensembl_geneName.tsv'

adata_raw =  smqpp.read_in_files(Indir, ftable_loc)

prelim_meta = pd.read_csv('meta_20048_fixed.csv',index_col=0)

prelim_meta = prelim_meta.loc[adata_raw.obs.index,:]


prelim_meta['QC_Unassigned_Unmapped'] = adata_raw.obs['QC_Unassigned_Unmapped'].values
prelim_meta['QC_Unassigned_NoFeatures'] = adata_raw.obs['QC_Unassigned_NoFeatures'].values
prelim_meta['QC_Unassigned_Ambiguity'] = adata_raw.obs['QC_Unassigned_Ambiguity'].values

adata_raw.obs = prelim_meta


##########################
#  defining qc parameters and tunning the qc function
###########################

cutoff = {}
cutoff['nMapped (log10)'] = np.log10(2*(10**5))
cutoff['nNuclear (log10)'] = 0
cutoff['fGenes:nTotal'] = 0.2
cutoff['nHCGenes'] = 4000
cutoff['mito:nGenes'] = 0.2
cutoff['nERCC:nMapped'] = 0.2
cutoff['QC_Unmapped:nTotal'] = 0.25    ### has to match what you choose down below
cutoff['nNuclear:nMapped'] = 0.89    ### has to match what you choose down below


threshold_qc  = {'nMapped (log10)': 0,
           'nNuclear (log10)': np.log10(120000),
           'fGenes:nTotal': 0.4,
           'nHCGenes': 5000,
           'mito:nGenes': 0.2,
           'nERCC:nMapped': 0.2,
           'QC_Unmapped:nTotal': 0.25,
            'nNuclear:nMapped':0.89
          }

adata_qc, QCdata , failed_idx, nTotal  = my_smartseq_qc2(adata_raw, cutoff=threshold_qc, 
                          MTpattern = 'MT-', 
                          title='SLX_20048')


adata_qc

######################
# combining qc results to a meta table
######################

meta_raw = adata_raw.obs

Counter(meta_raw.index == nTotal.index)
meta_raw['nTotal'] = nTotal.values
meta_raw['failed_qc'] = 'passed'
meta_raw['failed_qc']= np.where( failed_idx, 'failed', meta_raw['failed_qc'] ) 

for i in range(len(QCdata.keys())):
    
    meta_raw[list(QCdata.keys())[i]] = QCdata[list(QCdata.keys())[i]]

    print(list(QCdata.keys())[i])

meta_raw['comment'] = meta_raw['comment'].replace(np.nan, '-')


######################### 
#  saving objects
#########################

meta_raw.to_csv('meta_after_qc_for_plotting.csv')


sc.write('slx_20048_cells_225_genes_65988_adata_qc_meta_fixed',adata_qc)

sc.write('slx_20048_cells_276_genes_65988_adata_raw_meta_fixed',adata_raw)






#########################
#  this should be run first before the function is used
#########################
def my_smartseq_qc2(adata, cutoff=cutoff,
            MTpattern = 'MT-', ncols=4, figsize=(10,7), s=10, title=None, save=None):

    ## don't know why make_unique did not work sometimes, so do this again
    adata.var_names_make_unique()
    mito_genes = [name for name in adata.var_names if name.startswith(MTpattern)]
    if not mito_genes:
        raise ValueError('Pls enter correct MTpattern')
    else:
        print('mito_genes: '+str(mito_genes))
    mitoCNT = np.sum(adata[:,mito_genes].X, axis=1).copy()
    nuclearCNT = np.sum(adata[:,~np.in1d(adata.var_names,mito_genes)].X, axis=1).copy()
    if 'ERCC' not in adata.obsm_keys():
        erccCNT = np.zeros(adata.shape[0])
    else:
        erccCNT = np.sum(adata.obsm['ERCC'], axis=1)
    qcNames = [x for x in adata.obs_keys() if 'QC' in x]
    if not qcNames:
        qcCNT = np.zeros(adata.shape[0])
    else:
        qcCNT = np.sum(adata.obs[qcNames], axis=1).values
    nTotal = mitoCNT + nuclearCNT + erccCNT + qcCNT
    nMapped = mitoCNT + nuclearCNT + erccCNT
    nGenes = mitoCNT + nuclearCNT
    nHCGenes = np.sum(adata[:,~np.in1d(adata.var_names,mito_genes)].X.T*(10**6)/(nuclearCNT+1) > 10, axis=0)

    QCdata = {}
    QCdata['nMapped (log10)'] = np.log10(nMapped+1)
    QCdata['nNuclear (log10)'] = np.log10(nuclearCNT+1)
    QCdata['fGenes:nTotal'] = nGenes/(nTotal+1)
    QCdata['nHCGenes'] = nHCGenes
    QCdata['mito:nGenes'] = mitoCNT/(nGenes+1)
    QCdata['nERCC:nMapped'] = erccCNT/(nMapped+1)
    QCdata['nNuclear:nMapped'] = nuclearCNT/(nMapped+1)
    for qcIndex in qcNames:
        QCdata[qcIndex.replace('_Unassigned', '')+':nTotal'] = adata.obs[qcIndex]/(nTotal+1)
    
    # cells failed QC
    compara = {}
    compara['nMapped (log10)'] = '<'
    compara['nNuclear (log10)'] = '<'
    compara['fGenes:nTotal'] = '<'
    compara['nHCGenes'] = '<'
    compara['mito:nGenes'] = '>'
    compara['nERCC:nMapped'] = '>'
    compara['QC_Unmapped:nTotal'] = '>'     
    compara['nNuclear:nMapped'] = '<'    
    
    failed = []
    for k in cutoff.keys():
        if compara[k] == '<':
            failed.append(QCdata[k] < cutoff[k])
        else:
            failed.append(QCdata[k] > cutoff[k])
    failed = np.vstack(failed)
    failed_idx = np.sum(failed, axis=0)>0
    print('Number of passed cells: '+str(sum(np.sum(failed, axis=0)==0)))
    print('Number of failed cells: '+str(sum(failed_idx)))
    
    # plotting
    nrows = int(np.ceil(len(QCdata.keys())/ncols))
#     print(nrows)
    fig, ax = plt.subplots(nrows,ncols, figsize=figsize)
    for i in range(len(QCdata.keys())):
        colidx = i%ncols
        rowidx = np.floor(i/ncols).astype(int)
        ax[rowidx, colidx].scatter(nTotal, QCdata[list(QCdata.keys())[i]], s=s, color='black')
        ax[rowidx, colidx].scatter(nTotal[failed_idx], QCdata[list(QCdata.keys())[i]][failed_idx], s=s, color='red')
        if list(QCdata.keys())[i] in cutoff.keys():
            ax[rowidx, colidx].axhline(y=cutoff[list(QCdata.keys())[i]], color='orange', linestyle='dashed')
        #ax[rowidx, colidx].set_yscale('log',basey=10)
        ax[rowidx, colidx].set_ylabel(list(QCdata.keys())[i])
        ax[rowidx, colidx].grid(False)
    fig.text(0.5, -0.03, 'nTotal', ha='center')
    fig.suptitle(title)
    plt.tight_layout(pad=1)
    fig.subplots_adjust(top=0.88)
    
    if save is not None:
        plt.savefig(save)
        
    adata.obs['n_counts'] = nGenes
    adata.obs['percent_mito'] = mitoCNT/(nGenes+1)
    adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1)
    
    return adata[~failed_idx,:].copy() , QCdata , failed_idx, nTotal


