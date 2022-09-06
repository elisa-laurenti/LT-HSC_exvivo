#!/usr/bin/env python3

from scipy.stats import chi2
import smqpp
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import os
import smqpp
import re
import anndata
from collections import Counter
import seaborn as sns
from scipy import stats
import scipy
import statsmodels.api as sm
import statsmodels.stats.multitest as multi


def my_tech_var(adata, useERCC=True, cvThresh=.3, quant=.8, minBiolDisp=.5**2, 
           fdr=.1, meanForFit=None, copy=False):
    
    if 'sf_gene' not in adata.obs_keys():
        print("sf_gene is not found, redoing normalisation")
        smqpp.normalise_data(adata)
    
    data = np.exp(adata.X)-1
    aMean = np.mean(data, axis=0)
    aStd = np.std(data, axis=0)
    cv2a = (aStd/aMean)**2
    if useERCC:
        if 'ERCC' not in adata.obsm_keys():
            raise ValueError('ERCC does not exist, check data.')
        ercc_data = np.exp(adata.obsm['ERCC_norm'])-1
        sMean = np.mean(ercc_data, axis=0)
        sStd = np.std(ercc_data, axis=0)
        cv2s = (sStd/sMean)**2
    else:
        sMean = aMean
        cv2s = cv2a
    if meanForFit is None:
        meanForFit = np.quantile(sMean[cv2s>cvThresh], quant)
    print('MeanForFit: ', str(meanForFit))
    useForFit = (sMean>=meanForFit)
    print(np.sum(useForFit))
    
    a1tilde = 1/sMean[useForFit]
    x = sm.add_constant(a1tilde, prepend=False)
    y = cv2s[useForFit]
    link_func = sm.genmod.families.links.identity
    fit = sm.GLM(y, x, family=sm.families.Gamma(link=link_func)).fit()

    a0 = fit.params['const']
    a1t = fit.params[0]
    df = data.shape[0]-1
    m = data.shape[0]
    xi = None
    
    if useERCC:
        xi = np.mean(1/adata.obs['sf_ercc'])
        psi = xi + (a1t - xi)*np.mean(adata.obs['sf_ercc']/adata.obs['sf_gene'])
        cv2th = a0 + minBiolDisp + a0*minBiolDisp
        testDenom = (aMean*psi + cv2th*aMean**2)/(1+cv2th/m)
        pA = 1 - chi2.cdf((aStd**2)*df/testDenom, df=df)
    else:
        psi = a1t
        chi2_values = df * cv2s / (psi / sMean + a0)
        pA = 1 - chi2.cdf(chi2_values ,df=df)
    
    pA[np.isnan(pA)] = 1
    _, padj, _, _ = multi.multipletests(pA, method='fdr_bh')
    
    highVarGenes = adata.var_names[padj < fdr]
    print('Length of HVGs: '+ str(len(highVarGenes)))
    
    adata.uns['varGenes'] = {}
    adata.uns['varGenes']['parameters'] = {}
    adata.uns['varGenes']['genes'] = {}
    adata.uns['varGenes']['ercc'] = {}
    
    adata.uns['varGenes']['parameters']['minBiolDisp'] = minBiolDisp
    adata.uns['varGenes']['parameters']['a1tilde'] = a1t
    adata.uns['varGenes']['parameters']['a0'] = a0
    adata.uns['varGenes']['parameters']['psi'] = psi
    adata.uns['varGenes']['parameters']['xi'] = xi
    adata.uns['varGenes']['parameters']['df'] = df
    adata.uns['varGenes']['parameters']['useERCC'] = useERCC
    adata.uns['varGenes']['parameters']['meanForFit'] = meanForFit
    adata.uns['varGenes']['parameters']['useForFit'] = np.sum(useForFit)
    adata.uns['varGenes']['parameters']['cvThresh'] = cvThresh
    adata.uns['varGenes']['parameters']['quant'] = quant
    
    adata.uns['varGenes']['genes']['mean'] = aMean
    adata.uns['varGenes']['genes']['cv2'] = cv2a
    adata.uns['varGenes']['genes']['highVar'] = (padj < fdr)
    adata.uns['varGenes']['ercc']['mean'] = sMean
    adata.uns['varGenes']['ercc']['cv2'] = cv2s
    
    if copy:
        return adata.copy()



def my_normalise_data(adata, reCalSF=True, method='ExpAllC', copy=False):

    if reCalSF:
        print('Calculate SF for genes:')
        sf_genes = smqpp.est_size_factor(adata.X, method=method)
        adata.obs['sf_gene'] = sf_genes
        if 'ERCC' in adata.obsm_keys():
            if adata.obsm['ERCC'].size !=0:
                print('Calculate SF for erccs:')
                sf_ercc = smqpp.est_size_factor(adata.obsm['ERCC'].values, method=method)## original adata.obsm['ERCC']
                adata.obs['sf_ercc'] = sf_ercc
    else:
        if 'sf_gene' not in adata.obs_keys():
            raise ValueError('sf_gene is not found in .obs, please set reCalSF=True.')

    adata.X = np.log1p(adata.X/adata.obs['sf_gene'][:,None])
    if 'ERCC' in adata.obsm_keys():
        if adata.obsm['ERCC'].size !=0:
            adata.obsm['ERCC_norm'] = np.log1p(adata.obsm['ERCC']/adata.obs['sf_ercc'][:,None])
    if copy:
        return adata.copy()


