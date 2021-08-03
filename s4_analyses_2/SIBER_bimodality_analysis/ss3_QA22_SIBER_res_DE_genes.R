#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(SIBERG))

load('siber_QA22_by_timepoint.RData')

deg_genetable <- readRDS('DE_unique_genes.rds')


new_ln_list <- list()

for ( timepoint in names(df_ln_list)){
    
    temp_ln <- as.data.frame(df_ln_list[[timepoint]])
    
#     temp_ln$timepoint <- rep(timepoint, dim(temp_ln)[1])
#     temp_ln$model <- rep('LN', dim(temp_ln)[1])
    
    temp_ln <- temp_ln[which(rownames(temp_ln) %in% rownames(deg_geneTable)),]  ####  only in DEG used genes
    
    new_ln_list[[timepoint]] <- temp_ln
    
    print(dim(temp_ln))

}

saveRDS(new_ln_list, 'qa22___df_ln_list__DEG_genes.rds')


