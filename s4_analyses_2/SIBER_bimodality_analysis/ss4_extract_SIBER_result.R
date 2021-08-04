#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(SIBERG))


ln_22 <- readRDS('qa22___df_ln_list__DEG_genes.rds')

ln_81 <- readRDS('qa81___df_ln_list__DEG_genes.rds')


qa22_zero_inflated <- readRDS('qa22__results_zero_inflated_checks__DEG_genes.rds')

qa81_zero_inflated <- readRDS('qa81__results_zero_inflated_checks__DEG_genes.rds')

for (timepoint in names(ln_22)){
    
    ln_22[[timepoint]]$timepoint <- rep(timepoint, dim(ln_22[[timepoint]])[1])
    ln_22[[timepoint]]$model <- rep('LN', dim(ln_22[[timepoint]])[1])
    ln_22[[timepoint]]$batch <- rep('QA22', dim(ln_22[[timepoint]])[1])
    
    print(identical(rownames(ln_22[[timepoint]]), rownames(qa22_zero_inflated[[timepoint]])))
    
    ln_22[[timepoint]]$zero_inflated <- qa22_zero_inflated[[timepoint]]$zero_inflated_list
    ln_22[[timepoint]]$percentage_0 <- qa22_zero_inflated[[timepoint]]$percent0_list
    
    
#     print(dim(ln_22[[timepoint]]))
#     print(dim(qa22_zero_inflated[[timepoint]]))
}

for (timepoint in names(ln_81)){
    
    ln_81[[timepoint]]$timepoint <- rep(timepoint, dim(ln_81[[timepoint]])[1])
    ln_81[[timepoint]]$model <- rep('LN', dim(ln_81[[timepoint]])[1])
    ln_81[[timepoint]]$batch <- rep('QA81', dim(ln_81[[timepoint]])[1])
    
    print(identical(rownames(ln_81[[timepoint]]), rownames(qa81_zero_inflated[[timepoint]])))
    
    ln_81[[timepoint]]$zero_inflated <- qa81_zero_inflated[[timepoint]]$zero_inflated_list
    ln_81[[timepoint]]$percentage_0 <- qa81_zero_inflated[[timepoint]]$percent0_list
    
    
    
}


####################
#
####################

geneTable <- readRDS('geneTable.rds')

rownames(geneTable) <- geneTable$Ensembl_Gene_ID


for (timepoint in names(ln_22)){
    temp_genetable <- geneTable[rownames(ln_22[[timepoint]]),]
    
    print(identical(rownames(temp_genetable),rownames(ln_22[[timepoint]])))
    
    ln_22[[timepoint]]$Ensembl_Gene_ID <- rownames(ln_22[[timepoint]])
    ln_22[[timepoint]]$geneName <- temp_genetable$Associated_Gene_Name
    
}

for (timepoint in names(ln_81)){
        temp_genetable <- geneTable[rownames(ln_81[[timepoint]]),]
    
    print(identical(rownames(temp_genetable),rownames(ln_81[[timepoint]])))
    
    ln_81[[timepoint]]$Ensembl_Gene_ID <- rownames(ln_81[[timepoint]])
    ln_81[[timepoint]]$geneName <- temp_genetable$Associated_Gene_Name
    
}


dataframe_22 <- rbind(ln_22[['LT_0h']],
                     ln_22[['LT_6h']],
                     ln_22[['LT_24h_UNTR']],
                     ln_22[['LT_72h_PD']],
                     ln_22[['LT_72h_UNTR']])


dataframe_81 <- rbind(ln_81[['LT_0h']],
                     ln_81[['LT_6h']],
                     ln_81[['LT_24h_PD']],
                     ln_81[['LT_24h_UNTR']],
                     ln_81[['LT_72h_UNTR']])



not_zero_22 <- dataframe_22[(dataframe_22$zero_inflated == 'not_zero_inflated') ,]

not_zero_81 <- dataframe_81[(dataframe_81$zero_inflated == 'not_zero_inflated') ,]


no_pd_not_zero_22 <- not_zero_22[which(not_zero_22$timepoint != 'LT_72h_PD'),]

no_pd_not_zero_81 <- not_zero_81[which(not_zero_81$timepoint != 'LT_24h_PD'),]


combined_df <- rbind(no_pd_not_zero_22, no_pd_not_zero_81)


write.csv(combined_df,'QA81_22__BI_values_not_zero_inflated_noPD__dataframe.csv')



   