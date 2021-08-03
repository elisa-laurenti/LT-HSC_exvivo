#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(SIBERG))

##############################
#  loading in data
##############################

print('load in data')

meta_sub <- read.csv('QA22_meta_after_qc.csv',row.names=1)

count_175 <- read.csv('QA22_raw_counts_after_qc.csv',row.names=1)
 
 

#############################
#  filtering gene expression
##############################

print('filtering gene expression')

express_3cells <- rowSums(count_175!=0)>=3
print(table(express_3cells))

count_all_filtered <- count_175[express_3cells,]

print(dim(count_all_filtered))

print(colnames(count_all_filtered)[1:5])

print(rownames(count_all_filtered)[1:5])


#############################
#  estimating size factor
##############################

TMM <- calcNormFactors(count_all_filtered, method='TMM')

print(length(TMM))

#############################
#  
##############################

meta_sub <- droplevels(meta_sub)

table(meta_sub$Details)

meta_sub$Details <- factor(meta_sub$Details, levels = c('LT_0h','LT_6h','LT_24h_UNTR','LT_72h_PD','LT_72h_UNTR'))


unique(meta_sub$Details)


#############################
#  
##############################

base=exp(1)
eps=10
size_factor = 1/TMM

results_zero_inflated <- list()

for (timepoint in unique(meta_sub$Details)){
    
    temp_meta <- meta_sub[meta_sub$Details == timepoint,]
    
    temp_count <- count_all_filtered[,rownames(temp_meta)]
    
    temp_tmm <- size_factor[rownames(temp_meta)]
    
    
    gene_list <- c()
    percent0_list <- c()
    zero_inflated_list <- c()
    
    for (gene in rownames(temp_count)){  ### the original dataframe
    
        temp_gene_vals = temp_count[gene,]
        
        percent0 <- mean(temp_gene_vals==0, na.rm=TRUE)
        
        if (percent0> 0.2){
            zero_inflated <- 'zero_inflated'
        }else{
            zero_inflated <- 'not_zero_inflated'
        }

        gene_list <- c(gene_list, gene)
        
        percent0_list <- c(percent0_list, percent0)
        
        zero_inflated_list <- c(zero_inflated_list, zero_inflated)
        
    }
    
    temp_df <- data.frame(gene_list,percent0_list,zero_inflated_list)
    
    
    results_zero_inflated[[timepoint]] <- temp_df
    
    print(timepoint)
    print(dim(temp_meta))
    print(dim(temp_count))
}




#############################
#  
##############################

for(timepoint in names(results_zero_inflated)){
    temp_df <- results_zero_inflated[[timepoint]]
    
    sub_inflated <- temp_df[(temp_df$zero_inflated_list == 'zero_inflated'),]
    
    print(timepoint)
    print(dim(temp_df))
    print(dim(sub_inflated))
}



saveRDS(results_zero_inflated,'qa22__results_zero_inflated_checks.rds')


#############################
#  
##############################

de_genes <- readRDS('DE_unique_genes.rds')
 

new_ln_list <- list()

for ( timepoint in names(results_zero_inflated)){
    
    temp_ln <- results_zero_inflated[[timepoint]]
    rownames(temp_ln) <- temp_ln$gene_list
    
#     temp_ln$timepoint <- rep(timepoint, dim(temp_ln)[1])
#     temp_ln$model <- rep('LN', dim(temp_ln)[1])
    
    temp_ln <- temp_ln[which(rownames(temp_ln) %in% rownames(de_genes)),]  ####  only in DEG used genes
    
    new_ln_list[[timepoint]] <- temp_ln
    
    print(dim(temp_ln))

}

saveRDS(new_ln_list,'qa22__results_zero_inflated_checks__DEG_genes.rds')

