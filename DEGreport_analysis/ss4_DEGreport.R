#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))


load('ss3_DEGreport.RData')


sub_0 <- table[(table$timepoint == 'LT_0h'),]

sub_6 <- table[(table$timepoint == 'LT_6h'),]

sub_24 <- table[(table$timepoint == 'LT_24h_UNTR'),]

sub_72 <- table[(table$timepoint == 'LT_72h_UNTR'),]



gene_for_dataframe <- c()
cluster_for_dataframe <- c()

for (gene in sub_0$genes){
    
    gene_for_dataframe <- c(gene_for_dataframe, gene)
    
    cluster <- sub_0[(sub_0$genes == gene),'cluster']
    
    if(cluster != 2 && cluster !=7){
        
        cluster_for_dataframe <- c(cluster_for_dataframe, cluster)
        
    }
    if (cluster == 2){
        temp_sub_0 <- sub_0[(sub_0$genes == gene),'value']
        temp_sub_6 <- sub_6[(sub_6$genes == gene),'value']

        if(temp_sub_0 < temp_sub_6 ){
            
            cluster_for_dataframe <- c(cluster_for_dataframe,'2a')
            
        } else if (temp_sub_0 > temp_sub_6){
            cluster_for_dataframe <- c(cluster_for_dataframe,'2b')
        }
    }
    
    if (cluster == 7){
        temp_sub_0 <- sub_0[(sub_0$genes == gene),'value']
        temp_sub_6 <- sub_6[(sub_6$genes == gene),'value']

        if(temp_sub_0 < temp_sub_6 ){
            
            cluster_for_dataframe <- c(cluster_for_dataframe,'7a')
            
        } else if (temp_sub_0 > temp_sub_6){
            
            cluster_for_dataframe <- c(cluster_for_dataframe,'7b')
        }
    }

    
}



sub_0$old_cluster <- sub_0$cluster

sub_6$old_cluster <- sub_6$cluster

sub_24$old_cluster <- sub_24$cluster

sub_72$old_cluster <- sub_72$cluster


sub_0$cluster <- cluster_for_dataframe
sub_6$cluster <- cluster_for_dataframe

sub_24$cluster <- cluster_for_dataframe
sub_72$cluster <- cluster_for_dataframe



dictionary <-  table(sub_0$cluster)


new_group <- c()

for ( row in rownames(sub_0)){
    
    cluster <- sub_0[row,'cluster']
#     old_title <- sub_0[row,'title']
    
    num_gene <- dictionary[[cluster]]
    
    new_title <- paste0('Group:  ',cluster,' - genes: ',num_gene)
    
    new_group <- c(new_group, new_title)


}


sub_0$old_title <- sub_0$title
sub_6$old_title <- sub_6$title
sub_24$old_title <- sub_24$title
sub_72$old_title <- sub_72$title


sub_0$title <- new_group
sub_6$title <- new_group
sub_24$title <- new_group
sub_72$title <- new_group

combined_new_group <- rbind(sub_0,sub_6,sub_24,sub_72)


saveRDS(combined_new_group, 'DEGreport_result_table.rds')









