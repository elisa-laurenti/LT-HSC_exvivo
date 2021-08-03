#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(DESeq2))

result_list <- readRDS('deseq2_results_noPD.rds') ### deseq2 results between all timepoints

geneTable <- readRDS('genetable.rds',row.names=1)

fdr<-0.05
gene_list <- c()

for (i in result_list){
    i_Ordered <- i[!is.na(i$padj),]
    filtered <- i_Ordered[i_Ordered$padj < fdr,]
    gene_list <- c(gene_list, rownames(filtered))
}
unique_genes <- unique(gene_list)



rownames(geneTable) <- geneTable$Ensembl_Gene_ID

unique_table <- geneTable[unique_genes,]


saveRDS(unique_table,'DE_unique_genes.rds')






