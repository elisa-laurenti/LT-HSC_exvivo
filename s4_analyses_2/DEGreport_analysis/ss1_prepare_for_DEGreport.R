#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS('dds_noPD.rds')  ### deseq results from QA22 and qa81 combined with no PD cells

result_list <- readRDS('deseq2_results_noPD.rds') ### deseq2 results between all timepoints

######################
#
#####################


rld <- vst(dds,blind=FALSE)

ma <- assay(rld)


######################
#
#####################

ma <- limma::removeBatchEffect(ma,rld$batch)

assay(rld) <- ma


######################
#
#####################

fdr<-0.05
gene_list <- c()

for (i in result_list){
    i_Ordered <- i[!is.na(i$padj),]
    filtered <- i_Ordered[i_Ordered$padj < fdr,]
    gene_list <- c(gene_list, rownames(filtered))
}
unique_genes <- unique(gene_list)


######################
#
#####################

ma_padj05 <- ma[which(rownames(ma)%in%unique_genes),]

design <- as.data.frame(colData(dds))


saveRDS(ma_padj05,'ma_padj05.rds')

saveRDS(design,'design.rds')



