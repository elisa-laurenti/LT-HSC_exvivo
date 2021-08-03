#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS('dds_noPD.rds')  ### deseq results from QA22 and qa81 combined with no PD cells

unique_genes <- readRDS('DE_unique_genes.rds') 

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

ma_padj05 <- ma[which(rownames(ma)%in%unique_genes),]

design <- as.data.frame(colData(dds))


saveRDS(ma_padj05,'ma_padj05.rds')

saveRDS(design,'design.rds')



