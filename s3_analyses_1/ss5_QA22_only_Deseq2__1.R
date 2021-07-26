suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))


meta_22 <- read.csv('QA22_meta_after_qc.csv',row.names=1)
counts_22 <- read.csv('QA22_raw_counts_after_qc.csv',row.names=1)


 
count_175 <- t(counts_22)

identical(colnames(count_175),rownames(meta_22))
 
 

######################### 
# filtering
#########################

express_3cells <- rowSums(count_175!=0)>=3

table(express_3cells)

count_all_filtered <- count_175[express_3cells,]
dim(count_all_filtered)

meta_sub_2col <- meta_22[,c('Details','Cell_type_subtype')]

colnames(meta_sub_2col) <- c('timepoint','batch')

identical(colnames(count_all_filtered),rownames(meta_sub_2col))


######################### 
# dds
#########################
dds <- DESeq2::DESeqDataSetFromMatrix(count_all_filtered, meta_sub_2col, design = ~ timepoint)


dds$timepoint <- factor(dds$timepoint, 
                        levels = c('LT_0h','LT_6h','LT_24h_UNTR','LT_72h_PD','LT_72h_UNTR'))

saveRDS(dds, 'dds_qa22_175cells_no48.rds')

dds <- DESeq2::DESeq(dds, parallel=TRUE)


######################### 
# extract result
#########################
resultsNames(dds)

LT_72h_PD_vs_LT_0h <- results(dds,
                             name='timepoint_LT_72h_PD_vs_LT_0h')

LT_72h_UNTR_vs_LT_0h <- results(dds,
                               name='timepoint_LT_72h_UNTR_vs_LT_0h')

res_72pd_72untr <- results(dds, alpha=0.05, parallel=TRUE, contrast=c('timepoint','LT_72h_PD','LT_72h_UNTR'))

