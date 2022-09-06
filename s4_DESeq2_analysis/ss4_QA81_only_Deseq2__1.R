suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))


meta_81 <- read.csv('QA81_meta_after_qc.csv',row.names=1)
counts_81 <- read.csv('QA81_raw_counts_after_qc.csv',row.names=1)

 
count_all <- t(counts_81)

identical(colnames(count_all),rownames(meta_81))

  
 
######################### 
# filtering
#########################

express_3cells <- rowSums(count_all!=0)>=3
table(express_3cells)

count_all_filtered <- count_all[express_3cells,]
dim(count_all_filtered)

meta_81_2col <- meta_81[,c('Details','Cell_type_subtype')]

colnames(meta_81_2col) <- c('timepoint','batch')

identical(colnames(count_all_filtered),rownames(meta_81_2col))

dds <- DESeq2::DESeqDataSetFromMatrix(count_all_filtered, meta_81_2col, design = ~ timepoint)

dds$timepoint <- factor(dds$timepoint, levels = c('LT_0h','LT_6h', 'LT_24h_PD', 'LT_24h_UNTR','LT_72h_UNTR' ))

saveRDS(dds, 'dds_qa81.rds')


######################### 
# doing dds
#########################

dds <- DESeq2::DESeq(dds, parallel=TRUE)

######################### 
# extract result
#########################

resultsNames(dds)

LT_24h_PD_vs_LT_0h <- results(dds,
                              name= 'timepoint_LT_24h_PD_vs_LT_0h')

LT_24h_UNTR_vs_LT_0h <- results(dds,
                              name='timepoint_LT_24h_UNTR_vs_LT_0h' )

res_24pd_24untr <- results(dds, alpha=0.05, parallel=TRUE, contrast=c('timepoint','LT_24h_PD','LT_24h_UNTR'))






