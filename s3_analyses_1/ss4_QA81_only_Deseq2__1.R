suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))

scd453 <- readRDS('scd453_scran_bglab.rds')

meta453 <- pData(scd453)

meta453 <- subset(meta453 , Details !='LT_72h_PD')
meta453 <- droplevels(meta453)
head(meta453)

count_qa81 <- scd453@counts    ### 65988  570

dim(count_qa81)

count_qa81_no72PD <- count_qa81[,rownames(meta453)]

identical(colnames(count_qa81_no72PD), rownames(meta453))

######################### 
# filtering
#########################

express_3cells <- rowSums(count_qa81_no72PD!=0)>=3
table(express_3cells)

count_all_filtered <- count_qa81_no72PD[express_3cells,]
dim(count_all_filtered)

meta453_2col <- meta453[,c('Details','Cell_type_subtype')]

colnames(meta453_2col) <- c('timepoint','batch')

identical(colnames(count_all_filtered),rownames(meta453_2col))

dds <- DESeq2::DESeqDataSetFromMatrix(count_all_filtered, meta453_2col, design = ~ timepoint)

dds$timepoint <- factor(dds$timepoint, levels = c('LT_0h','LT_6h', 'LT_24h_PD', 'LT_24h_UNTR','LT_72h_UNTR' ))

saveRDS(dds, 'dds_qa81_no72PD.rds')


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






