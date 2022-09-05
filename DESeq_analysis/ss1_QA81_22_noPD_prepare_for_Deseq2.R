suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))


meta_22 <- read.csv('QA22_meta_after_qc.csv',row.names=1)
counts_22 <- read.csv('QA22_raw_counts_after_qc.csv',row.names=1)

meta_81 <- read.csv('QA81_meta_after_qc.csv',row.names=1)
counts_81 <- read.csv('QA81_raw_counts_after_qc.csv',row.names=1)


meta_all <- rbind(meta_22,meta_81)

meta_all <- subset(meta_all , Details !='LT_72h_PD' & Details !='LT_24h_PD')

counts_combined <- rbind(counts_22,counts_81)

count_all <- t(counts_combined)

identical(colnames(count_all),rownames(meta_all))



###############################
# filtering genes with expression in more than 3 cells
###############################

express_3cells <- rowSums(count_all!=0)>=3

table(express_3cells)

count_all_filtered <- count_all[express_3cells,]
dim(count_all_filtered)


###############################
# rearranging
###############################

meta_all_2col <- meta_all[,c('Details','Cell_type_subtype')]

head(meta_all_2col)

colnames(meta_all_2col) <- c('timepoint','batch')

identical(colnames(count_all),rownames(meta_all_2col))

dds <- DESeq2::DESeqDataSetFromMatrix(count_all_filtered, meta_all_2col, design = ~ batch + timepoint)


dds$timepoint <- factor(dds$timepoint, levels = c('LT_0h','LT_6h', 'LT_24h_UNTR','LT_72h_UNTR' ))


saveRDS(dds, 'dds_noPD.rds')


