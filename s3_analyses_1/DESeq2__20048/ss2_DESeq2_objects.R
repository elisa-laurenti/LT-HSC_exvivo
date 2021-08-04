suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))


###################################
#
###################################


meta <- read.csv(paste0(root_dir,'meta_20048.csv'),row.names=1)

counts <- fread(paste0(root_dir,'counts_20048.txt'))

var_df <- read.csv(paste0(root_dir,'adata_var_20048.csv'),row.names=1)

counts_df = setDF(counts)


colnames(counts_df) <- rownames(meta)

rownames(counts_df) <- rownames(var_df)


###################################
#
###################################

express_3cells <- rowSums(counts_df!=0)>=5

count_all_filtered <- counts_df[express_3cells,]
dim(count_all_filtered)


###################################
#
###################################

meta_sub <- meta[,c('time','condition_fixed')]

meta_sub$within_control <- as.factor(paste(meta_sub$time, meta_sub$condition_fixed, sep="_"))

colnames(meta_sub) <- c('time','condition','time_condition')

###################################
#
###################################


dds <- DESeq2::DESeqDataSetFromMatrix(count_all_filtered, meta_sub, design = ~ 0 + time_condition)


saveRDS(dds, 'dds_20048.rds')

dds <- DESeq2::DESeq(dds, parallel=TRUE)

save.image('ss2_DESeq2_objects.RData')

