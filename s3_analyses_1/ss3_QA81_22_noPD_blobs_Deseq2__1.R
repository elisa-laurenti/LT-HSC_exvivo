suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))


scd_qc <- readRDS('scd_qc.rds')
scd453 <- readRDS('scd453_scran_bglab.rds')

meta_qc <- pData(scd_qc)
meta453 <- pData(scd453)


# make sure they have the same columns and then combine them
meta_qc <- meta_qc[,which(colnames(meta_qc)%in%colnames(meta453))]
meta453 <- meta453[,which(colnames(meta453)%in%colnames(meta_qc))]
meta_all <- rbind(meta_qc,meta453)
meta_all <- droplevels(meta_all)

meta_all$Details <- mapvalues(meta_all$Details, 
                              from = c('LT-HSCs_T0','LT-HSCs_T24','LT-HSCs_T48','LT-HSCs_T6','LT-HSCs_T72_UN','LT-HSCs_T72_PD'),
                              to=c('LT_0h','LT_24h_UNTR','LT_48h','LT_6h','LT_72h_UNTR','LT_72h_PD'))

meta_all <- subset(meta_all , Details !='LT_72h_PD' & Details !='LT_48h' & Details !='LT_24h_PD')
meta_all <- droplevels(meta_all)
head(meta_all,3)

meta_qa22 <- meta_all[which(meta_all$Cell_type_subtype == 'QA22'),]
meta_qa81 <- meta_all[which(meta_all$Cell_type_subtype == 'QA81'),]


count_qa22 <- scd_qc@counts      ### 65988  576
count_qa81 <- scd453@counts    ### 65988  570
count_qa22_195 <- count_qa22[,which(colnames(count_qa22) %in% rownames(meta_qa22))]
dim(count_qa22_195)
count_qa81_361 <- count_qa81[,which(colnames(count_qa81) %in% rownames(meta_qa81))]
dim(count_qa81_361)
identical(colnames(count_qa22_195), rownames(meta_qa22))

identical(colnames(count_qa81_361), rownames(meta_qa81))


count_all<- cbind(count_qa22_195,count_qa81_361)
dim(count_all)

identical(colnames(count_all),rownames(meta_all))


#########################
#  load in meta table
#########################

first_blob = read.csv('first_blob.csv',row.names=1)
second_blob = read.csv('second_blob.csv', row.names=1)


meta_all$Detail_blob <- as.character(meta_all$Details)

meta_all$Detail_blob[which(rownames(meta_all)%in% rownames(first_blob))] <- '0h_A'
meta_all$Detail_blob[which(rownames(meta_all)%in% rownames(second_blob))] <- '0h_B'


##########################
#  filtering counts for Deseq
##########################

### lets set it to half of the smaller blob cluster
nMinCells = 18

numCellsGeneExpr <- apply(count_all, 1, function(x) sum(x>0))

length(numCellsGeneExpr)

gene <- names(numCellsGeneExpr)[numCellsGeneExpr > nMinCells]

length(gene)

count_all = count_all[gene,]
dim(count_all)


##########################
# make meta table
##########################

meta_all_2col <- meta_all[,c('Detail_blob','Cell_type_subtype')]

head(meta_all_2col,2)

colnames(meta_all_2col) <- c('timepoint','batch')

identical(colnames(count_all),rownames(meta_all_2col))

dim(count_all)

########################## 
#  make dds object
##########################

dds <- DESeq2::DESeqDataSetFromMatrix(count_all, meta_all_2col, design = ~ batch + timepoint)

dds

dds$timepoint <- factor(dds$timepoint, levels = c('0h_A',
                                                  '0h_B',
                                                  'LT_0h',
                                                  'LT_6h',
                                                  'LT_24h_UNTR',
                                                  'LT_72h_UNTR' ))

saveRDS(dds, 'dds_noPD_2blobs_new_filter.rds')

###########################
# doing deseq2 
###########################

dds <- DESeq2::DESeq(dds, parallel=TRUE)


###########################
# extract result
###########################

resultsNames(dds)

blobs <- results(dds, name='timepoint_0h_B_vs_0h_A')


###########################
#  load in geneTable
###########################
geneTable<-readRDS('geneTable.rds')
rownames(geneTable) <- geneTable$Ensembl_Gene_ID


###########################
# make dataframe
###########################
df_2nd_blob_vs_1st_blob <- as.data.frame(blobs)

table_2_1 <- geneTable[rownames(df_2nd_blob_vs_1st_blob),]

identical(rownames(table_2_1),rownames(df_2nd_blob_vs_1st_blob))

df_2nd_blob_vs_1st_blob$geneName <- table_2_1$Associated_Gene_Name


###########################
#  filtering
###########################

log2_df <- df_2nd_blob_vs_1st_blob[which(df_2nd_blob_vs_1st_blob$log2FoldChange >2 | df_2nd_blob_vs_1st_blob$log2FoldChange < -2 ),]

dim(log2_df)


basemean_df <- log2_df[which(log2_df$baseMean >= 10),]


write.csv(basemean_df,'df_blobs_logfold2_basemean10_new_filter.csv')
