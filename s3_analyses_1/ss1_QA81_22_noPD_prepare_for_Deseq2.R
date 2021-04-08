suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))

scd_qc <- readRDS('scd_qc.rds')
scd453 <- readRDS('scd453_scran_bglab.rds')

meta_qc <- pData(scd_qc)
meta453 <- pData(scd453)


meta_qc <- meta_qc[,which(colnames(meta_qc)%in%colnames(meta453))]
meta453 <- meta453[,which(colnames(meta453)%in%colnames(meta_qc))]
meta_all <- rbind(meta_qc,meta453)
meta_all <- droplevels(meta_all)

meta_all$Details <- mapvalues(meta_all$Details, 
                              from = c('LT-HSCs_T0','LT-HSCs_T24','LT-HSCs_T48','LT-HSCs_T6','LT-HSCs_T72_UN','LT-HSCs_T72_PD'),
                              to=c('LT_0h','LT_24h_UNTR','LT_48h','LT_6h','LT_72h_UNTR','LT_72h_PD'))



meta_all <- subset(meta_all , Details !='LT_72h_PD' & Details !='LT_48h' & Details !='LT_24h_PD')
meta_all <- droplevels(meta_all)
head(meta_all)



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


