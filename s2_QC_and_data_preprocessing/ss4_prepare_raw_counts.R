
library(bglab)

scd_qc <- readRDS('scd_qc.rds')
scd453 <- readRDS('scd453_scran_bglab.rds')

meta_qc <- pData(scd_qc)
meta453 <- pData(scd453)

temp_meta <- meta453[which(meta453$Details != 'LT_72h_PD'),]  # 361  30
meta361 <- temp_meta

# meta_all <- read.csv('meta_noOutlier_555cells.csv', row.names = 1)

# make sure they have the same columns and then combine them
meta_qc <- meta_qc[,which(colnames(meta_qc)%in%colnames(meta361))]
meta361 <- meta361[,which(colnames(meta361)%in%colnames(meta_qc))]
meta_all <- rbind(meta_qc,meta361)
meta_all <- droplevels(meta_all)

meta_all$Details <- mapvalues(meta_all$Details, 
                              from = c('LT-HSCs_T0','LT-HSCs_T24','LT-HSCs_T48','LT-HSCs_T6','LT-HSCs_T72_UN','LT-HSCs_T72_PD'),
                              to=c('LT_0h','LT_24h_UNTR','LT_48h','LT_6h','LT_72h_UNTR','LT_72h_PD'))



count_qa22 <- scd_qc@counts      ### 65988  576
count_qa81 <- scd453@counts    ### 65988  570

count_qa22_195 <- count_qa22[,which(colnames(count_qa22) %in% rownames(meta_qa22))]
dim(count_qa22_195)

count_qa81_361 <- count_qa81[,which(colnames(count_qa81) %in% rownames(meta_qa81))]
dim(count_qa81_361)

identical(colnames(count_qa22_195), rownames(meta_qa22))

identical(colnames(count_qa81_361), rownames(meta_qa81))

t_count_qa22_195 <- t(count_qa22_195)
t_count_qa81_361 <- t(count_qa81_361)

identical(rownames(t_count_qa22_195), rownames(meta_qa22))
identical(rownames(t_count_qa81_361), rownames(meta_qa81))

write.csv(t_count_qa81_361, 't_count_qa81_361.csv')
write.csv(t_count_qa22_195, 't_count_qa22_195.csv')


