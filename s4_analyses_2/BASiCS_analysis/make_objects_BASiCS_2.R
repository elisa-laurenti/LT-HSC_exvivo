suppressPackageStartupMessages(library(BASiCS))

counts <- readRDS('CountsFilter_8464_429.rds')

meta <- read.csv('meta_blobs_for_total_counts.csv', row.names=1)


identical(rownames(meta), colnames(counts))   ### TRUE

dim(counts) ## 8464  429

meta$BatchInfo <- meta$Cell_type_subtype


###############################
#  prepare for SingleCellExperiment, subsetting  meta table
###############################

m_0 <- meta[meta$Details == 'LT_0h',]
dim(m_0)

table(m_0$Details)


m_24 <- meta[meta$Details == 'LT_24h_UNTR',]

dim(m_24)

table(m_24$Details)


m_6 <- meta[meta$Details == 'LT_6h',]

dim(m_6)

table(m_6$Details)


m_72 <- meta[meta$Details == 'LT_72h_UNTR',]

dim(m_72)

table(m_72$Details)

######################################

m_0_A <- meta[meta$new_details == '0h_A',]

dim(m_0_A)

m_0_B <- meta[meta$new_details == '0h_B',]

dim(m_0_B)

m_0_11 <- meta[meta$new_details == 'LT_0h',]

dim(m_0_11)



###############################
#  prepare for SingleCellExperiment, subsetting  counts
###############################

c_0 <- counts[,rownames(m_0)]

dim(c_0)

c_6 <- counts[,rownames(m_6)]

dim(c_6)

c_24 <- counts[,rownames(m_24)]

dim(c_24)

c_72 <- counts[,rownames(m_72)]

dim(c_72)


######################################

c_0_A <- counts[,rownames(m_0_A)]

dim(c_0_A)

c_0_B <- counts[,rownames(m_0_B)]

dim(c_0_B)

c_0_11 <- counts[,rownames(m_0_11)]

dim(c_0_11)


###############################
#  making SingleCellExperiment object
###############################

Data_0 <- SingleCellExperiment(assays = list(counts = c_0), colData = m_0)

Data_6 <- SingleCellExperiment(assays = list(counts = c_6), colData = m_6)

Data_24 <- SingleCellExperiment(assays = list(counts = c_24), colData = m_24)

Data_72 <- SingleCellExperiment(assays = list(counts = c_72), colData = m_72)

Data_0_A <- SingleCellExperiment(assays = list(counts = c_0_A), colData = m_0_A)

Data_0_B <- SingleCellExperiment(assays = list(counts = c_0_B), colData = m_0_B)

Data_0_11 <- SingleCellExperiment(assays = list(counts = c_0_11), colData = m_0_11)


saveRDS(Data_0, 'Data_0__8464_85.rds')

saveRDS(Data_6, 'Data_6__8464_134.rds')

saveRDS(Data_24, 'Data_24__8464_86.rds')

saveRDS(Data_72, 'Data_72__8464_124.rds')

saveRDS(Data_0_A,'Data_0_A__8464_38.rds')

saveRDS(Data_0_B,'Data_0_B__8464_36.rds')

saveRDS(Data_0_11,'Data_0_11__8464_11.rds')

