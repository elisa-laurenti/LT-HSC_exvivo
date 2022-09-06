suppressPackageStartupMessages(library(bglab))
suppressPackageStartupMessages(library(DESeq2))

meta <- read.csv('meta_QA81_22_noPD.csv', row.names=1)

counts_22 <- read.csv('QA22_raw_counts_after_qc.csv',row.names=1)
counts_81 <- read.csv('QA81_raw_counts_after_qc.csv',row.names=1)

counts <- rbind(counts_22,counts_81)

count_429 <- counts[,rownames(meta)]

############################################################

CountsRPM=1000000*count_429/colSums(count_429)

Include =which(rowMeans(CountsRPM)>20)

CountsFilter <- count_429[ Include, ]


identical(rownames(meta), colnames(CountsFilter))


saveRDS(CountsFilter, 'CountsFilter_QA81_22_noPD.rds')

