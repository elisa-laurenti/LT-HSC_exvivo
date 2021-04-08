suppressPackageStartupMessages(library(bglab))


meta <- read.csv('meta_blobs_for_total_counts.csv', row.names=1)


scd_qc <- readRDS('scd_qc.rds')
scd453 <- readRDS('scd453_scran_bglab.rds')


count_qa22 <- scd_qc@counts      ### 65988  576
count_qa81 <- scd453@counts    ### 65988  570


identical(rownames(count_qa22),rownames(count_qa81))

counts<- cbind(count_qa22,count_qa81)

dim(counts)


count_429 <- counts[,rownames(meta)]

identical(rownames(meta), colnames(count_429))


dim(count_429)  ### 65988   429


######################################
#  filtering counts, using the filter === More than 20 RPM (on average), across all cells
###################################### 


CountsRPM=1000000*count_429/colSums(count_429)

dim(CountsRPM)

Include =which(rowMeans(CountsRPM)>20)

length(Include)


CountsFilter <- count_429[ Include, ]

dim(CountsFilter)  ### 8464  429

identical(rownames(meta), colnames(CountsFilter))


#####################################
# save object
#####################################
saveRDS(CountsFilter, 'CountsFilter_8464_429.rds')

