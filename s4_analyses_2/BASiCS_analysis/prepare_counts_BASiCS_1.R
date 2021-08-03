suppressPackageStartupMessages(library(bglab))


meta <- read.csv('meta_blobs_for_total_counts.csv', row.names=1)

counts_22 <- read.csv('QA22_raw_counts_after_qc.csv',row.names=1)

counts_81 <- read.csv('QA81_raw_counts_after_qc.csv',row.names=1)



counts_combined <- rbind(counts_22,counts_81)

 
count_all <- t(counts_combined)

count_429 <- count_all[,rownames(meta)]

identical(rownames(meta), colnames(count_429))


 
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

