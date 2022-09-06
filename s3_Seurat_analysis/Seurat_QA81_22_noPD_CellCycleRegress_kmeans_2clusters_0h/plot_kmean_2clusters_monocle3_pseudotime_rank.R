suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))

pseudot_meta = read.csv('meta_noPD_cycle_regress_pseudotime_UMAP.csv'),row.names=1)

blob_meta = read.csv('meta_noPD_cycle_regress_pseudotime_UMAP_KMeans_0h_2clusters_for_deseq.csv',row.names=1)

identical(rownames(pseudot_meta),rownames(blob_meta))

pseudot_meta$kmean = blob_meta$kmean

pseudot_meta$pseudotime_2D_rank = rank(pseudot_meta$pseudotime_2D)

only_kmean_clusters_meta = pseudot_meta[which(pseudot_meta$kmean %in% c(0,1)),]
only_kmean_clusters_meta$kmean = mapvalues(only_kmean_clusters_meta$kmean, c(0, 1), c("LT_0h_B", "LT_0h_A"))

ggplot(only_kmean_clusters_meta, aes(x=pseudotime_2D_rank, color=kmean)) +
  geom_density()+ggtitle('QA81_22 noPD CycleRegress pseudotime 2D t_0h kmean2 clusters')+
theme_bw()


