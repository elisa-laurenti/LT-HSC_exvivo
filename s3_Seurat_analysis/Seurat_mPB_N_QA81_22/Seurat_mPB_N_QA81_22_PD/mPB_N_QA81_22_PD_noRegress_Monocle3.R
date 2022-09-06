suppressMessages(library(monocle3))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(htmlwidgets))

seurat.combined.sct = readRDS('seurat.combined.QA81_22_PD_N_mPB_pca_umap.rds')

cds <- as.cell_data_set(seurat.combined.sct)

data_cds = cluster_cells(cds = cds, reduction_method = "UMAP")

data_cds <- learn_graph(data_cds, use_partition = FALSE)

data_cds <- order_cells(data_cds, 
                        reduction_method = "UMAP", 
                        root_cells = c('SLX.16064.i701_i506-0'))

saveRDS(data_cds,'PD_QA81_22_N_mPB__noRegress_pseudotime_monocle3.rds')

#####################################

cds_2D = readRDS('PD_QA81_22_N_mPB__noRegress_pseudotime_monocle3.rds')

pseudotime_2d = cds_2D@principal_graph_aux@listData$UMAP$pseudotime
identical(rownames(seurat.combined.sct@meta.data), names(pseudotime_2d))

seurat.combined.sct <- AddMetaData(
  object = seurat.combined.sct,
  metadata = cds_2D@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime_2D"
)


meta = seurat.combined.sct@meta.data

meta$pseudotime_2D_rank = rank(meta$pseudotime_2D)

write.csv(meta,'meta_PD_QA81_22_N_mPB__noRegress_pseudotime_monocle3.csv')
