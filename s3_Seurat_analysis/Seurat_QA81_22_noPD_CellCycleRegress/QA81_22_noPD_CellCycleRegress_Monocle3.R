
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(htmlwidgets)


seurat.combined.sct = readRDS('seurat.combined.QA81_22_pca_umap__CellCycleRegress.rds')

seurat.combined.sct

DimPlot(seurat.combined.sct, reduction = "umap", group.by = "Details")

##########
#  choose root cell
##########

umap = as.data.frame(seurat.combined.sct[["umap"]]@cell.embeddings)
umap$number = seq(dim(umap)[1])
umap[which(umap$UMAP_2 == max(umap$UMAP_2)),]

###############
#
###############

cds <- as.cell_data_set(seurat.combined.sct)

data_cds = cluster_cells(cds = cds, reduction_method = "UMAP")

data_cds <- learn_graph(data_cds, use_partition = TRUE)

data_cds <- order_cells(data_cds, 
                        reduction_method = "UMAP", 
                        root_cells = c('SLX.16064.i705_i502-0'))

plot_cells(
  cds = data_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)+ggtitle('QA81_22 noPD cycleRegress')



saveRDS(data_cds,'QA81_22_noPD__cycleRegress_monocle3_pseudotime.rds')


