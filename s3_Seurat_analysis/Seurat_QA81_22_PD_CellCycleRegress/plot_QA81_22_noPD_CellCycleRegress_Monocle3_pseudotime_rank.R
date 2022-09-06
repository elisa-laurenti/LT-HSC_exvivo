library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(htmlwidgets)


########
#
######


cds_2D = readRDS('QA81_22_PD__cycleRegress_monocle3_pseudotime.rds')

length(cds_2D@principal_graph_aux@listData$UMAP$pseudotime)

pseudotime_2d = cds_2D@principal_graph_aux@listData$UMAP$pseudotime

seurat.combined.sct <- AddMetaData(
  object = seurat.combined.sct,
  metadata = cds_2D@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime_2D"
)

meta = seurat.combined.sct@meta.data

ggplot(meta, aes(x=pseudotime_2D, color=Details)) +
  geom_density()

summary(meta$pseudotime_2D)

meta$pseudotime_2D_rank = rank(meta$pseudotime_2D)

ggplot(meta, aes(x=pseudotime_2D_rank, color=Details)) +
  geom_density()+ggtitle('QA81_22 PD CycleRegress pseudotime 2D')
