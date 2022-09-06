suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))

cds_2D = readRDS('QA81_22_noPD__cycleRegress_monocle3_pseudotime.rds')

length(cds_2D@principal_graph_aux@listData$UMAP$pseudotime)

seurat.combined.sct = readRDS('seurat.combined.QA81_22_pca_umap__CellCycleRegress.rds')


pseudotime_2d = cds_2D@principal_graph_aux@listData$UMAP$pseudotime

seurat.combined.sct <- AddMetaData(
  object = seurat.combined.sct,
  metadata = cds_2D@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime_2D"
)


meta = seurat.combined.sct@meta.data
umap = as.data.frame(seurat.combined.sct[["umap"]]@cell.embeddings)
meta$UMAP1 = umap$UMAP_1
meta$UMAP2 = umap$UMAP_2

write.csv(meta,'meta_noPD_cycle_regress_pseudotime_UMAP.csv')

