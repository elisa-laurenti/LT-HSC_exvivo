
library(Seurat)
library(sctransform)



seurat.combined.sct <- readRDS('seurat.combined.QA81_22__CellCycleRegress.rds')


seurat.combined.sct  <- RunPCA(seurat.combined.sct , verbose = FALSE)


seurat.combined.sct <- RunUMAP(seurat.combined.sct, 
                               reduction = "pca", dims = 1:30, return.model = TRUE)



saveRDS(seurat.combined.sct, 'seurat.combined.QA81_22_pca_umap__CellCycleRegress.rds')
