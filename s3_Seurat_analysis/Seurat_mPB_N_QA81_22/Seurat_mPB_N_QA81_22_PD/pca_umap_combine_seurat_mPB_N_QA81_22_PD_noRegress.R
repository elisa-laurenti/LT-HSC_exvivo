library(Seurat)
library(sctransform)

seurat.combined.sct <- readRDS('seurat.combined.QA81_22_PD_N_mPB.rds')

seurat.combined.sct  <- RunPCA(seurat.combined.sct , verbose = FALSE)

seurat.combined.sct <- RunUMAP(seurat.combined.sct, 
                               reduction = "pca", dims = 1:30, return.model = TRUE)

saveRDS(seurat.combined.sct, 'seurat.combined.QA81_22_PD_N_mPB_pca_umap.rds')

meta = seurat.combined.sct@meta.data

umap = as.data.frame(seurat.combined.sct[["umap"]]@cell.embeddings)

meta$UMAP1 = umap$UMAP_1
meta$UMAP2 = umap$UMAP_2

write.csv(meta,'meta_PD_QA81_22_N_mPB__UMAP.csv')
