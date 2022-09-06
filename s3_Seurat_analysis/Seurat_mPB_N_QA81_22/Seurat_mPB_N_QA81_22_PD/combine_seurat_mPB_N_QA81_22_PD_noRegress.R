library(Seurat)
library(sctransform)

qa81 <-readRDS('seurat_QA_81_cells_361_PD.rds')
qa22 <-readRDS('seurat_QA_22_cells_175_PD.rds')
mPB <- readRDS('seurat_mPB_418cells_noRegress.rds')


obj_list <- c(qa81,qa22,mPB)
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list =  obj_list, anchor.features = features)

seurat.anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
    anchor.features = features, k.filter =100)

seurat.combined.sct <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT")

saveRDS(seurat.combined.sct, 'seurat.combined.QA81_22_PD_N_mPB.rds')

