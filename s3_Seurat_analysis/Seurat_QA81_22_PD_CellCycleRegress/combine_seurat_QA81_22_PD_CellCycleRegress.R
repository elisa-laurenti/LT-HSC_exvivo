library(Seurat)
library(sctransform)

qa81 <- readRDS('seurat_QA_81_361cells__CycleRegress_PD.rds')

qa22 <- readRDS('seurat_QA_22_175cells__CycleRegress_PD.rds')

obj_list <- c(qa81,qa22)

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

obj_list <- PrepSCTIntegration(object.list =  obj_list, anchor.features = features)


seurat.anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
    anchor.features = features, k.filter =100)


seurat.combined.sct <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT")


print(table(seurat.combined.sct@meta.data$orig.ident))


saveRDS(seurat.combined.sct, 'seurat.combined.QA81_22_PD_CycleRegress.rds')

