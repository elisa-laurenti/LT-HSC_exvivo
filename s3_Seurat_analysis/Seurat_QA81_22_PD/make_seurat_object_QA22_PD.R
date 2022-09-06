library(Seurat)
library(sctransform)

meta <- read.csv('qa_22_meta_175cells_PD.csv',row.names=1)
counts <- read.csv('qa_22_raw_counts_175cells_PD.csv',row.names=1)

seurat <- CreateSeuratObject(counts =counts, project ='QA_22',
                                   min.cells = 3, min.features = 0,
                                   meta.data = meta)

Idents(seurat)<- 'QA_22'

seurat@meta.data$orig.ident <- 'QA_22'

seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt")

saveRDS(seurat, 'seurat_QA_22_cells_175_PD.rds')

