
library(Seurat)
library(sctransform)

meta <- read.csv('qa_22_meta_148cells.csv',row.names=1)
counts <- read.csv('qa22_raw_counts_148cells.csv',row.names=1)

seurat <- CreateSeuratObject(counts =counts, project ='QA_22',
                                   min.cells = 3, min.features = 0,
                                   meta.data = meta)

Idents(seurat)<- 'QA_22'
seurat@meta.data$orig.ident <- 'QA_22'

seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt")

saveRDS(seurat,'seurat_QA_22_cells_148.rds')

