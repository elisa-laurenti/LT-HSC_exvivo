
library(Seurat)
library(sctransform)

meta <- read.csv('qa_81_meta_281cells.csv',row.names=1)
counts <- read.csv('qa_81_raw_counts_281cells.csv',row.names=1)

seurat <- CreateSeuratObject(counts =counts, project ='QA_81',
                                   min.cells = 3, min.features = 0,
                                   meta.data = meta)

Idents(seurat)<- 'QA_81'
seurat@meta.data$orig.ident <- 'QA_81'

seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt")

saveRDS(seurat,'seurat_QA_81_cells_281.rds')


