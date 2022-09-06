library(Seurat)
library(sctransform)

meta <- read.csv('meta_mPB_418cells.csv',row.names=1)
counts <- read.csv('raw_counts_mPB_418cells.csv',row.names=1)

seurat <- CreateSeuratObject(counts =counts, project ='mPB',
                                   min.cells = 3, min.features = 0,
                                   meta.data = meta)

Idents(seurat)<- 'mPB'

seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt")

saveRDS(seurat, 'seurat_mPB_418cells_noRegress.rds')
