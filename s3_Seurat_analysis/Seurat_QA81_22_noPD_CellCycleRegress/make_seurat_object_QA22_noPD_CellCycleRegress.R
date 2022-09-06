
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

seurat_cell_cycle <- SCTransform(seurat, assay='RNA', new.assay.name='SCT', vars.to.regress = "percent.mt")


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


seurat_cell_cycle <- CellCycleScoring(seurat_cell_cycle,s.features = s.genes,g2m.features = g2m.genes,assay = 'SCT',set.ident = TRUE)



seurat_cell_cycle$CC.Difference <- seurat_cell_cycle$S.Score - seurat_cell_cycle$G2M.Score



seurat <- SCTransform(seurat_cell_cycle, assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'CC.Difference')
)


######################################################
#
#####################################################

saveRDS(seurat, 'seurat_QA_22_148cells__CellCycleRegress.rds')




