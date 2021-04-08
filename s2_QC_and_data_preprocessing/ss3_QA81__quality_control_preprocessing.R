library(bglab)
library(geneplotter)
library(gplots)
require(stringr) 
library(gatepoints)
library(plotly)
library(plyr)

########################################################################################################################
# Load counts and metadata 
########################################################################################################################
# 
counts_14878 <- read.csv('SLX-14878_HTSeq_counts.txt',
                         sep='\t',row.names =1)

counts_14930 <- read.csv('SLX-14930_HTSeq_counts.txt',
                         sep='\t',row.names =1)

counts_16064 <- read.csv('SLX-16064_HTSeq_counts.txt',
                         sep='\t',row.names =1)

counts <- cbind(counts_14878, counts_14930, counts_16064)
counts <- as.matrix(counts)

meta <- read.csv('QA81_metadata_correct.csv')

counts <- counts[,meta$ID]

########################################################################################################################
#Download the appropriate gene and ensembl IDs from ensembl.
########################################################################################################################

library(biomaRt)
host <- "jul2015.archive.ensembl.org"

listMarts() #tells you which BioMart databases are available; we are chosing Ensembl
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host=host)

listDatasets(ensembl) #tells you which datasets are available in the chosen database; we are chosing Ensembl
ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)

#values we are interested in to retrieve from Biomart
attributes <- listAttributes(ensembl)
geneTable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),mart = ensembl)
colnames(geneTable) <- c("Ensembl_Gene_ID", "Associated_Gene_Name")
geneTable <- geneTable[!duplicated(geneTable$Ensembl_Gene_ID),]


#######################################################################
# read gene table and clear ERCC
#####################################################################

qcInd <- grep("^__", rownames(counts))
erccInd <- grep("^ERCC", rownames(counts))

ercc <- counts[erccInd,]
qc <- counts[qcInd,]
counts <- counts[-c(qcInd, erccInd),]

scd <- newSCD("RNAseq", counts = counts, genoData = geneTable, 
              spike = ercc, qc = qc, phenoData = meta)

# metalaneID is the CRI identifier
scd2 <- performQC(scd, pdf = "QA81_comrade_align_qc_2x10e5", metaLaneID ='CRI_identifier')

fdr <- 0.05

hvg_all <- read.csv('var.out.nospike.csv',row.names = 1)

hvgs <- hvg_all[which(hvg_all$FDR<fdr),] # 8522 

hvgs_list <- rownames(hvgs)

origFilter <- bglab:::saveFilters(scd2)
filterGene(scd2) <- FALSE
scd2 <- bglab::selectVariableGenes(scd2, includeGenes = hvgs_list,
                              reset = TRUE)
scd2 <- bglab:::restoreFilters(scd2, origFilter)

filterGene(scd2) <- TRUE

dim(scd2) # 8522 456

ColorsOrdered <- c('blue','cornflowerblue', 'khaki2','lavenderblush3','hotpink1','red')

scd3<- runPCA(scd2, scale. = TRUE)

plot(scd3, reduceMethod = "pca", colorBy = "Details", plotCols=ColorsOrdered, outline = "black", plotLegend=FALSE, main="PCA")

X11()
plot(scd3, reduceMethod = "pca", colorBy = "Details", plotCols=ColorsOrdered, outline = "black", plotLegend=FALSE, main="PCA")
outlier1 <-fhs(eigenvecs(getPCA(scd3))[,1:2]) 

outlier <- outlier1[1:3]  #"SLX.16064.i711_i517" "SLX.16064.i726_i515" "SLX.16064.i716_i518"

scd453 <- excludeCells(scd3, cellNames = outlier) # 8522, 453

meta453<-pData(scd453)

#Perform PCA (on technically variable genes only)
scd453 <- runPCA(scd453)


scd453 <- runTSNE(scd453, scale. = TRUE)

saveRDS(scd453, 'scd453_scran_bglab.rds')

