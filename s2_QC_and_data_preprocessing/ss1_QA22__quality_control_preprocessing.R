library(bglab)
library(geneplotter)
library(gplots)
library(gatepoints)

########################################################################################################################
# Load counts and metadata 
########################################################################################################################

counts_0_6_h <- read.table("SLX-12561_HTSeq_counts_21Mar17.txt", header = TRUE, row.names = 1)
counts_24_48_h <- read.table("SLX-12562_HTSeq_counts_22Mar17.txt", header = TRUE, row.names = 1)
counts_72_h_PD <- read.table("SLX-12563_HTSeq_counts.txt", header = TRUE, row.names = 1)
counts <- cbind(counts_0_6_h, counts_24_48_h, counts_72_h_PD)
counts <- as.matrix(counts)


meta <- read.csv("metadata_final_SB_230317b.csv")

#also need to add a column called "cluster" which will specify color for each group to be used in graphs
ClusterCol <- as.character(meta$Details)
ClusterCol <-replace(ClusterCol, ClusterCol=="GMP_T0", "grey82")
ClusterCol <-replace(ClusterCol, ClusterCol=="LT-HSCs_T0", "blue")
ClusterCol <-replace(ClusterCol, ClusterCol=="GMP_T6", "grey62")
ClusterCol <-replace(ClusterCol, ClusterCol=="LT-HSCs_T6", "cornflowerblue")
ClusterCol <-replace(ClusterCol, ClusterCol=="GMP_T24", "grey42")
ClusterCol <-replace(ClusterCol, ClusterCol=="LT-HSCs_T24", "khaki2")
ClusterCol <-replace(ClusterCol, ClusterCol=="GMP_T48", "grey22")
ClusterCol <-replace(ClusterCol, ClusterCol=="LT-HSCs_T48", "orange")
ClusterCol <-replace(ClusterCol, ClusterCol=="GMP_T72_UN", "black")
ClusterCol <-replace(ClusterCol, ClusterCol=="GMP_T72_PD", "darkgreen")
ClusterCol <-replace(ClusterCol, ClusterCol=="LT-HSCs_T72_UN", "red")
ClusterCol <-replace(ClusterCol, ClusterCol=="LT-HSCs_T72_PD", "hotpink1")
ClusterCol <-replace(ClusterCol, ClusterCol=="NO_cDNA", "black")
ClusterCol <-replace(ClusterCol, ClusterCol=="NO_CDNA", "black")


meta <- cbind(meta, ClusterCol)

########################################################################################################################
# read in geneTable
########################################################################################################################

geneTable <- readRDS('geneTable.rds')

########################################################################################################################
#Extract QC produced by HTSeq and the ERCCs from the previously downloaded counts table.
########################################################################################################################

qcInd <- grep("^__", rownames(counts))
erccInd <- grep("^ERCC", rownames(counts))

ercc <- counts[erccInd,]
qc <- counts[qcInd,]
counts <- counts[-c(qcInd, erccInd),]

########################################################################################################################
#single cell dataset (SCD) object pipeline for the bglab package.
########################################################################################################################


scd <- newSCD("RNAseq", counts = counts, genoData = geneTable, spike = ercc, qc = qc, phenoData = meta)


########################################################################################################################
# 1- Perform QC on all cells
########################################################################################################################


scd2 <- performQC(scd, pdf = "qc_2x10e5", metaLaneID = 4)
PassedCells_2e5 <- colnames(counts(scd2))
PassedCells_pos <- match(PassedCells_2e5, meta$ID)
PassedCells_2e5_anno <- data.frame (Cell= PassedCells_2e5, Condition= meta$Details[PassedCells_pos])
table(PassedCells_2e5_anno$Condition) #322 cells passed                                   

scd1 <- performQC(scd, pdf = "qc_1x10e5", metaLaneID = 4, cutoffs = c(100000, 0.1, 0.2, 0,0,1,0))
PassedCells_1e5 <- colnames(counts(scd1))
PassedCells_pos1 <- match(PassedCells_1e5, meta$ID)
PassedCells_1e5_anno <- data.frame (Cell= PassedCells_1e5, Condition= meta$Details[PassedCells_pos1])
table(PassedCells_1e5_anno$Condition)  #348 cells passed

PassedCellsTable <- rbind(table(PassedCells_2e5_anno$Condition), table(PassedCells_1e5_anno$Condition))
rownames(PassedCellsTable) <- c("Normal", "LowCutOff")
PassedCellsTable_Ordered <- PassedCellsTable[,c(1,4,2,3,6,5,7,10,8,9,12,11)]


PassedCellsTable_Percent <- PassedCellsTable_Ordered
PassedCellsTable_Percent[,1:6] <- PassedCellsTable_Ordered[,1:6]/24*100
PassedCellsTable_Percent[,7:12] <- PassedCellsTable_Ordered[,7:12]/72*100
 
########################################################################################################################
# 2- Eliminate GMP samples and select technically variable genes as per Brennecke et al.
########################################################################################################################
ColorsOrdered <- c("grey82", "grey42", "grey22", "grey62", "darkgreen", "black", "blue", "khaki2", "orange", "cornflowerblue", "hotpink1", "red", "white", "white")

GMP_pos <- grep("GMP", PassedCells_2e5_anno$Condition)
GMPcellNames <- as.character(PassedCells_2e5_anno[GMP_pos,1])
scd_LT <- excludeCells(scd2, cellNames = GMPcellNames)
scd_LT <- techVar(scd_LT, useERCC = TRUE, meanForFit = 10)
plotTechVar(scd_LT@technicalNoise, pdf="QA22_LT_TechnicalNoisePlot_wERCC.pdf", main="with ERCC")
#238 LT-HSC passed

techVarGenes_wERCC <- fData(scd_LT)
dim(scd_LT)
#2105 HVGs


PassedLT <- colnames(exprs(scd_LT))
PassedLTpos <- match(PassedLT, meta$ID)
meta_LT <- meta[PassedLTpos,]

########################################################################################################################
# 3- Perform dimensionality reduction techniques
########################################################################################################################

#Perform PCA (on technically variable genes only)
scd_LT <- runPCA(scd_LT, scale. = TRUE)
PopsShown <- levels(meta_LT$Details)[c(7,10,8,9,12,11)]
plot(scd_LT, reduceMethod = "pca", colorBy = "Details", plotCols=ColorsOrdered, outline = "black", plotLegend=FALSE, main="PCA")
legend("bottomleft",  PopsShown, fill=ColorsOrdered[c(7,10,8,9,12,11)], cex = 0.9, bty = "n")

#eliminating 6h outlier by restricting axes (not eliminated from analysis)
plot(scd_LT, reduceMethod = "pca", colorBy = "Details", plotCols=ColorsOrdered, outline = "black", plotLegend=FALSE, main="PCA", ylim=c(-20,20))
legend("bottomright",  PopsShown, fill=ColorsOrdered[c(7,10,8,9,12,11)], cex = 0.7, bty = "n")


#this cell is an outlier also on pca3 and pca4 -> decided to eliminate altogether

X11()
plot(scd_LT, reduceMethod = "pca", colorBy = "Details", plotCols=ColorsOrdered, outline = "black", plotLegend=FALSE, main="PCA")
outlier <-fhs(eigenvecs(getPCA(scd_LT))[,1:2]) #SLX.12561.N716_S506
dev.off()

########################################################################################################################
# 4- Eliminate GMP and LT6h outlier samples and select technically variable genes as per Brennecke et al.
########################################################################################################################

GMP_pos <- grep("GMP", PassedCells_2e5_anno$Condition)
GMPcellNames <- as.character(PassedCells_2e5_anno[GMP_pos,1])
CellsToExclude <- c(GMPcellNames, outlier)
scd_LT <- excludeCells(scd2, cellNames = CellsToExclude)
scd_LT <- techVar(scd_LT, useERCC = TRUE, meanForFit = 10)
plotTechVar(scd_LT@technicalNoise, pdf="QA22_LT_TechnicalNoisePlot_wERCC.pdf", main="with ERCC")
#237 LT-HSC passed

techVarGenes_wERCC <- fData(scd_LT)
dim(scd_LT)
#2469 HVGs


PassedLT <- colnames(exprs(scd_LT))
PassedLTpos <- match(PassedLT, meta$ID)
meta_LT <- meta[PassedLTpos,]


saveRDS(scd_LT, 'scd_LT.rds')


