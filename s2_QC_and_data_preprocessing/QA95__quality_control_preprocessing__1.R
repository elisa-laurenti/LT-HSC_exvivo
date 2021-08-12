
suppressMessages(library(bglab))

#############################################
#
#############################################

count_19175 <- read.csv('SLX-19175_HTSeq_counts.txt',sep='\t',row.names =1)

count_19306 <- read.csv('SLX-19306_HTSeq_counts.txt',sep='\t',row.names =1)


counts <- cbind(count_19175 , count_19306)
counts <- as.matrix(counts)

meta <- read.csv("meta_qa95_original.csv",row.names=1)

meta <- meta[colnames(counts),]

#############################################
#
#############################################

geneTable <- readRDS('geneTable.rds')

#############################################
#
#############################################

qcInd <- grep("^__", rownames(counts))
erccInd <- grep("^ERCC", rownames(counts))

ercc <- counts[erccInd,]
qc <- counts[qcInd,]
counts <- counts[-c(qcInd, erccInd),]

#############################################
#
#############################################

scd <- newSCD("RNAseq", counts = counts, genoData = geneTable, spike = ercc, qc = qc, phenoData = meta)


scd2 <- performQC(scd, pdf = "qc_2x10e5", metaLaneID = 2)

scd2 <- techVar(scd2, useERCC = TRUE, meanForFit = 10)

#############################################
#          getting outliers
#############################################

scd2 <- runPCA(scd2, scale. = TRUE)


plot_pca <- scd2@pca@eigenvectors

plot_pca <- as.data.frame(plot_pca)
plot_pca$Sample_name <- pData(scd2)$Sample_name

ggplot(plot_pca, aes(x=PC1, y=PC2, color=Sample_name)) +
  geom_point()


outlier2 <- plot_pca[plot_pca[,1] < -80 | plot_pca[,2] > 40,]

#############################################
#
#############################################

scd3 <- excludeCells(scd2, cellNames = rownames(outlier2))

scd3 <- techVar(scd3, useERCC = TRUE, meanForFit = 10)

scd3 <- runPCA(scd3, scale. = TRUE)

plot_pca3 <- scd3@pca@eigenvectors
plot_pca3 <- as.data.frame(plot_pca3)
plot_pca3$Sample_name <- pData(scd3)$Sample_name

ggplot(plot_pca3, aes(x=PC1, y=PC2, color=Sample_name)) +
  geom_point()

plot_pca3$Plate.number <- pData(scd3)$Plate.number

ggplot(plot_pca3, aes(x=PC1, y=PC2, color=Plate.number)) +
  geom_point()


save.image('aa1_QA95__quality_control_preprocessing.RData')

