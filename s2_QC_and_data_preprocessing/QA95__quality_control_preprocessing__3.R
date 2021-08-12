
suppressMessages(library(bglab))

load('aa1_QA95__quality_control_preprocessing.RData')


meta_new_qc <- readRDS('meta_new_qc.rds')

##################################
#
##################################

failed_meta <- pData(scd)[-which(rownames(pData(scd)) %in% rownames(meta_new_qc)),]

scd4 <- excludeCells(scd2, cellNames = rownames(failed_meta))

scd4 <- techVar(scd4, useERCC = TRUE, meanForFit = 10)

scd4 <- runPCA(scd4, scale. = TRUE)

plot_pca4 <- scd4@pca@eigenvectors
plot_pca4 <- as.data.frame(plot_pca4)
plot_pca4$Sample_name <- pData(scd4)$Sample_name

ggplot(plot_pca4, aes(x=PC1, y=PC2, color=Sample_name)) +
  geom_point()


plot_pca4$Plate.number <- pData(scd4)$Plate.number


save.image('QA95__quality_control_preprocessing__3.RData')

write.csv(meta_new_qc,'meta_new_qc_220_21.csv')


