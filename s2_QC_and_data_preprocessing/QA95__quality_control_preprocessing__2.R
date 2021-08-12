
suppressMessages(library(bglab))

load('aa1_QA95__quality_control_preprocessing.RData')

###################################
#  more cell  filtering
###################################


meta <- pData(scd)
dim(meta)

qcoutput <- scd3@qcOutput

failqc <- c(qcoutput$failQC$SLX.19175 , qcoutput$failQC$SLX.19306)

meta$not_aligned <- qcoutput$qcPlots$`__not_aligned:Total`

meta$total_reads <- qcoutput$totalReads

meta$failqc <- 'passed'

meta$failqc[which(rownames(meta) %in% failqc)] <- 'failed_qc'

meta$failqc[which(rownames(meta) %in% rownames(outlier2))] <- 'outliers'

ggplot(meta, aes(x=total_reads,y=not_aligned, color=failqc)) +
  geom_point()

meta$new_qc <- 'passed'

meta$new_qc[which(meta$not_aligned > 0.3)] <- 'failed_new'

meta$new_qc[which(meta$failqc != 'passed')] <- 'failed_new'

ggplot(meta, aes(x=total_reads,y=not_aligned, color=new_qc)) +
  geom_point()

passed_meta <- meta[which(meta$new_qc == 'passed'),]
dim(passed_meta)


saveRDS(passed_meta, 'meta_new_qc.rds')

