
library(bglab)
library(geneplotter)
library(gplots)
# library(vioplot)
library(easyGgplot2)
#============================================================================================
# QCplots and outliers
#============================================================================================

scd_LT <- readRDS("scd_LT.rds")
 

scd_qc <- performQC(scd_LT, metaLaneID = "CRI_identifier" ,cutoffs = c(2e5, .3, .15, 0.75, 2000, 0.2, 0),pdf = NULL)

saveRDS(scd_qc, 'scd_qc.rds')
