suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS('dds_noPD.rds') 

############################
#  doing deseq
############################
dds <- DESeq2::DESeq(dds, parallel=TRUE)

###########################
# extracting results
############################
resultsNames(dds)

res_6_0 <- results(dds, alpha=0.05, tidy=TRUE, parallel=TRUE, name= 'timepoint_LT_6h_vs_LT_0h')

res_24_0 <- results(dds, alpha=0.05, tidy=TRUE, parallel=TRUE, name= 'timepoint_LT_24h_UNTR_vs_LT_0h')

res_72_0 <- results(dds, alpha=0.05, tidy=TRUE, parallel=TRUE, name= 'timepoint_LT_72h_UNTR_vs_LT_0h')


rownames(res_6_0) <- res_6_0$row
rownames(res_24_0) <- res_24_0$row
rownames(res_72_0) <- res_72_0$row

res_24_6 <- results(dds, alpha=0.05, parallel=TRUE, contrast=c('timepoint','LT_24h_UNTR','LT_6h'))

res_72_6 <- results(dds, alpha=0.05, parallel=TRUE, contrast=c('timepoint','LT_72h_UNTR','LT_6h'))

res_72_24 <- results(dds, alpha=0.05, parallel=TRUE, contrast=c('timepoint','LT_72h_UNTR','LT_24h_UNTR'))


########################################################################################################################
# read in geneTable
########################################################################################################################

geneTable <- readRDS('geneTable.rds')


rownames(geneTable) <- geneTable$Ensembl_Gene_ID


#############################################################
#  saving results in csv
#############################################################

differential_expressed_dir <- 'deseq2_results_noPD_NOfilter_oneList'
dir.create(differential_expressed_dir)

for (i in names(result_list)){
    temp_df <- result_list[[i]]

    new_geneTable <- geneTable[rownames(temp_df),]

    temp_df$geneName <- mapvalues(rownames(temp_df),from=new_geneTable$Ensembl_Gene_ID,to=new_geneTable$Associated_Gene_Name)

    print(dim(temp_df))
    print(i)
    dir_name <- paste0(differential_expressed_dir,'/',i)
    print(dir_name)
    dir.create(dir_name)

    write.csv(temp_df, paste0(dir_name,'/',i,'_NOfilter.csv'))
    
}

