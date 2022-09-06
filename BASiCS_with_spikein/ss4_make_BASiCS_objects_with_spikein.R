suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(BASiCS))


df_ercc = read.csv('ERCC_92_num_molecules_for_BASiCS.csv',row.names=1)

ercc_22 = read.csv('ercc_qa22_92.csv',row.names=1)

ercc_81 = read.csv('ercc_qa81_92.csv',row.names=1)

meta <- read.csv('meta_QA81_22_noPD.csv', row.names=1)

combine_ercc = cbind(ercc_22,ercc_81)


counts = readRDS('CountsFilter_QA81_22_noPD.rds')


ercc_429cells = combine_ercc[,colnames(counts)]

identical(colnames(ercc_429cells),colnames(counts))

rownames(df_ercc) = df_ercc$SpikeID


timepoints = c('LT_0h','LT_6h','LT_24h_UNTR','LT_72h_UNTR')


list_counts=list()
list_tech=list()
list_spikeinfo=list()
list_batchinfo=list()

for (time in timepoints){
    print(time)
    
    temp_meta = meta[meta$Details == time,]
    print(dim(temp_meta))
    
    temp_counts = counts[,rownames(temp_meta)]
    print(dim(temp_counts))
    
    temp_ercc = ercc_429cells[,rownames(temp_meta)]
    print(dim(temp_ercc))
    
    noZero_ercc = temp_ercc[rowSums(temp_ercc[])>0,]
    print(dim(noZero_ercc))
    
    sub_ercc = df_ercc[rownames(noZero_ercc),]
    print(dim(sub_ercc))
    
    counts_ercc = rbind(temp_counts,noZero_ercc)
    
    erccInd <- grepl("^ERCC", rownames(counts_ercc))
    
    list_counts[[time]] = as.matrix(counts_ercc)
    list_tech[[time]] = erccInd
    list_spikeinfo[[time]] = sub_ercc
    list_batchinfo[[time]] = temp_meta$Cell_type_subtype
    
    print('=========================')
}


data_0h = newBASiCS_Data(Counts= list_counts[["LT_0h"]],    ## as.matrix(c_0), 
                         Tech= list_tech[["LT_0h"]],    ##erccInd, 
                         SpikeInfo = list_spikeinfo[["LT_0h"]],   ##sub_df_ercc,
                         BatchInfo = list_batchinfo[["LT_0h"]])    ## m_0$Cell_type_subtype) ## Cell_type_subtype

data_6h = newBASiCS_Data(Counts= list_counts[["LT_6h"]],    ## as.matrix(c_0), 
                         Tech= list_tech[["LT_6h"]],    ##erccInd, 
                         SpikeInfo = list_spikeinfo[["LT_6h"]],   ##sub_df_ercc,
                         BatchInfo = list_batchinfo[["LT_6h"]])    ## m_0$Cell_type_subtype) ## Cell_type_subtype

data_72h = newBASiCS_Data(Counts= list_counts[["LT_72h_UNTR"]],    ## as.matrix(c_0), 
                         Tech= list_tech[["LT_72h_UNTR"]],    ##erccInd, 
                         SpikeInfo = list_spikeinfo[["LT_72h_UNTR"]],   ##sub_df_ercc,
                         BatchInfo = list_batchinfo[["LT_72h_UNTR"]])    ## m_0$Cell_type_subtype) ## Cell_type_subtype


saveRDS(data_0h,'BASiCS_LT_0h.rds')

saveRDS(data_6h,'BASiCS_LT_6h.rds')

saveRDS(data_72h,'BASiCS_LT_72h.rds')

############################################################

counts_0 = list_counts[['LT_24h_UNTR']]

meta_0 = meta[meta$Details == 'LT_24h_UNTR',]

identical(colnames(counts_0), rownames(meta_0))

meta_minus1 = meta_0[-which(rownames(meta_0) == 'SLX.14930.i701_i504'),]

counts_minus1 = counts_0[,rownames(meta_minus1)]

identical(colnames(counts_minus1), rownames(meta_minus1))

spike_in = list_spikeinfo[["LT_24h_UNTR"]]

batch_info = meta_minus1$Cell_type_subtype

tech_info  =list_tech[["LT_24h_UNTR"]]


data_24h = newBASiCS_Data(Counts= counts_minus1,    ## as.matrix(c_0), 
                         Tech= tech_info,    ##erccInd, 
                         SpikeInfo = spike_in,   ##sub_df_ercc,
                         BatchInfo = batch_info)    ## m_0$Cell_type_subtype) ## Cell_type_subtype

saveRDS(data_24h,'BASiCS_LT_24h.rds')


