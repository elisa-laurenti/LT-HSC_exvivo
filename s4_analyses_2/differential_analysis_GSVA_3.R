suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(GSVAdata))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(limma))


load('gsva_qa81_22.RData')

meta_all <- read.csv('meta_qa81_22__429cells.csv',row.names=1)

meta <- meta_all[,c('Details','Cell_type_subtype')]

colnames(meta) <- c('timepoint','batch')

identical(rownames(meta), colnames(res))

######################
#  making design matrix
######################
design <- model.matrix(~ 0+ meta$timepoint)

colnames(design) <- c('LT_0h',  'LT_24h_UNTR','LT_6h','LT_72h_UNTR' )

fit <- lmFit(res, design)

contrast.matrix <- makeContrasts(LT_72h_UNTR-LT_24h_UNTR, LT_72h_UNTR-LT_6h, LT_72h_UNTR-LT_0h,
                                 LT_24h_UNTR-LT_6h, LT_24h_UNTR-LT_0h,
                                 LT_6h-LT_0h,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

adjPvalueCutoff <- 0.05

allgenes_table <- list()

degenes_table <- list()

for (i in comparisons){
    comp_name <- gsub(' ','',i, fixed=TRUE)
    print(comp_name)
    
    temp_allGenes <- topTable(fit2, coef=i, number=Inf)
    print(dim(temp_allGenes))
#     print(head(temp_table, 1))
    temp_DEgeneSets <- topTable(fit2, coef=i, number=Inf, p.value=adjPvalueCutoff, adjust="BH")
    
    allgenes_table[[comp_name]] <- temp_allGenes
    
    degenes_table[[comp_name]] <- temp_DEgeneSets
    
}

new_colnames <- c()
for ( i in colnames(gsva_df)){
    print(i)
    comp_name <- gsub(' ','',i, fixed=TRUE)
    print(comp_name)
    new_colnames <- c(new_colnames,comp_name)
}

colnames(gsva_df) <- new_colnames

for (name in names(allgenes_table)){
    df <- allgenes_table[[name]]
    df$neg_log_pvalue <- -log10(df$P.Value)

    df$adj_P_005 <- 'nonsignificant'
    
    df$adj_P_005[which(df$adj.P.Val < 0.05)] <- 'significant'
    
    
    df<- df[rownames(gsva_df),]
    if (identical(rownames(gsva_df), rownames(df))){
        df$gsva_score <- gsva_df[[name]]
        
    }
    print(dim(df))
    
    df_list[[name]] <- df
}

for ( name in names(df_list)){

    file_name <- paste0('gsva_csv/' , name ,'_gsva.csv')
    print(file_name)
    
    print(dim(df_list[[name]]))
    write.csv(df_list[[name]], file_name)
}

