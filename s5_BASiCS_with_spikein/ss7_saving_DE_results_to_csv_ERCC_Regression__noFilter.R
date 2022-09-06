suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(BASiCS))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(ggplot2))

load('BASiCS_TestDE__all_comparisons_ERCC_Regression.RData')
ls()

geneTable <- readRDS('geneTable.rds')

de_list <-list('de_0h__24h'=de_0h__24h,
               'de_0h__6h'=de_0h__6h,
               'de_0h__72h'=de_0h__72h,
               'de_24h__72h'=de_24h__72h,
               'de_6h__24h'=de_6h__24h,
               'de_6h__72h'=de_6h__72h)

separate_folder <- 'BASiCs_DE_results_separate_full_csv_ERCC_Regression__noFilter'
dir.create(separate_folder)

for (name in names(de_list)){
    
    mean_name <- paste0(separate_folder,'/',name, '__mean_ERCC_Regression__noFilter.csv')
    dispersion_name <- paste0(separate_folder,'/',name, '__over-dispersion_ERCC_Regression__noFilter.csv')
    residual_dispersion_name <- paste0(separate_folder,'/',name, '__residual_over-dispersion_ERCC_Regression__noFilter.csv')
    
    print(name)
    print(mean_name)
    print(dispersion_name)
    print(residual_dispersion_name)

    temp_resdisp = de_list[[name]]@Results$ResDisp@Table
    temp_disp = de_list[[name]]@Results$Disp@Table
    temp_mean = de_list[[name]]@Results$Mean@Table
    
    temp_resdisp$ensembl_id = temp_resdisp$GeneName
    temp_disp$ensembl_id = temp_disp$GeneName
    temp_mean$ensembl_id = temp_mean$GeneName
    
    temp_resdisp$GeneName=mapvalues(temp_resdisp$GeneName,
                                    from=geneTable$Ensembl_Gene_ID,
                                    to=geneTable$Associated_Gene_Name,
                                    warn_missing = FALSE)
    
    temp_disp$GeneName=mapvalues(temp_disp$GeneName,
                                    from=geneTable$Ensembl_Gene_ID,
                                    to=geneTable$Associated_Gene_Name,
                                    warn_missing = FALSE)
    
    
    temp_mean$GeneName=mapvalues(temp_mean$GeneName,
                                    from=geneTable$Ensembl_Gene_ID,
                                    to=geneTable$Associated_Gene_Name,
                                    warn_missing = FALSE)
    
    
    write.csv(temp_mean, mean_name)
    write.csv(temp_disp, dispersion_name)
    write.csv(temp_resdisp, residual_dispersion_name)
    
}






