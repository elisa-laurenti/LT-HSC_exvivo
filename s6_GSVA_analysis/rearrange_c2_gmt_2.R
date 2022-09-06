
save_dir<- 'gene_set_gmt/c2_all_v7_1_symbols_csv'

file_list <- list.files(save_dir,'.csv',full.names = TRUE)

temp_df <- read.csv(paste0(save_dir,'/',file_list[1]))[2]

gene_set_list<- list()

for (file in file_list){
    temp_df <- read.csv(file)[2]
    pathway_name <- names(temp_df)
    
    gene_list <- as.character(temp_df[,pathway_name])
    
    gene_set_list[[pathway_name]] = gene_list
}


saveRDS(gene_set_list,'gene_set_list_c2_all_v7_1_symbols.rds')

