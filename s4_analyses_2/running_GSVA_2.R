#!/usr/bin/env Rscript

suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(GSVAdata))

print('loading gene set object')
gene_set_list<- readRDS('gene_set_list_c2_all_v7_1_symbols.rds')


print('loading counts')
counts <- readRDS('qa81_22__33411_429_symbols.rds')



print('making gene set object')
new_gene_set_list <- c()

for (i in names(gene_set_list)){
    # print(i)
    # print(length(gene_set_list[[i]]))
    
    temp_geneset <- GeneSet(gene_set_list[[i]], setName=i)
    
    new_gene_set_list <- c(new_gene_set_list, temp_geneset)
}

print('making GeneSetCollection')
gsc <- GeneSetCollection(new_gene_set_list)

print('converting counts to matix')
mtx_counts <- as.matrix(counts)


print('doing gsva')
res <- gsva(expr = mtx_counts, gset.idx.list=gsc, 
			min.sz=10, max.sz=500, parallel.sz=12, verbose=TRUE)

print('finished gsva')

print('saving image')

save.image('gsva_qa81_22.RData')

print('done')

############################
# saving
############################

res_df <- as.data.frame(res)

write.csv(res_df, 'GSVA_result_matrix.csv')

