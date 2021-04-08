########################################################################################################################
#Download the appropriate gene and ensembl IDs from ensembl.
########################################################################################################################

library(biomaRt)
host <- "jul2015.archive.ensembl.org"

listMarts() #tells you which BioMart databases are available; we are chosing Ensembl
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host=host)

listDatasets(ensembl) #tells you which datasets are available in the chosen database; we are chosing Ensembl
ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)

#values we are interested in to retrieve from Biomart
attributes <- listAttributes(ensembl)
geneTable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),mart = ensembl)
colnames(geneTable) <- c("Ensembl_Gene_ID", "Associated_Gene_Name")
geneTable <- geneTable[!duplicated(geneTable$Ensembl_Gene_ID),]

saveRDS(geneTable, 'geneTable.rds')


