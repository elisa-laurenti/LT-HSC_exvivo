#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))




source('my_cluster_functions.R')

# ######################
#
# ##########################

ma_padj05<- readRDS('ma_padj05.rds')
design <- readRDS('design.rds')

print(ls())

print(dim(ma_padj05))

print(dim(design))


ma <- ma_padj05
metadata <- design
minc=15
summarize="merge"
time="timepoint"
col=NULL
consensusCluster = FALSE
reduce=TRUE
cutoff=0.70
scale=TRUE
pattern = NULL
groupDifference = NULL
eachStep = TRUE
plot=TRUE
fixy=NULL


benchmarking <- NULL
metadata <- as.data.frame(metadata)
ma = ma[, row.names(metadata)]
rownames(ma) = make.names(rownames(ma))
if (is.null(col)){
    col = "colored"
    metadata[,col] = rep("one_group", nrow(metadata))
}

if (!summarize %in% names(metadata))
    metadata[,summarize] = as.factor(paste0(metadata[,col],
                                            metadata[,time]))

# ensure there are no missing levels in summarize
metadata[,summarize] = droplevels(metadata[,summarize])

stopifnot(class(metadata)[1] == "data.frame")
stopifnot(class(ma)[1] == "matrix" | class(ma)[1] == "data.frame")
stopifnot(summarize %in% names(metadata))
stopifnot(time %in% names(metadata))

ma <- as.matrix(ma)


if (!is.null(fixy))
        stopifnot(length(fixy) == 2)
if (nrow(ma)>3000 & is.null(pattern))
        message("A large number of genes was given-- please, ",
                "make sure this is not an error. Normally, ",
                "only DE genes will be useful for this function.")
    message("Working with ", nrow(ma), " genes.")

counts_group <- .summarize_scale(ma,
                                 metadata[[summarize]],
                                 FALSE)

if (!is.null(groupDifference))
    counts_group <- .remove_low_difference(counts_group,
                                           groupDifference,
                                           eachStep)
if (scale){
    norm_sign <- t(apply(counts_group, 1, .scale))
}else {
    norm_sign <- counts_group
}

colnames(norm_sign) <- colnames(counts_group)

metadata_groups <- metadata %>%
        dplyr::distinct_(summarize, .keep_all = TRUE)
    rownames(metadata_groups) = metadata_groups[,summarize]

norm_sign <- norm_sign[, row.names(metadata_groups), drop = FALSE]


print('doing the clutering')

if (!consensusCluster & is.null(pattern)){
    print('yes')
    cluster_genes = .make_clusters(counts_group)
#     groups <- .select_genes(cluster_genes, norm_sign, minc,
#                            reduce = reduce,
#                            cutoff = cutoff)
#     benchmarking <- .benckmark_cutoff(cluster_genes, norm_sign, minc)
}


print('done and saving image')

save.image('ss2_DEGreport.RData')










