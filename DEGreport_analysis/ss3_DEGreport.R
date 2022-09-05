#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))


source('my_cluster_functions.R')

load('ss2_DEGreport.RData')

##########################
#   change the number cluster you want here
##########################
num_of_cluster <- 13  #(numeric)
select <- cutree(as.hclust(cluster_genes), k=num_of_cluster)
#******************************************


########################################
select <- select[select %in% names(table(select))[table(select) > minc]]

if (reduce)
    select <- .reduce(select, counts_group)
select <- select[select %in% names(table(select))[table(select) > minc]]
message("Working with ", length(select), " genes after filtering: minc > ",minc)


groups <- select

df <- data.frame(genes = names(groups), 
                cluster = groups, stringsAsFactors = FALSE)

raw <- counts_group %>% as.data.frame %>% 
    rownames_to_column("genes") %>%
    gather(!!sym(summarize), "value", -genes) %>%
    inner_join(metadata_groups %>%
                   mutate_if(is.factor, as.character)) %>%
    inner_join(df, by = "genes")


summarise <- raw %>%
    group_by(!!sym(summarize), !!sym("cluster"),
             !!sym(time), !!sym(col)) %>%
    summarise(abundance = median(value),
              n_genes = n()) %>% 
    ungroup()

normalized <- norm_sign %>% as.data.frame() %>% 
rownames_to_column("genes") %>%
gather(!!sym(summarize), "value", -genes) %>%
inner_join(metadata_groups %>%
               mutate_if(is.factor, as.character)) %>%
inner_join(df, by = "genes") 



table <- normalized
# time <-
color = NULL
min_genes = 10
process = FALSE
points = TRUE
boxes = TRUE
smooth = TRUE
lines = TRUE
facet = TRUE
cluster_column = "cluster"
prefix_title = "Group: "

table$timepoint <- factor(table$timepoint, levels = c('LT_0h',
                                                'LT_6h',
                                                'LT_24h_UNTR',
                                                'LT_72h_UNTR'))


stopifnot(class(table)[1] == "data.frame")

if (cluster_column  %in% colnames(table)){
    table[["cluster"]] = table[[cluster_column]]
}
if (process){
    table <- .process(table, time, color)
}

if ("cluster"  %in% colnames(table)){
    counts <- table(distinct(table, genes, cluster)[["cluster"]])
    counts <- counts[counts>=min_genes]
    if (length(counts)==0)
        stop("No clusters with min_genes > ", min_genes)
    table <- inner_join(table,
                        data.frame(cluster = as.integer(names(counts)),
                                   title = paste(prefix_title,
                                                  names(counts),
                                                  "- genes:" ,
                                                  counts),
                                   stringsAsFactors = FALSE),
                        by = "cluster")
}

if (is.null(color)){
    color = "dummy"
    table[[color]] = ""
    lines = FALSE
}
table[["line_group"]] = paste(table[["genes"]],
                              table[[color]])


splan <- length(unique(table[[time]])) - 1L



save.image('ss3_DEGreport.RData')






