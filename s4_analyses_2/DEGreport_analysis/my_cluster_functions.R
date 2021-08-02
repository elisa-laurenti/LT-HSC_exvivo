#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scde))
suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(limma))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))


#############################
#  functions
#################################

.summarize_scale <- function(ma, group, scale = TRUE){
    counts_group = t(sapply(rownames(ma), function(g){
        sapply(levels(group), function(i){
            idx = which(group == i)
            mean(ma[g, idx], na.rm = TRUE)
        })
    }))
    colnames(counts_group) = levels(group)
#     .logger(head(counts_group), "summarize_scale::counts_group")
#     .logger(head(group), "summarize_scale::group")
    if (scale) {
        norm_sign <- t(apply(counts_group, 1, .scale))
    }else{
        norm_sign <- counts_group
    }
    colnames(norm_sign) = colnames(counts_group)
    norm_sign
}

.logger = function(toout,msg=""){
    logdebug(paste("\n\nchecking data:" , msg, "\n\n"))
    if (getLogger()[['level']]!=20)
        print(toout)
}

.scale <- function(e){
    scale(e)
}

.make_clusters <- function(counts_group){
    m <- (1 - cor(t(counts_group), method = "kendall"))
    d <- as.dist(m^2)
    c <- diana(d, diss = TRUE, stand = FALSE)
    c
}

.reduce <- function(groups, counts_group){
    lapply(unique(groups), function(g){
        ma <- counts_group[names(groups)[groups==g],]
        nokeep <- apply(ma, 2L, function(x){
            out <- boxplot(x, plot = FALSE)$out
            rownames(ma)[which(x %in% out)]
        }) %>% as.vector() %>% unlist() %>% unique()
        keep <- setdiff(rownames(ma), nokeep)
        group <- rep(g, length(keep))
        names(group) <- keep
        group
    }) %>%  unlist()
}

.benckmark_cutoff <- function(tree, counts, minc = 15){
    # browser()
    series <- unique(round(tree$height, digits = 3))
    series <- sort(series[2:length(series)])
    list_clusters <- lapply(series, function(s) {
        select <- cutree(as.hclust(tree), h = s)
        select <- select[select %in% names(table(select))[table(select) > minc]]
    })
    list_pct <- lapply(list_clusters, function(c) {
        .pct_var(counts, c)
    })
    names(list_clusters) <- paste0("cutoff", series)
    names(list_pct) <- paste0("cutoff", series)
    df_clusters <- lapply(names(list_clusters), function(s){
        c = list_clusters[[s]]
        if (length(c) == 0) {return(data.frame())}
        data.frame(genes = make.names(names(c)),
                   cluster = c,
                   cutoff = s,
                   stringsAsFactors = F)
    }) %>% bind_rows() %>% 
        distinct() %>% 
        spread(cutoff, cluster)
    list(genes=df_clusters, pcts=list_pct)
}

.pct_var <- function(counts, clusters){
    clusters_sd <- lapply(unique(clusters), function(c){
        counts[names(clusters[clusters == c]),, drop = FALSE] %>% 
            apply(., 2, sd) %>% 
            median(.)
    }) %>% unlist() %>%  median
    (1 - clusters_sd / median(apply(counts, 2, sd))) * 100
}


