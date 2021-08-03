#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(SIBERG))

##############################
#  loading in data
##############################

print('load in data')


count_qa81 <- read.csv('QA81_raw_counts_after_qc.csv',row.names=1)

meta453 <- read.csv('QA81_meta_after_qc.csv',row.names=1)



#############################
#  filtering gene expression
##############################

print('filtering gene expression')

express_3cells <- rowSums(count_qa81!=0)>=3
print(table(express_3cells))

count_all_filtered <- count_qa81[express_3cells,]

print(dim(count_all_filtered))

print(colnames(count_all_filtered)[1:5])

print(rownames(count_all_filtered)[1:5])



#############################
#  estimating size factor
##############################

TMM <- calcNormFactors(count_all_filtered, method='TMM')

print(length(TMM))



#############################
#  
##############################

print('subsetting meta and counts by timepoint')

meta_sub <- droplevels(meta453)

meta_list <- list()
tmm_list <- list()
count_list <- list()

for (timepoint in unique(meta_sub$Details)){
    temp_meta <- meta_sub[meta_sub$Details == timepoint,]
    
    temp_tmm <- TMM[rownames(temp_meta)]
    
    temp_count <- count_all_filtered[,rownames(temp_meta)]
    
    meta_list[[timepoint]] <- temp_meta
    
    tmm_list[[timepoint]] <- temp_tmm
    
    count_list[[timepoint]] <- temp_count
    
    print(timepoint)
    print(dim(temp_meta))
    print(length(temp_tmm))
    print(dim(temp_count))
}




#############################
#  looping
##############################

print('looping right now')



df_ln_list <- list()
df_nb_list <- list()
df_gp_list <- list()

for (timepoint in unique(meta_sub$Details)){
    
    temp_count <- count_list[[timepoint]]
    
    temp_tmm <- tmm_list[[timepoint]]
    
    ln_list <- list()
    nb_list <- list()
    gp_list <- list()

    for(i in 1:nrow(temp_count)){   #####nrow(temp_count)
    #     print(i)
        ln_res <- SIBER(y=temp_count[i, ]+1, d=1/temp_tmm, model='LN')  ### if you dont add 1 will get an errorr
        nb_res <- SIBER(y=temp_count[i, ], d=1/temp_tmm, model='NB')
        gp_res <- SIBER(y=temp_count[i, ], d=1/temp_tmm, model='GP')


        temp_ln_list <- c()
        temp_nb_list <- c()
        temp_gp_list <- c()

        for (j in 1:length(ln_res)){
            temp_ln_list[j] <- ln_res[[j]]
        }
        for (j in 1:length(nb_res)){
            temp_nb_list[j] <- nb_res[[j]]
        }
        for (j in 1:length(gp_res)){
            temp_gp_list[j] <- gp_res[[j]]
        }

        ln_list[[i]] <- temp_ln_list
        nb_list[[i]] <- temp_nb_list
        gp_list[[i]] <- temp_gp_list


    }



    df_ln <- do.call('rbind',ln_list)
    df_nb <- do.call('rbind',nb_list)
    df_gp <- do.call('rbind',gp_list)


    colnames(df_ln) <- c('mu1','mu2','sigma1','sigma2','pi1','delta','BI')
    colnames(df_nb) <- c('mu1','mu2','sigma1','sigma2','pi1','delta','BI')
    colnames(df_gp) <- c('mu1','mu2','sigma1','sigma2','pi1','delta','BI')


    print('do they have the same length')
    print(dim(df_ln)[1] == dim(temp_count)[1])
    print(dim(df_nb)[1] == dim(temp_count)[1])
    print(dim(df_gp)[1] == dim(temp_count)[1])



    rownames(df_ln) <- rownames(temp_count)
    rownames(df_nb) <- rownames(temp_count)
    rownames(df_gp) <- rownames(temp_count)
    
    df_ln_list[[timepoint]] <- df_ln
    df_nb_list[[timepoint]] <- df_nb
    df_gp_list[[timepoint]] <- df_gp

}



###################################
# saving iamge
###################################

save.image('siber_QA81_by_timepoint.RData')



print('all done')


