suppressPackageStartupMessages(library(BASiCS))



print('loading objects')

Data_0 <- readRDS('BASiCS_LT_0h.rds')

Data_6 <- readRDS('BASiCS_LT_6h.rds')

Data_24 <- readRDS('BASiCS_LT_24h.rds')

Data_72 <- readRDS('BASiCS_LT_72h.rds')


N <- 10000
Thin <- 10
Burn <- 1000


print('doing 0h')

Chain_0h <- BASiCS_MCMC(Data = Data_0, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = TRUE,  Regression = TRUE,
                             PrintProgress = TRUE, Threads=6)

print('done 0h')

print('saving 0h')

saveRDS(Chain_0h, 'MCMC_Chain_0h_morethan_20RPM_ERCC_Regression.rds')


print('doing 6h')
Chain_6h <- BASiCS_MCMC(Data = Data_6, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = TRUE,  Regression = TRUE,
                             PrintProgress = TRUE, Threads=6)
print('doing 6h')

print('saving 6h')

saveRDS(Chain_6h, 'MCMC_Chain_6h_morethan_20RPM_ERCC_Regression.rds')



print('doing 24h')
Chain_24h <- BASiCS_MCMC(Data = Data_24, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = TRUE,  Regression = TRUE,
                             PrintProgress = TRUE, Threads=6)

print('doing 24h')

print('saving 24h')

saveRDS(Chain_24h, 'MCMC_Chain_24h_morethan_20RPM_ERCC_Regression.rds')



print('doing 72h')

Chain_72h <- BASiCS_MCMC(Data = Data_72, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = TRUE,  Regression = TRUE,
                             PrintProgress = TRUE, Threads=6)

print('doing 72h')

print('saving 72h')

saveRDS(Chain_72h, 'MCMC_Chain_72h_morethan_20RPM_ERCC_Regression.rds')



print('saving image')

save.image('calculate_MCMC_morethan_20RPM_ERCC_Regression.RData')

print('all done')

