#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BASiCS))


print('loading objects')

Data_0 <- readRDS('Data_0__8464_85.rds')

Data_6 <- readRDS('Data_6__8464_134.rds')

Data_24 <- readRDS('Data_24__8464_86.rds')

Data_72 <- readRDS('Data_72__8464_124.rds')


N <- 20000
Thin <- 20
Burn <- 10000

print('doing 0h')

Chain_0h <- BASiCS_MCMC(Data = Data_0, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = FALSE,  Regression = TRUE,
                             PrintProgress = FALSE)

print('done 0h')

print('saving 0h')

saveRDS(Chain_0h, 'MCMC_Chain_0h_morethan_20RPM.rds')


print('doing 6h')
Chain_6h <- BASiCS_MCMC(Data = Data_6, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = FALSE,  Regression = TRUE,
                             PrintProgress = FALSE)
print('doing 6h')

print('saving 6h')

saveRDS(Chain_6h, 'MCMC_Chain_6h_morethan_20RPM.rds')



print('doing 24h')
Chain_24h <- BASiCS_MCMC(Data = Data_24, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = FALSE,  Regression = TRUE,
                             PrintProgress = FALSE)

print('doing 24h')

print('saving 24h')

saveRDS(Chain_24h, 'MCMC_Chain_24h_morethan_20RPM.rds')



print('doing 72h')

Chain_72h <- BASiCS_MCMC(Data = Data_72, N =N, 
                             Thin = Thin, Burn =Burn, 
                             WithSpikes = FALSE,  Regression = TRUE,
                             PrintProgress = FALSE)

print('doing 72h')

print('saving 72h')

saveRDS(Chain_72h, 'MCMC_Chain_72h_morethan_20RPM.rds')



print('saving image')

save.image('calculate_MCMC_morethan_20RPM.RData')

print('all done')



