suppressPackageStartupMessages(library(BASiCS))
suppressPackageStartupMessages(library(coda))

load('calculate_MCMC_morethan_20RPM.RData')

################################
#  MU
#################################

par(mfrow=c(2,2))
plot(apply(Chain_0h@parameters$mu,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('0h'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)

plot(apply(Chain_6h@parameters$mu,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('6h'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)

plot(apply(Chain_24h@parameters$mu,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('24h'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)

plot(apply(Chain_72h@parameters$mu,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('72h'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)


################################
#  delta
#################################

par(mfrow=c(2,2))
plot(apply(Chain_0h@parameters$delta,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('0h_delta'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)

plot(apply(Chain_6h@parameters$delta,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('6h_delta'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)

plot(apply(Chain_24h@parameters$delta,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('24h_delta'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)

plot(apply(Chain_72h@parameters$delta,1,median), type='l', btyt='n', xlab='Iterations',
     ylab='Draws', main=expression('72h_delta'),
    cex.main=2.5, cex.lab=1.5, cex.axis=1.5)


#################################
# coda package
################################# 

mu_mcmc_0h = mcmc(Chain_0h@parameters$mu)
mu_mcmc_6h = mcmc(Chain_6h@parameters$mu)
mu_mcmc_24h = mcmc(Chain_24h@parameters$mu)
mu_mcmc_72h = mcmc(Chain_72h@parameters$mu)

delta_mcmc_0h = mcmc(Chain_0h@parameters$delta)
delta_mcmc_6h = mcmc(Chain_6h@parameters$delta)
delta_mcmc_24h = mcmc(Chain_24h@parameters$delta)
delta_mcmc_72h = mcmc(Chain_72h@parameters$delta)

geweke_mu_0h = geweke.diag(mu_mcmc_0h)
geweke_mu_6h = geweke.diag(mu_mcmc_6h)
geweke_mu_24h = geweke.diag(mu_mcmc_24h)
geweke_mu_72h = geweke.diag(mu_mcmc_72h)

geweke_delta_0h = geweke.diag(delta_mcmc_0h)
geweke_delta_6h = geweke.diag(delta_mcmc_6h)
geweke_delta_24h = geweke.diag(delta_mcmc_24h)
geweke_delta_72h = geweke.diag(delta_mcmc_72h)

summary(cbind(geweke_mu_0h$z, geweke_mu_6h$z, 
              geweke_mu_24h$z, geweke_mu_72h$z, 
              geweke_delta_0h$z, geweke_delta_6h$z,
             geweke_delta_24h$z, geweke_delta_72h$z))


#####################################
# saving
####################################

save.image('convergence_check_all_timepoints_morethan_20RPM.RData')


