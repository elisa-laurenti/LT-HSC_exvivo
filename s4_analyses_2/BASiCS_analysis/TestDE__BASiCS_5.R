suppressPackageStartupMessages(library(BASiCS))

load('calculate_MCMC_morethan_20RPM.RData')

###############################
# 0h and 6h
###############################
de_0h__6h <- BASiCS_TestDE(Chain1 = Chain_0h, Chain2 = Chain_6h,
                      GroupLabel1 = "0h", GroupLabel2 = "6h",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

###############################
# 0h and 24h
###############################
de_0h__24h <- BASiCS_TestDE(Chain1 = Chain_0h, Chain2 = Chain_24h,
                      GroupLabel1 = "0h", GroupLabel2 = "24h",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

###############################
# 0h and 72h
###############################
de_0h__72h <- BASiCS_TestDE(Chain1 = Chain_0h, Chain2 = Chain_72h,
                      GroupLabel1 = "0h", GroupLabel2 = "72h",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

###############################
# 6h and 24h
###############################
de_6h__24h <- BASiCS_TestDE(Chain1 = Chain_6h, Chain2 = Chain_24h,
                      GroupLabel1 = "6h", GroupLabel2 = "24h",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

###############################
# 6h and 72h
###############################
de_6h__72h <- BASiCS_TestDE(Chain1 = Chain_6h, Chain2 = Chain_72h,
                      GroupLabel1 = "6h", GroupLabel2 = "72h",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

###############################
# 24h and 72h
###############################
de_24h__72h <- BASiCS_TestDE(Chain1 = Chain_24h, Chain2 = Chain_72h,
                      GroupLabel1 = "24h", GroupLabel2 = "72h",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)


###############################
# saving
###############################
save.image('BASiCS_TestDE__all_comparisons__morethan_20RPM.RData')


