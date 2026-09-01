rm(list=ls()) # clean

require(RTMB)
library(rstan)
library(StanHeaders)
library(tmbstan)
library(shinystan)

setwd("C:\\Research\\Iwc26\\NP_Humpbacks\\")

load(file="Bayes.RData")
launch_shinystan(mcmcout)

# Could save this as a data.frame
post2 <- rstan::extract(mcmcout)
Nsim <- length(post2$lp__)
print(attributes(post2))
