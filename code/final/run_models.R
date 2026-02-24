library(RTMB)
library(tidyverse)
library(here)
# ========================================================================================

# all the functions we need
source(here('code','final','read_data.R'))
source(here('code','final',"write_outputs.R"))
source(here('code','final',"plotting_functions.R"))

# Fitting Function
DoRun <- function(Code,SensCase,StrayBase=0,Nage=11,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),WithMirror=0,Yr1=1970,Yr2=2023,YrSDevs=2000, 
                  UseKPrior=1, Kmax=60000,
                  rvars=NULL, # which parameters are random
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",
                  envOpt="rS", splineK=5, envlag=0, envVars=c("sst","mld"),
                  AllPlots=T,DoBoot=F,BootUse,
                  ByCatchFile="BycatchActual_2024_04_24.csv",
                  DataFileName= here('data',"Hump.dat"),
                  FullDiag=T,
                  # PlotsDir= paste0(here('plots'),"/"),
                  WghtTotal=1,Idirichlet=1,MaxN=100,seed=19101,
                  DoBayes=F,Init=NULL,subdir="")
{
  
  cat("Doing ",paste0(Code,"-",SensCase)," using data file ",DataFileName,"\n")
  
  # =================================================================================================================================
  # Make Model Data
  # =================================================================================================================================
  # Pass the correct parameters to the data-making function
  run_args <- as.list(environment())
  
  # make the dataset
  dat <- do.call(MakeDataScenario,run_args)
  
  # make the data objects available to this function
  list2env(dat, envir = environment())
  # =================================================================================================================================
  # Load Model Version
  # =================================================================================================================================
  # source the correct model code, depending on "envOpt"
  # i.e., choice about how to drive survival with environmental variables
  if(envOpt== 'rS') source(here('code','final','rS.R'))
  if(envOpt== 'env-survival') source(here('code','final','env-survival.R'))
  if(envOpt== 'ddOnly') source(here('code','final','ddOnly.R'))
  if(envOpt== 'env-K') source(here('code','final','env-K.R'))
  # =================================================================================================================================
  # Build TMB Model
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling MakeADM")
  model <- MakeADFun(cmb(f,dat), parameters,random=rvars, map=map,DLL="Hump",silent=F)

  # =================================================================================================================================
  # Optimize
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling Fit")
  FitBest  <- 1.0e20
  fit <- nlminb(model$par, model$fn, model$gr,verbose=T)
  while(abs(FitBest-fit$objective)>0.01){
      FitBest <- fit$objective
      fit <- nlminb(fit$par, model$fn, model$gr,verbose=T)
      if (FullDiag==T) print(fit$objective)
  }

  # Save model
  fildir <- here("Diags",subdir)
  if(!dir.exists(fildir)) dir.create(fildir,recursive = T)
  readr::write_rds(model,paste0(fildir,"/",Code,SensCase,"_TMB.rds"))
  
  # =================================================================================================================================
  # Extract Results
  # =================================================================================================================================
  
  rept <- model$report()
  best <- model$env$last.par.best
  adr <-sdreport(model)
  sdr <- summary(adr)
  sdrf <- summary(adr, "fixed", p.value = TRUE)
  
  
  cat("final nll:",rept$neglogL,"\nTotal data likelihood:",rept$datalike,"\nTotal penalty:",rept$Penal,
      "\nSurvey likelihood:",rept$LogLike1,"\nMixing BtF likelihood:",rept$LogLike2a,"\nMixing FtB likelihood",rept$LogLike2b,"\n")
  #print(exp(LogNT))
  
 
  
  # =================================================================================================================================
  # Write and Plot
  # =================================================================================================================================
  
  WriteOut(Code,Abbrev=SensCase,rept=rept,sdr=sdr,sdrf=sdrf,data=dat,subdir = subdir)
  
  if(AllPlots){
    obj <- list(Code=Code,SensCase=SensCase,input=dat,report=rept,sdreport=sdr,sdfixed=sdrf)
    fildir <- here("plots",subdir)
    if(!dir.exists(fildir)) dir.create(fildir,recursive = T)
    
    suppressWarnings({
      plot_title <- paste0(fildir,"/",Code,"-",SensCase)
      # basic fixed effect parameters (minus environmental covariates)
      p0 <- plot_fixed_p(obj = obj)
      # abundance plots
      p1 <- plot_abundance(obj=obj,opt = "total")
      p2 <- plot_abundance(obj=obj,opt= "breed")
      p3 <- plot_abundance(obj=obj,opt="feed")
      # Mixing parameters plots
      p4 <- plot_proportions(obj,direction="B-F")
      p5 <- plot_proportions(obj,direction="F-B")
      p13 <- plot_mixing(obj)
      # Survival plots
      p6 <- plot_survival(obj,opt = "F")
      p7 <- plot_survival_curve(obj,opt="F")
      # Mortality plots
      p8 <- plot_mortdiff(obj,opt="feed",type="raw")
      p9 <- plot_compare_mort(obj,opt="feed",type="rate")
      p10 <- plot_compare_mort(obj,opt="total",type="rate")
      p11 <- plot_compare_mort(obj,opt="total",type="cumulative")

      ggsave(paste(plot_title,"fixed effects.png"),p0,height=6,width=9)
      ggsave(paste(plot_title,"total abundance.png"),p1,height=4,width=6)
      ggsave(paste(plot_title,"breeding ground abundance.png"),p2,height=6,width=9)
      ggsave(paste(plot_title,"feeding ground abundance.png"),p3,height=6,width=9)
      ggsave(paste(plot_title,"proportions breed to feed.png"),p4,height=5,width=10)
      ggsave(paste(plot_title,"proportions feed to breed.png"),p5,height=5,width=10)
      ggsave(paste(plot_title,"feeding ground survival.png"),p6,height=6,width=8)
      ggsave(paste(plot_title,"environment vs survival.png"),p7,height=6,width=8)
      ggsave(paste(plot_title,"mortality feeding grounds.png"),p8,height=6,width=8)
      ggsave(paste(plot_title,"mortality rate feeding grounds.png"),p9,height=6,width=8)
      ggsave(paste(plot_title,"mortality rate total.png"),p10,height=4,width=6)
      ggsave(paste(plot_title,"cumulative total mortality.png"),p11,height=4,width=6)
      ggsave(paste(plot_title,"mixing.png"),p13,width=10,height=8)
      if(envOpt=="env-K"){
        p12 <- plot_varK(obj)
        ggsave(paste(plot_title,"varying K.png"),p12,height=6,width=8)
      }
      })
    
  }
  
  # =================================================================================================================================
  # Full Bayesian sampling
  # =================================================================================================================================
  if (DoBayes==T){
    require(tmbstan)
    # Now consider mcmc sampling (stan)
    print("trying Bayes")
    options(mc.cores=4)
    mcmcout <- tmbstan(obj=model,iter=1000,refresh=100,warmup = 10, chains=4,cores=1,
                       init = list(best,best,best,best),
                       seed = 1916, thin = 1)
    ## Key information from run. Including the two recommended
    ## convergence diagnostics:
    # print(summary(mcmcout))
    # create file directory if it does not exist yet
    fildir <- here("Diags",subdir)
    if(!dir.exists(fildir)) dir.create(fildir,recursive = T)
    # save Bayes out
    readr::write_rds(mcmcout,paste0(fildir,"/",Code,SensCase,"_Bayes.rds"))
  }
  
  return(best)
}

# ==========================================================================================
# THIS IS WHERE ALL THE ACTION HAPPENS- ACTUALLY RUN THE MODEL
###################################################################################################

# Survival as random effect
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="final/rS test",
            envOpt="rS",
            rvars=c("SFdev"),
            SF=c(0,1,1,1,1,1), WithMirror=0)
# environmentally driven survival
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="final/env-survival test",
            envOpt="env-survival", envVars=c("sst","mld"),
            rvars=c("epsEnv"),
            SF=c(0,1,1,1,1,1), WithMirror=0)
# density-dependent survival only, but with no environmental drivers and no variation in K
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="final/ddOnly test",
            envOpt="ddOnly", UseKPrior = 1, Kmax=60000,
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL)
# environmental index with variable K
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="final/env-K test",
            rvars=c("Kdev"),
            envOpt="env-K",envVars=c("sst","mld"), UseKPrior = 1, Kmax=60000,
            SF=c(0,1,1,1,0,0),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL)

###################################################################################################
# Bayesian Versions
###################################################################################################

# Survival as random effect (with Bayesian sampling)
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="final/rS Bayes test",
            envOpt="rS",
            rvars=c("SFdev"),
            SF=c(0,1,1,1,1,1), WithMirror=0,
            DoBayes = T)
# environment direct (with Bayesian sampling)
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FAenvDirect",envOpt="direct",
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = T,Init=NULL)