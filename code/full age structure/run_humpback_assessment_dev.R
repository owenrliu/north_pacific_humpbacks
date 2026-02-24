library(RTMB)
library(tidyverse)
library(here)
# ========================================================================================

# all the functions we need
source(here('code','full age structure',"data_reading_functions_dev.R"))
source(here('code','full age structure',"output_writing_functions.R"))
source(here('code','full age structure',"plotting_functions.R"))

# Fitting Function
DoRun <- function(Code,SensCase,StochSopt=1,StrayBase=0,Nage=11,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),SigmaDevS=6,SigmaDevF=0.01,WithMirror=0,Yr1=1970,Yr2=2023,YrSDevs=2000, 
                  UseKPrior=1, Kmax=60000,
                  rvars=NULL, # which parameters are random
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",
                  envOpt="none", splineK=5, EF= c(0,1,0,1,1,1), envlag=0,
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
  if(envOpt== 'none') source(here('code','full age structure','base.R'))
  if(envOpt== 'direct') source(here('code','full age structure','envDirect.R'))
  # if(envOpt== 'index') source(here('code','full age structure','envIndex.R'))
  if(envOpt== 'ddOnly') source(here('code','full age structure','ddOnly.R'))
  if(envOpt== 'varK') source(here('code','full age structure','varK.R'))
  if(envOpt== 'randomS') source(here('code','full age structure','base_randomSDevs.R'))
  if(envOpt== 'index_randomS') source(here('code','full age structure','envIndex_randomSDevs.R'))
  # =================================================================================================================================
  # Build TMB Model
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling MakeADM")
  model <- MakeADFun(cmb(f,dat), parameters,random=rvars, map=map,DLL="Hump",silent=F)
  # if(envOpt=='randomS'|envOpt=="index_randomS") model <-MakeADFun(cmb(f,dat), parameters, random=c("SFdev"), map=map,DLL="Hump",silent=F)
  # if(envOpt=='direct') model <-MakeADFun(cmb(f,dat), parameters, random=c("epsEnv","SFdev"), map=map,DLL="direct",silent=F)
  # if(envOpt=='varK') model <-MakeADFun(cmb(f,dat), parameters, random=c("Kdev"), map=map,DLL="Hump",silent=F)
  # # if(envOpt=='varK') model <-MakeADFun(cmb(f,dat), parameters, map=map,DLL="Hump",silent=F)
  # else model <- MakeADFun(cmb(f,dat), parameters, map=map,DLL="Hump",silent=T)
  
  # testing
  # model <- MakeADFun(cmb(f,dat),parameters,random=c("SFdev"),map=map,silent=F)
  
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
  stdreport <-sdreport(model)
  rep <- summary(stdreport)
  rep2<- summary(stdreport, "fixed", p.value = TRUE)
  
  # Extract N vectors
  Index <- which(rownames(rep)=="LogNT")
  LogNT<- rep[Index,]
  Index <- which(rownames(rep)=="LogNb")
  LogNb<- rep[Index,]
  Index <- which(rownames(rep)=="LogNf")
  LogNf<- rep[Index,]
  Index <- which(rownames(rep)=="SurvOutB")
  SurvOutB<- rep[Index,]
  Index <- which(rownames(rep)=="SurvOutF")
  SurvOutF<- rep[Index,]
  cat("final nll:",rept$neglogL,"\nTotal data likelihood:",rept$datalike,"\nTotal penalty:",rept$Penal,
      "\nSurvey likelihood:",rept$LogLike1,"\nMixing BtF likelihood:",rept$LogLike2a,"\nMixing FtB likelihood",rept$LogLike2b,"\n")
  #print(exp(LogNT))
  
  # Stuff to pass to other routines
  StockDef <- NULL
  StockDef$StochSopt <- StochSopt
  StockDef$StrayBase <- StrayBase
  StockDef$IAmat <- IAmat
  StockDef$SA <- SA
  StockDef$TimeLag <- TimeLag
  StockDef$DensDepOpt <- DensDepOpt
  StockDef$SF <- SF
  StockDef$SigmaDevS <- SigmaDevS
  StockDef$SigmaDevF <- SigmaDevF
  StockDef$WithMirror <- WithMirror
  StockDef$AddCV <- AddCV
  StockDef$MixWeights <- MixWeights
  
  # =================================================================================================================================
  # Write and Plot
  # =================================================================================================================================
  
  WriteOut(Code,Abbrev=SensCase,Yr1=Yr1,Yr2=Yr2,BreedNames=BreedNames,FeedNames=FeedNames,
           rept=rept,rep=rep,rep2=rep2,StockDef=StockDef,data=dat,subdir = subdir)
  
  if(AllPlots){
    obj <- list(Code=Code,SensCase=SensCase,input=dat,report=rept,sdreport=rep,sdfixed=rep2)
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
      if(envOpt=="varK"){
        p12 <- plot_varK(obj)
        ggsave(paste(plot_title,"varying K.png"),p12,height=6,width=8)
      }
      })
    
  }
  
  # =================================================================================================================================
  # Bootstrap
  # =================================================================================================================================
  
  # Bootstrap
  BootUse <- matrix(1,nrow=4,ncol=500)
  if (DoBoot==T) Bootstrap(Code=paste0(Code,SensCase),dat=dat,parameters,map,rept,Yr1,Yr2,
                           BreedNames,FeedNames,Nboot=10,seed,BootUse,bestOrig=best,envOpt=envOpt,subdir=subdir)
  
  # =================================================================================================================================
  # Full Bayesian sampling
  # =================================================================================================================================
  if (DoBayes==T){
    require(tmbstan)
    # Now consider mcmc sampling (stan)
    print("trying Bayes")
    options(mc.cores=4)
    mcmcout <- tmbstan(obj=model,iter=10000,refresh=100,warmup = 10000/10, chains=4,cores=1,
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
# BOOTSTRAPPING FUNCTION
# ==========================================================================================

Bootstrap <- function(Code,dat,parameters,map,rept,Yr1,Yr2,BreedNames,FeedNames,Nboot=2,seed=19101,BootUse,bestOrig,envOpt,subdir=""){
  
  fildir <- here("Diags",subdir)
  if(!dir.exists(fildir)) dir.create(fildir,recursive=T)
  FileNameBoot <- paste0(fildir,"/",Code,".Boot")
  
  write("",FileNameBoot)
  TrueSurv <- dat$SurveyR
  SurveyI <- dat$SurveyI
  Nbreed <- length(BreedNames)
  Nfeed <- length(FeedNames)
  Nyear <- Yr2-Yr1+1
  PredSurv <- rept$PredSurv
  PredMixOut <- rept$PredMixOut
  AddV <-rept$AddV
  set.seed(seed)
  seeds <- floor(runif(Nboot+1,1,10000))
  Nyr <- (Yr2+1)-Yr1+1
  
  # Do bootstrap
  for (Iboot in 1:Nboot){
    # Bootstrapped dataset
    datboot <- dat
    set.seed(seeds[Iboot+1])
    #cat("Boot",Iboot,seeds[Iboot+1],"\n")
    if (Iboot > 1)
      for (Irep in 1:BootUse[Iboot])
      {
        # Generate data (Abundance)
        cat(Iboot,Irep,"\n")
        for (II in 1:length(TrueSurv[,1]))
        {
          SEL <- TrueSurv[II,2]
          if (SurveyI[II,7]>0) SEL <- sqrt(SEL^2+AddV[SurveyI[II,7]])
          datboot$SurveyR[II,1] <- as.numeric(PredSurv[II]*exp(rnorm(1,0,SEL)-SEL^2.0/2.0))
          if (datboot$SurveyR[II,1]<=0) print("oops1")
        }
        # Generate data (mixing)
        for (II in 1:2)
          for (Ibreed in 1:Nbreed)
          {
            xx <- rmultinom(1:Nfeed,round(datboot$ObsMixBtoFO[II,Ibreed]),prob=0.00001+PredMixOut[1,Ibreed,])  
            datboot$ObsMixBtoFE[II,Ibreed,] <- xx/sum(xx)
            if (any(datboot$ObsMixBtoFE[II,Ibreed,]<0)) print("oops2")
          }
        for (II in 1:2)
          for (Ifeed in 1:Nfeed)
          {
            xx <- rmultinom(1:Nbreed,round(datboot$ObsMixFtoBO[II,Ifeed]),prob=0.00001+PredMixOut[2,,Ifeed])  
            datboot$ObsMixFtoBE[II,,Ifeed] <- xx/sum(xx)
            if (any(datboot$ObsMixFtoBE[II,,Ifeed]<0)) print("oops3")
          }
        datboot$NmixData <- 1000
      } # Irep and Iboot
    
    # Apply RTMB
    
    if(envOpt=='randomS'|envOpt=="index_randomS") model <-MakeADFun(cmb(f,datboot), parameters, random=c("SFdev"), map=map,DLL="Hump",silent=F)
    else model <- MakeADFun(cmb(f,datboot), parameters, map=map,DLL="Hump",silent=T)
    
    #model$par <- bestOrig
    FitBest  <- 1.0e20
    fit <- nlminb(model$par, model$fn, model$gr)
    while(abs(FitBest-fit$objective)>0.01)
    {
      FitBest <- fit$objective
      if (!is.na(model$fn(fit$par)))
        fit <- nlminb(fit$par, model$fn, model$gr) 
      else
        cat("Quitting",Iboot,"\n")
      #print(fit$objective)
    }
    best <- model$env$last.par.best
    stdreport <-sdreport(model)
    rep <- summary(stdreport)
    rep2<- summary(stdreport, "fixed", p.value = TRUE)
    rept <- model$report()
    cat("final",Iboot,rept$neglogL,rept$datalike,rept$Penal,rept$LogLike1,rept$LogLike2a,rept$LogLike2b,"\n")
    
    Index <- which(rownames(rep)=="LogNT") 
    LogNT<- rep[Index,]
    Index <- which(rownames(rep)=="LogNb")
    LogNb<- rep[Index,]
    Index <- which(rownames(rep)=="LogNf")
    LogNf<- rep[Index,]
    BreedK <- rept$BreedK
    FeedK <- rept$FeedK
    
    write("\nParameters",FileNameBoot,append=T)
    for (II in 1:length(best))
      write(best[II],FileNameBoot,append=T,ncol=4)
    
    write("\nTotal abundance: Estimate SD_log",FileNameBoot,append=T)
    for (Icomp in 1:3)
      for (iyr in Yr1:(Yr2+1))
      { xx <- c(Iboot,Icomp,iyr,exp(LogNT[3*(iyr-Yr1+1-1)+Icomp,1]),LogNT[3*(iyr-Yr1+1-1)+Icomp,2]); write(xx,FileNameBoot,append=T,ncol=5); }
    
    write("\nBreeding stock abundance",FileNameBoot,append=T)
    write(paste0(unlist(BreedNames)," ",unlist(BreedNames)," ",unlist(BreedNames)),FileNameBoot,append=T,ncol=Nbreed)  
    BreedOut <- cbind(rep(Iboot,Nyear),Yr1:Yr2)
    for (Ibreed in 1:Nbreed)
      for (Icomp in 1:3)
      {
        xx <- NULL
        for (iyr in Yr1:Yr2) 
        { Jpnt <- Nbreed*3*(iyr-Yr1+1-1)+(Ibreed-1)*3+Icomp; xx <- c(xx,exp(LogNb[Jpnt,1])); }
        BreedOut <-cbind(BreedOut,xx)
      }
    write(t(BreedOut),FileNameBoot,append=T,ncol=3*Nbreed+2); 
    
    write("\nFeeding stock abundance",FileNameBoot,append=T)
    write(paste0(unlist(FeedNames)," ",unlist(FeedNames)," ",unlist(FeedNames)),FileNameBoot,append=T,ncol=Nfeed)  
    FeedOut <- cbind(rep(Iboot,Nyear),Yr1:Yr2)
    for (Ifeed in 1:Nfeed)
      for (Icomp in 1:3)
      {
        xx <- NULL
        for (iyr in Yr1:Yr2) 
        { Jpnt <- Nfeed*3*(iyr-Yr1+1-1)+(Ifeed-1)*3+Icomp; xx <- c(xx,exp(LogNf[Jpnt,1])); }
        FeedOut <-cbind(FeedOut,xx)
      }
    write(t(FeedOut),FileNameBoot,append=T,ncol=3*Nfeed+2); 
    
    write("\nSurvival by feeding ground",FileNameBoot,append=T)
    write(t(unlist(FeedNames)),FileNameBoot,append=T,sep=',',ncol=Nfeed)
    for (iyr in Yr1:Yr2)
    { xx <- c(Iboot,iyr,rept$SurvOutF[,iyr-Yr1+1]); write(xx,FileNameBoot,append=T,ncol=Nfeed+2); }
    
  } # NBoot loop
  
}

# ==========================================================================================
# THIS IS WHERE ALL THE ACTION HAPPENS- ACTUALLY RUN THE MODEL
###################################################################################################
### TESTING ZONE
# Survival as random effect
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 rS/test",
            envOpt="randomS",
            rvars=c("SFdev"),
            SF=c(0,1,1,1,1,1), WithMirror=0,
            EF=c(0,0,0,0,0,0))
# environment direct
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 direct/no RUS_WAL",
            envOpt="direct",
            rvars=c("epsEnv"),
            SF=c(0,1,1,1,1,1), WithMirror=0,
            EF=c(0,1,1,1,1,1))
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 direct/no RUS_WAL lag1",
            envOpt="direct",
            rvars=c("epsEnv"), envlag = 1, # environment vars lagged 1 year
            SF=c(0,1,1,1,1,1), WithMirror=0,
            EF=c(0,1,1,1,1,1))
# density-dependent survival only, but with no environmental drivers and no variation in K
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FAdd",
            envOpt="ddOnly", UseKPrior = 1, Kmax=60000,
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL)
# environmental index with variable K
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 varK",
            rvars=c("Kdev"),
            envOpt="varK", UseKPrior = 1, Kmax=60000,
            SF=c(0,1,1,1,0,0),
            EF=c(0,1,1,1,0,0),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL)
###################################################################################################
# Bayesian Versions
###################################################################################################

# Survival as random effect (with Bayesian sampling)
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FArS",envOpt="randomS",
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = T,Init=NULL)
# environment direct (with Bayesian sampling)
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FAenvDirect",envOpt="direct",
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = T,Init=NULL)

# Boostraps (BootUse is used to trap horrible bootstraps)
# if (Case==4)
# {
#   Codes <- c("B1F1","B1F2","B2F1","B2F2")
#   BootUse <- matrix(1,nrow=4,ncol=500)
#   for (Icode in 2:2)
#   {
#     Code <- Codes[Icode]  
#     xx <- DoRun(Code,SensCase="BC",AllPlots=T,DoBoot=T,seed=1234,BootUse=BootUse[Icode,])
#   }
# }