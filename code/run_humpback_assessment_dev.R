library(RTMB)
library(tidyverse)
library(here)
# ========================================================================================

# all the functions we need
source(here('code',"data_reading_functions_dev.R"))
source(here('code',"output_writing_functions.R"))
source(here('code',"plotting_functions.R"))

# Fitting Function
DoRun <- function(Code,SensCase,StochSopt=1,StrayBase=0,Nage=11,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),SigmaDevS=6,SigmaDevF=0.01,WithMirror=0,Yr1=1970,Yr2=2023,YrSDevs=1995, 
                  UseKPrior=1, Kmax=60000,
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",
                  envOpt="none",
                  AllPlots=F,DoBoot=F,BootUse,
                  ByCatchFile="BycatchActual_2024_04_24.csv",
                  DataFileName= here('data',"Hump.dat"),
                  FullDiag=T,
                  # PlotsDir= paste0(here('plots'),"/"),
                  WghtTotal=1,Idirichlet=1,MaxN=100,seed=19101,
                  DoBayes=F,SetNew=0,Init=NULL,subdir="")
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
  if(envOpt== 'index') source(here('code','full age structure','envIndex.R'))
  if(envOpt== 'varK') source(here('code','full age structure','varK.R'))
  if(envOpt== 'randomS') source(here('code','full age structure','base_randomSDevs.R'))
  if(envOpt== 'index_randomS') source(here('code','full age structure','envIndex_randomSDevs.R'))
  # =================================================================================================================================
  # Build TMB Model
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling MakeADM")
  if(envOpt=='randomS'|envOpt=="index_randomS") model <-MakeADFun(cmb(f,dat), parameters, random=c("SFdev"), map=map,DLL="Hump",silent=F)
  if(envOpt=='varK') model <-MakeADFun(cmb(f,dat), parameters, random=c("Kdev"), map=map,DLL="Hump",silent=F)
  else model <- MakeADFun(cmb(f,dat), parameters, map=map,DLL="Hump",silent=T)
  
  # =================================================================================================================================
  # Optimize
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling Fit")
  FitBest  <- 1.0e20
  if (SetNew==1) model$par <- Init
  #rept <<- model$report(model$par)
  #AAA
  if (SetNew==0)
  {
    fit <- nlminb(model$par, model$fn, model$gr,verbose=T)
    while(abs(FitBest-fit$objective)>0.01)
    {
      FitBest <- fit$objective
      fit <- nlminb(fit$par, model$fn, model$gr,verbose=T)
      if (FullDiag==T) print(fit$objective)
    }
  }
  if (SetNew==3)
  {
    fit <- optim(model$par, model$fn, model$gr)
    fit$objective <- fit$value
    while(abs(FitBest-fit$objective)>0.01)
    {
      FitBest <- fit$objective
      fit <- nlminb(fit$par, model$fn, model$gr,verbose=T)
      if (FullDiag==T) print(fit$objective)
    }
  }
  if (SetNew==2)
  {
    fit <- optim(model$par, model$fn,model$gr,method="CG")
    print(fit)
    vv <-model$report(fit$par)
    while(abs(FitBest-fit$value)>0.01)
    {
      FitBest <- fit$value
      print(model$fn(fit$par))
      fit <- optim(fit$par, model$fn, model$gr,method="CG")
      vv <-model$report(fit$par)
      print(model$fn(fit$par))
      fit <- nlminb(fit$par, model$fn, model$gr,verbose=T)
      vv <-model$report(fit$par)
      print(model$fn(fit$par))
      fit$value <- fit$objective
      if (FullDiag==T) print(fit$value)
    }
  }
  
  # Save model
  fildir <- here("Diags",subdir)
  if(!dir.exists(fildir)) dir.create(fildir)
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
  cat("final",rept$neglogL,rept$datalike,rept$Penal,rept$LogLike1,rept$LogLike2a,rept$LogLike2b,"\n")
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
    if(!dir.exists(fildir)) dir.create(fildir)
    plot_title <- paste0(fildir,"/",Code,"-",SensCase)
    p1 <- plot_abundance(obj=obj,opt = "total")
    p2 <- plot_abundance(obj=obj,opt= "breed")
    p3 <- plot_abundance(obj=obj,opt="feed")
    p4 <- plot_proportions(obj,direction="B-F")
    p5 <- plot_proportions(obj,direction="F-B")
    p6 <- plot_survival(obj,opt = "F")
    if(envOpt=="varK"){
      p7 <- plot_varK(obj)
      ggsave(paste(plot_title,"varying K.png"),p7,height=6,width=8)
    }
    ggsave(paste(plot_title,"total abundance.png"),p1,height=4,width=6)
    ggsave(paste(plot_title,"breeding ground abundance.png"),p2,height=6,width=9)
    ggsave(paste(plot_title,"feeding ground abundance.png"),p3,height=6,width=9)
    ggsave(paste(plot_title,"proportions breed to feed.png"),p4,height=5,width=10)
    ggsave(paste(plot_title,"proportions feed to breed.png"),p5,height=5,width=10)
    ggsave(paste(plot_title,"feeding ground survival.png"),p6,height=6,width=8)

    
    png(paste(plot_title,"mixing.png"),width = 1000,height=1000,units='px')
    plot_mixing(obj=obj)
    dev.off()
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
    if(!dir.exists(fildir)) dir.create(fildir)
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
  if(!dir.exists(fildir)) dir.create(fildir)
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
### TESTING ZONES
# Base model (no environmental drivers, with mirroring of survival
xx <- DoRun(Code="B2F1",SensCase="BC",subdir="B2F1 FAbase",envOpt="none",SF=c(0,1,0,1,1,1),YrSDevs=2000,WithMirror = 1,AllPlots=T,DoBoot=F,Init=NULL,SetNew=0)
# Base model, direct driving of survival for all feeding grounds
xx <- DoRun(Code="B2F1",SensCase="BC",subdir="B2F1 FAenvDirect",envOpt="direct",SF=c(1,1,1,1,1,1),YrSDevs=1995,WithMirror = 0,AllPlots=T,DoBoot=F,Init=NULL,SetNew=0)
# Base model, environmental index of survival
xx <- DoRun(Code="B2F1",SensCase="BC",subdir="B2F1 FAenvIndex",envOpt="index",SF=c(0,1,0,1,1,1),YrSDevs=2000,WithMirror = 1, AllPlots=T,DoBoot=F,Init=NULL,SetNew=0)
# Variable K model
xx <- DoRun(Code="B2F1",SensCase="BC",subdir="B2F1 FAvarK",envOpt="varK",SF=c(1,1,1,1,1,1),YrSDevs=2000,WithMirror = 0, AllPlots=T,DoBoot=F,Init=NULL,SetNew=0)
# With survival as random effects
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FArS",envOpt="randomS",
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL,SetNew=0)
# (with Bayesian sampling)
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FArS",envOpt="randomS",
            SF=c(1,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
            DoBayes = T,Init=NULL,SetNew=0)
# environmental index with survival as random effects
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FAenvIndex_rS",envOpt="index_randomS",
            SF=c(0,1,0,1,1,1),WithMirror = 1,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL,SetNew=0)

# environmental index with random K
xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
            SensCase="BC",subdir="B2F1 FAvarK",envOpt="varK",
            SF=c(0,1,0,1,1,1),WithMirror = 1,AllPlots=T,DoBoot=F,
            DoBayes = F,Init=NULL,SetNew=0)
###################################################################################################

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