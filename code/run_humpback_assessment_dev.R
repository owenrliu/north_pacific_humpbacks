library(RTMB)
library(tidyverse)
library(here)
# ========================================================================================

# all the functions we need
source(here('code',"data_reading_functions_dev.R"))
# source(here('code',"data_reading_functions.R"))
source(here('code',"output_writing_functions.R"))
source(here('code',"plotting_functions.R"))

FullDiag <- F
PlotsDir <- paste0(here('plots'),"/")

DataFileName <- here('data',"Hump.dat")

Kmax <- 60000
Ks <- seq(from=0,to=3*Kmax,length=100)

join1 <- 1/(1+exp(30*(Ks-Kmax)))
join2 <- 1/(1+exp(-30*(Ks-Kmax)))

K1 <- 2*Kmax; Y1 <- log(1/0.0001-1)
K2 <-   Kmax; Y2 <- log(1/0.9999-1)
Prior_slope <- (Y2-Y1)/(K2-K1)
Prior_int <- Y2 - Prior_slope*K2
#cat(Prior_slope,Prior_int,"\n")
Priors <- 1.0/(1+exp(Prior_int+Prior_slope*Ks))

UseKPrior <- 0

# FUNCTION THAT ACTUALLY FITS THE ASSESSMENT MODEL
DoRun <- function(Code,SensCase,StochSopt=1,StrayBase=0,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),SigmaDevS=6,SigmaDevF=0.01,WithMirror=0,Yr1=1970,Yr2=2023,YrSDevs=1995,
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",envOpt="none",
                  AllPlots=F,DoBoot=F,
                  ByCatchFile="BycatchActual_2024_04_24.csv",
                  WghtTotal=1,Idirichlet=1,MaxN=100,seed=19101,
                  BootUse,SetNew=0,Init=NULL,subdir="")
{
  
  cat("Doing ",paste0(Code,"-",SensCase)," using data file ",DataFileName,"\n")
  
  dat <- MakeDataScenario(Code=Code,SensCase=SensCase,StochSopt=StochSopt,StrayBase=StrayBase,
                          DataFileName=DataFileName,Yr1=Yr1,Yr2=Yr2,YrSDevs=YrSDevs,CatchSer=CatchSer,
                          envOpt=envOpt,ByCatchFile=ByCatchFile,AddCV=AddCV,MixWeights=MixWeights,
                          MaxN=MaxN,SF=SF,WithMirror=WithMirror,IAmat=IAmat,SA=SA,SC=SC,
                          TimeLag=TimeLag,DensDepOpt=DensDepOpt,WghtTotal=WghtTotal)
  list2env(dat, envir = environment()) # make the data objects available to this function
  #print(str(data))
  
  # source the correct model code, depending on "envOpt"
  # i.e., choice about how to drive survival with environmental variables
  if(envOpt== 'none') source(here('code','humpback_assessment_model.R'))
  if(envOpt== 'direct') source(here('code','humpback_assessment_model_envDirect.R'))
  if(envOpt== 'index') source(here('code','humpback_assessment_model_envIndex.R'))
  if(envOpt== 'varK') source(here('code','humpback_assessment_model_varK.R'))

  
  ################################################################################
  if (FullDiag==T) print("Calling MakeADM")
  model <- MakeADFun(cmb(f,dat), parameters, map=map,DLL="Hump",silent=T)
  prefit <- model$report()
  #print(rept$neglogL)
  
  # Iterate to get convergence (usually not needed)
  if (FullDiag==T) print("Calling Fit")
  FitBest  <- 1.0e20
  if (SetNew==1) model$par <- Init
  #rept <<- model$report(model$par)
  #AAA
  if (SetNew==0)
  {
    fit <- nlminb(model$par, model$fn, model$gr)
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
  
  # Extract output
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
  
  # Print outs and plots
  
  WriteOut(Code,Abbrev=SensCase,Yr1=Yr1,Yr2=Yr2,BreedNames=BreedNames,FeedNames=FeedNames,
           rept=rept,rep=rep,rep2=rep2,StockDef=StockDef,data=dat,subdir = subdir)
  
  if(AllPlots){
    obj <- list(input=dat,report=rept,sdreport=rep,sdfixed=rep2)
    fildir <- here("plots",subdir)
    if(!dir.exists(fildir)) dir.create(fildir)
    plot_title <- paste0(fildir,"/",Code,"-",SensCase)
    p1 <- plot_abundance(obj=obj,opt = "total")
    p2 <- plot_abundance(obj=obj,opt= "breed")
    p3 <- plot_abundance(obj=obj,opt="feed")
    p4 <- plot_proportions(obj,direction="B-F")
    p5 <- plot_proportions(obj,direction="F-B")
    p6 <- plot_survival(obj,opt="F")
    ggsave(paste(plot_title,"total abundance.png"),p1,height=4,width=6)
    ggsave(paste(plot_title,"breeding ground abundance.png"),p2,height=6,width=8)
    ggsave(paste(plot_title,"feeding ground abundance.png"),p3,height=6,width=8)
    ggsave(paste(plot_title,"proportions breed to feed.png"),p4,height=5,width=10)
    ggsave(paste(plot_title,"proportions feed to breed.png"),p5,height=5,width=10)
    ggsave(paste(plot_title,"feeding ground survival.png"),p6,height=6,width=8)
    png(paste(plot_title,"mixing.png"),width = 1000,height=1000,units='px')
    plot_mixing(obj=obj)
    dev.off()
  }
  
  # Bootstrap
  if (DoBoot==T) Bootstrap(Code=paste0(Code,SensCase),data=dat,parameters,map,rept,Yr1,Yr2,
                           BreedNames,FeedNames,Nboot=500,seed,BootUse,bestOrig=best)
  
  return(best)
}

# ==========================================================================================
# BOOTSTRAPPING FUNCTION

Bootstrap <- function(Code,data,parameters,map,rept,Yr1,Yr2,BreedNames,FeedNames,Nboot=2,seed=19101,BootUse,bestOrig)
{
  FileNameBoot <- paste0("Diags/",Code,".Boot")
  write("",FileNameBoot)
  TrueSurv <- data$SurveyR
  SurveyI <- data$SurveyI
  Nbreed <- length(BreedNames)
  Nfeed <- length(FeedNames)
  Nyear <- Yr2-Yr1+1
  PredSurv <- rept$PredSurv
  PredMixOut <- rept$PredMixOut
  AddV <-rept$AddV
  set.seed(seed)
  seeds <- floor(runif(Nboot+1,1,10000))
  
  # Do bootstrap
  for (Iboot in 1:Nboot)
  {
    set.seed(seeds[Iboot+1])
    #cat("Boot",Iboot,seeds[Iboot+1],"\n")
    if (Iboot > 0)
      for (Irep in 1:BootUse[Iboot])
      {
        # Generate data (Abundance)
        #cat(Iboot,Irep,"\n")
        for (II in 1:length(TrueSurv[,1]))
        {
          SEL <- TrueSurv[II,2]
          if (SurveyI[II,7]>0) SEL <- sqrt(SEL^2+AddV[SurveyI[II,7]])
          dat$SurveyR[II,1] <<- as.numeric(PredSurv[II]*exp(rnorm(1,0,SEL)-SEL^2.0/2.0))
          #new <- PredSurv[II]*exp(rnorm(1,0,SEL)-SEL^2.0/2.0)
          #new <- dat$SurveyR[II,1]
          #print(c(dat$SurveyR[II,1],PredSurv[II],new))
          if (dat$SurveyR[II,1]<=0) print("oops1")
        }
        # Generate data (mixing)
        for (II in 1:2)
          for (Ibreed in 1:Nbreed)
          {
            xx <- rmultinom(1:Nfeed,round(data$ObsMixBtoFO[II,Ibreed]),prob=0.00001+PredMixOut[1,Ibreed,])  
            dat$ObsMixBtoFE[II,Ibreed,] <<- xx/sum(xx)
            if (any(dat$ObsMixBtoFE[II,Ibreed,]<0)) print("oops2")
          }
        for (II in 1:2)
          for (Ifeed in 1:Nfeed)
          {
            xx <- rmultinom(1:Nbreed,round(data$ObsMixFtoBO[II,Ifeed]),prob=0.00001+PredMixOut[2,,Ifeed])  
            dat$ObsMixFtoBE[II,,Ifeed] <<- xx/sum(xx)
            if (any(dat$ObsMixFtoBE[II,,Ifeed]<0)) print("oops3")
          }
        dat$NmixData <- 1000
      } # Irep and Iboot
    
    # Apply RTMB
    model <- MakeADFun(f, parameters, map=map,DLL="Hump",silent=T)
    if (FullDiag==T) print("Calling Fit")
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
    for (iyr in Yr1:(Yr2+1))
    { xx <- c(Iboot,iyr,exp(LogNT[iyr-Yr1+1,1]),LogNT[iyr-Yr1+1,2]); write(xx,FileNameBoot,append=T,ncol=4); }
    
    write("\nBreeding stock abundance",FileNameBoot,append=T)
    write(paste0(unlist(BreedNames)),FileNameBoot,append=T,ncol=Nbreed)  
    BreedOut <- cbind(rep(Iboot,Nyear),Yr1:Yr2)
    for (Ibreed in 1:Nbreed)
    {
      Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
      xx <- NULL
      for (iyr in Yr1:Yr2) xx <- c(xx,exp(LogNb[Ipnt[iyr-Yr1+1],1]));
      BreedOut <-cbind(BreedOut,xx)
    }
    write(t(BreedOut),FileNameBoot,append=T,ncol=Nbreed+2); 
    
    write("\nFeeding stock abundance",FileNameBoot,append=T)
    write(paste0(unlist(FeedNames)),FileNameBoot,append=T,ncol=Nfeed)  
    FeedOut <- cbind(rep(Iboot,Nyear),Yr1:Yr2)
    for (Ifeed in 1:Nfeed)
    {
      Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
      xx <- NULL
      for (iyr in Yr1:Yr2)xx <- c(xx,exp(LogNf[Ipnt[iyr-Yr1+1],1])); 
      FeedOut <-cbind(FeedOut,xx)
    }
    write(t(FeedOut),FileNameBoot,append=T,ncol=Nfeed+2); 
    
    write("\nSurvival by feeding ground",FileNameBoot,append=T)
    write(t(unlist(FeedNames)),FileNameBoot,append=T,sep=',',ncol=Nfeed)
    for (iyr in Yr1:Yr2)
    { xx <- c(Iboot,iyr,rept$SurvOutF[,iyr-Yr1+1]); write(xx,FileNameBoot,append=T,ncol=Nfeed+2); }
    
  } # Boot
  
}
# ==========================================================================================
# THIS IS WHERE ALL THE ACTION HAPPENS- ACTUALLY RUN THE MODEL
###################################################################################################
### TESTING ZONE
# Base model (no environmental drivers, with mirroring of survival
xx <- DoRun("B2F2",SF=c(0,1,0,1,1,1),YrSDevs=2000,WithMirror = 1,SensCase="BC",AllPlots=T,DoBoot=F,Init=NULL,SetNew=0,envOpt="none",subdir="B2F2 base")
# Base model, direct driving of survival
xx <- DoRun("B2F2",SF=c(1,1,1,1,1,1),YrSDevs = 1993,SensCase="BC",AllPlots=T,DoBoot=F,Init=NULL,SetNew=0,envOpt='direct',subdir="B2F2 direct")
# Base model, environmental index of survival
xx <- DoRun("B2F2",SF=c(1,1,1,1,1,1),YrSDevs = 1993,SensCase="BC",AllPlots=T,DoBoot=F,Init=NULL,SetNew=0,envOpt='index',subdir="B2F2 index")

# Variable K model?
xx <- DoRun("B2F2",SF=c(1,1,1,1,1,1),YrSDevs = 2000,SensCase="BC",AllPlots=T,DoBoot=F,Init=NULL,SetNew=0,envOpt='varK',subdir="B2F2 varK")

###################################################################################################
# Select the case for the analysis
Case <- 1

# Base run (for testing)
if (Case==1)
{
  #xx <- DoRun("B1F1",SensCase="BC",AllPlots=F,DoBoot=F,Init=best,SetNew=1)
  xx <- DoRun("B2F2",SF=c(1,1,1,1,1,1),SensCase="BC",AllPlots=T,DoBoot=F,Init=NULL,SetNew=0)
  #save(xx,file="D:\\sav2.sav")
  #DoRun("B2F2",SensCase="BC",AllPlots=T,DoBoot=T)
}

# All base-case models
if (Case==2)
{
  Codes <- c("B1F1","B1F2","B2F1","B2F2")
  for (Icode in 1:4)
  { Code <- Codes[Icode]; xx <- DoRun(Code,SensCase="BC",AllPlots=T,DoBoot=F,SetNew=0);}
}

# Sensitivity tests
if (Case==3)
{
  Codes <- c("B1F1","B1F2","B2F1","B2F2")
  for (Icode in 1:4)
  {
    Code <- Codes[Icode]  
    BestEst <-DoRun(Code,SensCase="BC")
    DoRun(Code,SensCase="S1")
    DoRun(Code,SensCase="S2",MixWeights=c(1,0))
    DoRun(Code,SensCase="S3",MixWeights=c(0,1))
    DoRun(Code,SensCase="S4",SA=0.94)
    DoRun(Code,SensCase="S5",SA=0.98)
    DoRun(Code,SensCase="S6",StochSopt=0)
    DoRun(Code,SensCase="S7",SF=c(0,1,0,1,0,0))
    DoRun(Code,SensCase="S8",StrayBase=0.01)
    DoRun(Code,SensCase="S9",SigmaDevS=4)
    DoRun(Code,SensCase="S10",SigmaDevS=1)
    DoRun(Code,SensCase="S11",WithMirror=0)
    DoRun(Code,SensCase="S12",WithMirror=2,SF=c(0,1,0,1,1,0))
    DoRun(Code,SensCase="S13",CatchSer="L")
    DoRun(Code,SensCase="S14",CatchSer="H")
    DoRun(Code,SensCase="S15",AddCV=F)
    DoRun(Code,SensCase="S16",Yr1=1965)
    DoRun(Code,SensCase="S17",Yr1=1975)
    DoRun(Code,SensCase="S18",WghtTotal=10)
    DoRun(Code,SensCase="S19",DensDepOpt=1)
    DoRun(Code,SensCase="S20",DensDepOpt=2)
    DoRun(Code,SensCase="S21",MaxN=50)
    DoRun(Code,SensCase="S22",MaxN=200)
  }
}

# Boostraps (BootUse is used to trap horrible bootstraps)
if (Case==4)
{
  Codes <- c("B1F1","B1F2","B2F1","B2F2")
  BootUse <- matrix(1,nrow=4,ncol=500)
  for (Icode in 2:2)
  {
    Code <- Codes[Icode]  
    xx <- DoRun(Code,SensCase="BC",AllPlots=T,DoBoot=T,seed=1234,BootUse=BootUse[Icode,])
  }
}