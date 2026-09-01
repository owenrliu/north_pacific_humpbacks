library(RTMB)

# ========================================================================================

setwd("C:\\Research\\Iwc25\\NP_Humpbacks\\")
source("HumpFigs2.R")
source("ReadData.R")
source("WriteOut.R")

FullDiag <- F
PlotsDir <- "Plots/"

DataFileName <- "Hump.dat"

sumSQ <- function(x){
  TheMean <- sum(x)/length(x);
  Output <- 0;
  for (ii in 1:length(x)) Output <- Output + (x[ii]-TheMean)^2.0;
  return(Output);
}

# ========================================================================================
TheSLA <- function(AreaType,AreaNum,Need,Catches,SurveyData)
 {
  Use <- SurveyData[,5] == 1
  SurveyData <- SurveyData[Use,]
  Use <- which(SurveyData[,2]==max(SurveyData[,2]))
  Est <- SurveyData[Use,8]
  MaxStrikes <- Est*0.01
  cat(AreaType,AreaNum,Need,Est,MaxStrikes,"\n")
  SLA <- min(MaxStrikes,Need)
  return(SLA)  
 }

# ========================================================================================
RunModelProj <- function(data,modelpars,SLASpecs,DataSpecs,Nproj=20)
 {
  # extract the model parameters
  rval <- modelpars$rval
  logK <- modelpars$logK
  logBK <- modelpars$logBK
  MixPars <-modelpars$MixPars
  AddV <- modelpars$AddV
  SFdev <- modelpars$SFdev
  InfluxP <- modelpars$InfluxP
  Sigma_SBdev <- modelpars$Sigma_SBdev
  Sigma_SFdev <- modelpars$Sigma_SFdev
  Sigma_FBdev <- modelpars$Sigma_FBdev
  Sigma_FFdev <- modelpars$Sigma_FFdev
  inert_par <- modelpars$inert_par

  # Extract all data that are needed
  Yr1 <- data$Yr1
  Yr2 <- data$Yr2
  Nbreed <- data$Nbreed
  Nfeed <- data$Nfeed
  IAmat <- data$IAmat 
  SA <- data$SA
  SC <- data$SC
  TimeLag <- data$TimeLag
  MixI <- data$MixI
  SBdevEst <- data$SBdevEst*0
  SBdevMat <- data$SBdevMat 
  FBdevEst <- data$FBdevEst
  FBdevMat <- data$FBdevMat 
  SFdevEst <- data$SFdevEst
  SFdevMat <- data$SFdevMat 
  FFdevEst <- data$FFdevEst
  FFdevMat <- data$FFdevMat 
  Nmirror <- data$Nmirror
  Mirror <- data$Mirror
  StrayBase <- data$StrayBase 
  DensDepOpt <- data$DensDepOpt
  StochSopt <- data$StochSopt
   
  # Now find limits 
  NyrHist <- Yr2-Yr1+1;                                             # Number of years of historical data
  Nyr <- Yr2-Yr1+1+Nproj;                                           # Total number of years
  MixPar <- length(MixPars);                                        # Number of mixing parameters
   
  ThreshPop = 1.0e-20;                                              # Minimum population size
  Influx = 1.0/(1+exp(InfluxP));                                    # Transfer influx
  
  # Create the data file
  TheSurveyData <- cbind(data$SurveyI,data$SurveyR)
  #print(summary(TheSurveyData))
  
  # ========================================================================================================================
  # Create a mixing matrix
  # ========================================================================================================================
  Mix <-matrix(0,nrow=Nbreed,ncol=Nfeed);
  Icnt <- 1; 
  for (Ibreed in 1:Nbreed)
  {
    Total <- 0;
    for (Ifeed in 1:Nfeed)
      if (MixI[Ibreed,Ifeed] > 0)
      { Mix[Ibreed,Ifeed] <- exp(MixPars[Icnt]); Icnt <- Icnt + 1; Total <- Total+ Mix[Ibreed,Ifeed]; }
     else
      { Mix[Ibreed,Ifeed] <- abs(1.0*MixI[Ibreed,Ifeed]); Total <- Total + Mix[Ibreed,Ifeed]; }
    for (Ifeed in 1:Nfeed) Mix[Ibreed,Ifeed] <- Mix[Ibreed,Ifeed]/Total;
  }
  
  # ========================================================================================================================
  # Population projecton
  # ========================================================================================================================
    
  
  # Basic demographics
  Amat <- IAmat*1.0;
  # Zerbini et al. (2010): Mar Biol. 157: 1225-1236. Note that Amat is age-at-maturity (age-at-first parturition less 1)
  #fmax <- 2*((1.0+rval)^(Amat+1.0)-SA*(1.0+rval)^Amat)
  #f0 <- 2*(1.0-SA)
  #fmax <- 2*((1.0+rval)^(Amat+1.0)-SA*(1.0+rval)^Amat)/(SC*SA^(Amat))
  #f0 <- 2*(1.0-SA)/(SC*SA^(Amat))
  # Corrected
  fmax <- 2*(exp(rval*(Amat+1.0))-SA*exp(rval*Amat))/(SC*SA^(Amat))
  #fmax <- fmax/(1.0+exp(30.0*(fmax-0.99)))+0.99/(1.0+exp(-30.0*(fmax-0.99))) 
  f0 <- 2*(1.0-SA)/(SC*SA^(Amat))
  # Fecundity at carrying capacity
  ParA <- (fmax-f0)/f0;                                          # Resilience parameter
  #print(c(rval,fmax,f0,ParA))
  
  # Extract the SBdevs
  SBdevPnt <- 1;
  SBdevYr <- matrix(0,nrow=Nbreed,ncol=Nyr);
  for (Ibreed in 1:Nbreed)
   if (SBdevEst[Ibreed] > 0)
    for (Year in SBdevMat[Ibreed,1]:SBdevMat[Ibreed,2])
      { SBdevYr[Ibreed,Year+1] <- SBdev[SBdevPnt]; SBdevPnt <- SBdevPnt + 1; }
  
  # Extract the SFdevs
  SFdevPnt <- 1;
  SFdevYr <- matrix(0,nrow=Nfeed,ncol=Nyr);
  for (Ifeed in 1:Nfeed)
   if (SFdevEst[Ifeed] > 0)
    for (Year in SFdevMat[Ifeed,1]:SFdevMat[Ifeed,2])
      { SFdevYr[Ifeed,Year+1] <- SFdev[SFdevPnt]; SFdevPnt <- SFdevPnt + 1; }
  
  # Mirror as needed
  if (Nmirror>0)
   for (Im in 1:Nmirror)
    {
     Imi = Mirror[Im,1]+1; Omi = Mirror[Im,2]+1;
     for (Year in 1:Nyr) SFdevYr[Omi,Year] <- SFdevYr[Imi,Year];
    }
  
  # Extract the FBdevs
  FBdevPnt <- 1;
  FBdevYr <- matrix(0,nrow=Nbreed,ncol=Nyr+1);
  for (Ibreed in 1:Nbreed)
   if (FBdevEst[Ibreed] > 0)
    for (Year in FBdevMat[Ibreed,1]:FBdevMat[Ibreed,2])
     { FBdevYr[Ibreed,Year] <- FBdev[FBdevPnt]; FBdevPnt <- FBdevPnt + 1; }
  
  # Set up basis for A and K (not used in the base mode)
  ParAV <- matrix(0,nrow=Nfeed,ncol=Nyr+1);
  MultK <- matrix(0,nrow=Nfeed,ncol=Nyr+1);
  for (Year in 1:(Nyr+1))
    for (Ifeed in 1:Nfeed)
    { ParAV[Ifeed,Year] <- ParA; MultK[Ifeed,Year] <- 1; }
  
  # Set up K by herd, breeding stock, and feeding gound
  BreedK <- rep(0,Nbreed); FeedK <- rep(0,Nfeed); 
  NN <- array(0,dim=c(Nbreed,Nfeed,Nyr+1)); NNK <- matrix(0,nrow=Nbreed,ncol=Nfeed)
  Nb <- matrix(0,nrow=Nbreed,ncol=Nyr+1);  Nf <- matrix(0,nrow=Nfeed,ncol=Nyr+1)
  for (Ibreed in 1:Nbreed)
   { BreedK[Ibreed] <- exp(logK[Ibreed]); for (Ifeed in 1:Nfeed) NNK[Ibreed,Ifeed] <- BreedK[Ibreed]*Mix[Ibreed,Ifeed]; }
  for (Ifeed in 1:Nfeed)
   { for (Ibreed in 1:Nbreed) FeedK[Ifeed] <- FeedK[Ifeed] + NNK[Ibreed,Ifeed]; }
  
  # Temporary varibales
  NbS <- matrix(0,nrow=Nbreed,ncol=Nyr+1); NfS <- matrix(0,nrow=Nfeed,ncol=Nyr+1)
  NNS <- array(0,dim=c(Nbreed,Nfeed,Nyr+1)); FecOut <-matrix(0,nrow=Nfeed,ncol=Nyr)
  PropMixEst <- array(0,dim=c(Nfeed,Nbreed,Nyr+1))
  SurvOutB <- matrix(0,nrow=Nbreed,ncol=Nyr)
  SurvOutF <- matrix(0,nrow=Nfeed,ncol=Nyr)
  SurvTot <- array(0,dim=c(Nbreed,Nfeed,Nyr));
  MortDiff <- array(0,dim=c(Nbreed,Nfeed,Nyr));
  Temp1 <- rep(0,Nbreed)
  
  # Set up Initial year (note that the same depletion applied to all feeding grounds within a breeding stock)
  for (Ibreed in 1:Nbreed)
   {
    Nb[Ibreed,1] <- exp(logK[Ibreed])*1.0/(1+exp(logBK[Ibreed]));
    for (Ifeed in 1:Nfeed) NN[Ibreed,Ifeed,1] <- Nb[Ibreed,1]*Mix[Ibreed,Ifeed];
   }
  
  # Now project forward
  CatchB2 <- matrix(0,Nyr,Nbreed)
  CatchB2[1:NyrHist,] <- CatchB
  CatchF2 <- matrix(0,Nyr,Nfeed)
  CatchF2[1:NyrHist,] <- CatchF
  for (Year in 1:Nyr)
   {
    # Call the SLA (As needed)
    if (Year>NyrHist)
     {
      for (II in 1:length(SLASpecs[,1]))
       {
        SLA_area <-SLASpecs[II,2]; SLA_need <- SLASpecs[II,3]
        if (SLASpecs[II,1]==1)  #  Breeding group SLA
         {
          CatchPass <- CatchB2[1:(Year-1),SLA_area]
          Use <- TheSurveyData[,3] == 2 & TheSurveyData[,4] == SLA_area
          SurveyPass <- TheSurveyData[Use,]
          CatchB2[Year,SLA_area] <- TheSLA("Breed",SLA_area,SLA_need,CatchPass,SurveyPass)
         }
        if (SLASpecs[II,1]==2) # Feeding ground SLA
         {
          CatchPass <- CatchF2[1:(Year-1),SLA_area]
          Use <- TheSurveyData[,3] == 3 & TheSurveyData[,4] == SLA_area
          SurveyPass <- TheSurveyData[Use,]
          CatchF2[Year,SLA_area] <- TheSLA("Feed",SLA_area,SLA_need,CatchPass,SurveyPass)
        }
       }
     }
     
    #cat("Year = ",Year,Nyr,"\n")
    # Save breeding ground numbers (start of the year)
    for (Ibreed in 1:Nbreed)
      for (Ifeed in 1:Nfeed)
        NNS[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year];
    
    # Save the breeding ground numbers
    for (Ibreed in 1:Nbreed)
      for (Ifeed in 1:Nfeed)
        NbS[Ibreed,Year] <- NbS[Ibreed,Year] + NN[Ibreed,Ifeed,Year];
    
    # Remove breeding ground catches
    for (Ibreed in 1:Nbreed)
    {
      Temp1[Ibreed] <- Nb[Ibreed,Year] - CatchB2[Year,Ibreed];
      # Trick to avoid negative population sizes
      MultC <- 1.0/(1.0+exp(-30.0*(Temp1[Ibreed]-ThreshPop)))
      Temp1[Ibreed] <- ThreshPop + MultC * Temp1[Ibreed];
      # Allocate breeding ground catches to "herd" (Eqn B.3)
      for (Ifeed in 1:Nfeed) NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year]*Temp1[Ibreed]/Nb[Ibreed,Year];
    }
    
    # Allow for straying before feeding ground catches (Eqn B.7)
    for (Ibreed in 1:Nbreed)
    {
      Nstray <- rep(0,Nfeed)
      for (Ifeed in 1:Nfeed)
        if (MixI[Ibreed,Ifeed] != 0)
        {
          if (Ifeed > 1)
            if(MixI[Ibreed,Ifeed-1] !=0)
            {
              Depl1 <- NN[Ibreed,Ifeed,Year]/NNK[Ibreed,Ifeed];
              Depl2 <- NN[Ibreed,Ifeed-1,Year]/NNK[Ibreed,Ifeed-1];
              DeplRat <- Depl1/Depl2;
              
              Ver1 <- 1.0/(1.0+exp(-100.0*(DeplRat-1.0)));
              Ver2 <- 1.0/(1.0+exp(100.0*(DeplRat-2.0)));
              Ver3 <- 1.0/(1.0+exp(-100.0*(DeplRat-2.0)));
              Stray <- StrayBase*(DeplRat-1.0)*Ver1*Ver2+StrayBase*Ver3;
              Nstray[Ifeed] <- Nstray[Ifeed] - Stray*NN[Ibreed,Ifeed,Year];
              Nstray[Ifeed-1] <- Nstray[Ifeed-1] + Stray*NN[Ibreed,Ifeed,Year];
            }
          if (Ifeed < Nfeed)
            if(MixI[Ibreed,Ifeed+1] !=0)
            {
              Depl1 <- NN[Ibreed,Ifeed,Year]/NNK[Ibreed,Ifeed];
              Depl2 <- NN[Ibreed,Ifeed+1,Year]/NNK[Ibreed,Ifeed+1];
              DeplRat <- Depl1/Depl2;
              
              Ver1 = 1.0/(1.0+exp(-100.0*(DeplRat-1.0)));
              Ver2 = 1.0/(1.0+exp(100.0*(DeplRat-2.0)));
              Ver3 = 1.0/(1.0+exp(-100.0*(DeplRat-2.0)));
              Stray = StrayBase*(DeplRat-1.0)*Ver1*Ver2+StrayBase*Ver3;
              Nstray[Ifeed] <- Nstray[Ifeed] - Stray*NN[Ibreed,Ifeed,Year];
              Nstray[Ifeed+1] <- Nstray[Ifeed+1] + Stray*NN[Ibreed,Ifeed,Year];
            }
        }
      for (Ifeed in 1:Nfeed) NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year] + Nstray[Ifeed];
    }  # Each breeding stock
    
    # Save the feeding ground numbers
    for (Ifeed in 1:Nfeed)
      for (Ibreed in 1:Nbreed) NfS[Ifeed,Year] <- NfS[Ifeed,Year] + NN[Ibreed,Ifeed,Year];

    # Remove feeding ground numbers (Eqn B.2)
    CBreed <- rep(0,Nbreed);
    for (Ifeed in 1:Nfeed)
    {
      Nf[Ifeed,Year] <- 0; 
      for (Ibreed in 1:Nbreed) Nf[Ifeed,Year] <- Nf[Ifeed,Year] + NN[Ibreed,Ifeed,Year];
      for (Ibreed in 1:Nbreed) PropMixEst[Ifeed,Ibreed,Year] <- NN[Ibreed,Ifeed,Year]/Nf[Ifeed,Year];
      for (Ibreed in 1:Nbreed) CBreed[Ibreed] <- CBreed[Ibreed] + PropMixEst[Ifeed,Ibreed,Year]*CatchF2[Year,Ifeed];
      for (Ibreed in 1:Nbreed)
        if (MixI[Ibreed,Ifeed] != 0)
        {
          NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year] - PropMixEst[Ifeed,Ibreed,Year]*CatchF2[Year,Ifeed];
          # Trick to avoid negative population sizes
          MultC <- 1.0/(1.0+exp(-30.0*(NN[Ibreed,Ifeed,Year]-ThreshPop)))
          NN[Ibreed,Ifeed,Year] <- ThreshPop + MultC * NN[Ibreed,Ifeed,Year];
          NN[Ibreed,Ifeed,Year] <- sqrt(NN[Ibreed,Ifeed,Year]*NN[Ibreed,Ifeed,Year])
        }
      Nf[Ifeed,Year] <- 0;  for (Ibreed in 1:Nbreed) Nf[Ifeed,Year] <- Nf[Ifeed,Year] + NN[Ibreed,Ifeed,Year];
    }
    
    # Influx (part 1) [Note used for this assessment]
    for (Ifeed in 1:Nfeed) MultK[Ifeed,Year+1] <- MultK[Ifeed,Year];
    if (Year==2014-Yr1+1)
    {
      MultK[4,Year+1] <- MultK[4,Year]*(1.0-Influx);
      MultK[5,Year+1] <- MultK[5,Year]*(1.0-Influx);
      MultK[6,Year+1] <- (MultK[6,Year]*Nf[6,1] + MultK[4,Year]*Influx*Nf[4,1]+MultK[5,Year]*Influx*Nf[5,1])/(MultK[6,Year]*Nf[6,1]);
      for (Ibreed in 1:Nbreed)
      {
        NN[Ibreed,6,Year] <- NN[Ibreed,6,Year] + NN[Ibreed,4,Year]*Influx + NN[Ibreed,5,Year]*Influx;
        NN[Ibreed,4,Year] <- NN[Ibreed,4,Year]*(1.0-Influx);
        NN[Ibreed,5,Year] <- NN[Ibreed,5,Year]*(1.0-Influx);
      }
    }
    
    # Update feeding ground numbers
    for (Ifeed in 1:Nfeed)
    { Nf[Ifeed,Year] <- 0; for (Ibreed in 1:Nbreed) Nf[Ifeed,Year] <- Nf[Ifeed,Year] + NN[Ibreed,Ifeed,Year]; }
    
    # Density-dependence is on feeding ground numbers
    for (Ibreed in 1:Nbreed)
    {
      # Old
      #Ioffset <- Year-IAmat+1-TimeLag; if (Ioffset < 1) Ioffset <- 1;
      # NEW
      Ioffset <- Year-(IAmat+1)+1-TimeLag; if (Ioffset < 1) Ioffset <- 1;
      
      # Feeding ground density-dependence but recruitment related to recent (Eqn B.1)
      if (DensDepOpt==0)
        for (Ifeed in 1:Nfeed)
        {
          # Survival (Eqn B.6)
          LogitSA <- log(1.0/SA-1.0);
          SurvOutB[Ibreed,Year] <- 1.0/(1.0+exp(LogitSA+SBdevYr[Ibreed,Year]*Sigma_SBdev));
          SurvOutF[Ifeed,Year] <- 1.0/(1.0+exp(LogitSA+SFdevYr[Ifeed,Year]*Sigma_SFdev));
          if (StochSopt==0) SAuse <- SurvOutB[Ibreed,Year];
          if (StochSopt==1) SAuse <- SurvOutF[Ifeed,Year];
          # Eqn B.4
          Depl <- Nf[Ifeed,Ioffset]/(FeedK[Ifeed]*MultK[Ifeed,Ioffset]);
          #Depl <- Depl/(1.0+exp(30.0*(Depl-0.99)))+0.99/(1.0+exp(-30.0*(Depl-0.99))) 
          Term1 <- 1.0 - Depl;
          Term1 <- f0*(1.0+ParAV[Ifeed,Ioffset]*Term1);
          if (Term1 < 0.0001) Term1 <- 0.0001
          Term1 <- 0.0001+(Term1-0.0001)/(1+exp(-30.0*Term1)) 
          LogitFec <- log(1.0/Term1-1.0);
          Term1 <- 1.0/(1.0+exp(LogitFec+FBdevYr[Ibreed,Year]*Sigma_FBdev));
          FecOut[Ifeed,Year] <- Term1;
          SurvTot[Ibreed,Ifeed,Year] <- SAuse*NN[Ibreed,Ifeed,Year]/NNS[Ibreed,Ifeed,Year];
          MortDiff[Ibreed,Ifeed,Year] <- (SAuse-SA)*NN[Ibreed,Ifeed,Year];
          
          # Eqn B.1
          SurvJuv <- SC/SA
          for (Iage in 0:(Amat))
          {
            if (StochSopt==0) SurvJuv <- SurvJuv*SurvOutB[Ibreed,max(1,Year-Iage)]
            if (StochSopt==1) SurvJuv <- SurvJuv*SurvOutF[Ifeed,max(1,Year-Iage)]
            #SurvJuv <- SurvJuv * SA
          }
          NN[Ibreed,Ifeed,Year+1] <- SAuse*NN[Ibreed,Ifeed,Year] + 0.5*Term1*SurvJuv*NN[Ibreed,Ifeed,Ioffset];
          ParAV[Ifeed,Year+1] <- ParAV[Ifeed,Ioffset]*exp(inert_par * log(FeedK[Ifeed]/Nf[Ifeed,Ioffset]));
        }
      
      # Feeding ground density-dependence but recruitment related to unfished (Eqn B9.a)
      if (DensDepOpt==1)
        for (Ifeed in 1:Nfeed)
        {
          # Survival (Eqn B.6)
          LogitSA <- log(1.0/SA-1.0);
          SurvOutB[Ibreed,Year] <- 1.0/(1.0+exp(LogitSA+SBdevYr[Ibreed,Year]*Sigma_SBdev));
          SurvOutF[Ifeed,Year] <- 1.0/(1.0+exp(LogitSA+SFdevYr[Ifeed,Year]*Sigma_SFdev));
          if (StochSopt==0) SAuse <- SurvOutB[Ibreed,Year];
          if (StochSopt==1) SAuse <- SurvOutF[Ifeed,Year];
          # Eqn B.4
          Term1 <- 1.0 - Nf[Ifeed,Ioffset]/(FeedK[Ifeed]*MultK[Ifeed,Ioffset]);
          Term1 <- f0*(1.0+ParAV[Ifeed,Ioffset]*Term1);
          Term1 <- 0.000001+(Term1-0.000001)/(1+exp(-30*Term1)) 
          LogitFec <- log(1.0/Term1-1.0);
          Term1 <- 1.0/(1.0+exp(LogitFec+FBdevYr[Ibreed,Year]*Sigma_FBdev));
          FecOut[Ifeed,Year] <- Term1;
          SurvTot[Ibreed,Ifeed,Year] <- SAuse*NN[Ibreed,Ifeed,Year]/NNS[Ibreed,Ifeed,Year];
          MortDiff[Ibreed,Ifeed,Year] <- (SAuse-SA)*NN[Ibreed,Ifeed,Year];
          # Eqn B.9a
          SurvJuv <- SC/SA
          for (Iage in 0:(Amat))
          {
            if (StochSopt==0) SurvJuv <- SurvJuv*SurvOutB[Ibreed,max(1,Year-Iage)]
            if (StochSopt==1) SurvJuv <- SurvJuv*SurvOutF[Ifeed,max(1,Year-Iage)]
          }
          NN[Ibreed,Ifeed,Year+1] <- SAuse*NN[Ibreed,Ifeed,Year] + 0.5*Term1*SurvJuv*Nb[Ibreed,Ioffset]*Mix[Ibreed,Ifeed];
          ParAV[Ifeed,Year+1] <- ParAV[Ifeed,Ioffset]*exp(inert_par * log(FeedK[Ifeed]/Nf[Ifeed,Ioffset]));
        }
      
      # Feeding ground density-dependence but recruitment related to unfished & recent (Eqn B.9b)
      if (DensDepOpt==2)
        for (Ifeed in 1:Nfeed)
        {
          # Survival (Eqn B.6)
          LogitSA <- log(1.0/SA-1.0);
          SurvOutB[Ibreed,Year] <- 1.0/(1.0+exp(LogitSA+SBdevYr[Ibreed,Year]*Sigma_SBdev));
          SurvOutF[Ifeed,Year] <- 1.0/(1.0+exp(LogitSA+SFdevYr[Ifeed,Year]*Sigma_SFdev));
          if (StochSopt==0) SAuse = SurvOutB[Ibreed,Year];
          if (StochSopt==1) SAuse = SurvOutF[Ifeed,Year];
          # Eqn B.4
          Term1 <- 1.0 - Nf[Ifeed,Ioffset]/(FeedK[Ifeed]*MultK[Ifeed,Ioffset]);
          Term1 <- f0*(1.0+ParAV[Ifeed,Ioffset]*Term1);
          Term1 <- 0.000001+(Term1-0.000001)/(1+exp(-30*Term1)) 
          LogitFec <- log(1.0/Term1-1.0);
          Term1 <- 1.0/(1.0+exp(LogitFec+FBdevYr[Ibreed,Year]*Sigma_FBdev));
          FecOut[Ifeed,Year] <- Term1;
          SurvTot[Ibreed,Ifeed,Year] <- SAuse*NN[Ibreed,Ifeed,Year]/NNS[Ibreed,Ifeed,Year];
          MortDiff[Ibreed,Ifeed,Year] <- (SAuse-SA)*NN[Ibreed,Ifeed,Year];
          # Eqn B.9b
          SurvJuv <- SC/SA
          for (Iage in 0:(Amat))
          {
            if (StochSopt==0) SurvJuv <- SurvJuv*SurvOutB[Ibreed,max(1,Year-Iage)]
            if (StochSopt==1) SurvJuv <- SurvJuv*SurvOutF[Ifeed,max(1,Year-Iage)]
          }
          NN[Ibreed,Ifeed,Year+1] <- SAuse*NN[Ibreed,Ifeed,Year] + 0.5*0.5*Term1*SurvJuv*(Nb[Ibreed,Ioffset]*Mix[Ibreed,Ifeed]+NN[Ibreed,Ifeed,Ioffset]);
          ParAV[Ifeed,Year+1] <- ParAV[Ifeed,Ioffset]*exp(inert_par * log(FeedK[Ifeed]/Nf[Ifeed,Ioffset]));
        }
      
      # Update breeding ground numbers
      Nb[Ibreed,Year+1] <- 0;  for (Ifeed in 1:Nfeed) Nb[Ibreed,Year+1] <- Nb[Ibreed,Year+1] + NN[Ibreed,Ifeed,Year+1];
    }

    # Generate new surveys
    for (Isurv in 1:length(DataSpecs[,1]))
    {
      SurvType <- DataSpecs[Isurv,1]
      # Feeding ground survey
      if (SurvType==2 & Year > NyrHist & ((Year+1-DataSpecs[Isurv,3]) %% DataSpecs[Isurv,4])==0)
      {
        SurvArea <- DataSpecs[Isurv,2]
        SurvCV <- DataSpecs[Isurv,5]
        Est <- sum(NN[,SurvArea,Year])*exp(rnorm(1,0,SurvCV)-SurvCV^2/2)
        cat(Year+1,SurvType,SurvArea,SurvCV,Est,"\n")
        Vec <- c(Year,Year,3,SurvArea,1,1,1,Est,SurvCV)
        TheSurveyData <- rbind(TheSurveyData,Vec)
      }  
      # Breeding ground survey
      if (SurvType==1 & Year > NyrHist & ((Year+1-DataSpecs[Isurv,3]) %% DataSpecs[Isurv,4])==0)
      {
        SurvArea <- DataSpecs[Isurv,2]
        SurvCV <- DataSpecs[Isurv,5]
        Est <- sum(NN[SurvArea,,Year+1])*exp(rnorm(1,0,SurvCV)-SurvCV^2/2)
        cat(Year+1,SurvType,SurvArea,SurvCV,Est,"\n")
        Vec <- c(Year,Year,2,SurvArea,1,1,1,Est,SurvCV)
        TheSurveyData <- rbind(TheSurveyData,Vec)
      }
    } # Gen survey
    
    
  }
  
  # Save herd numbers (start of the year year beyond the end of the simulation
  for (Ibreed in 1:Nbreed)
   for (Ifeed in 1:Nfeed) NNS[Ibreed,Ifeed,Nyr+1] <- NN[Ibreed,Ifeed,Nyr+1];
  # Save the breeding ground numbers
  for (Ibreed in 1:Nbreed)
   for (Ifeed in 1:Nfeed) NbS[Ibreed,Nyr+1] <- NbS[Ibreed,Nyr+1] + NN[Ibreed,Ifeed,Nyr+1];
  # Save the breeding ground numbers
  for (Ifeed in 1:Nfeed)
   for (Ibreed in 1:Nbreed) NfS[Ifeed,Nyr+1] <- NfS[Ifeed,Nyr+1] + NN[Ibreed,Ifeed,Nyr+1];
  
  # Summary outputs
  NT<- rep(0,Nyr+1);                                             # Total numbers
  Nb <- matrix(0,nrow=Nbreed,ncol=Nyr+1);                        # Breeding numbers
  Nf <- matrix(0,nrow=Nfeed,ncol=Nyr+1);                         # Feeding numbers
  for (Ibreed in 1:Nbreed) for (Iyr in 1:Nyr) Nb[Ibreed,Iyr] <- NbS[Ibreed,Iyr];
  for (Ifeed in 1:Nfeed)
    for (Iyr in 1:(Nyr+1)) Nf[Ifeed,Iyr] = NfS[Ifeed,Iyr];
  for (Iyr in 1:(Nyr+1))
   {
    TotalNT <- 0;
    for (Ibreed in 1:Nbreed) TotalNT <- TotalNT + NbS[Ibreed,Iyr];
    NT[Iyr] = TotalNT;
   }
  print(NT)

  Outs <- NULL
  Outs$NT <- NT
  Outs$Nbs <- Nb
  Outs$Nfs <- Nf
  Outs$CatchB <- CatchB2
  Outs$CatchF <- CatchF2
  Outs$TheSurveyData <- TheSurveyData
  
  return(Outs);
 }  
  
# ====================================================================================================================

DoRun <- function(Code,SensCase,StochSopt=1,StrayBase=0,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),SigmaDevS=6,SigmaDevF=0.01,WithMirror=1,Yr1=1970,Yr2=2023,
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",AllPlots=F,DoBoot=F,
                  ByCatchFile="BycatchActual_2024_04_24.csv",
                  WghtTotal=1,Idirichlet=1,MaxN=100,seed=19101,
                  BootUse)
{
  
  cat("Doing ",paste0(Code,"-",SensCase)," using data file ",DataFileName,"\n")
  
  # Specific model variants
  BreedingOpt <- substr(Code,1,2)
  FeedingOpt <- substr(Code,3,4)

  # -----------------------------------------------------------------------------------------------------------------------------------
  
  # Read in the structure data file
  DatFile <- read.table(DataFileName,fill=T,col.names=c(1:200),comment.char="?")
  
  # Details of the breeding and feeding grounds
  Index <- which(DatFile[,1]=="Number_of_breeding_and_feeding_grounds:"); Nbreed <- as.numeric(DatFile[Index,2]); Nfeed <- as.numeric(DatFile[Index,3])
  Index <- which(DatFile[,1]=="Breeding_grounds" & DatFile[,2]==BreedingOpt); Nbreed <- as.numeric(DatFile[Index,3]); BreedNames <- DatFile[Index+1,1:Nbreed]
  Index <- which(DatFile[,1]=="Feeding_grounds" & DatFile[,2]==FeedingOpt); Nfeed <- as.numeric(DatFile[Index,3]); FeedNames <- DatFile[Index+1,1:Nfeed]
  
  # Mixing stuff
  Index <- which(DatFile[,1]=="Mixing_matrix" & DatFile[,2]==Code)+1
  MixI <- matrix(0,nrow=Nbreed,ncol=Nfeed)
  rownames(MixI) <- BreedNames; colnames(MixI) <- FeedNames; NmixPar <- 0
  for (Ibreed in 1:Nbreed)
    for (Ifeed in 1:Nfeed)  
    {
      MixI[Ibreed,Ifeed] <- as.numeric(DatFile[Index+Ibreed,Ifeed])
      if (MixI[Ibreed,Ifeed] >0) NmixPar <- NmixPar + 1
    }
  
  Years <- Yr1:(Yr2+1); Nyear <- length(Years)
  Index <- which(DatFile[,2]=="Year_for_feeding_to_breeding"); YearFeedBreed <- as.numeric(DatFile[Index+1,1])-Yr1+1
  Index <- which(DatFile[,2]=="Dirichlet_(1)_or_normal_(0)"); Idirichlet <- as.numeric(DatFile[Index+1,1])

  # -----------------------------------------------------------------------------------------------------------------------------------
  # Read the catch data
  Outs <- ReadCatches(DatFile,BreedingOpt,FeedingOpt,Nbreed,Nfeed,BreedNames,FeedNames,Yr1,Yr2,CatchSer,ByCatchFile)
  CatchF <<- Outs$CatchF 
  CatchB <<- Outs$CAtchB
  
  # Read the survey data (note adjustments to indices to reflect C++ to R)
  Outs <- ReadSurveyData(Code,BreedingOpt,FeedingOpt,BreedNames,FeedNames,Yr1,Yr2,SensTest=SensCase)
  SurveyI <- Outs$SurveyI
  SurveyI[,1] <- SurveyI[,1] + 1
  SurveyI[,2] <- SurveyI[,2] + 1
  SurveyI[,5] <- SurveyI[,5] + 1
  Index <- which(SurveyI[,4]!=-1)
  SurveyI[Index,4] <- SurveyI[Index,4] + 1
  SurveyR <- Outs$SurveyR
  NsurveyData <- Outs$NsurveyData
  NsurveySeries <- Outs$NsurveySeries
  SurveySeries <- Outs$SurveySeries
  NextraCV1 <- Outs$NextraCV1
  
  Outs <- ReadMixingData(Code,DatFile,Yr1,Yr2,Nbreed,Nfeed,BreedNames,FeedNames,MixWeights,SensTest=SensCase,MaxN=MaxN)
  ObsMixBtoFE <- Outs$ObsMixBtoFE
  ObsMixBtoFP <- Outs$ObsMixBtoFP
  ObsMixBtoFO <- Outs$ObsMixBtoFO
  ObsMixFtoBE <- Outs$ObsMixFtoBE
  ObsMixFtoBP <- Outs$ObsMixFtoBP
  ObsMixFtoBO <- Outs$ObsMixFtoBO
  NmixData <- Outs$NmixData
  
  # ============================================================================================================================== 
  
  # Specify devs for survival and fecundity
  if (BreedingOpt=="B1")
   {
    SBdevEst <- c(0,1,1,0)
    SBdevMat <- matrix(c(0,0,2000,Yr2,2000,Yr2,0,0),ncol=2,byrow=T)-Yr1
    FBdevEst <- c(0,0,0,0)
    FBdevMat <- matrix(c(0,0,2000,Yr2,2000,Yr2,0,0),ncol=2,byrow=T)-Yr1
   }
  if (BreedingOpt=="B2")
   {
    SBdevEst <- c(0,1,0,1,0)
    SBdevMat <- matrix(c(0,0,2000,Yr2,0,0,2000,Yr2,0,0),ncol=2,byrow=T)-Yr1
    FBdevEst <- c(0,0,0,0,0)
    FBdevMat <- matrix(c(0,0,2000,Yr2,0,0,2000,Yr2,0,0),ncol=2,byrow=T)-Yr1
   }
  if (FeedingOpt=="F1" || FeedingOpt=="F2")
   {
    Yrs <- c(2000,Yr2)
    SFdevEst <- NULL; SFdevMat <- matrix(0,nrow=Nfeed,ncol=2,byrow=T)
    for (Ifeed in 1:Nfeed)
     { 
      SFdevEst <- c(SFdevEst,SF[Ifeed])
      if (SF[Ifeed]==1) SFdevMat[Ifeed,] <- Yrs-Yr1
     }
    FFdevEst <- c(0,0,0,0,0,0)
    FFdevMat <- matrix(c(0,0,rep(c(2000,Yr2),5)),ncol=2,byrow=T)-Yr1
   }
  
  # How many devs to estimate
  nBsdevs <- 0
  for (Ibreed in 1:Nbreed)
    if (SBdevEst[Ibreed]==1)
      for (Iyr in SBdevMat[Ibreed,1]:SBdevMat[Ibreed,2]) nBsdevs <- nBsdevs+1
  nBfdevs <- 0
  for (Ibreed in 1:Nbreed)
    if (FBdevEst[Ibreed]==1)
      for (Iyr in FBdevMat[Ibreed,1]:FBdevMat[Ibreed,2]) nfdevs <- nBfdevs+1
  nFsdevs <- 0
  for (Ifeed in 1:Nfeed)
    if (SFdevEst[Ifeed]==1)
      for (Iyr in SFdevMat[Ifeed,1]:SFdevMat[Ifeed,2]) nFsdevs <- nFsdevs+1
  nFfdevs <- 0
  for (Ibreed in 1:Nfeed)
    if (FFdevEst[Ifeed]==1)
      for (Iyr in FFdevMat[Ifeed,1]:FFdevMat[Ifeed,2]) nFfdevs <- nFfdevs+1
  
  # Mirror devs between feeding grounds
  Nmirror <- 0
  Mirror <- matrix(0,nrow=2,ncol=2)
  if (WithMirror==1)
   {
    Nmirror <- 1
    Mirror[1,1] <- 4; Mirror[1,2] <- 3;
    Mirror <- Mirror-1
   } 
  if (WithMirror==2)
   {
    Nmirror <- 2
    Mirror[1,1] <- 4; Mirror[1,2] <- 3;
    Mirror[2,1] <- 6; Mirror[2,2] <- 5;
    Mirror <- Mirror-1
   }
  
  data2 <<- list(Nbreed=Nbreed, Nfeed=Nfeed,Yr1=Yr1,Yr2=Yr2,IAmat=IAmat,SA=SA,SC=SC,TimeLag=TimeLag,
               SurveyI=SurveyI,SurveyR=SurveyR,CatchB=CatchB,CatchF=CatchF,
               SBdevEst=SBdevEst,SBdevMat=SBdevMat,FBdevEst=FBdevEst,FBdevMat=FBdevMat,
               SFdevEst=SFdevEst,SFdevMat=SFdevMat,FFdevEst=FFdevEst,FFdevMat=FFdevMat,
               MixI=MixI,
               NsurveyData=NsurveyData,NsurveySeries=NsurveySeries,SurveySeries=SurveySeries,
               NmixData=NmixData,YearFeedBreed=YearFeedBreed,Idirichlet=Idirichlet,
               ObsMixBtoFE=ObsMixBtoFE,ObsMixBtoFP=ObsMixBtoFP,ObsMixBtoFO=ObsMixBtoFO,
               ObsMixFtoBE=ObsMixFtoBE,ObsMixFtoBP=ObsMixFtoBP,ObsMixFtoBO=ObsMixFtoBO,
               StochSopt=StochSopt,StrayBase=StrayBase,DensDepOpt=DensDepOpt,
               Nmirror=Nmirror,Mirror=Mirror,WghtTotal=WghtTotal)
  #print(str(data))
  
  # Additional variances
  NextraCV <- NextraCV1
  if (NextraCV==0) NextraCV1 <- 1
  AddV=rep(0,NextraCV1)
  
  # Mixing parameters
  MixPars <- rep(-1,NmixPar)
  MixPars <- NULL
  for (Ibreed in 1:Nbreed)
   {
    Iref <- which(MixI[Ibreed,]==-1)
    Temp <- rep(0,Nfeed)
    for (IdataS in 1:2) 
      for (Ifeed in 1:Nfeed) Temp[Ifeed] <- Temp[Ifeed] + (0.001+ObsMixBtoFE[IdataS,Ibreed,Ifeed])/2.0
    for (Ifeed in 1:Nfeed) if (MixI[Ibreed,Ifeed]!=0 && Ifeed!=Iref) MixPars <- c(MixPars,log(Temp[Ifeed]/Temp[Iref]))
   }
  
  # Survival and fecundity devs
  if (nBfdevs==0) nBfdevs <- 1
  if (nFfdevs==0) nFfdevs <- 1
  SBdev <- rep(0,nBsdevs);FBdev <- rep(0,nBfdevs);
  SFdev <- rep(0,nFsdevs);FFdev <- rep(0,nFfdevs);
  
  parameters <- list(rval=0.07,logK=rep(log(20000),Nbreed),logBK=rep(3,Nbreed),
                     InfluxP=10,inert_par=0,MixPars=MixPars,AddV=AddV,
                     Sigma_SBdev=SigmaDevS, SBdev=SBdev,
                     Sigma_FBdev=SigmaDevF, FBdev=FBdev,
                     Sigma_SFdev=SigmaDevS, SFdev=SFdev,
                     Sigma_FFdev=SigmaDevF, FFdev=FFdev)
  

  # Bootstrap
  xx <- Bootstrap(Code=paste0(Code,SensCase),data2,parameters,map,rept,Yr1,Yr2,
                           BreedNames,FeedNames,Nboot=3,seed,BootUse)

  return(xx)
 }
# ==========================================================================================

Bootstrap <- function(Code,data,parameters,map,rept,Yr1,Yr2,BreedNames,FeedNames,Nboot=2,seed=19101,BootUse)
{
  FileNameBoot <- paste0("Diags/",Code,".Boot")
  TablePars <- read.table(FileNameBoot,col.names=paste0("C",1:20),fill=T)
  Nbreed <- length(BreedNames)
  Nfeed <- length(FeedNames)
  Nyear <- Yr2-Yr1+1
  set.seed(seed)
  seeds <- floor(runif(Nboot+1,1,10000))
  #print(parameters)
 
  SLASpecs <- matrix(0,nrow=2,ncol=3)
  SLASpecs[1,] <- c(1,2,50)                 # Breeding ground 2
  SLASpecs[2,] <- c(2,6,20)                 # Feeding ground 6 
  DataSpecs <- matrix(0,nrow=2,ncol=6)
  DataSpecs[1,] <- c(1,2,Nyear+1,10,0.13,1)     # Breeding ground
  DataSpecs[2,] <- c(2,6,Nyear+2, 1,0.20,1)       # Feeding ground
  
  Results <- vector(mode="list",length=Nboot)
  
  # Do bootstrap
  IbootPnt <- 1
  for (Iboot in 1:Nboot)
   {
    set.seed(seeds[Iboot+1])
    #cat("Boot",Iboot,seeds[Iboot+1],"\n")
    logK <- rep(0,Nbreed)
    
    # Read from the bootstrap file
    IbootPnt <- IbootPnt
    rval <- as.numeric(TablePars[IbootPnt+1,1]); IbootPnt <- IbootPnt+1
    logK =as.numeric(TablePars[IbootPnt+c(1:Nbreed),1]); IbootPnt <- IbootPnt+Nbreed
    logBK =as.numeric(TablePars[IbootPnt+c(1:Nbreed),1]); IbootPnt <- IbootPnt+Nbreed
    Nmix <- length(parameters$MixPars)
    MixPars =as.numeric(TablePars[IbootPnt+c(1:Nmix),1]); IbootPnt <- IbootPnt+Nmix
    Nadvar <- length(parameters$AddV)
    AddV =as.numeric(TablePars[IbootPnt+c(1:Nadvar),1]); IbootPnt <- IbootPnt+Nadvar
    NsurvPar <- length(parameters$SFdev)
    SFdev = as.numeric(TablePars[IbootPnt+c(1:NsurvPar),1]); IbootPnt <- IbootPnt+NsurvPar
    
    # Skip Total abundance
    IbootPnt <- IbootPnt + (Nyear+2) + (Nyear+2) + (Nyear+2) + (Nyear+2)+1

    #print(str(data))
    modelpars <- NULL
    modelpars$rval <- rval
    modelpars$logK <- logK
    modelpars$logBK <- logBK
    modelpars$MixPars <- MixPars
    modelpars$AddV <- AddV
    modelpars$SFdev <- SFdev
    modelpars$InfluxP <- parameters$InfluxP
    modelpars$Sigma_SBdev <- parameters$Sigma_SBdev
    modelpars$Sigma_SFdev <- parameters$Sigma_SFdev
    modelpars$Sigma_FBdev <- parameters$Sigma_FBdev
    modelpars$Sigma_FFdev <- parameters$Sigma_FFdev
    modelpars$inert_par <- parameters$inert_par
    print("here")
    #print(str(modelpars))
    xx <- RunModelProj(data,modelpars,SLASpecs,DataSpecs,Nproj=20)
    Results[[Iboot]] <- xx
 
  } # Boot
  print(str(Results))
  return(Results)
}
# ==========================================================================================

Codes <- c("B1F1","B1F2","B2F1","B2F2")
BootUse <- matrix(1,nrow=4,ncol=500)
BootUse[2,145] <- 2
for (Icode in 3:3)
 {
  Code <- Codes[Icode]  
  xx <- DoRun(Code,SensCase="BC",AllPlots=T,DoBoot=T,seed=1234,BootUse=BootUse[Icode,])
 }

