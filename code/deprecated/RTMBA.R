library(RTMB)

# ========================================================================================

setwd("C:\\Research\\Iwc25\\NP_Humpbacks\\")
source("HumpFigs2.R")
source("ReadData.R")
source("WriteOut.R")

FullDiag <- F
PlotsDir <- "Plots/"

DataFileName <- "Hump.dat"

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

UseKPrior <- 1
# ========================================================================================
f <- function(parms)
 {
  sumSQ <- function(x){
    TheMean <- sum(x)/length(x);
    Output <- 0;
    for (ii in 1:length(x)) Output <- Output + (x[ii]-TheMean)^2.0;
    return(Output);
  }
  # 
  getAll(data2, parms, warn=FALSE)                                  # Make the data and parameters global to this function
  Nyr <- Yr2-Yr1+1;                                                 # Number of years
  MixPar <- length(MixPars);                                        # Number of mixing parameters
  
  ThreshPop = 1.0e-20;                                              # Minimum population size
  Influx = 1.0/(1+exp(InfluxP));                                    # Transfer influx
  
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
  # Population projection
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

  Penal2 <- 100.0/(1+exp(-1000.0*(fmax-0.98)));
  
  
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
  
  # Temporary variables
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
  for (Year in 1:Nyr)
   {
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
      Temp1[Ibreed] <- Nb[Ibreed,Year] - CatchB[Year,Ibreed];
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
      for (Ibreed in 1:Nbreed) CBreed[Ibreed] <- CBreed[Ibreed] + PropMixEst[Ifeed,Ibreed,Year]*CatchF[Year,Ifeed];
      for (Ibreed in 1:Nbreed)
       if (MixI[Ibreed,Ifeed] != 0)
        {
         NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year] - PropMixEst[Ifeed,Ibreed,Year]*CatchF[Year,Ifeed];
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

  # =================================================================================================================================
  # Likelihood
  # =================================================================================================================================

  LogLike1 = 0;
  LikeCompSurv <- rep(0,NsurveyData); 
  LikeMixComp <- rep(0,NmixData);
  LogLike2a <- rep(0,2); LogLike2b <- rep(0,2)
  PredSurv <- rep(0,NsurveyData); SurveySD <- rep(0,NsurveyData); Qest <- rep(0,NsurveySeries)
  PredMixOut <- array(0,dim=c(2,Nbreed,Nfeed))
  PredMix <- rep(0,NmixData)
  ObsMixPropI <- matrix(0,nrow=NmixData,ncol=4);                           # Mixing observations (integers)
  ObsMixProp <- matrix(0,nrow=NmixData,ncol=2);                            # Mixing Observations (reals)
  
  # Survey data (Absolute abundance estimates) [Eqn B.1]
  for (Iclass in 1:NsurveySeries)
   {
    if (SurveySeries[Iclass]==1)
     {
      for (Isurv in 1:NsurveyData)
       if (SurveyI[Isurv,5] == Iclass)
        {
         Pred <- 0;
         for (Iyr in SurveyI[Isurv,1]:SurveyI[Isurv,2])
          {
           if (SurveyI[Isurv,3]==1) for (Ibreed in 1:Nbreed) Pred <- Pred + NbS[Ibreed,Iyr];
           if (SurveyI[Isurv,3]==2) Pred <- Pred + NbS[SurveyI[Isurv,4],Iyr];
           if (SurveyI[Isurv,3]==3) Pred <- Pred + NfS[SurveyI[Isurv,4],Iyr];
          }
         Pred <- Pred / (SurveyI[Isurv,2]-SurveyI[Isurv,1]+1);
         PredSurv[Isurv] <- Pred;
         SD2 <- SurveyR[Isurv,2]*SurveyR[Isurv,2];
         if (SurveyI[Isurv,7]>0) SD2 <- SD2 + AddV[SurveyI[Isurv,7]];
         SurveySD[Isurv] <- sqrt(SD2);
         LikeCompSurv[Isurv] <- -dnorm(log(SurveyR[Isurv,1]), log(PredSurv[Isurv]), sqrt(SD2), log=T);
         if (SurveyI[Isurv,3]==1) LikeCompSurv[Isurv] <- LikeCompSurv[Isurv] * WghtTotal;
         if (SurveyI[Isurv,6]==1) LogLike1 <- LogLike1 + LikeCompSurv[Isurv];
        }
    } # Absolute indicies
  
  # Survey data (Relative abundance estimates)
  if (SurveySeries[Iclass]==2)
   {
    # Find ML estimare of survey Q
    logq = 0; ndat = 0;
    for (Isurv in 1:NsurveyData)
     if (SurveyI[Isurv,5]==Iclass)
      {
       Pred = 0;
       for (Iyr in SurveyI[Isurv,1]:SurveyI[Isurv,2])
        {
         if (SurveyI[Isurv,3]==1) for (Ibreed in 1:Nbreed) Pred <- Pred + NbS[Ibreed,Iyr];
         if (SurveyI[Isurv,3]==2) Pred <- Pred + NbS[SurveyI[Isurv,4],Iyr];
         if (SurveyI[Isurv,3]==3) Pred <- Pred + NfS[SurveyI[Isurv,4],Iyr];
        }
       Pred = Pred / (SurveyI[Isurv,2]-SurveyI[Isurv,1]+1);
       PredSurv[Isurv] <- Pred;
       SD2 = SurveyR[Isurv,2]^2;
       if (SurveyI[Isurv,7]>0) SD2 <- SD2 + AddV[SurveyI[Isurv,7]];
       if (SurveyI[Isurv,6]==1) ndat <- ndat + 1.0/SD2;
       if (SurveyI[Isurv,6]==1) logq <- logq + log(SurveyR[Isurv,1]/Pred)/SD2;
      }
    q <- exp(logq/ndat);
    Qest[Iclass] <- q;
    # Actual likelihood (Eqn C.2)
    for (Isurv in 1:NsurveyData)
     if (SurveyI[Isurv,5] == Iclass)
      {
       PredSurv[Isurv] <- q*PredSurv[Isurv];
       SD2 <- SurveyR[Isurv,2]^2;
       if (SurveyI[Isurv,7]>0) SD2 <- SD2 + AddV[SurveyI[Isurv,7]];
       SurveySD[Isurv] <- sqrt(SD2);
       LikeCompSurv[Isurv] <- -dnorm(log(SurveyR[Isurv,1]), log(PredSurv[Isurv]), sqrt(SD2), log=T);
       if (SurveyI[Isurv,6]==1) LogLike1 <- LogLike1 + LikeCompSurv[Isurv];
      }
    } 
   } # Relative indices

  # Mixing proportions (Breeding to Feeding) (Eqn C.3)
  Icnt <- 1; Jcnt <- 1;
  for (IdataS in 1:2)
   {
    for (Ibreed in 1:Nbreed)
     {
      for (Ifeed in 1:Nfeed) PredMixOut[1,Ibreed,Ifeed] = NNS[Ibreed,Ifeed,YearFeedBreed]/NbS[Ibreed,YearFeedBreed];
      for (Ifeed in 1:Nfeed)
       if (ObsMixBtoFP[IdataS,Ibreed,Ifeed] > 0)
        {
         Pred <- NNS[Ibreed,Ifeed,YearFeedBreed]/NbS[Ibreed,YearFeedBreed];
         Obs <- ObsMixBtoFE[IdataS,Ibreed,Ifeed];
         SD <- ObsMixBtoFP[IdataS,Ibreed,Ifeed];
         PredMix[Icnt] = Pred;
         ObsMixPropI[Icnt,1] = 1;
         ObsMixPropI[Icnt,2] = IdataS;
         ObsMixPropI[Icnt,3] = Ibreed;
         ObsMixPropI[Icnt,4] = Ifeed;
         ObsMixProp[Icnt,1] = Obs;
         ObsMixProp[Icnt,2] = SD;
         if (Obs > 0)
          {
           if (Idirichlet == 0) LikeMixComp[Icnt] = -dnorm(Obs, Pred, SD2, log=T);
           if (Idirichlet == 1) LikeMixComp[Icnt] = -ObsMixBtoFO[IdataS,Ibreed]*Obs*log((Pred+0.0001)/Obs);
          }
         LogLike2a[Jcnt] <- LogLike2a[Jcnt] + LikeMixComp[Icnt];
         Icnt <- Icnt + 1;
        }
      }
    Jcnt <- Jcnt + 1;
   }
  #  Mixing proportions (Feeding to Breeding)
  Jcnt <- 1;
  for (IdataS in 1:2)
   {
    for (Ifeed in 1:Nfeed)
     {
      for (Ibreed in 1:Nbreed)
       PredMixOut[2,Ibreed,Ifeed] <- NNS[Ibreed,Ifeed,YearFeedBreed]/NfS[Ifeed,YearFeedBreed];
      for (Ibreed in 1:Nbreed)
       if (ObsMixFtoBP[IdataS,Ibreed,Ifeed] > 0)
        {
          Pred <- NNS[Ibreed,Ifeed,YearFeedBreed]/NfS[Ifeed,YearFeedBreed];
          Obs <- ObsMixFtoBE[IdataS,Ibreed,Ifeed];
          SD <- ObsMixFtoBP[IdataS,Ibreed,Ifeed];
          PredMix[Icnt] <- Pred;
          ObsMixPropI[Icnt,1] <- 2;
          ObsMixPropI[Icnt,2] <- IdataS;
          ObsMixPropI[Icnt,3] <- Ibreed;
          ObsMixPropI[Icnt,4] <- Ifeed;
          ObsMixProp[Icnt,1] = Obs
          ObsMixProp[Icnt,2] = SD;
          if (Obs > 0)
           {
            if (Idirichlet == 0) LikeMixComp[Icnt] <- -dnorm(Obs, Pred, SD2, log=T);
            if (Idirichlet == 1) LikeMixComp[Icnt] <- -ObsMixFtoBO[IdataS,Ibreed]*Obs*log((Pred+0.0001)/Obs);
           }
          LogLike2b[Jcnt] <- LogLike2b[Jcnt] + LikeMixComp[Icnt];
          Icnt <- Icnt + 1;
        }
    }
    Jcnt <- Jcnt + 1;
  }

  # Weak penalties for stability
  Penal <- 0;
  Penal <- Penal2;
  for (Ibreed in 1:Nbreed) Penal <- Penal + 0.0001*logK[Ibreed]*logK[Ibreed];
  for (Imix in 1:MixPar) Penal <- Penal + MixPars[Imix]*MixPars[Imix];
  Penal <- Penal + sumSQ(logBK)
  Penal <- Penal - sum(dnorm(SBdev,0.0,1.0,log=T));
  Penal <- Penal - sum(dnorm(FBdev,0.0,1.0,log=T));
  Penal <- Penal - sum(dnorm(SFdev,0.0,1.0,log=T));
  TotalK <- sum(BreedK)
  if (UseKPrior > 0)
   Penal <- Penal - UseKPrior*log(1.0/(1.0+exp(Prior_int+Prior_slope*TotalK)))
  datalike <- LogLike1 + sum(LogLike2a) + sum(LogLike2b);
  neglogL <-  Penal + datalike;
  #print(neglogL)

  # Summary outputs
  LogNT<- rep(0,Nyr+1);                                             # Log total numbers
  LogNb <- matrix(0,nrow=Nbreed,ncol=Nyr+1);                        # Breeding numbers
  LogNf <- matrix(0,nrow=Nfeed,ncol=Nyr+1);                         # Feeding numbers
  for (Ibreed in 1:Nbreed) for (Iyr in 1:(Nyr+1)) LogNb[Ibreed,Iyr] <- log(NbS[Ibreed,Iyr]);
  for (Ifeed in 1:Nfeed)
    for (Iyr in 1:(Nyr+1)) LogNf[Ifeed,Iyr] = log(NfS[Ifeed,Iyr]);
  for (Iyr in 1:(Nyr+1))
   {
    TotalNT <- 0;
    for (Ibreed in 1:Nbreed) TotalNT <- TotalNT + NbS[Ibreed,Iyr];
    LogNT[Iyr] = log(TotalNT);
   }
  #print(exp(LogNT))

  ADREPORT(LogNT);
  ADREPORT(LogNb);
  ADREPORT(LogNf);
  ADREPORT(SBdevYr);
  ADREPORT(SurvOutB);
  ADREPORT(SurvOutF);
  
  REPORT(Nb);
  REPORT(Nf);
  REPORT(NbS);
  REPORT(NfS);
  REPORT(LikeCompSurv);
  REPORT(LikeMixComp);
  REPORT(PredSurv);
  REPORT(Qest);
  REPORT(f0);
  REPORT(fmax)
  REPORT(neglogL);
  REPORT(Penal);
  REPORT(LogLike1);
  REPORT(LogLike2a);
  REPORT(LogLike2b);
  REPORT(BreedK);
  REPORT(FeedK);
  REPORT(Mix);
  REPORT(ObsMixPropI);
  REPORT(ObsMixProp);
  REPORT(PredMix);
  REPORT(BreedK);
  REPORT(FeedK);
  REPORT(NN);
  REPORT(NNS);
  REPORT(FecOut);
  REPORT(SBdevYr);
  REPORT(FBdevYr);
  REPORT(Influx);
  REPORT(AddV);
  REPORT(datalike);
  REPORT(SurvOutB);
  REPORT(SurvOutF);
  REPORT(ParAV);
  REPORT(MultK);
  REPORT(MortDiff);
  REPORT(PredMixOut);
  REPORT(SurvTot);
  REPORT(AddV);
  REPORT(SurveySD);
  
  return(neglogL);
 }  
  
# ====================================================================================================================

DoRun <- function(Code,SensCase,StochSopt=1,StrayBase=0,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),SigmaDevS=6,SigmaDevF=0.01,WithMirror=1,Yr1=1970,Yr2=2023,
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",AllPlots=F,DoBoot=F,
                  ByCatchFile="BycatchActual_2024_04_24.csv",
                  WghtTotal=1,Idirichlet=1,MaxN=100,seed=19101,
                  BootUse,SetNew=0,Init=NULL)
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
  
  map1=list(#rval=factor(NA),#logBK=rep(factor(NA),Nbreed),
    InfluxP=factor(NA),
    inert_par=factor(NA),
    #Sigma_SBdev=factor(NA),
    Sigma_SBdev=factor(NA),SBdev=rep(factor(NA),length(SBdev)),
    Sigma_FBdev=factor(NA),FBdev=rep(factor(NA),length(FBdev)),
    Sigma_SFdev=factor(NA),
    #Sigma_SFdev=factor(NA),SFdev=rep(factor(NA),length(SFdev)),
    Sigma_FFdev=factor(NA),FFdev=rep(factor(NA),length(FFdev))
  )
  map2=list(#rval=factor(NA),#logBK=rep(factor(NA),Nbreed),
    InfluxP=factor(NA),
    inert_par=factor(NA),
    Sigma_SBdev=factor(NA),
    #Sigma_SBdev=factor(NA),SBdev=rep(factor(NA),length(SBdev)),
    Sigma_FBdev=factor(NA),FBdev=rep(factor(NA),length(FBdev)),
    #Sigma_SFdev=factor(NA),
    Sigma_SFdev=factor(NA),SFdev=rep(factor(NA),length(SFdev)),
    Sigma_FFdev=factor(NA),FFdev=rep(factor(NA),length(FFdev))
  )
  if (StochSopt==0) map <- map2
  if (StochSopt==1) map <- map1
  if (AddCV==F) map$AddV <- rep(factor(NA),length(AddV))
  
  ################################################################################
  if (FullDiag==T) print("Calling MakeADM")
  model <- MakeADFun(f, parameters, map=map,DLL="Hump",silent=T)
  rept <<- model$report()
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
  best <- model$env$last.par.best
  stdreport <-sdreport(model)
  rep <- summary(stdreport)
  rep2<- summary(stdreport, "fixed", p.value = TRUE)
  #print(rep2)
  #print(best)
  rept <<- model$report()
  #print(best)
  #print(rep2)

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
  
  # Print outs and graphs
  WriteOut(Code,Abbrev=SensCase,Yr1,Yr2,BreedNames,FeedNames,rept,rep,rep2,StockDef,data2)
  
  PDF <- F
  if (AllPlots==T) PDF <- T
  Code2 <- paste0(Code,"-",SensCase)
  if (AllPlots==T) PlotPropns(Nprop,rept$ObsMixPropI,rept$ObsMixProp,rept$PredMix,PDF=PDF,Code=Code2,BreedNames,FeedNames,Nyear)
  PlotIndex(cbind(SurveyI,SurveyR),LogNT,LogNb,LogNf,rept$Qest,PDF=PDF,Code=Code2,BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed,rept,AllPlots=AllPlots)
  if (AllPlots==T) PlotTraj(LogNT,LogNb,LogNf,PDF=PDF,Code=Code2,BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed)
  PlotSurvF(SurvOutF,baseS=0.96,PDF=PDF,Code=Code2,BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed)
  if (AllPlots==T) PlotSurvF(SurvOutF,baseS=0.96,PDF=PDF,Code=Code2,BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed)
  
  # Bootstrap
  if (DoBoot==T) Bootstrap(Code=paste0(Code,SensCase),data2,parameters,map,rept,Yr1,Yr2,
                           BreedNames,FeedNames,Nboot=500,seed,BootUse,bestOrig=best)

  return(best)
 }
# ==========================================================================================

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
         data2$SurveyR[II,1] <<- as.numeric(PredSurv[II]*exp(rnorm(1,0,SEL)-SEL^2.0/2.0))
         #new <- PredSurv[II]*exp(rnorm(1,0,SEL)-SEL^2.0/2.0)
         #new <- data2$SurveyR[II,1]
         #print(c(data2$SurveyR[II,1],PredSurv[II],new))
         if (data2$SurveyR[II,1]<=0) print("oops1")
       }
       # Generate data (mixing)
       for (II in 1:2)
        for (Ibreed in 1:Nbreed)
         {
          xx <- rmultinom(1:Nfeed,round(data$ObsMixBtoFO[II,Ibreed]),prob=0.00001+PredMixOut[1,Ibreed,])  
          data2$ObsMixBtoFE[II,Ibreed,] <<- xx/sum(xx)
          if (any(data2$ObsMixBtoFE[II,Ibreed,]<0)) print("oops2")
         }
       for (II in 1:2)
        for (Ifeed in 1:Nfeed)
         {
          xx <- rmultinom(1:Nbreed,round(data$ObsMixFtoBO[II,Ifeed]),prob=0.00001+PredMixOut[2,,Ifeed])  
          data2$ObsMixFtoBE[II,,Ifeed] <<- xx/sum(xx)
          if (any(data2$ObsMixFtoBE[II,,Ifeed]<0)) print("oops3")
        }
       data2$NmixData <- 1000
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

# Select the case for the analysis
Case <- 3

# Base run (for testing)
if (Case==1)
 {
  #xx <- DoRun("B1F1",SensCase="BC",AllPlots=F,DoBoot=F,Init=best,SetNew=1)
  xx <- DoRun("B1F1",SensCase="BC",AllPlots=F,DoBoot=F,Init=NULL,SetNew=0)
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

