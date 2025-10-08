# FUNCTION THAT DEFINES THE HUMPBACK ASSESSMENT MODEL
cmb <- function(f, d) function(p) f(p, d) # combine objective function with specific data
f <- function(parms,dat){
  
  getAll(dat, parms, warn=FALSE)                                  # Make the data and parameters global to this function
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
  logBK_sumSQ <- sum((logBK-mean(logBK))^2)
  Penal <- Penal + logBK_sumSQ
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
