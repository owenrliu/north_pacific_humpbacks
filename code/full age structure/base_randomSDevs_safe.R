calcKPrior <- function(Kmax){
  
  Ks <- seq(from=0,to=3*Kmax,length=100)
  
  join1 <- 1/(1+exp(30*(Ks-Kmax)))
  join2 <- 1/(1+exp(-30*(Ks-Kmax)))
  
  K1 <- 2*Kmax; Y1 <- log(1/0.0001-1)
  K2 <-   Kmax; Y2 <- log(1/0.9999-1)
  Prior_slope <- (Y2-Y1)/(K2-K1)
  Prior_int <- Y2 - Prior_slope*K2
  #cat(Prior_slope,Prior_int,"\n")
  Priors <- 1.0/(1+exp(Prior_int+Prior_slope*Ks))
  KPriors <- list(Prior_int=Prior_int,Prior_slope=Prior_slope)
  
  KPriors
}

# BASELINE HUMPBACK ASSESSMENT MODEL - NO ENVIRONMENTAL COVARIATES
# ========================================================================================
cmb <- function(f, d) function(p) f(p, d) # combine objective function with specific data
f <- function(parms,dat)
 {
  # ========================================================================================================================
  # Initial parameters
  # ========================================================================================================================
  # Get the parameters from the input data
  getAll(dat, parms, warn=FALSE)
  # Number of years
  Nyr <- length(Yr1:Yr2)
  # Minimum population size
  ThreshPop = 1.0e-20;
  # Transfer influx
  Influx = 1.0/(1+exp(InfluxP));
  
  # =================================================================================================================================
  # Use Prior for K?
  # =================================================================================================================================
  if(UseKPrior==1){
    KPriors <- calcKPrior(Kmax=Kmax)
    Prior_int <- KPriors$Prior_int
    Prior_slope <- KPriors$Prior_slope
  }
  # ========================================================================================================================
  # Set up mixing matrix
  # ========================================================================================================================
  # tMixI <- t(MixI)
  # Mix <- as.vector(tMixI)
  # # Fill in the estimated mixing parameters in the right spot
  # Mix[tMixI>0] <- exp(MixPars)
  # Mix[tMixI<=0] <- abs(tMixI)[tMixI<=0]
  # # Re-form the matrix and normalize across feeding grounds
  # Mix <- matrix(Mix,nrow=Nbreed,ncol=Nfeed,byrow=T)
  # Mix <- Mix/rowSums(Mix)
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
    # SAFETY: Protect against zero total
    Total <- max(Total, 1e-10)
    for (Ifeed in 1:Nfeed) Mix[Ibreed,Ifeed] <- Mix[Ibreed,Ifeed]/Total;
  }

  # ========================================================================================================================
  # Set up population projection
  # ========================================================================================================================

  # Basic demographics
  Amat <- IAmat*1.0
  # Minimum maturity age; min age for inclusion in density depl; min age for use in catch removal
  MatAgeT <- Amat+2
  # Minimum age for inclusion in density-dependence
  MinAgeD <- 2
  # Minimum age vulnerable to catch removal
  MinAgeC <- 1
  
  # SAFETY: Cap survival to prevent division by zero
  SA <- min(SA, 0.9999)
  SC <- min(SC, 0.9999)
  
  # Equlibrium numbers
  Neqn <- rep(0,Nage); Neqn[1] <- 1; Neqn[2] <- SC;
  for (Iage in 3:Nage) Neqn[Iage] <- Neqn[Iage-1]*SA;
  Neqn[Nage] <- Neqn[Nage]/(1-SA);
  Neqn <- Neqn/sum(Neqn)
  NfecEqn <- sum(Neqn[MinAgeD:Nage]) # Proportional 1+ abundance
  # SAFETY: Protect against zero denominator
  NfecEqn <- max(NfecEqn, 1e-10)
 
  # Zerbini et al. (2010): Mar Biol. 157: 1225-1236. Note that Amat is age-at-maturity (age-at-first parturition less 1)
  # SAFETY: Cap exponential calculations to prevent overflow
  exp_term1 <- min(rval*(Amat+1.0), 700)
  exp_term2 <- min(rval*Amat, 700)
  fmax <- 2*(exp(exp_term1)-SA*exp(exp_term2))/(SC*SA^(Amat))
  f0 <- 2*(1.0-SA)/(SC*SA^(Amat))
  # Fecundity at carrying capacity
  ParA <- (fmax-f0)/f0

  # print(c(rval,fmax,f0,ParA))

  # Penalty on fmax
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
  
  # Set up K (1+) by herd, breeding stock, and feeding ground
  # K by breeding ground
  BreedK <- exp(logK)
  # K by herd
  NNK <- BreedK*Mix
  # K by feeding ground
  FeedK <- colSums(NNK)
  
  # ========================================================================================================================
  # Abundance Tracking Matrices
  # ========================================================================================================================
  
  # Numbers at age in each herd (updated after each process)
  NNN <- array(0,dim=c(Nbreed,Nfeed,Nyr+1,Nage))
  # Numbers in each herd, all ages, start of year
  NNS <- array(0,dim=c(Nbreed,Nfeed,Nyr+1));
  # breeding ground numbers, start of year
  NbS <- matrix(0,nrow=Nbreed,ncol=Nyr+1)
  
  # breeding ground numbers before removals
  Nb <- matrix(0,nrow=Nbreed,ncol=Nyr+1); 
  # feeding ground numbers, before removals
  NfS <- matrix(0,nrow=Nfeed,ncol=Nyr+1)
  # feeding grounds numbers after removals
  Nf <- matrix(0,nrow=Nfeed,ncol=Nyr+1)
  
  # Total numbers in each herd, after removals
  NN <- array(0,dim=c(Nbreed,Nfeed,Nyr+1));
  
  # Fecundity
  FecOut <-matrix(0,nrow=Nfeed,ncol=Nyr)
  # 
  PropMixEst <- array(0,dim=c(Nfeed,Nbreed,Nyr+1))
  # Time varying survival, breeding grounds (0-1)
  SurvOutB <- matrix(0,nrow=Nbreed,ncol=Nyr)
  # Time varying survival, feeding grounds (0-1)
  SurvOutF <- matrix(0,nrow=Nfeed,ncol=Nyr)
  # Total proportional survival, including breeding ground removals and survival
  SurvTot <- array(0,dim=c(Nbreed,Nfeed,Nyr));
  # Difference in mortality on feeding grounds, relative to baseline survival
  MortDiff <- array(0,dim=c(Nbreed,Nfeed,Nyr));
  
  # Breeding and feeding output (1=0+; 2=1+; 3=mature)
  NfitBreed <- array(0,dim=c(3,Nbreed,Nyr+1))
  NfitFeed <- array(0,dim=c(3,Nfeed,Nyr+1))
  
  # ========================================================================================================================
  # Initial year
  # ========================================================================================================================
  
  for (Ibreed in 1:Nbreed)
   {
    Nb[Ibreed,1] <- BreedK[Ibreed]*1.0/(1+exp(logBK[Ibreed]));
    for (Ifeed in 1:Nfeed) 
     {
      NN[Ibreed,Ifeed,1] <- Nb[Ibreed,1]*Mix[Ibreed,Ifeed];
      NNN[Ibreed,Ifeed,1,] <- NN[Ibreed,Ifeed,1]*Neqn/NfecEqn
     }
   }
  
  # ========================================================================================================================
  # Projection
  # ========================================================================================================================
  
  for (Year in 1:Nyr)
   {
    #cat("Year = ",Year,Nyr,"\n")
    # Save breeding ground numbers (1+) (start of the year)
    for (Ibreed in 1:Nbreed)
    {
      for (Ifeed in 1:Nfeed){
        NNS[Ibreed,Ifeed,Year] <- sum(NNN[Ibreed,Ifeed,Year,2:Nage])
      }
      NbS[Ibreed,Year] <- sum(NNS[Ibreed,,Year])
      for (Icomp in 1:3)
      {
        if (Icomp==1) NfitBreed[Icomp,Ibreed,Year] <- sum(NNN[Ibreed,,Year,])
        if (Icomp==2) NfitBreed[Icomp,Ibreed,Year] <- sum(NNN[Ibreed,,Year,2:Nage])
        if (Icomp==3) NfitBreed[Icomp,Ibreed,Year] <- sum(NNN[Ibreed,,Year,MatAgeT:Nage])
        # SAFETY: Protect against zero abundance
        NfitBreed[Icomp,Ibreed,Year] <- max(NfitBreed[Icomp,Ibreed,Year], ThreshPop)
      }
    }

    # Save feeding ground numbers (1+) (start of the year)
    for (Ifeed in 1:Nfeed)
    {
      NfS[Ifeed,Year] <- sum(NNS[,Ifeed,Year])
      for (Icomp in 1:3)
      {
        if (Icomp==1) NfitFeed[Icomp,Ifeed,Year] <- sum(NNN[,Ifeed,Year,])
        if (Icomp==2) NfitFeed[Icomp,Ifeed,Year] <- sum(NNN[,Ifeed,Year,2:Nage])
        if (Icomp==3) NfitFeed[Icomp,Ifeed,Year] <- sum(NNN[,Ifeed,Year,MatAgeT:Nage])
        # SAFETY: Protect against zero abundance
        NfitFeed[Icomp,Ifeed,Year] <- max(NfitFeed[Icomp,Ifeed,Year], ThreshPop)
      }
    }

    # Breeding ground survival and removals
    for (Ibreed in 1:Nbreed)
    {
      Frac <- 0; DevSpec <- 0;
      if (SBdevEst[Ibreed]>0) DevSpec <- exp(SBdevYr[Ibreed,Year])
      Nb[Ibreed,Year] <- NbS[Ibreed,Year]
      # SAFETY: Protect against division by zero
      Nb[Ibreed,Year] <- max(Nb[Ibreed,Year], ThreshPop)
      BreedK[Ibreed] <- max(BreedK[Ibreed], ThreshPop)
      
      Frac <- Nb[Ibreed,Year]/BreedK[Ibreed]
      if (Frac > 1.05) Frac <- 1.05;
      SurvOutB[Ibreed,Year] <- (1.0+exp(-10.0*(Frac-1.05)))/(1.0+exp(-10.0*(Frac-0.05)))
      SurvOutB[Ibreed,Year] <- SurvOutB[Ibreed,Year]*DevSpec
      for (Ifeed in 1:Nfeed)
      {
        Temp <- NNS[Ibreed,Ifeed,Year]
        NNS[Ibreed,Ifeed,Year] <- NNS[Ibreed,Ifeed,Year]*SurvOutB[Ibreed,Year]
        if (Year <= NremoveY)
          for (Icatch in 1:Nremove)
            if (RemoveB[Icatch]==Ibreed & RemoveI[Icatch,1]==Year)
            {
              Remove = RemoveI[Icatch,2]
              NNS[Ibreed,Ifeed,Year] <- NNS[Ibreed,Ifeed,Year]*(1.0-Remove/Temp)
            }
        SurvTot[Ibreed,Ifeed,Year] <- NNS[Ibreed,Ifeed,Year]/Temp
      }
    }

    # Feeding ground survival and removals
    for (Ifeed in 1:Nfeed)
    {
      NowNf <- NfS[Ifeed,Year];
      # SAFETY: Protect against division by zero
      NowNf <- max(NowNf, ThreshPop)
      FeedK[Ifeed] <- max(FeedK[Ifeed], ThreshPop)
      
      Frac <- NowNf/FeedK[Ifeed];
      DevSpec <- 0;
      if (SFdevEst[Ifeed]>0) DevSpec <- exp(SFdevYr[Ifeed,Year])
      SurvOutF[Ifeed,Year] <- (1.0+exp(-10.0*(Frac-1.05)))/(1.0+exp(-10.0*(Frac-0.05)))
      SurvOutF[Ifeed,Year] <- SurvOutF[Ifeed,Year]*DevSpec
      for (Ibreed in 1:Nbreed)
      {
        Temp <- NNS[Ibreed,Ifeed,Year]
        NNS[Ibreed,Ifeed,Year] <- NNS[Ibreed,Ifeed,Year]*SurvOutF[Ifeed,Year]
        if (Year <= NremoveY)
          for (Icatch in 1:Nremove)
            if (RemoveF[Icatch]==Ifeed & RemoveI[Icatch,1]==Year)
            {
              Remove = RemoveI[Icatch,2]
              NNS[Ibreed,Ifeed,Year] <- NNS[Ibreed,Ifeed,Year]*(1.0-Remove/Temp)
            }
        SurvTot[Ibreed,Ifeed,Year] <- SurvTot[Ibreed,Ifeed,Year]*NNS[Ibreed,Ifeed,Year]/Temp
        MortDiff[Ibreed,Ifeed,Year] <- Temp - NNS[Ibreed,Ifeed,Year]
        Nf[Ifeed,Year] <- Nf[Ifeed,Year] + NNS[Ibreed,Ifeed,Year]
      }
    }

    # Recruitment
    for (Ibreed in 1:Nbreed)
    {
      for (Ifeed in 1:Nfeed)
      {
        NN[Ibreed,Ifeed,Year] <- 0;
        for (Iage in 1:Nage)
        {
          if (Iage==1)
          {
            # SAFETY: Protect against division by zero in fecundity calculation
            NowNf <- max(Nf[Ifeed,Year], ThreshPop)
            FeedK_safe <- max(FeedK[Ifeed]*MultK[Ifeed,Year], ThreshPop)
            
            Frac <- NowNf/FeedK_safe
            DevSpec <- 0;
            if (FBdevEst[Ibreed]>0) DevSpec <- exp(FBdevYr[Ibreed,Year])
            Fec <- f0*(1.0-Frac/(1.0+ParAV[Ifeed,Year]))*DevSpec;
            FecOut[Ifeed,Year] <- FecOut[Ifeed,Year] + Fec;
            Recruits <- Fec*(SC/2.0)*sum(NNN[Ibreed,Ifeed,Year,MatAgeT:Nage])
            NNN[Ibreed,Ifeed,Year+1,Iage] <- Recruits
            NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year] + Recruits
          }
          else
          {
            # Age the population
            Temp <- NNN[Ibreed,Ifeed,Year,Iage-1]
            if (Iage < Nage)
            {
              NNN[Ibreed,Ifeed,Year+1,Iage] <- Temp*SA;
              NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year] + NNN[Ibreed,Ifeed,Year+1,Iage];
            }
            if (Iage==Nage)
            {
              NNN[Ibreed,Ifeed,Year+1,Nage] <- NNN[Ibreed,Ifeed,Year+1,Nage] + Temp*SA;
              NN[Ibreed,Ifeed,Year] <- NN[Ibreed,Ifeed,Year] + Temp*SA;
            }
          }
        }
      }
    }

    # Transfer feeding to breeding
    for (Ibreed in 1:Nbreed)
      for (Ifeed in 1:Nfeed)
      {
        Frac <- NN[Ibreed,Ifeed,Year]/Nf[Ifeed,Year]
        if (is.nan(Frac)) Frac <- 0;
        if (Frac==0) Frac <- 1.0e-10;
        PropMixEst[Ifeed,Ibreed,Year] <- Frac
      }

    # Influx from other populations
    if (Influx > 0)
      for (Ibreed in 1:Nbreed)
        for (Ifeed in 1:Nfeed)
          for (Iage in 1:Nage)
            NNN[Ibreed,Ifeed,Year+1,Iage] <- (1.0-Influx)*NNN[Ibreed,Ifeed,Year+1,Iage]
   }

  # Fill in the final year
  for (Ibreed in 1:Nbreed)
  {
    for (Ifeed in 1:Nfeed) NNS[Ibreed,Ifeed,Nyr+1] <- sum(NNN[Ibreed,Ifeed,Nyr+1,2:Nage])
    NbS[Ibreed,Nyr+1] <- sum(NNS[Ibreed,,Nyr+1])
    for (Icomp in 1:3)
    {
      if (Icomp==1) NfitBreed[Icomp,Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,])
      if (Icomp==2) NfitBreed[Icomp,Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,2:Nage])
      if (Icomp==3) NfitBreed[Icomp,Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,MatAgeT:Nage])
      # SAFETY: Protect against zero abundance
      NfitBreed[Icomp,Ibreed,Nyr+1] <- max(NfitBreed[Icomp,Ibreed,Nyr+1], ThreshPop)
    }
  }

  for (Ifeed in 1:Nfeed)
  {
    NfS[Ifeed,Nyr+1] <- sum(NNS[,Ifeed,Nyr+1])
    for (Icomp in 1:3)
    {
      if (Icomp==1) NfitFeed[Icomp,Ifeed,Nyr+1] <- sum(NNN[,Ifeed,Nyr+1,])
      if (Icomp==2) NfitFeed[Icomp,Ifeed,Nyr+1] <- sum(NNN[,Ifeed,Nyr+1,2:Nage])
      if (Icomp==3) NfitFeed[Icomp,Ifeed,Nyr+1] <- sum(NNN[,Ifeed,Nyr+1,MatAgeT:Nage])
      # SAFETY: Protect against zero abundance
      NfitFeed[Icomp,Ifeed,Nyr+1] <- max(NfitFeed[Icomp,Ifeed,Nyr+1], ThreshPop)
    }
  }
  
  # ========================================================================================================================
  # Likelihood evaluation
  # ========================================================================================================================

  # Setup
  PredSurv <- rep(0,NsurveyData);
  SurveySD <- rep(0,NsurveyData);
  LikeCompSurv <- rep(0,NsurveyData);
  LikeMixComp <- rep(0,NmixData);
  Qest <- rep(0,Nsurvey);
  PredMix <- rep(0,NmixData)
  PredMixOut <- array(0,dim=c(2,Nbreed,Nfeed))
  ObsMixPropI <- matrix(0,nrow=NmixData,ncol=4)
  ObsMixProp <- matrix(0,nrow=NmixData,ncol=2)

  LogLike1 <- 0;
  LogLike2a <- rep(0,2);
  LogLike2b <- rep(0,2);

  for (Iclass in 1:Nsurvey)
  {
    # Survey data (Absolute abundance estimates)
    if (SurveySeries[Iclass]==1)
    {
      for (Isurv in 1:NsurveyData)
        if (SurveyI[Isurv,5]==Iclass)
        {
          Pred = 0;
          for (Iyr in SurveyI[Isurv,1]:SurveyI[Isurv,2])
          {
            Icomp <- SurveyI[Isurv,8]+1
            if (SurveyI[Isurv,3]==1) for (Ibreed in 1:Nbreed) Pred <- Pred + NfitBreed[Icomp,Ibreed,Iyr];
            if (SurveyI[Isurv,3]==2) Pred <- Pred + NfitBreed[Icomp,SurveyI[Isurv,4],Iyr];
            if (SurveyI[Isurv,3]==3) Pred <- Pred + NfitFeed[Icomp,SurveyI[Isurv,4],Iyr];
          }
          Pred <- Pred / (SurveyI[Isurv,2]-SurveyI[Isurv,1]+1);
          # SAFETY: Protect against zero/negative predictions
          Pred <- max(Pred, ThreshPop)
          PredSurv[Isurv] <- Pred;
          SD2 <- SurveyR[Isurv,2]*SurveyR[Isurv,2];
          if (SurveyI[Isurv,7]>0) SD2 <- SD2 + AddV[SurveyI[Isurv,7]];
          SurveySD[Isurv] <- sqrt(SD2);
          LikeCompSurv[Isurv] <- -dnorm(log(SurveyR[Isurv,1]), log(Pred), sqrt(SD2), log=T);
          if (SurveyI[Isurv,3]==1) LikeCompSurv[Isurv] <- LikeCompSurv[Isurv] * WghtTotal;
          if (SurveyI[Isurv,6]==1) LogLike1 <- LogLike1 + LikeCompSurv[Isurv];
        }
    } # Absolute indices

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
            Icomp <- SurveyI[Isurv,8]+1
            if (SurveyI[Isurv,3]==1) for (Ibreed in 1:Nbreed) Pred <- Pred + NfitBreed[Icomp,Ibreed,Iyr];
            if (SurveyI[Isurv,3]==2) Pred <- Pred + NfitBreed[Icomp,SurveyI[Isurv,4],Iyr];
            if (SurveyI[Isurv,3]==3) Pred <- Pred + NfitFeed[Icomp,SurveyI[Isurv,4],Iyr];
          }
          Pred = Pred / (SurveyI[Isurv,2]-SurveyI[Isurv,1]+1);
          # SAFETY: Protect against zero/negative predictions
          Pred <- max(Pred, ThreshPop)
          PredSurv[Isurv] <- Pred;
          SD2 = SurveyR[Isurv,2]^2;
          if (SurveyI[Isurv,7]>0) SD2 <- SD2 + AddV[SurveyI[Isurv,7]];
          if (SurveyI[Isurv,6]==1) ndat <- ndat + 1.0/SD2;
          if (SurveyI[Isurv,6]==1) logq <- logq + log(SurveyR[Isurv,1]/Pred)/SD2;
        }
      # SAFETY: Protect against division by zero
      ndat <- max(ndat, 1e-10)
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
      # SAFETY: Protect against division by zero
      NbS_safe <- max(NbS[Ibreed,YearFeedBreed], ThreshPop)
      for (Ifeed in 1:Nfeed) PredMixOut[1,Ibreed,Ifeed] = NNS[Ibreed,Ifeed,YearFeedBreed]/NbS_safe;
      for (Ifeed in 1:Nfeed)
        if (ObsMixBtoFP[IdataS,Ibreed,Ifeed] > 0)
        {
          Pred <- PredMixOut[1,Ibreed,Ifeed]
          Obs <- ObsMixBtoFE[IdataS,Ibreed,Ifeed];
          SD <- ObsMixBtoFP[IdataS,Ibreed,Ifeed];
          # SAFETY: Define SD2
          SD2 <- SD^2
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
  
  # Mixing proportions (Feeding to Breeding)
  Jcnt <- 1;
  for (IdataS in 1:2)
  {
    for (Ifeed in 1:Nfeed)
    {
      # SAFETY: Protect against division by zero
      NfS_safe <- max(NfS[Ifeed,YearFeedBreed], ThreshPop)
      for (Ibreed in 1:Nbreed)
        PredMixOut[2,Ibreed,Ifeed] <- NNS[Ibreed,Ifeed,YearFeedBreed]/NfS_safe;
      for (Ibreed in 1:Nbreed)
        if (ObsMixFtoBP[IdataS,Ibreed,Ifeed] > 0)
        {
          Pred <- PredMixOut[2,Ibreed,Ifeed]
          Obs <- ObsMixFtoBE[IdataS,Ibreed,Ifeed];
          SD <- ObsMixFtoBP[IdataS,Ibreed,Ifeed];
          # SAFETY: Define SD2
          SD2 <- SD^2
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
  
  # SF devs as REs
  SFsigma <- exp(log_SFsigma)
  LogLikeSDevs <- -sum(dnorm(SFdev,mean=0,sd=SFsigma,log=T))

  # Weak penalties for stability
  Penal <- 0;
  Penal <- Penal2;
  for (Ibreed in 1:Nbreed) Penal <- Penal + 0.0001*logK[Ibreed]*logK[Ibreed];
  for (Imix in 1:NmixPar) Penal <- Penal + MixPars[Imix]*MixPars[Imix];
  logBK_sumSQ <- sum((logBK-mean(logBK))^2)
  Penal <- Penal + logBK_sumSQ
  Penal <- Penal - sum(dnorm(SBdev,0.0,1.0,log=T));
  Penal <- Penal - sum(dnorm(FBdev,0.0,1.0,log=T));
  # Penal <- Penal - sum(dnorm(SFdev,0.0,1.0,log=T));
  
  TotalK <- sum(BreedK)
  if (UseKPrior > 0)  Penal <- Penal - UseKPrior*log(1.0/(1.0+exp(Prior_int+Prior_slope*TotalK)))
  datalike <- LogLike1 + sum(LogLike2a) + sum(LogLike2b)+LogLikeSDevs;
  neglogL <-  Penal + datalike;
  # print likelihood?
  # cat(paste("Survey likelihood: ",LogLike1,"\nMixingBtF likelihood: ",sum(LogLike2a),"\nMixingFtB likelihood: ",sum(LogLike2b),
  #       "\nSurvival REs likelihood: ",LogLikeSDevs))
  #print(neglogL)
  
  # =================================================================================================================================
  # Report Outputs
  # =================================================================================================================================
  
  # Summary outputs
  LogNT<- matrix(0,nrow=3,ncol=Nyr+1);                                     # Log total numbers
  LogNb <- array(0,dim=c(3,Nbreed,Nyr+1));                                 # Breeding numbers
  LogNf <- array(0,dim=c(3,Nfeed,Nyr+1));                                  # Feeding numbers
  for (Icomp in 1:3)
  {
    for (Ibreed in 1:Nbreed) for (Iyr in 1:(Nyr+1)) LogNb[Icomp,Ibreed,Iyr] <- log(NfitBreed[Icomp,Ibreed,Iyr]);
    for (Ifeed in 1:Nfeed)  for (Iyr in 1:(Nyr+1)) LogNf[Icomp,Ifeed,Iyr] = log(NfitFeed[Icomp,Ifeed,Iyr]);
    for (Iyr in 1:(Nyr+1)) { TotalNT <- sum(NfitBreed[Icomp,,Iyr]); LogNT[Icomp,Iyr] = log(TotalNT); }
  }

  ADREPORT(LogNT);
  ADREPORT(LogNb);
  ADREPORT(LogNf);
  ADREPORT(SBdevYr);
  ADREPORT(SFdevYr);
  ADREPORT(SurvOutB);
  ADREPORT(SurvOutF);
  
  REPORT(Nb);
  REPORT(Nf);
  REPORT(NbS);
  REPORT(NfS);
  REPORT(NfitBreed)
  REPORT(NfitFeed)
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
  REPORT(NNN);
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
