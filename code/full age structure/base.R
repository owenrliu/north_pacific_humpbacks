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
  
  # ========================================================================================================================
  # Set up mixing matrix
  # ========================================================================================================================
  tMixI <- t(MixI)
  Mix <- as.vector(tMixI)
  # Fill in the estimated mixing parameters in the right spot
  Mix[tMixI>0] <- exp(MixPars)
  Mix[tMixI<=0] <- abs(tMixI)[tMixI<=0]
  # Re-form the matrix and normalize across feeding grounds
  Mix <- matrix(Mix,nrow=Nbreed,ncol=Nfeed,byrow=T)
  Mix <- Mix/rowSums(Mix)

  # ========================================================================================================================
  # Set up population projection
  # ========================================================================================================================

  # Basic demographics
  Amat <- IAmat*1.0;
  Nage <- 11
  
  # Equlibrium numbers
  Neqn <- rep(0,11); Neqn[1] <- 1; Neqn[2] <- SC;
  for (Iage in 3:Nage) Neqn[Iage] <- Neqn[Iage-1]*SA;
  Neqn[Nage] <- Neqn[Nage]/(1-SA);
  Neqn <- Neqn/sum(Neqn)
  NfecEqn <- sum(Neqn[((Amat+2):Nage)])
 
  
  # Zerbini et al. (2010): Mar Biol. 157: 1225-1236. Note that Amat is age-at-maturity (age-at-first parturition less 1)
  #fmax <- 2*((1.0+rval)^(Amat+1.0)-SA*(1.0+rval)^Amat)
  #f0 <- 2*(1.0-SA)
  #fmax <- 2*((1.0+rval)^(Amat+1.0)-SA*(1.0+rval)^Amat)/(SC*SA^(Amat))
  #f0 <- 2*(1.0-SA)/(SC*SA^(Amat))
  
  # Fecundity in the limit of zero population size
  fmax <- 2*(exp(rval*(Amat+1.0))-SA*exp(rval*Amat))/(SC*SA^(Amat))
  #fmax <- fmax/(1.0+exp(30.0*(fmax-0.99)))+0.99/(1.0+exp(-30.0*(fmax-0.99))) 
  
  # Fecundity at carrying capacity
  f0 <- 2*(1.0-SA)/(SC*SA^(Amat))
  # print(f0)
  # print(2/(NfecEqn/Neqn[1]))
  
  # Resilience parameter
  ParA <- (fmax-f0)/f0;
  #print(c(rval,fmax,f0,ParA))

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
  
  # Set up K by herd, breeding stock, and feeding ground
  # K by breeding ground and feeding ground
  BreedK <- rep(0,Nbreed); FeedK <- rep(0,Nfeed);
  # K by herd
  NNK <- matrix(0,nrow=Nbreed,ncol=Nfeed)
  for (Ibreed in 1:Nbreed)
  { BreedK[Ibreed] <- exp(logK[Ibreed]); for (Ifeed in 1:Nfeed) NNK[Ibreed,Ifeed] <- BreedK[Ibreed]*Mix[Ibreed,Ifeed]; }
  for (Ifeed in 1:Nfeed)
  { for (Ibreed in 1:Nbreed) FeedK[Ifeed] <- FeedK[Ifeed] + NNK[Ibreed,Ifeed]; }
  
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
  
  Temp1 <- rep(0,Nbreed)
  Temp2 <- rep(0,Nbreed)
  
  MinAgeT <- Amat+2
  
  # ========================================================================================================================
  # Initial year
  # ========================================================================================================================
  
  for (Ibreed in 1:Nbreed)
   {
    Nb[Ibreed,1] <- exp(logK[Ibreed])*1.0/(1+exp(logBK[Ibreed]));
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
    # Save breeding/feeding ground numbers (start of the year)
    for (Ibreed in 1:Nbreed)
     for (Ifeed in 1:Nfeed) NNS[Ibreed,Ifeed,Year] <- sum(NNN[Ibreed,Ifeed,Year,MinAgeT:Nage]);
    
    # Save the breeding ground numbers
    for (Ibreed in 1:Nbreed) NbS[Ibreed,Year] <-sum(NNS[Ibreed,,Year]);
      
    # ========================================================================================================================
    # Removals
    # ========================================================================================================================
    
    for (Ibreed in 1:Nbreed)
     {
      Nb0P <- sum(NNN[Ibreed,,1,MinAgeT:Nage])
      Temp1[Ibreed] <- Nb0P - CatchB[Year,Ibreed];
      # Trick to avoid negative population sizes
      MultC <- 1.0/(1.0+exp(-30.0*(Temp1[Ibreed]-ThreshPop)))
      Temp1[Ibreed] <- ThreshPop + MultC * Temp1[Ibreed];
      # Allocate breeding ground catches to "herd" (Eqn B.3)
      for (Ifeed in 1:Nfeed) NNN[Ibreed,Ifeed,Year,MinAgeT:Nage] <- NNN[Ibreed,Ifeed,Year,MinAgeT:Nage]*Temp1[Ibreed]/Nb0P;
     }
  
    # Allow for straying before feeding ground catches (Eqn B.7)
    if (6==9)
    for (Ibreed in 1:Nbreed)
     {
      Nstray <- rep(0,Nfeed)
      for (Ifeed in 1:Nfeed)
       if (MixI[Ibreed,Ifeed] != 0)
        {
         print("do strays")
         if (Ifeed > 1)
          if(MixI[Ibreed,Ifeed-1] !=0)
           {
            Depl1 <- sum(NNN[Ibreed,Ifeed,Year,MinAgeT:Nage])/NNK[Ibreed,Ifeed];
            Depl2 <- sum(NNN[Ibreed,Ifeed-1,Year,MinAgeT:Nage])/NNK[Ibreed,Ifeed-1];
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
            Depl1 <- sum(NNN[Ibreed,Ifeed,Year,MinAgeT:Nage])/NNK[Ibreed,Ifeed];
            Depl2 <- sum(NNN[Ibreed,Ifeed+1,Year,MinAgeT:Nage])/NNK[Ibreed,Ifeed+1];
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
    for (Ifeed in 1:Nfeed) NfS[Ifeed,Year] <- sum(NNN[,Ifeed,Year,MinAgeT:Nage]);
  
    # Remove feeding ground numbers (Eqn B.2)
    CBreed <- rep(0,Nbreed);
    for (Ifeed in 1:Nfeed)
     {
      for (Ibreed in 1:Nbreed) Temp2[Ibreed] <- sum(NNN[Ibreed,Ifeed,Year,MinAgeT:Nage])
      Nf[Ifeed,Year] <- sum(Temp2)
      PropMixEst[Ifeed,,Year] <- Temp2/Nf[Ifeed,Year];
      CBreed[Ibreed] <- CBreed[Ibreed] + PropMixEst[Ifeed,Ibreed,Year]*CatchF[Year,Ifeed];
      for (Ibreed in 1:Nbreed)
       if (MixI[Ibreed,Ifeed] != 0)
        {
         Temp3 <- Temp2[Ibreed] - PropMixEst[Ifeed,Ibreed,Year]*CatchF[Year,Ifeed];
         # Trick to avoid negative population sizes
         MultC <- 1.0/(1.0+exp(-30.0*(Temp3-ThreshPop)))
         Temp3 <- ThreshPop + MultC * Temp3;
         Temp3 <- sqrt(Temp3*Temp3)
         NNN[Ibreed,Ifeed,Year,MinAgeT:Nage] <- NNN[Ibreed,Ifeed,Year,MinAgeT:Nage]*Temp3/Temp2[Ibreed];
        }
      Nf[Ifeed,Year] <-  sum(NNN[,Ifeed,Year,MinAgeT:Nage]);
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
        NNN[Ibreed,6,Year,MinAgeT:Nage] <- NNN[Ibreed,6,Year,MinAgeT:Nage] + NNN[Ibreed,4,Year,MinAgeT:Nage]*Influx + NNN[Ibreed,5,Year,MinAgeT:Nage]*Influx;
        NNN[Ibreed,4,Year,MinAgeT:Nage] <- NNN[Ibreed,4,Year,MinAgeT:Nage]*(1.0-Influx);
        NNN[Ibreed,5,Year,MinAgeT:Nage] <- NNN[Ibreed,5,Year,MinAgeT:Nage]*(1.0-Influx);
       }
     }
    
    # Update feeding ground numbers
    for (Ifeed in 1:Nfeed)
    { Nf[Ifeed,Year] <- 0; for (Ibreed in 1:Nbreed) Nf[Ifeed,Year] <- sum(NNN[,Ifeed,Year,MinAgeT:Nage]); }

    
    # ========================================================================================================================
    # Survival and Fecundity
    # ========================================================================================================================
    
    # Density-dependence is on feeding ground numbers
    for (Ibreed in 1:Nbreed)
     {
      # Old
      #Ioffset <- Year-IAmat+1-TimeLag; if (Ioffset < 1) Ioffset <- 1;
      # NEW
      Ioffset <- Year-(IAmat+1)+1-TimeLag; if (Ioffset < 1) Ioffset <- 1;

      for (Ifeed in 1:Nfeed)
       {
        # Survival (Eqn B.6)
        LogitSA <- log(1.0/SA-1.0);
        SurvOutB[Ibreed,Year] <- 1.0/(1.0+exp(LogitSA+SBdevYr[Ibreed,Year]*Sigma_SBdev));
        SurvOutF[Ifeed,Year] <- 1.0/(1.0+exp(LogitSA+SFdevYr[Ifeed,Year]*Sigma_SFdev));
        if (StochSopt==0) SAuse <- SurvOutB[Ibreed,Year];
        if (StochSopt==1) SAuse <- SurvOutF[Ifeed,Year];
        if (StochSopt==0) SurvJuv <-  SurvOutB[Ibreed,Year]*SC/SA
        if (StochSopt==1) SurvJuv <-  SurvOutF[Ifeed,Year]*SC/SA
        # Fecundity (Eqn B.4)
        Depl <- Nf[Ifeed,Ioffset]/(FeedK[Ifeed]*MultK[Ifeed,Ioffset]);
        Term1 <- 1.0 - Depl;
        Term1 <- f0*(1.0+ParAV[Ifeed,Ioffset]*Term1);
        Term1 <- 0.0001+(Term1-0.0001)/(1+exp(-30.0*Term1)) 
        LogitFec <- log(1.0/Term1-1.0);
        Term1 <- 1.0/(1.0+exp(LogitFec+FBdevYr[Ibreed,Year]*Sigma_FBdev));
        FecOut[Ifeed,Year] <- Term1;
        # Check this
        SurvTot[Ibreed,Ifeed,Year] <- SAuse*NN[Ibreed,Ifeed,Year]/NNS[Ibreed,Ifeed,Year];
        MortDiff[Ibreed,Ifeed,Year] <- (SAuse-SA)*NN[Ibreed,Ifeed,Year];
        
        # ========================================================================================================================
        # Next Year's Starting Numbers
        # ========================================================================================================================
        
        # Juveniles and Adults (2+)
        NNN[Ibreed,Ifeed,Year+1,2] <- NNN[Ibreed,Ifeed,Year,1] * SurvJuv
        for (Iage in 3:Nage)
         NNN[Ibreed,Ifeed,Year+1,Iage] <- NNN[Ibreed,Ifeed,Year,Iage-1] * SAuse
        NNN[Ibreed,Ifeed,Year+1,Nage] <-  NNN[Ibreed,Ifeed,Year+1,Nage] + NNN[Ibreed,Ifeed,Year,Nage] * SAuse
        
        # Recruits - Feeding ground density-dependence but recruitment related to recent (Eqn B.1)
        if (DensDepOpt==0)
         {
          NNN[Ibreed,Ifeed,Year+1,1] <-  0.5*Term1*sum(NNN[Ibreed,Ifeed,Year,MinAgeT:Nage]);
         }

        # Recruits - Feeding ground density-dependence but recruitment related to unfished (Eqn B9.a)
        if (DensDepOpt==1)
         {
          NNN[Ibreed,Ifeed,Year+1,1] <-  0.5*Term1*Nb[Ibreed,Ioffset]*Mix[Ibreed,Ifeed];
         }

        # Recruits - Feeding ground density-dependence but recruitment related to unfished & recent (Eqn B.9b)
        if (DensDepOpt==2)
         {
          NNN[Ibreed,Ifeed,Year+1,1] <-  0.5*0.5*Term1*(Nb[Ibreed,Ioffset]*Mix[Ibreed,Ifeed]+sum(NNN[Ibreed,Ifeed,Year,MinAgeT:Nage]));
        }
        ParAV[Ifeed,Year+1] <- ParAV[Ifeed,Ioffset]*exp(inert_par * log(FeedK[Ifeed]/Nf[Ifeed,Ioffset]));
      } # close Ifeed
      # Update breeding ground numbers after feeding ground processes
      Nb[Ibreed,Year+1] <-  sum(NNN[Ibreed,,Year+1,MinAgeT:Nage]);
     } # close Ibreed
    } # close Year   
  
  # Save herd numbers for the start of the final year (Nyr+1)
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
  for (Imix in 1:NmixPar) Penal <- Penal + MixPars[Imix]*MixPars[Imix];
  logBK_sumSQ <- sum((logBK-mean(logBK))^2);
  Penal <- Penal + logBK_sumSQ;
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