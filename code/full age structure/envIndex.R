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
  
  # Equlibrium numbers
  Neqn <- rep(0,Nage); Neqn[1] <- 1; Neqn[2] <- SC;
  for (Iage in 3:Nage) Neqn[Iage] <- Neqn[Iage-1]*SA;
  Neqn[Nage] <- Neqn[Nage]/(1-SA);
  Neqn <- Neqn/sum(Neqn)
  NfecEqn <- sum(Neqn[MinAgeD:Nage]) # Proportional 1+ abundance
 
  # Zerbini et al. (2010): Mar Biol. 157: 1225-1236. Note that Amat is age-at-maturity (age-at-first parturition less 1)
  fmax <- 2*(exp(rval*(Amat+1.0))-SA*exp(rval*Amat))/(SC*SA^(Amat))
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
  
  # ========================================================================================================================
  # Environmental Index of Survival
  # ========================================================================================================================
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
  
  # calculate index
  env_index <- omega_sst*sst+omega_mld*mld
  # Difference between SFdevs and environmental index
  IenvStart=Nyr-(Yr2-YrSDevs)
  envPenal <- c(SFdevYr[,IenvStart:Nyr]-env_index)
  # ========================================================================================================================
  
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
      NbS[Ibreed,Year] <-sum(NNS[Ibreed,,Year]);
      NfitBreed[1,Ibreed,Year] <- sum(NNN[Ibreed,,Year,1:Nage])
      NfitBreed[2,Ibreed,Year] <- sum(NNN[Ibreed,,Year,2:Nage])
      NfitBreed[3,Ibreed,Year] <- sum(NNN[Ibreed,,Year,MatAgeT:Nage])
    }
      
    # ========================================================================================================================
    # Removals
    # ========================================================================================================================
    
    # Remove breeding ground catches (Eqn B.3)
    for (Ibreed in 1:Nbreed)
    {
      Nb0P <- sum(NNN[Ibreed,,Year,MinAgeC:Nage])
      Temp1 <- Nb0P - CatchB[Year,Ibreed];
      # Trick to avoid negative population sizes
      MultC <- 1.0/(1.0+exp(-30.0*(Temp1-ThreshPop)))
      Temp1 <- ThreshPop + MultC * Temp1;
      # Allocate breeding ground catches to "herd" 
      for (Ifeed in 1:Nfeed) 
        if (MixI[Ibreed,Ifeed]!=0) NNN[Ibreed,Ifeed,Year,MinAgeC:Nage] <- NNN[Ibreed,Ifeed,Year,MinAgeC:Nage]*Temp1/Nb0P;
    }
    
    # Allow for straying before feeding ground catches (Eqn B.9)
    for (Ibreed in 1:Nbreed)
    {
      Nstray <- matrix(0,nrow=Nfeed,ncol=Nage)
      for (Ifeed in 1:Nfeed)
        if (MixI[Ibreed,Ifeed] != 0)
        {
          if (Ifeed > 1)
            if(MixI[Ibreed,Ifeed-1] !=0)
            {
              Depl1 <- sum(NNN[Ibreed,Ifeed,Year,2:Nage])/NNK[Ibreed,Ifeed];
              Depl2 <- sum(NNN[Ibreed,Ifeed-1,Year,2:Nage])/NNK[Ibreed,Ifeed-1];
              DeplRat <- Depl1/Depl2;
              
              Ver1 <- 1.0/(1.0+exp(-100.0*(DeplRat-1.0)));
              Ver2 <- 1.0/(1.0+exp(100.0*(DeplRat-2.0)));
              Ver3 <- 1.0/(1.0+exp(-100.0*(DeplRat-2.0)));
              Stray <- StrayBase*(DeplRat-1.0)*Ver1*Ver2+StrayBase*Ver3;
              Nstray[Ifeed,] <- Nstray[Ifeed] - Stray*NNN[Ibreed,Ifeed,Year,];
              Nstray[Ifeed-1,] <- Nstray[Ifeed-1] + Stray*NNN[Ibreed,Ifeed,Year,];
            }
          if (Ifeed < Nfeed)
            if(MixI[Ibreed,Ifeed+1] !=0)
            {
              Depl1 <- sum(NNN[Ibreed,Ifeed,Year,2:Nage])/NNK[Ibreed,Ifeed];
              Depl2 <- sum(NNN[Ibreed,Ifeed+1,Year,2:Nage])/NNK[Ibreed,Ifeed+1];
              DeplRat <- Depl1/Depl2;
              
              Ver1 = 1.0/(1.0+exp(-100.0*(DeplRat-1.0)));
              Ver2 = 1.0/(1.0+exp(100.0*(DeplRat-2.0)));
              Ver3 = 1.0/(1.0+exp(-100.0*(DeplRat-2.0)));
              Stray = StrayBase*(DeplRat-1.0)*Ver1*Ver2+StrayBase*Ver3;
              Nstray[Ifeed,] <- Nstray[Ifeed,] - Stray*NNN[Ibreed,Ifeed,Year,];
              Nstray[Ifeed+1,] <- Nstray[Ifeed+1,] + Stray*NNN[Ibreed,Ifeed,Year,];
            }
        }
      for (Ifeed in 1:Nfeed) NNN[Ibreed,Ifeed,Year,] <- NNN[Ibreed,Ifeed,Year,] + Nstray[Ifeed,];
    }  # Each breeding stock
    
    # Save the feeding ground numbers
    for (Ifeed in 1:Nfeed) NfS[Ifeed,Year] <- sum(NNN[,Ifeed,Year,2:Nage]);
    
    # Remove feeding ground numbers (Eqn B.2)
    Temp2 <- rep(0,Nbreed)
    for (Ifeed in 1:Nfeed)
    {
      for (Ibreed in 1:Nbreed) Temp2[Ibreed] <- sum(NNN[Ibreed,Ifeed,Year,MinAgeC:Nage])
      Nf0P <- sum(Temp2)
      Temp1 <- Nf0P - CatchF[Year,Ifeed];
      # Trick to avoid negative population sizes
      MultC <- 1.0/(1.0+exp(-30.0*(Temp1-ThreshPop)))
      Temp1 <- ThreshPop + MultC * Temp1;
      # Allocate feeding ground catches to "herd" 
      for (Ibreed in 1:Nbreed) 
        if (MixI[Ibreed,Ifeed]!=0) NNN[Ibreed,Ifeed,Year,MinAgeC:Nage] <- NNN[Ibreed,Ifeed,Year,MinAgeC:Nage]*Temp1/Nf0P;
    }
 
    # Influx (part 1) [Not used for this assessment]
    for (Ifeed in 1:Nfeed) MultK[Ifeed,Year+1] <- MultK[Ifeed,Year];
    if (Year==2014-Yr1+1)
    {
      MultK[4,Year+1] <- MultK[4,Year]*(1.0-Influx);
      MultK[5,Year+1] <- MultK[5,Year]*(1.0-Influx);
      MultK[6,Year+1] <- (MultK[6,Year]*Nf[6,1] + MultK[4,Year]*Influx*Nf[4,1]+MultK[5,Year]*Influx*Nf[5,1])/(MultK[6,Year]*Nf[6,1]);
      for (Ibreed in 1:Nbreed)
      {
        NNN[Ibreed,6,Year,MinAgeC:Nage] <- NNN[Ibreed,6,Year,MinAgeC:Nage] + NNN[Ibreed,4,Year,MinAgeC:Nage]*Influx + NNN[Ibreed,5,Year,MinAgeC:Nage]*Influx;
        NNN[Ibreed,4,Year,MinAgeC:Nage] <- NNN[Ibreed,4,Year,MinAgeC:Nage]*(1.0-Influx);
        NNN[Ibreed,5,Year,MinAgeC:Nage] <- NNN[Ibreed,5,Year,MinAgeC:Nage]*(1.0-Influx);
      }
    }
    
    # Update feeding ground numbers
    for (Ifeed in 1:Nfeed) 
    {
      Nf[Ifeed,Year] <- sum(NNN[,Ifeed,Year,2:Nage]);
      NfitFeed[1,Ifeed,Year] <- sum(NNN[,Ifeed,Year,1:Nage])
      NfitFeed[2,Ifeed,Year] <- sum(NNN[,Ifeed,Year,2:Nage])
      NfitFeed[3,Ifeed,Year] <- sum(NNN[,Ifeed,Year,MatAgeT:Nage])
    }
    
    # ========================================================================================================================
    # Survival and Fecundity
    # ========================================================================================================================
    
    # Density-dependence is on feeding ground numbers
    for (Ibreed in 1:Nbreed){
      for (Ifeed in 1:Nfeed){
        if (MixI[Ibreed,Ifeed] != 0){
          
          # Survival (Eqn B.8)
          LogitSA <- log(1.0/SA-1.0);
          SurvOutB[Ibreed,Year] <- 1.0/(1.0+exp(LogitSA+SBdevYr[Ibreed,Year]*Sigma_SBdev));
          # SurvOutF[Ifeed,Year] <- 1.0/(1.0+exp(LogitSA+SFdevYr[Ifeed,Year]*Sigma_SFdev));
          SurvOutF[Ifeed,Year] <- 1.0/(1.0+exp(LogitSA+SFdevYr[Ifeed,Year])); # no sigma in this version
          if (StochSopt==0) SAuse <- SurvOutB[Ibreed,Year];
          if (StochSopt==1) SAuse <- SurvOutF[Ifeed,Year];
          LogitSC <- log(1.0/SC-1.0);
          SurvSCB <- 1.0/(1.0+exp(LogitSC+SBdevYr[Ibreed,Year]*Sigma_SBdev));
          # SurvSCF <- 1.0/(1.0+exp(LogitSC+SFdevYr[Ifeed,Year]*Sigma_SFdev));
          SurvSCF <- 1.0/(1.0+exp(LogitSC+SFdevYr[Ifeed,Year])); # no sigma in this version
          if (StochSopt==0) SurvJuv <-  SurvSCB
          if (StochSopt==1) SurvJuv <-  SurvSCF
          
          # Fecundity (Density-dependent)
          Depl <- Nf[Ifeed,Year]/(FeedK[Ifeed]*MultK[Ifeed,Year]);
          Term1 <- f0*(1.0+ParAV[Ifeed,Year]*(1.0-Depl));
          Term1 <- 0.0001+(Term1-0.0001)/(1+exp(-30.0*Term1)) 
          LogitFec <- log(1.0/Term1-1.0);
          Term1 <- 1.0/(1.0+exp(LogitFec+FBdevYr[Ibreed,Year]*Sigma_FBdev));
          FecOut[Ifeed,Year] <- Term1;
          # Check this
          SurvTot[Ibreed,Ifeed,Year] <- SAuse*sum(NNN[Ibreed,Ifeed,Year,1:Nage])/NNS[Ibreed,Ifeed,Year];
          MortDiff[Ibreed,Ifeed,Year] <- (SAuse-SA)*sum(NNN[Ibreed,Ifeed,Year,1:Nage]);
          
          # ========================================================================================================================
          # Next Year's Starting Numbers
          # ========================================================================================================================
          # Ages 2+
          NNN[Ibreed,Ifeed,Year+1,2] <- NNN[Ibreed,Ifeed,Year,1] * SurvJuv
          for (Iage in 3:Nage)
            NNN[Ibreed,Ifeed,Year+1,Iage] <- NNN[Ibreed,Ifeed,Year,Iage-1] * SAuse
          NNN[Ibreed,Ifeed,Year+1,Nage]  <- NNN[Ibreed,Ifeed,Year+1,Nage] + NNN[Ibreed,Ifeed,Year,Nage] * SAuse
          
          # Recruits
          # Feeding ground density-dependence but recruitment related to recent (Eqn B.7a)
          if (DensDepOpt==0) NNN[Ibreed,Ifeed,Year+1,1] <-  0.5*Term1*sum(NNN[Ibreed,Ifeed,Year,MatAgeT:Nage]);
          # Feeding ground density-dependence but recruitment related to unfished (Eqn B7.b)
          if (DensDepOpt==1)  NNN[Ibreed,Ifeed,Year+1,1] <-  0.5*Term1*Nb[Ibreed,Year+1]*Mix[Ibreed,Ifeed];
          # Feeding ground density-dependence but recruitment related to unfished & recent (Eqn B.7c)
          if (DensDepOpt==2) NNN[Ibreed,Ifeed,Year+1,1] <-  0.5*0.5*Term1*(Nb[Ibreed,Year+1]*Mix[Ibreed,Ifeed]+sum(NNN[Ibreed,Ifeed,Year,MatAgeT:Nage]));
          ParAV[Ifeed,Year+1] <- ParAV[Ifeed,Year]*exp(inert_par * log(FeedK[Ifeed]/Nf[Ifeed,Year]));
          
        }} # close Ifeed
      
      # Update breeding ground numbers after feeding ground processes
      Nb[Ibreed,Year+1] <-  sum(NNN[Ibreed,,Year+1,MatAgeT:Nage])
      
     } # close Ibreed
    } # close Year   
  
  # Save herd numbers for the start of the final year (Nyr+1)
  for (Ibreed in 1:Nbreed)
    for (Ifeed in 1:Nfeed) NNS[Ibreed,Ifeed,Nyr+1] <- sum(NNN[Ibreed,Ifeed,Nyr+1,2:Nage]);
  # Save the breeding ground numbers
  for (Ibreed in 1:Nbreed) NbS[Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,2:Nage]);
  # Save the breeding ground numbers
  for (Ifeed in 1:Nfeed) NfS[Ifeed,Nyr+1] <-sum(NNN[,Ifeed,Nyr+1,2:Nage]);
  
  # Save breeding ground summaries
  for (Ibreed in 1:Nbreed) 
  {
    NfitBreed[1,Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,1:Nage])
    NfitBreed[2,Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,2:Nage])
    NfitBreed[3,Ibreed,Nyr+1] <- sum(NNN[Ibreed,,Nyr+1,MatAgeT:Nage])
  } 
  # Save feeding ground summaries
  for (Ifeed in 1:Nfeed) 
  {
    NfitFeed[1,Ifeed,Nyr+1] <- sum(NNN[,Ifeed,Nyr+1,1:Nage])
    NfitFeed[2,Ifeed,Nyr+1] <- sum(NNN[,Ifeed,Nyr+1,2:Nage])
    NfitFeed[3,Ifeed,Nyr+1] <- sum(NNN[,Ifeed,Nyr+1,MatAgeT:Nage])
  }
  
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
            Icomp <- SurveyI[Isurv,8]+1
            # Icomp <- SurveyI[Isurv,8]
            if (SurveyI[Isurv,3]==1) for (Ibreed in 1:Nbreed) Pred <- Pred + NfitBreed[Icomp,Ibreed,Iyr];
            if (SurveyI[Isurv,3]==2) Pred <- Pred + NfitBreed[Icomp,SurveyI[Isurv,4],Iyr];
            if (SurveyI[Isurv,3]==3) Pred <- Pred + NfitFeed[Icomp,SurveyI[Isurv,4],Iyr];
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
          Pred <- PredMixOut[1,Ibreed,Ifeed]
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
  
  # Mixing proportions (Feeding to Breeding)
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
          Pred <- PredMixOut[2,Ibreed,Ifeed]
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
  logBK_sumSQ <- sum((logBK-mean(logBK))^2)
  Penal <- Penal + logBK_sumSQ
  Penal <- Penal - sum(dnorm(SBdev,0.0,1.0,log=T));
  Penal <- Penal - sum(dnorm(FBdev,0.0,1.0,log=T));
  # Penal <- Penal - sum(dnorm(SFdev,0.0,1.0,log=T)); # Taking this out because driving directly with environment
  Penal <- Penal - sum(dnorm(envPenal,0.0,1.0,log=T));
  TotalK <- sum(BreedK)
  if (UseKPrior > 0)  Penal <- Penal - UseKPrior*log(1.0/(1.0+exp(Prior_int+Prior_slope*TotalK)))
  datalike <- LogLike1 + sum(LogLike2a) + sum(LogLike2b);
  neglogL <-  Penal + datalike;
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