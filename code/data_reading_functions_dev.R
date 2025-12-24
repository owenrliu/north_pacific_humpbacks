# Function to read the humpback bycatch data
ReadCatches <- function(DatFile,BreedingOpt,FeedingOpt,Nbreed,Nfeed,BreedNames,FeedNames,Yr1,Yr2,CatchSer,ByCatchFile="Null.csv")
 {
  # Read in the catches
  #=====================
  Index <- which(DatFile[,1]=="Total_areas_in_catch_files:"); MaxAreas <- as.numeric(DatFile[Index,2])
  Index <- which(DatFile[,1]=="First_year_with_catches:"); Y1Cat <- as.numeric(DatFile[Index,2])
  if (CatchSer=="L") Offset <- -2
  if (CatchSer=="B") Offset <- -1
  if (CatchSer=="H") Offset <-  0
  
  Index <- which(DatFile[,1]=="Columns_with_breeding_gound_catches:" & DatFile[,2]==BreedingOpt)+1; 
  CatchLinkB <- matrix(0,nrow=Nbreed,ncol=MaxAreas)
  rownames(CatchLinkB) <- BreedNames
  for (Ibreed in 1:Nbreed)
    for (Iarea in 1:MaxAreas) CatchLinkB[Ibreed,Iarea] <- as.numeric(DatFile[Index+Ibreed,Iarea])
  Index <- which(DatFile[,1]=="Columns_with_feeding_gound_catches:" & DatFile[,2]==FeedingOpt)+1; 
  CatchLinkF <- matrix(0,nrow=Nfeed,ncol=MaxAreas)
  rownames(CatchLinkF) <- FeedNames
  for (Ifeed in 1:Nfeed)
    for (Iarea in 1:MaxAreas) CatchLinkF[Ifeed,Iarea] <- as.numeric(DatFile[Index+Ifeed,Iarea])
  Tnames <- c("Asia","Hawaii","Mexico","MX_AR","MX_ML","Central Am","RUS+WAL","EAL+BER","WGOA","NGOA","SEA-NBC","SBC+WA","OR+CA")
  colnames(CatchLinkF) <- Tnames
  colnames(CatchLinkB) <- Tnames

  #Commercial, bycatch and ship trikes
  CatchCom <- read.csv(here('data',"CatchCom.csv"),skip=1,head=F)[,-1]
  CatchBy <- read.csv(paste0(here('data'),'/',ByCatchFile),skip=1,head=F)[,-1]

  # Now allocate catches by breeding and feeding groups
  CatchB <- matrix(0,nrow=(Yr2-Yr1+1),ncol=Nbreed); CatchF <- matrix(0,nrow=(Yr2-Yr1+1),ncol=Nfeed);
  for (Iyr in Yr1:Yr2)
  {
    Jyr <- Iyr-Yr1+1; Kyr <- Iyr-Y1Cat+1
    for (Ibreed in 1:Nbreed)
      for (Iarea in 1:MaxAreas)
        if (CatchLinkB[Ibreed,Iarea]==1) CatchB[Jyr,Ibreed] = CatchB[Jyr,Ibreed]+ CatchCom[Kyr,3*Iarea+Offset]+CatchBy[Kyr,3*Iarea+Offset]
    for (Ifeed in 1:Nfeed)
      for (Iarea in 1:MaxAreas)
        if (CatchLinkF[Ifeed,Iarea]==1) CatchF[Jyr,Ifeed] = CatchF[Jyr,Ifeed]+ CatchCom[Kyr,3*Iarea+Offset]+CatchBy[Kyr,3*Iarea+Offset]
  }
  Outs <- NULL
  Outs$CatchF <- CatchF
  Outs$CatchB <- CatchB
  return(Outs)
  
 }

# ===========================================================================================
# Function to read the humpback survey data
ReadSurveyData <- function(Code,BreedingOpt,FeedingOpt,BreedNames,FeedNames,Yr1,Yr2,SensTest="")
 {

  Surveys <- read.csv(here('data',"SurveyAll.duringWorkshop1.csv"),fill=T,comment.char="?",header=T,row.names=NULL)[,1:12]
  colnames <- c("Year1","Year2","Estimate","CV","Area","Rel","Use","Add.cv","Hypothesis","Class","SensUse","Reference")
  colnames(Surveys) <- colnames
  #print(head(Surveys))
  if (FullDiag==T) cat("Initial number of lines",length(Surveys[,1]),"\n")

  # Extract the surveys to use
  #Index <- which (Surveys$Use == "Yes")
  Index <- which ((Surveys$Use == "Yes" | Surveys$Use == SensTest) & Surveys$SensUse != SensTest)
  Surveys <- Surveys[Index,]
  SearchString <- c("All",paste0(BreedingOpt,FeedingOpt),paste0(BreedingOpt," only"),paste0(FeedingOpt," only"))
  if (FullDiag==T) print(SearchString)
  Index <- which (Surveys$Hypothesis %in% SearchString)
  Surveys <- Surveys[Index,]
  if (FullDiag==T) cat("Final number of lines",length(Surveys[,1]),"\n")
  
  # Deal with incorrect names
  if (BreedingOpt=="B1")
  { Index <- which(Surveys$Area=="MX_ML"); Surveys$Area[Index] <- "Mexico"; }
  
  Index <- (Surveys$Add.cv == "No" | Surveys$Add.cv == "Maybe")
  Surveys$Add.cv[Index] <- 0
  
  # check if the filepath exists; if not, create it
  if(!dir.exists(here('Diags'))) {dir.create(here('Diags'))}
  if(!dir.exists(here('Diags','SurveyUse'))) {dir.create(here('Diags','SurveyUse'))}

  write.csv(Surveys,file=paste0(here('Diags','SurveyUse'),'/',Code,SensTest,".csv"),row.names=F)
  Surveys$Year1 <- as.numeric(Surveys$Year1)
  Surveys$Year2 <- as.numeric(Surveys$Year2)
  Surveys$Estimate <- as.numeric(Surveys$Estimate)
  Surveys$CV <- as.numeric(Surveys$CV)
  Surveys$Add.cv <- as.numeric(Surveys$Add.cv)
  
  # Now process the data
  NextraCV1 <- 0
  SurveyI<- NULL; SurveyR <- NULL
  for (Iline in 1:length(Surveys[,1]))
  {
    Type <- -1
    if (Surveys$Area[Iline]=="Total") { Type <- 1; Area <- -1 }
    if (Surveys$Area[Iline] %in% BreedNames) {Type <- 2; Area <- which(Surveys$Area[Iline]==BreedNames)-1 }
    if (Surveys$Area[Iline] %in% FeedNames)  {Type <- 3; Area <- which(Surveys$Area[Iline]==FeedNames)-1 }
    if (Type != -1)
    {
      if (Surveys$Add.cv[Iline] > NextraCV1) NextraCV1 <- Surveys$Add.cv[Iline]
      Vec1 <- c(Surveys$Year1[Iline]-Yr1,Surveys$Year2[Iline]-Yr1,Type,Area,Surveys$Class[Iline]-1,1,Surveys$Add.cv[Iline])
      Vec2 <- c(Surveys$Estimate[Iline],Surveys$CV[Iline])
      SurveyI <- rbind(SurveyI,Vec1);SurveyR <- rbind(SurveyR,Vec2);
    }
    else
    {
      print(Surveys[Iline,])
      cat("Error, line:",Iline,"\n")  
    }
  }
  colnames(SurveyI) <- c("Year1","Year2","Type","Area","Class","?","Add.cv")
  
  NsurveyData <- length (SurveyI[,1])
  # Number of survey series
  NsurveySeries <- length(unique(SurveyI[,5]))
  SurveySeries <- c(1,rep(2,NsurveySeries-1))
  
  Outs <- NULL
  Outs$SurveyI <- SurveyI
  Outs$SurveyR <- SurveyR
  Outs$NsurveyData <- NsurveyData
  Outs$NsurveySeries <- NsurveySeries
  Outs$SurveySeries <- SurveySeries
  Outs$NextraCV1 <- NextraCV1
  return(Outs)
 }

# ===========================================================================================
# Function to read the breeding-feeding mixing data
ReadMixingData <- function(Code,DatFile,Yr1,Yr2,Nbreed,Nfeed,BreedNames,FeedNames,MixWeights,SensTest="",MaxN=-1000)
{
  # Maximum effective sample size
  #MaxN <- 100
  
  # Read in mixing data
  # =================== 
  MixFile1 <- read.csv(here('data',"Genetics_mixing_data_allscenarios_table_Long.csv"))
  MixFile2 <- read.csv(here('data',"Mark-Recapture_mixing_data_allscenarios_table_Long.csv"))
  MixFile <- rbind(MixFile1,MixFile2)
  Index <- which(MixFile$Hypothesis==Code)
  MixFile <- MixFile[Index,]
  
  Type <- c("Mark-Recapture","Genetics")
  NmixData <- 0; Nprop <-rep(0,4)
  Index <- which(DatFile[,2]=="Minimum_CV_for_the_mixing_data"); MinMixCV <- as.numeric(DatFile[Index+1,1])
  Index <- which(DatFile[,2]=="Minimum_SD_for_the_mixing_data"); MinMixSD <- as.numeric(DatFile[Index+1,1])
  
  # check if the filepath exists; if not, create it
  if(!dir.exists(here('Diags'))) {dir.create(here('Diags'))}
  if(!dir.exists(here('Diags','MixUse'))) {dir.create(here('Diags','MixUse'))}
  FileName <- paste0(here('Diags'),"/MixUse/",Code,SensTest,".csv")
  
  # Now pull out results (breeding to feeding)
  ObsMixBtoFE <- array(0,dim=c(2,Nbreed,Nfeed))
  ObsMixBtoFP <- array(0,dim=c(2,Nbreed,Nfeed))
  ObsMixBtoFO <- matrix(0,nrow=2,ncol=Nbreed)
  write("Breeding to feeding",FileName)
  if (FullDiag==T) print("Breeding to feeding mixing")
  for (Itype in 1:2)
  {
    write(Type[Itype],FileName,append=T)
    Index <- MixFile$Method == Type[Itype] & MixFile$Direction == "BreedingtoFeeding"
    MixFile2 <- MixFile[Index,]
    for (Ibreed in 1:length(BreedNames))
    {
      Total <- 0; Top <- 0; Bot <- 0
      for (Ifeed in 1:length(FeedNames))
      {
        Index <- MixFile2$Feeding %in% FeedNames[Ifeed] & MixFile2$Breeding %in% BreedNames[Ibreed]
        Est <- as.numeric(MixFile2$Estimate[Index])
        CV <- as.numeric(MixFile2$CV[Index])
        SD <- CV*Est
        if (is.na(CV)) 
        {
          Lower95 <- as.numeric(MixFile2$lowerCI[Index])
          Upper95 <- as.numeric(MixFile2$upperCI[Index])
          CV <- (Upper95-Lower95)/(2.0*1.96)/Est
          if (Est<=0) CV <- 0
          SD <- CV*Est
        }
        Top <- Top + Est*(1-Est)
        Bot <- Bot + SD^2
        if (CV < MinMixCV) CV <- MinMixCV
        SD <- Est*CV
        if (Est >0 & SD < MinMixSD) SD <- MinMixSD
        if (Est > 0 || SD > 0) { NmixData <- NmixData + 1; Nprop[Itype] <- Nprop[Itype] + 1}
        Total <- Total + Est
        
        #cat(Itype,Ibreed,Ifeed,Est,CV,SD,"\n")
        ObsMixBtoFE[Itype,Ibreed,Ifeed] <- Est
        ObsMixBtoFP[Itype,Ibreed,Ifeed] <- SD
      }# Feed
      ObsMixBtoFO[Itype,Ibreed] <- min(MaxN,Top/Bot*MixWeights[Itype])
      #print(Total)
    } # Breed 
    write("Estimates",FileName,append=T)
    write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixBtoFE[Itype,Ibreed,]))
      xx <- c(xx,",#",unlist(BreedNames[Ibreed]))
      write(xx,FileName,append=T,ncol=Nfeed+2,sep=",")
    }   
    write("SD",FileName,append=T)
    write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixBtoFP[Itype,Ibreed,]))
      xx <- c(xx,",#",unlist(BreedNames[Ibreed]))
      write(xx,FileName,append=T,ncol=Nfeed+2,sep=",")
    }    
  }  # Type
  write("EffN",FileName,append=T)
  write(t(unlist(BreedNames)),FileName,append=T,sep=',',ncol=Nfeed)
  write(t(ObsMixBtoFO),FileName,append=T,ncol=Nbreed,sep=",")
  
  # Now pull out results (feeding to breeding)
  ObsMixFtoBE <- array(0,dim=c(2,Nbreed,Nfeed))
  ObsMixFtoBP <- array(0,dim=c(2,Nbreed,Nfeed))
  ObsMixFtoBO <- matrix(0,nrow=2,ncol=Nfeed)
  write("Feeding to Breeding",FileName,append=T)
  if (FullDiag==T) print("Feeding to breeding mixing")
  for (Itype in 1:2)
  {
    write(Type[Itype],FileName,append=T)
    Index <- MixFile$Method == Type[Itype] & MixFile$Direction == "FeedingtoBreeding"
    MixFile2 <- MixFile[Index,]
    for (Ifeed in 1:length(FeedNames))
    {
      Total <- 0; Top <- 0; Bot <- 0
      for (Ibreed in 1:length(BreedNames))
      {
        Index <- MixFile2$Feeding %in% FeedNames[Ifeed] & MixFile2$Breeding %in% BreedNames[Ibreed]
        Est <- as.numeric(MixFile2$Estimate[Index])
        CV <- as.numeric(MixFile2$CV[Index])
        SD <- CV*Est
        if (is.na(CV)) 
        {
          Lower95 <- as.numeric(MixFile2$lowerCI[Index])
          Upper95 <- as.numeric(MixFile2$upperCI[Index])
          CV <- (Upper95-Lower95)/(2.0*1.96)/Est
          if (Est<=0) CV <- 0
          SD <- CV*Est
        }
        Top <- Top + Est*(1-Est)
        Bot <- Bot + SD^2
        if (CV < MinMixCV) CV <- MinMixCV
        SD <- Est*CV
        if (Est > 0 & SD < MinMixSD) SD <- MinMixSD
        if (Est > 0 || SD > 0) { NmixData <- NmixData + 1; Nprop[Itype+2] <- Nprop[Itype+2] + 1}
        Total <- Total+ Est
        
        #cat(Ibreed,Ifeed,Est,CV,SD,"\n")
        ObsMixFtoBE[Itype,Ibreed,Ifeed] <- Est
        ObsMixFtoBP[Itype,Ibreed,Ifeed] <- SD
      } # Breed
      ObsMixFtoBO[Itype,Ifeed] <- min(MaxN,Top/Bot*MixWeights[Itype])
      #print(Total)
    } # Feed 
    write("Estimates",FileName,append=T)
    write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixFtoBE[Itype,Ibreed,]))
      xx <- c(xx,",#",unlist(BreedNames[Ibreed]))
      write(xx,FileName,append=T,ncol=Nfeed+2,sep=",")
    }   
    write("SD",FileName,append=T)
    write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixFtoBP[Itype,Ibreed,]))
      xx <- c(xx,",#",unlist(BreedNames[Ibreed]))
      write(xx,FileName,append=T,ncol=Nfeed+2,sep=",")
    }   
  }  # Type
  write("EffN",FileName,append=T)
  write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
  write(t(ObsMixFtoBO),FileName,append=T,ncol=Nfeed,sep=",")
  #print(NmixData)
  #print(Nprop)
  
  Outs <- NULL
  Outs$ObsMixBtoFE <- ObsMixBtoFE
  Outs$ObsMixBtoFP <- ObsMixBtoFP
  Outs$ObsMixBtoFO <- ObsMixBtoFO
  Outs$ObsMixFtoBE <- ObsMixFtoBE
  Outs$ObsMixFtoBP <- ObsMixFtoBP
  Outs$ObsMixFtoBO <- ObsMixFtoBO
  Outs$NmixData <- NmixData
  
  return(Outs)
}

# ===========================================================================================
ReadGLORYS <- function(DatFile,FeedingOpt,YrStart=2000,YrEnd=2023){
  envdat <- read_rds(here('data','processed covariates','glorys_anomalies_by_zone.rds'))
  # pull the right timeseries for feeding grounds
  which_feed <- which(DatFile[,1]=="Feeding_grounds" & DatFile[,2]==FeedingOpt)
  Nfeed <- as.numeric(DatFile[which_feed,3])
  FeedNames <- as.character(DatFile[which_feed+1,1:Nfeed])
  sst <- envdat$sst |> 
    mutate(zone=str_remove_all(zoneID,"_\\d*"))
  sstL <- map(FeedNames,\(x){
    ssti <- sst |> filter(zone==x,year>=YrStart,year<=YrEnd) |> 
      pull(anom)
    ssti
  })
  sstMat <- do.call(rbind,sstL)
  
  chl <- envdat$chl |> 
    mutate(zone=str_remove_all(zoneID,"_\\d*"))
  chlL <- map(FeedNames,\(x){
    chli <- chl |> filter(zone==x,year>=YrStart,year<=YrEnd) |> 
      pull(anom)
    chli
  })
  chlMat <- do.call(rbind,chlL)
  # fill zeroes
  Outs <- list(
    sst=sstMat,
    chl=chlMat
  )
  Outs
}

# ===========================================================================================
# Function to make a data list for a specific user-designated scenario
MakeDataScenario <- function(Code,SensCase,StochSopt=1,StrayBase=0,DataFileName,Yr1=1970,Yr2=2023,
                             YrSDevs=1993, CatchSer="B",envOpt="none",
                             ByCatchFile="BycatchActual_2024_04_24.csv",AddCV=T,MixWeights=c(1,1),
                             MaxN=100,SigmaDevS=6,SigmaDevF=0.01,WithMirror=1,
                             SF=c(0,1,0,1,1,1),IAmat=8,SA=0.96,SC=0.8,
                             TimeLag=0,DensDepOpt=0,WghtTotal=1){
  
  # Specific model variants
  BreedingOpt <- substr(Code,1,2)
  FeedingOpt <- substr(Code,3,4)
  
  # ==============================================================================================================================
  # Read in the data file
  DatFile <- read.table(DataFileName,fill=T,col.names=c(1:200),comment.char="?")
  
  # Set up breeding and feeding grounds
  which_breed <- which(DatFile[,1]=="Breeding_grounds" & DatFile[,2]==BreedingOpt)
  Nbreed <- as.numeric(DatFile[which_breed,3])
  BreedNames <- DatFile[which_breed+1,1:Nbreed]
  which_feed <- which(DatFile[,1]=="Feeding_grounds" & DatFile[,2]==FeedingOpt)
  Nfeed <- as.numeric(DatFile[which_feed,3])
  FeedNames <- DatFile[which_feed+1,1:Nfeed]
  
  # Find the correct mixing data in the data file
  mix_start <- which(DatFile[,1]=="Mixing_matrix" & DatFile[,2]==Code)+1
  # make the mixing data matrix
  MixI <- map(1:Nbreed,\(x) as.numeric(DatFile[mix_start+x,1:Nfeed]))
  MixI <- do.call(rbind,MixI)
  rownames(MixI) <- BreedNames; colnames(MixI) <- FeedNames
  # Number of mixing parameters
  NmixPar <- sum(MixI>0)
  
  Years <- Yr1:(Yr2+1); Nyear <- length(Years)
  which_year_feedbreed <- which(DatFile[,2]=="Year_for_feeding_to_breeding")
  YearFeedBreed <- as.numeric(DatFile[which_year_feedbreed+1,1])-Yr1+1 # note: where does this get used??
  which_dirichlet <- which(DatFile[,2]=="Dirichlet_(1)_or_normal_(0)")
  Idirichlet <- as.numeric(DatFile[which_dirichlet+1,1])
  
  # ==============================================================================================================================
  # Read the catch data
  OutsCatch <- ReadCatches(DatFile,BreedingOpt,FeedingOpt,Nbreed,Nfeed,BreedNames,FeedNames,Yr1,Yr2,CatchSer,ByCatchFile)
  CatchF <- OutsCatch$CatchF 
  CatchB <- OutsCatch$CatchB
  
  # Read the survey data (note adjustments to indices to reflect C++ to R)
  OutsSurv <- ReadSurveyData(Code,BreedingOpt,FeedingOpt,BreedNames,FeedNames,Yr1,Yr2,SensTest=SensCase)
  SurveyI <- OutsSurv$SurveyI
  SurveyI[,1] <- SurveyI[,1] + 1
  SurveyI[,2] <- SurveyI[,2] + 1
  SurveyI[,5] <- SurveyI[,5] + 1
  which_surveys <- which(SurveyI[,4]!=-1)
  SurveyI[which_surveys,4] <- SurveyI[which_surveys,4] + 1
  SurveyR <- OutsSurv$SurveyR
  NsurveyData <- OutsSurv$NsurveyData
  NsurveySeries <- OutsSurv$NsurveySeries
  SurveySeries <- OutsSurv$SurveySeries
  NextraCV1 <- OutsSurv$NextraCV1
  
  # Read the mixing data
  OutsMix <- ReadMixingData(Code,DatFile,Yr1,Yr2,Nbreed,Nfeed,BreedNames,FeedNames,MixWeights,SensTest=SensCase,MaxN=MaxN)
  ObsMixBtoFE <- OutsMix$ObsMixBtoFE
  ObsMixBtoFP <- OutsMix$ObsMixBtoFP
  ObsMixBtoFO <- OutsMix$ObsMixBtoFO
  ObsMixFtoBE <- OutsMix$ObsMixFtoBE
  ObsMixFtoBP <- OutsMix$ObsMixFtoBP
  ObsMixFtoBO <- OutsMix$ObsMixFtoBO
  NmixData <- OutsMix$NmixData
  
  
  # ============================================================================================================================== 
  # Specify devs for survival and fecundity
  # Make devs for survival and fecundity
  # The form and inclusion of these will differ based on option "envOpt"
  OutsEnv <- ReadGLORYS(DatFile,FeedingOpt,YrStart=YrSDevs,YrEnd = 2023) # list with sst and chl anomalies by year (YrSDevs to Yr2)
  sst <- OutsEnv$sst
  chl <- OutsEnv$chl
  omega_sst <- 0 # initial coefficient for sst
  omega_chl <- 0 # initial coefficient for chl
  
  # Set up the right survival matrices for the chosen scenario
  if (BreedingOpt=="B1")
  {
    SBdevEst <- c(0,1,1,0)
    SBdevMat <- matrix(c(0,0,YrSDevs,Yr2,YrSDevs,Yr2,0,0),ncol=2,byrow=T)-Yr1
    FBdevEst <- c(0,0,0,0)
    FBdevMat <- matrix(c(0,0,YrSDevs,Yr2,YrSDevs,Yr2,0,0),ncol=2,byrow=T)-Yr1
  }
  if (BreedingOpt=="B2")
  {
    SBdevEst <- c(0,1,0,1,0)
    SBdevMat <- matrix(c(0,0,YrSDevs,Yr2,0,0,YrSDevs,Yr2,0,0),ncol=2,byrow=T)-Yr1
    FBdevEst <- c(0,0,0,0,0)
    FBdevMat <- matrix(c(0,0,YrSDevs,Yr2,0,0,YrSDevs,Yr2,0,0),ncol=2,byrow=T)-Yr1
  }
  if (FeedingOpt=="F1" || FeedingOpt=="F2")
  {
    Yrs <- c(YrSDevs,Yr2)
    SFdevEst <- NULL
    SFdevMat <- matrix(0,nrow=Nfeed,ncol=2,byrow=T)
    for (Ifeed in 1:Nfeed)
    {
      SFdevEst <- c(SFdevEst,SF[Ifeed])
      if (SF[Ifeed]==1) SFdevMat[Ifeed,] <- Yrs-Yr1
    }
    FFdevEst <- c(0,0,0,0,0,0)
    FFdevMat <- matrix(c(0,0,rep(c(YrSDevs,Yr2),5)),ncol=2,byrow=T)-Yr1
  }
  # How many deviates to estimate (survival and fecundity, breeding and feeding grounds)
  nBsdevs <- sum(SBdevEst)*length(YrSDevs:Yr2) #survival, breeding grounds
  nBfdevs <- sum(FBdevEst)*length(YrSDevs:Yr2) #fecundity, breeding grounds
  nFsdevs <- sum(SFdevEst)*length(YrSDevs:Yr2) #survival, feeding grounds
  nFfdevs <- sum(FFdevEst)*length(YrSDevs:Yr2) #fecundity, feeding grounds
  
  # Carrying capacity deviates
  nKdevs <- Nfeed*length(YrSDevs:Yr2)
  
  if (nBfdevs==0) nBfdevs <- 1
  if (nFfdevs==0) nFfdevs <- 1
  SBdev <- rep(0,nBsdevs)
  FBdev <- rep(0,nBfdevs)
  SFdev <- rep(0,nFsdevs)
  FFdev <- rep(0,nFfdevs)
  Kdev <- rep(0.001,nKdevs)
  
  # Mirror devs between feeding grounds? (currently not used for environmentally-driven survival model)
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
  # Set up additional variances
  NextraCV <- ifelse(NextraCV1==0,1,NextraCV1)
  AddV=rep(0,NextraCV)
  
  # Set up mixing parameters
  # MixPars <- rep(-1,NmixPar)
  MixPars <- numeric()
  
  for (Ibreed in 1:Nbreed)
  {
    Iref <- which(MixI[Ibreed,]==-1)
    Temp <- rep(0,Nfeed)
    for (IdataS in 1:2) 
      for (Ifeed in 1:Nfeed) Temp[Ifeed] <- Temp[Ifeed] + (0.001+ObsMixBtoFE[IdataS,Ibreed,Ifeed])/2.0
    for (Ifeed in 1:Nfeed) if (MixI[Ibreed,Ifeed]!=0 && Ifeed!=Iref) MixPars <- c(MixPars,log(Temp[Ifeed]/Temp[Iref]))
  }
  
  # ==============================================================================================================================
  ## List of parameters for TMB to estimate
  parameters <- list(rval=0.07, # Maximum growth rate
                     logK=rep(log(20000),Nbreed), # Breeding ground K
                     logBK=rep(3,Nbreed), # Relative depletion of each breeding stock
                     InfluxP=10, # Transfer influx (?)
                     inert_par=0, # related to fecundity calculation
                     MixPars=MixPars, # parameters of the mixing matrix
                     AddV=AddV, # additional variance parameters for three time-series
                     omega_sst = omega_sst, # coefficient on the sst timeseries affecting survival
                     omega_chl = omega_chl, # coefficent on the chl timeseries affecting survival
                     Sigma_SBdev=SigmaDevS, SBdev=SBdev, # survival devs, breeding grounds
                     Sigma_FBdev=SigmaDevF, FBdev=FBdev, # fecundity devs, breeding grounds
                     Sigma_SFdev=SigmaDevS,
                     SFdev=SFdev, # survival devs, feeding grounds
                     Kdev = Kdev, # carrying capacity devs, feeding grounds
                     Sigma_FFdev=SigmaDevF, FFdev=FFdev) # fecundity devs, feeding grounds
  
  # If driving survival without environment, remove omegas from parameter list
  if(envOpt=="none"){
    parameters$omega_sst <- parameters$omega_chl <- parameters$Kdev <- NULL
  }
  
  # If driving survival directly with environment, remove SFdevs from parameter list
  if(envOpt=="direct"){
    parameters$SFdev <- NULL
    parameters$Sigma_SFdev <- NULL
    parameters$Kdev <- NULL
  }
  if(envOpt=="index"){
    parameters$Sigma_SFdev <- NULL
    parameters$Kdev <- NULL
  }
  
  
  # Fixing parameters
  # factor(NA) means DON'T ESTIMATE THIS, USE FIXED PARAM
  # These will go into estimation at their INITIAL VALUES and will not change
  mapv1=list(#rval=factor(NA),#logBK=rep(factor(NA),Nbreed),
    InfluxP=factor(NA),
    inert_par=factor(NA),
    #Sigma_SBdev=factor(NA),
    Sigma_SBdev=factor(NA),SBdev=rep(factor(NA),length(SBdev)),
    Sigma_FBdev=factor(NA),FBdev=rep(factor(NA),length(FBdev)),
    Sigma_SFdev=factor(NA),
    #Sigma_SFdev=factor(NA),SFdev=rep(factor(NA),length(SFdev)),
    Sigma_FFdev=factor(NA),FFdev=rep(factor(NA),length(FFdev))
  )
  mapv2=list(#rval=factor(NA),#logBK=rep(factor(NA),Nbreed),
    InfluxP=factor(NA),
    inert_par=factor(NA),
    Sigma_SBdev=factor(NA),
    Sigma_SBdev=factor(NA),SBdev=rep(factor(NA),length(SBdev)),
    Sigma_FBdev=factor(NA),FBdev=rep(factor(NA),length(FBdev)),
    Sigma_SFdev=factor(NA),SFdev=rep(factor(NA),length(SFdev)),
    Sigma_FFdev=factor(NA),FFdev=rep(factor(NA),length(FFdev))
  )
  mapv3=list(#rval=factor(NA),#logBK=rep(factor(NA),Nbreed),
    InfluxP=factor(NA),
    inert_par=factor(NA),
    Sigma_SBdev=factor(NA),SBdev=rep(factor(NA),length(SBdev)),
    Sigma_FBdev=factor(NA),FBdev=rep(factor(NA),length(FBdev)),
    #Sigma_SFdev=factor(NA),SFdev=rep(factor(NA),length(SFdev)),
    Sigma_FFdev=factor(NA),FFdev=rep(factor(NA),length(FFdev))
  )
  if (StochSopt==1) map <- mapv1
  if (StochSopt==0 | envOpt=="varK") map <- mapv2
  if (envOpt=="direct"|envOpt=="index") map <- mapv3
  if (AddCV==F) map$AddV <- rep(factor(NA),length(AddV))
  
  # ==============================================================================================================================
  # Return a big list of data for the model
  datout <- list(Nbreed=Nbreed, 
                 Nfeed=Nfeed,
                 BreedNames=BreedNames,
                 FeedNames=FeedNames,
                 Yr1=Yr1,
                 Yr2=Yr2,
                 YrSDevs=YrSDevs,
                 Years=Years,
                 IAmat=IAmat,
                 SA=SA,
                 SC=SC,
                 TimeLag=TimeLag,
                 SurveyI=SurveyI,
                 SurveyR=SurveyR,
                 NextraCV1 = NextraCV1,
                 CatchB=CatchB,
                 CatchF=CatchF,
                 sst=sst,
                 chl=chl,
                 omega_sst = omega_sst,
                 omega_chl = omega_chl,
                 SF=SF,
                 SBdevEst=SBdevEst,
                 SBdevMat=SBdevMat,
                 FBdevEst=FBdevEst,
                 FBdevMat=FBdevMat,
                 SFdevEst=SFdevEst,
                 SFdevMat=SFdevMat,
                 FFdevEst=FFdevEst,
                 FFdevMat=FFdevMat,
                 MixI=MixI,
                 NsurveyData=NsurveyData,
                 NsurveySeries=NsurveySeries,
                 SurveySeries=SurveySeries,
                 NmixData=NmixData,
                 NmixPar=NmixPar,
                 YearFeedBreed=YearFeedBreed,
                 Idirichlet=Idirichlet,
                 ObsMixBtoFE=ObsMixBtoFE,
                 ObsMixBtoFP=ObsMixBtoFP,
                 ObsMixBtoFO=ObsMixBtoFO,
                 ObsMixFtoBE=ObsMixFtoBE,
                 ObsMixFtoBP=ObsMixFtoBP,
                 ObsMixFtoBO=ObsMixFtoBO,
                 StochSopt=StochSopt,
                 StrayBase=StrayBase,
                 DensDepOpt=DensDepOpt,
                 Nmirror=Nmirror,
                 Mirror=Mirror,
                 WghtTotal=WghtTotal,
                 map=map,
                 parameters=parameters)
  return(datout)
}
