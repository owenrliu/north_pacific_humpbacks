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
  Outs$CAtchB <- CatchB
  return(Outs)
  
 }

# ===========================================================================================

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
  
  # Deal with incorrect names0
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
      AAAA
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
