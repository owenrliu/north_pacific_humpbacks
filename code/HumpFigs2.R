# ------------------------------------------------------------------------------------

PlotPropns <- function(Nprop,ObsMixPropI,ObsMixProp,PredMix,PDF=T,Code="No model",BreedNames,FeedNames,Nyear)
{
  FileName <- paste0(PlotsDir,"Propn",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  
  # Image size
  par(mfrow=c(2,1),mar=c(10,4,2,2)+0.1,oma=c(1,1,3,1))
  
  # Extract the breeding to feeding ground proportions
  VecOut <- matrix(0,ncol=9,nrow=100)
  Ipnt <- 0
  cols <- NULL
  for (Ibreed in 1:4)
    for (Ifeed in 1:6) 
    {
      IndexA <- which(ObsMixPropI[,1]==1 & ObsMixPropI[,3]==Ibreed & ObsMixPropI[,4]==Ifeed) 
      if (length(IndexA) > 0)
      {
        Pred <- PredMix[IndexA[1]]
        vec <- c(Ibreed,Ifeed,Pred)
        Ipnt <- Ipnt + 1
        for(IdataS in 1:2)
        {
          Jpnt <- IndexA[which(ObsMixPropI[IndexA,2]==IdataS)]
          if (IdataS %% 2 ==1) cols <- c(cols,"green")
          if (IdataS %% 2 ==1) cols <- c(cols,"red")
          ObsA <- c(-1,-1,-1)
          if (length(Jpnt) > 0)
            ObsA <- c(ObsMixProp[Jpnt,1],ObsMixProp[Jpnt,1]-1.96*ObsMixProp[Jpnt,2],ObsMixProp[Jpnt,1]+1.96*ObsMixProp[Jpnt,2])
          vec <- c(vec,ObsA)
        }
        #print(vec)
        VecOut[Ipnt,] <- vec
      }
    }
  VecOut <- VecOut[1:Ipnt,]
  NpropV <- length(VecOut[,1])
  #print(VecOut)
  #print(cols)
  
  # Now print them
  plot(0,0,xlab="",ylab="Breeding -> Feeding proportions",xlim=c(0.5,NpropV+0.5),ylim=c(0,1.05),axes=F,type="n")
  points(c(1:NpropV),VecOut[,3],pch=17,col=cols)
  box();  axis(2)
  XLAB <- NULL
  for (II in 1:NpropV)
    XLAB <- c(XLAB,paste(BreedNames[VecOut[II,1]],"-",FeedNames[VecOut[II,2]]))  
  for (II in 1:NpropV) 
  { points(II-0.2,VecOut[II,4],pch=1,col=cols[II]); lines(c(II-0.2,II-0.2),c(VecOut[II,5],VecOut[II,6]),col=cols[II]); }
  for (II in 1:NpropV) 
  { points(II+0.2,VecOut[II,7],pch=16,col=cols[II]); lines(c(II+0.2,II+0.2),c(VecOut[II,8],VecOut[II,9]),col=cols[II]); }
  axis(1,1:NpropV,XLAB,las=2,cex.axis=1) 
  abline(h=0,col="black")
  
  # Extract the feeding to breeding grounds proprtions
  VecOut <- matrix(0,ncol=9,nrow=100)
  Ipnt <- 0
  cols <- NULL
  for (Ifeed in 1:6) 
    for (Ibreed in 1:4)
    {
      IndexA <- which(ObsMixPropI[,1]==2 & ObsMixPropI[,3]==Ibreed & ObsMixPropI[,4]==Ifeed) 
      if (length(IndexA) > 0)
      {
        Pred <- PredMix[IndexA[1]]
        vec <- c(Ibreed,Ifeed,Pred)
        Ipnt <- Ipnt + 1
        for(IdataS in 1:2)
        {
          Jpnt <- IndexA[which(ObsMixPropI[IndexA,2]==IdataS)]
          if (IdataS %% 2 ==1) cols <- c(cols,"green")
          if (IdataS %% 2 ==1) cols <- c(cols,"red")
          ObsA <- c(-1,-1,-1)
          if (length(Jpnt) > 0)
            ObsA <- c(ObsMixProp[Jpnt,1],ObsMixProp[Jpnt,1]-1.96*ObsMixProp[Jpnt,2],ObsMixProp[Jpnt,1]+1.96*ObsMixProp[Jpnt,2])
          vec <- c(vec,ObsA)
        }
        #print(vec)
        VecOut[Ipnt,] <- vec
      }
    }
  VecOut <- VecOut[1:Ipnt,]
  NpropV <- length(VecOut[,1])

  # Now plot the graph
  plot(0,0,xlab="",ylab="Feeding -> Breeding proportions",xlim=c(0.5,NpropV+0.5),ylim=c(0,1.05),axes=F,type="n")
  points(c(1:NpropV),VecOut[,3],pch=17,col=cols)
  box();  axis(2)
  XLAB <- NULL
  for (II in 1:NpropV)
    XLAB <- c(XLAB,paste(BreedNames[VecOut[II,1]],"-",FeedNames[VecOut[II,2]]))  
  for (II in 1:NpropV) 
  { points(II-0.2,VecOut[II,4],pch=1,col=cols[II]); lines(c(II-0.2,II-0.2),c(VecOut[II,5],VecOut[II,6]),col=cols[II]); }
  for (II in 1:NpropV) 
  { points(II+0.2,VecOut[II,7],pch=16,col=cols[II]); lines(c(II+0.2,II+0.2),c(VecOut[II,8],VecOut[II,9]),col=cols[II]); }
  axis(1,1:NpropV,XLAB,las=2,cex.axis=1) 
  abline(h=0,col="red")
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
}

# ------------------------------------------------------------------------------------

PlotIndex <- function(SurveyD,LogNT,LogNb,LogNf,Qest,PDF=T,Code="No model",BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed,rept,AllPlots)
 {
  
  Lag <- 45
  LastYr <- Nyear-1
  
  PlotYrs <- (Nyear-Lag):LastYr
  Tyears <- Years[PlotYrs]

  FileName <- paste0(PlotsDir,"TotalAbund",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(1,1),oma=c(3,3,20,3),mar=c(4,4,2,2))
  
  #Total population absolute abundance data
  Index2 <- SurveyD[SurveyD[,3]== 1,] 
  SD1 <- Index2[,9];
  SDD <- SD1
  if (Index2[1,7]>0) SDD <- sqrt(SD1^2+rept$AddV[Index2[1,7]])
  Low1 <- Index2[,8]*exp(-1.96*SD1)
  Upp1 <- Index2[,8]*exp(1.96*SD1)
  Low2 <- Index2[,8]*exp(-1.96*SDD)
  Upp2 <- Index2[,8]*exp(1.96*SDD)
  ymax2 <- max(exp(LogNT[,1]),max(Upp2))*1.1
  plot(Tyears,exp(LogNT[PlotYrs,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
  for (II in 1:length(Index2[,1]))
  {
    Yr <- (Index2[II,1]+Index2[II,2])/2+Yr1
    #cat(Yr,Low[II],Upp[II],"\n")
    points(Yr,Index2[II,8],pch=16)  
    if (Index2[1,7]==0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=1)
    if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=2,col="red",lty-1)
    if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low2[II],Upp2[II]),lty=1,lwd=1)
  }  
  title("All")
  lines(Tyears,exp(LogNT[PlotYrs,1]),lwd=2,lty=1)
  lines(Tyears,exp(LogNT[PlotYrs,1]-1.96*LogNT[PlotYrs,2]),lty=2)
  lines(Tyears,exp(LogNT[PlotYrs,1]+1.96*LogNT[PlotYrs,2]),lty=2)
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
  #if (AllPlots==F) return()

  #Breeding area absolute abundance data
  FileName <- paste0(PlotsDir,"BreedAbund",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(3,2),oma=c(1,1,3,1),mar=c(4,4,2,2))
  if (FullDiag==T) print("Doing absolute abundance data")
  for (Ibreed in 1:Nbreed)
  {
    Index2 <- SurveyD[SurveyD[,3]== 2 & SurveyD[,4]== Ibreed & SurveyD[,5]== 1,]
    if (length(Index2)==9) Index2 <- t(Index2)
    Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
    Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
    Ipnt <- rev(rev(Ipnt[-c(Nyear-LastYr)]))

    SD1 <- Index2[,9];
    SDD <- SD1
    if (Index2[1,7]>0) SDD <- sqrt(SD1^2+rept$AddV[Index2[1,7]])
    Low1 <- Index2[,8]*exp(-1.96*SD1)
    Upp1 <- Index2[,8]*exp(1.96*SD1)
    Low2 <- Index2[,8]*exp(-1.96*SDD)
    Upp2 <- Index2[,8]*exp(1.96*SDD)
    ymax2 <- max(exp(LogNb[Ipnt,1]),max(Upp2))*1.1
    plot(Tyears,exp(LogNb[Ipnt,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
    for (II in 1:length(Index2[,1]))
    {
      Yr <- (Index2[II,1]+Index2[II,2])/2+Yr1
      #cat(Yr,Low[II],Upp[II],"\n")
      points(Yr,Index2[II,8],pch=16)  
      if (Index2[1,7]==0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=1)
      if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=2,col="red")
      if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low2[II],Upp2[II]),lty=1,lwd=1)
    }  
    title(BreedNames[Ibreed])
    lines(Tyears,exp(LogNb[Ipnt,1]),lwd=2,lty=1)
    lines(Tyears,exp(LogNb[Ipnt,1]-1.96*LogNb[Ipnt,2]),lty=2)
    lines(Tyears,exp(LogNb[Ipnt,1]+1.96*LogNb[Ipnt,2]),lty=2)
  } 
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()

  #Feeding area absolute abundance data
  FileName <- paste0(PlotsDir,"FeedAbund",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(3,2),oma=c(1,1,3,1),mar=c(4,4,2,2))
  for (Ifeed in 1:Nfeed)
  {
    Index2 <- SurveyD[SurveyD[,3]== 3 & SurveyD[,4]== Ifeed & SurveyD[,5]== 1,]
    if (length(Index2)==9) Index2 <- t(Index2)
    Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
    Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
    Ipnt <- rev(rev(Ipnt[-c(Nyear-LastYr)]))
    SD1 <- Index2[,9]
    SDD <- SD1
    if (Index2[1,7]>0) SDD <- sqrt(SD1^2+rept$AddV[Index2[1,7]])
    Low1 <- Index2[,8]*exp(-1.96*SD1)
    Upp1 <- Index2[,8]*exp(1.96*SD1)
    Low2 <- Index2[,8]*exp(-1.96*SDD)
    Upp2 <- Index2[,8]*exp(1.96*SDD)
    ymax2 <- max(exp(LogNf[Ipnt,1]),max(Upp2))*1.1
    plot(Tyears,exp(LogNf[Ipnt,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
    for (II in 1:length(Index2[,1]))
    {
      Yr <- (Index2[II,1]+Index2[II,2])/2+Yr1
      #cat(Yr,Low[II],Upp[II],"\n")
      points(Yr,Index2[II,8],pch=16)  
      if (Index2[1,7]==0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=1)
      if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=2,col="red")
      if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low2[II],Upp2[II]),lty=1,lwd=1)
    }  
    title(FeedNames[Ifeed])
    lines(Tyears,exp(LogNf[Ipnt,1]),lwd=2,lty=1)
    lines(Tyears,exp(LogNf[Ipnt,1]-1.96*LogNf[Ipnt,2]),lty=2)
    lines(Tyears,exp(LogNf[Ipnt,1]+1.96*LogNf[Ipnt,2]),lty=2)
  } 
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()

  ## Relative indices
  FileName <- paste0(PlotsDir,"RelAbund",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(3,3),oma=c(1,1,3,1),mar=c(4,4,2,2))
  if (FullDiag==T) print("Doing relative abundance data")
  for (Iseries in 2:max(SurveyD[,5]))
  {
    Index2 <- SurveyD[SurveyD[,5]== Iseries,]
    if (length(Index2)==9) Index2 <- t(Index2)
    if (Index2[1,3]==2)
     {
      Ibreed <- Index2[1,4]
      Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
      Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
      Ipnt <- rev(rev(Ipnt[-c(Nyear-LastYr)]))
      SD1 <- Index2[,9];
      SDD <- SD1
      if (Index2[1,7]>0) SDD <- sqrt(SD1^2+rept$AddV[Index2[1,7]])
      Low1 <- Index2[,8]*exp(-1.96*SD1)
      Upp1 <- Index2[,8]*exp(1.96*SD1)
      Low2 <- Index2[,8]*exp(-1.96*SDD)
      Upp2 <- Index2[,8]*exp(1.96*SDD)
      ymax2 <- max(Qest[Iseries]*exp(LogNb[Ipnt,1]),max(Upp2))*1.1
      plot(Tyears,Qest[Iseries]*exp(LogNb[Ipnt,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
      for (II in 1:length(Index2[,1]))
       {
        Yr <- (Index2[II,1]+Index2[II,2])/2+Yr1
        #cat(Yr,Low[II],Upp[II],"\n")
        points(Yr,Index2[II,8],pch=16)  
        if (Index2[1,7]==0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=1)
        if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=2,col="red")
        if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low2[II],Upp2[II]),lty=1,lwd=1)
       }  
      lines(Tyears,Qest[Iseries]*exp(LogNb[Ipnt,1]),lwd=2,lty=1)
      Ibreed <- Index2[1,4]
      title(BreedNames[Ibreed])
      lines(Tyears,Qest[Iseries]*exp(LogNb[Ipnt,1]),lwd=2,lty=1)
      lines(Tyears,Qest[Iseries]*exp(LogNb[Ipnt,1]-1.96*LogNb[Ipnt,2]),lty=2)
      lines(Tyears,Qest[Iseries]*exp(LogNb[Ipnt,1]+1.96*LogNb[Ipnt,2]),lty=2)
    }
    if (Index2[1,3]==3)
     {
      Ifeed <- Index2[1,4]
      Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
      Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
      Ipnt <- rev(rev(Ipnt[-c(Nyear-LastYr)]))
      SD1 <- Index2[,9];
      SDD <- SD1
      if (Index2[1,7]>0) SDD <- sqrt(SD1^2+rept$AddV[Index2[1,7]])
      Low1 <- Index2[,8]*exp(-1.96*SD1)
      Upp1 <- Index2[,8]*exp(1.96*SD1)
      Low2 <- Index2[,8]*exp(-1.96*SDD)
      Upp2 <- Index2[,8]*exp(1.96*SDD)
      ymax2 <- max(Qest[Iseries]*exp(LogNf[Ipnt,1]),max(Upp2))*1.1
      plot(Tyears,Qest[Iseries]*exp(LogNf[Ipnt,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
      for (II in 1:length(Index2[,1]))
      {
        Yr <- (Index2[II,1]+Index2[II,2])/2+Yr1
        #cat(Yr,Low[II],Upp[II],"\n")
        points(Yr,Index2[II,8],pch=16)  
        if (Index2[1,7]==0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=1)
        if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low1[II],Upp1[II]),lty=1,lwd=2,col="red")
        if (Index2[1,7]>0) lines(c(Yr,Yr),c(Low2[II],Upp2[II]),lty=1,lwd=1)
      }  
      Ibreed <- Index2[1,4]
      title(FeedNames[Ibreed])
      lines(Tyears,Qest[Iseries]*exp(LogNf[Ipnt,1]),lwd=2,lty=1)
      lines(Tyears,Qest[Iseries]*exp(LogNf[Ipnt,1]-1.96*LogNf[Ipnt,2]),lty=2)
      lines(Tyears,Qest[Iseries]*exp(LogNf[Ipnt,1]+1.96*LogNf[Ipnt,2]),lty=2)
    }
  
  }
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
  
}
  
  # ------------------------------------------------------------------------------------
  
ENvfun <- function()
{
  
  # Index #1
  Years <- 2000:(Yr2-1)
  
  print("CaliforniaMOCI.csv")
  Env <- read.csv("CaliforniaMOCI.csv",skip=2,header=F)
  Index <- which(Env[,2] %in% Years)
  Env <- Env[Index,]
  pairs(Env[,4:6])
  Seas <- unique(Env[,3])
  par(mfrow=c(3,3))
  for (Ibreed in 2:3)
    for (II in 1:length(Seas))
      for (JJ in 4:6) 
      {
        Index <- which(Env[,3] %in% Seas[II])
        Env2 <- Env[Index,JJ]
        plot(Env2,Sdevs[Ibreed,31:50])
        xx <- lm(Sdevs[Ibreed,31:50]~Env2)
        print(summary(xx)$r.squared)
      }
  
  print("ERSST PDO")
  Env <- read.csv("ERSST PDO.csv",skip=1,header=F)
  Index <- which(Env[,1] %in% Years)
  Env <- Env[Index,14]
  par(mfrow=c(2,2))
  for (Ibreed in 2:3)
  {
    Env2 <- Env
    plot(Env2,Sdevs[Ibreed,31:50])
    xx <- lm(Sdevs[Ibreed,31:50]~Env2)
    print(summary(xx)$r.squared)
  }
  
  print("Oceanic nino")
  Env <- read.csv("Oceanic nino index.csv",skip=1,header=F)
  Index <- which(Env[,2] %in% Years)
  Env <- Env[Index,]
  Seas <- unique(Env[,1])
  par(mfrow=c(3,3))
  for (Ibreed in 2:3)
    for (II in 1:length(Seas))
    {
      Index <- which(Env[,1] %in% Seas[II])
      Env2 <- Env[Index,3]
      plot(Env2,Sdevs[Ibreed,31:50])
      xx <- lm(Sdevs[Ibreed,31:50]~Env2)
      print(summary(xx)$r.squared)
    }
 }

# ====================================================================================

PlotTraj <- function(LogNT,LogNb,LogNf,PDF=T,Code="No model",BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed)
{
  Lag <- 45
  LastYr <- Nyear-1
  PlotYrs <- (Nyear-Lag):LastYr
  Tyears <- Years[PlotYrs]
  #Breeding area absolute abundance data
  FileName <- paste0(PlotsDir,"BreedTotal",Code,".png")
  if (PDF==T) png(filename=FileName,width=600,height=800)
  par(mfrow=c(3,2),oma=c(1,1,3,1),mar=c(4,4,2,2))
  if (FullDiag==T) print("Doing absolute abundance")
  for (Ibreed in 1:Nbreed)
  {
    Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
    Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
    Ipnt <- rev(rev(Ipnt[-c(Nyear-LastYr)]))
    ymax2 <- max(exp(LogNb[Ipnt,1]))*1.4
    plot(Tyears,exp(LogNb[Ipnt,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
    title(BreedNames[Ibreed])
    lines(Tyears,exp(LogNb[Ipnt,1]-1.96*LogNb[Ipnt,2]),lty=2)
    lines(Tyears,exp(LogNb[Ipnt,1]+1.96*LogNb[Ipnt,2]),lty=2)
  } 
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
  
  #Feeding area absolute abundance 
  FileName <- paste0(PlotsDir,"FeedTotal",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(3,2),oma=c(1,1,3,1),mar=c(4,4,2,2))
  for (Ifeed in 1:Nfeed)
  {
    Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
    Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
    Ipnt <- rev(rev(Ipnt[-c(Nyear-LastYr)]))
    ymax2 <- max(exp(LogNf[Ipnt,1]))*1.4
    plot(Tyears,exp(LogNf[Ipnt,1]),xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
    title(FeedNames[Ifeed])
    lines(Tyears,exp(LogNf[Ipnt,1]-1.96*LogNf[Ipnt,2]),lty=2)
    lines(Tyears,exp(LogNf[Ipnt,1]+1.96*LogNf[Ipnt,2]),lty=2)
  } 
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
}

# =======================================================================================

PlotSurvB <- function(SurvOutB,baseS=0.96,PDF=T,Code="No model",BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed)
{
  Lag <- 35
  PlotYrs <- (Nyear-Lag):(Nyear-1)
  Tyears <- Years[PlotYrs]
  #Survival estimated by breeding area
  FileName <- paste0(PlotsDir,"SurvivalB",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(3,2),oma=c(1,1,3,1),mar=c(4,4,2,2))
  if (FullDiag==T) print("Doing survival")
  for (Ibreed in 1:Nbreed)
  {
    Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
    Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
    Ipnt <- Ipnt[-length(Ipnt)]
    ymax2 <- max(SurvOutB[Ipnt,1])*1.4
    plot(Tyears,SurvOutB[Ipnt,1],xlab="Year",ylab="Survival rate (Breeding stock)",ylim=c(0.4,ymax2),type="l",col="red",lwd=2) 
    abline(h=baseS,col="blue")
    abline(h=1,col="blue")
    title(BreedNames[Ibreed])
    lines(Tyears,SurvOutB[Ipnt,1]-1.96*SurvOutB[Ipnt,2],lty=2)
    lines(Tyears,SurvOutB[Ipnt,1]+1.96*SurvOutB[Ipnt,2],lty=2)
    abline(v=2014,lwd=3,col="black",lty=2)
  } 
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
}

# =======================================================================================

PlotSurvF <- function(SurvOutF,baseS=0.96,PDF=T,Code="No model",BreedNames,FeedNames,Nyear,Years,Yr1,Yr2,Nbreed,Nfeed)
{
  Lag <- 35
  PlotYrs <- (Nyear-Lag):(Nyear-1)
  Tyears <- Years[PlotYrs]
  #Survival estimated by feeding ground
  FileName <- paste0(PlotsDir,"SurvivalF",Code,".png")
  if (PDF==T) png(filename=FileName,width=700,height=800)
  par(mfrow=c(3,2),oma=c(1,1,3,1),mar=c(4,4,2,2))
  if (FullDiag==T) print("Doing survival")
  for (Ifeed in 1:Nfeed)
  {
    Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
    Ipnt <- Ipnt[which(Years>=Yr2+1-Lag)]
    Ipnt <- Ipnt[-length(Ipnt)]
    ymax2 <- max(SurvOutF[Ipnt,1])*1.4
    plot(Tyears,SurvOutF[Ipnt,1],xlab="Year",ylab="Survival rate (Feeding ground)",ylim=c(0.4,ymax2),type="l",col="red",lwd=2) 
    abline(h=baseS,col="blue")
    abline(h=1,col="blue")
    title(FeedNames[Ifeed])
    lines(Tyears,SurvOutF[Ipnt,1]-1.96*SurvOutF[Ipnt,2],lty=2)
    lines(Tyears,SurvOutF[Ipnt,1]+1.96*SurvOutF[Ipnt,2],lty=2)
    abline(v=2014,lwd=3,col="black",lty=2)
  } 
  mtext(paste0("Model run = ",Code),side=3,line=0,outer=T,cex=1.5)
  if (PDF==T) dev.off()
}

