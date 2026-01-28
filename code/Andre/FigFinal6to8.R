setwd("C:\\Research\\Iwc26\\NP_humpbacks\\")

ReadData <- function(Code,Code2,Yr1=1965,Yr2=2023,Nboot,type="a",PDF=F)
{
  # Specific model variants
  BreedingOpt <- substr(Code,1,2)
  FeedingOpt <- substr(Code,3,4)
  
  # -----------------------------------------------------------------------------------------------------------------------------------
  
  # Read in the structure data file
  DatFile <- read.table("Hump.dat",fill=T,col.names=c(1:200),comment.char="?")
  
  # Details of the breeding and feeding grounds
  Index <- which(DatFile[,1]=="Number_of_breeding_and_feeding_grounds:"); Nbreed <- as.numeric(DatFile[Index,2]); Nfeed <- as.numeric(DatFile[Index,3])
  Index <- which(DatFile[,1]=="Breeding_grounds" & DatFile[,2]==BreedingOpt); Nbreed <- as.numeric(DatFile[Index,3]); BreedNames <- DatFile[Index+1,1:Nbreed]
  Index <- which(DatFile[,1]=="Feeding_grounds" & DatFile[,2]==FeedingOpt); Nfeed <- as.numeric(DatFile[Index,3]); FeedNames <- DatFile[Index+1,1:Nfeed]
  
  BaseS <- 0.96
  Nyear <- Yr2-Yr1+1 
  Years <- 1980:2021
  Use <- Years-Yr1+1
  FileName <- paste0("Diags/",Code,Code2,".Boot")
  DatFile <- read.table(FileName,fill=T,col.names=c(1:200),comment.char="?")
  Index <- which(DatFile[,1]=="Total" & DatFile[,2]=="abundance:")
  print(length(Index))
  Nboot <- length(Index)-1
  cat(FileName,Nboot,"\n")
  
  SummaryT <- matrix(0,nrow=Nyear,ncol=1+Nboot)  
  SummarySur <- array(0,dim=c(Nfeed,Nyear,1+Nboot))
  IndexS <- which(DatFile[,1]=="Survival" & DatFile[,3]=="feeding")
  SDsT <- rep(0,Nyear)
  Index1 <- Index[1]
  for (Iyr in 1:Nyear) SDsT[Iyr] <- as.numeric(DatFile[Index1+Iyr,4]) 
  for (II in 0:Nboot)
   {
    Index1 <- Index[II+1]
    for (Iyr in 1:Nyear)
      SummaryT[Iyr,II+1] <- as.numeric(DatFile[Index1+Iyr,3]) 
    Index1 <- IndexS[II+1]
    for (Ifeed in 1:Nfeed)
      for (Iyr in 1:Nyear)
        SummarySur[Ifeed,Iyr,II+1] <- as.numeric(DatFile[Index1+Iyr,2+Ifeed])  
   }
  
  # Remove years not of interest and extract the SDs for abundance
  SummaryT2 <- SummaryT
  SummaryT <- SummaryT[Use,]
  SummarySur <- SummarySur[,Use,]
  SDsT <- SDsT[Use]
  ActT <- SummaryT[,1]
  ActSur <- SummarySur[,,1]
  SummaryT <- SummaryT[,-1]
  SummarySur <- SummarySur[,,-1]

  Outs <- NULL
  Outs$Nboot <- Nboot
  Outs$BaseS <- BaseS
  Outs$SummaryT <- SummaryT
  Outs$SummaryT2 <- SummaryT2
  Outs$SummarySur <- SummarySur
  Outs$ActT <- ActT
  Outs$ActSur <- ActSur
  Outs$SDsT <- SDsT
  Outs$Code <- Code2
  Outs$Code2 <- Code2
  Outs$Yr1 <- Yr1
  Outs$Yr2 <- Yr2
  Outs$Nbreed <- Nbreed
  Outs$Nfeed <- Nfeed
  Outs$FeedNames <- FeedNames
  Outs$Years <- Yr1:Yr2
  return(Outs)
}

##===================================================================================================

Fig1a <- function(Resu,Years=1980:2021,type="a")
{
  # Extract what we need
  Code <- Resu$Code
  Code2 <- Resu$Code2
  Nfeed <- Resu$Nfeed
  Nbreed <- Resu$Nbreed
  BaseS <- Resu$BaseS
  Yr1 <- Resu$Yr1
  Yr2 <- Resu$Yr2
  Nyear <- Yr2-Yr1+1 
  Nboot <- Resu$Nboot
  SummaryT <- Resu$SummaryT
  SummarySur <- Resu$SummarySur
  SDsT <- Resu$SDs
  ActT <- Resu$ActT
  ActSur <- Resu$ActSur

  Years <- 1980:2021
  Use <- Years-Yr1+1
  
  Qants <- matrix(0,nrow=length(SummaryT[,1]),ncol=5)
  for (II in 1:length(SummaryT[,1]))
   Qants[II,] <- quantile(SummaryT[II,],prob=c(0.05,0.25,0.5,0.75,0.95))     
  ymax <- max(Qants)*1.1

  # Plot for total 
  plot(Years,ActT,xlab="",ylab="",type="l",ylim=c(0,ymax),yaxs="i",lwd=2,col="Red")
  lines(Years,ActT*exp(1.645*SDsT),lty=2)
  lines(Years,ActT*exp(-1.645*SDsT),lty=2)
  mtext(type,side=3,line=1)
  plot(Years,ActT,xlab="",ylab="",type="l",ylim=c(0,ymax),yaxs="i")
  xx <- c(Years,rev(Years))
  yy <- c(Qants[,1],rev(Qants[,5]))
  polygon(xx,yy,col="gray50")
  yy <- c(Qants[,2],rev(Qants[,4]))
  polygon(xx,yy,col="gray95")
  lines(Years,Qants[,3],col="blue",lwd=2)
  lines(Years,ActT,col="red",lwd=2)
}

# ======================================================================================

Fig1b <- function(Resu,Years=1980:2021,type="a")
{
  # Extract what we need
  Code <- Resu$Code
  Code2 <- Resu$Code2
  Nfeed <- Resu$Nfeed
  Nbreed <- Resu$Nbreed
  BaseS <- Resu$BaseS
  Yr1 <- Resu$Yr1
  Yr2 <- Resu$Yr2
  Nyear <- Yr2-Yr1+1 
  Nboot <- Resu$Nboot
  SummaryT <- Resu$SummaryT
  SummarySur <- Resu$SummarySur
  SDsT <- Resu$SDs
  ActT <- Resu$ActT
  ActSur <- Resu$ActSur
  FeedNames <- Resu$FeedNames
  
  Qants2 <- array(0,dim=c(Nfeed,length(SummaryT[,1]),5))
  for (Ifeed in 1:Nfeed)
    for (II in 1:length(SummaryT[,1]))
      Qants2[Ifeed,II,] <- quantile(SummarySur[Ifeed,II,],prob=c(0.05,0.25,0.5,0.75,0.95))     
  ymax2 <- max(Qants2)*1.1
 
  # Plot for survival
  par(mfrow=c(3,2),oma=c(5,4,3,1),mar=c(4,4,2,2))
  for (Ifeed in 1:Nfeed)
   {
    plot(Years,ActSur[Ifeed,],xlab="Year",ylab="",ylim=c(0.4,ymax2),type="l",col="red",lwd=2) 
    abline(h=BaseS,col="blue")
    abline(h=1,col="blue")
    title(FeedNames[Ifeed])
    xx <- c(Years,rev(Years))
    yy <- c(Qants2[Ifeed,,1],rev(Qants2[Ifeed,,5]))
    polygon(xx,yy,col="gray50")
    yy <- c(Qants2[Ifeed,,2],rev(Qants2[Ifeed,,4]))
    polygon(xx,yy,col="gray95")
    lines(Years,Qants2[Ifeed,,3],col="blue",lwd=2)
    lines(Years,ActSur[Ifeed,],col="red",lwd=2)
    abline(v=2014,lwd=3,col="black",lty=2)
   } 
  mtext("Year",side=1,line=1,outer=T) 
  mtext("Survival rate (Feeding ground)",side=2,line=1,outer=T) 
  mtext(type,side=3,line=1,outer=T) 
}

# ======================================================================================

Fig1c <- function(Resus,Years=1980:2021,labs,type="a")
{
  Nmodels <- length(Resus)
  print(Nmodels)
  
  par(mfrow=c(3,1),oma=c(2,2,2,2),mar=c(5,4,2,2))
  AllRatio <- NULL
  for (Imodel in 1:Nmodels)
   {
    SummaryT <- Resus[[Imodel]]$SummaryT2
    Use <- which(Resus[[Imodel]]$Years %in% c(2002,2021))
    Vals <- SummaryT[Use,]
    Ratio <- 100*(Vals[2,]/Vals[1,]-1.0)
    hist(Ratio,main=labs[Imodel])
    AllRatio <- c(AllRatio,Ratio)
   }  
  hist(AllRatio,main=type)
  cat("Ratio",round(c(quantile(AllRatio,prob=0.05),mean(AllRatio),quantile(AllRatio,prob=c(0.5,0.9))),2),"\n")
  
 }

# ======================================================================================

Resu <- vector(mode="list",length=4)
Resu[[1]] <- ReadData("B1F1","BC",Yr1=1970,Yr2=2023,Nboot=500,PDF=T)
Resu[[2]] <- ReadData("B1F2","BC",Yr1=1970,Yr2=2022,Nboot=500,PDF=T)
Resu[[3]] <- ReadData("B2F1","BC",Yr1=1970,Yr2=2022,Nboot=500,PDF=T)
Resu[[4]] <- ReadData("B2F2","BC",Yr1=1970,Yr2=2022,Nboot=500,PDF=T)


png(filename="E:/FigD7a1.png",width=700,height=800)
Fig1b(Resu[[1]],type="B1-F1")
dev.off()
png(filename="E:/FigD7a2.png",width=700,height=800)
Fig1b(Resu[[2]],type="B1-F2")
dev.off()
png(filename="E:/FigD7a3.png",width=700,height=800)
Fig1b(Resu[[3]],type="B2-F1")
dev.off()
png(filename="E:/FigD7a4.png",width=700,height=800)
Fig1b(Resu[[4]],type="B2-F2")
dev.off()

Resu2 <- vector(mode="list",length=2)
Resu2[[1]] <- ReadData("B2F1","BC",Yr1=1970,Yr2=2023,Nboot=500,PDF=T)
Resu2[[2]] <- ReadData("B2F2","BC",Yr1=1970,Yr2=2023,Nboot=500,PDF=T)

png(filename="E:/FigD8.png",width=700,height=800)
Fig1c(Resu2,labs=c("B2-F1","B2-F2"),type="B2 breeding stock")
dev.off()


