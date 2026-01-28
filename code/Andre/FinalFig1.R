setwd("C:\\Research\\Iwc25\\NP_Humpbacks\\")

ReadData <- function(Code,Code2,Yr1=1965,Yr2=2023,Breeds,Feeds)
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
  FileName <- paste0("Diags/",Code,Code2,".Out")
  
  DatFile <- read.table(FileName,fill=T,col.names=c(1:200),comment.char="?")
  
  Index <- which(DatFile[,1]=="Total" & DatFile[,2]=="abundance:")
  Total <- matrix(0,nrow=Nyear,ncol=2)  
  for (Iy in 1:Nyear)
  { Total[Iy,1] <-as.numeric(DatFile[Index+Iy,2]); Total[Iy,2] <-as.numeric(DatFile[Index+Iy,3]); }

  Index <- which(DatFile[,1]=="Breeding" & DatFile[,3]=="abundance:")
  Breed <- array(0,dim=c(Nbreed,Nyear,2))
  for (Ibreed in 1:Nbreed)
   {
    Ioffset <- Index+(Ibreed-1)*(Nyear+1)+1
    for (Iy in 1:Nyear)
    { Breed[Ibreed,Iy,1] <-as.numeric(DatFile[Ioffset+Iy,2]); Breed[Ibreed,Iy,2] <-as.numeric(DatFile[Ioffset+Iy,3]); }
   }
  Index <- which(DatFile[,1]=="Feeding" & DatFile[,3]=="abundance:")
  Feed <- array(0,dim=c(Nfeed,Nyear,2))
  for (Ifeed in 1:Nfeed)
   {
    Ioffset <- Index+(Ifeed-1)*(Nyear+1)+1
    for (Iy in 1:Nyear)
    { Feed[Ifeed,Iy,1] <-as.numeric(DatFile[Ioffset+Iy,2]); Feed[Ifeed,Iy,2] <-as.numeric(DatFile[Ioffset+Iy,3]); }
   }

  Outs <- NULL
  Outs$Code <- Code2
  Outs$Code2 <- Code2
  Outs$Yr1 <- Yr1
  Outs$Yr2 <- Yr2
  Outs$Yrs <- Yr2-Yr1+1
  Outs$Years <- Yr1:Yr2
  Outs$Nbreed <- Nbreed
  Outs$Nfeed <- Nfeed
  Outs$BreedNames <- BreedNames
  Outs$FeedNames <- FeedNames
  Outs$Breeds <- Breeds
  Outs$Feeds <- Feeds
  Outs$Breed <- Breed
  Outs$Feed <- Feed
  Outs$Total <- Total
  return(Outs)
}
## ==================================================================================================================
## ==================================================================================================================

Fig1A <- function(Resu,Years,BaseM=4)
 {
  Nmodel <-length(Resu)
  col <- c("red","red","blue","blue")
  lty <- c(1,3,1,3)
  
  # Find maximum
  ymax <- 0
  for (II in 1:Nmodel) 
    { Use <- which(Resu[[II]]$Years %in% Years); ymax <- max(ymax,Resu[[II]]$Total[Use,1]) }
  ymax <- ymax*1.2
  
  png("D:\\FigF1.png",width=600,height=800)
  par(mfrow=c(1,1),oma=c(2,4,20,4))
  Use <- which(Resu[[1]]$Years %in% Years)
  plot(Years,Resu[[1]]$Total[Use,1],ylim=c(0,ymax),yaxs="i",type="n",xlab="Year",ylab="Abundance")
  
  for (Imodel in 1:Nmodel)
   {
    Use <- which(Resu[[Imodel]]$Years %in% Years)
    Est <-Resu[[Imodel]]$Total[Use,1]
    Sd <- Resu[[Imodel]]$Total[Use,2]
    if (Imodel == 1) Lab <- cbind(Years,Est,Est*exp(-1.96*Sd),Est*exp(1.96*Sd))
    if (Imodel > 1) Lab <- cbind(Lab,Est,Est*exp(-1.96*Sd),Est*exp(1.96*Sd))
    lines(Years,Resu[[Imodel]]$Total[Use,1],lty=lty[Imodel],col=col[Imodel])
    if (Imodel==BaseM)
     {
      lines(Years,Est*exp(-1.96*Sd),lty=lty[Imodel],col=col[Imodel])
      lines(Years,Est*exp(1.96*Sd),lty=lty[Imodel],col=col[Imodel])
     }
   }
  write.csv(Lab,file="D:\\For_Debi.csv")
  legend("bottomright",leg=c("B1F1","B1F2","B2F1","B2F2"),lty=lty,col=col,cex=1,lwd=1)
  dev.off()
}

## ==================================================================================================================
  
Fig1B <- function(Resu,Years,BaseM=4,IsBreed,Names)
{
  Nmodel <-length(Resu)
  col <- c("red","red","blue","blue")
  lty <- c(1,3,1,3)
  Narea <- length(Names)
  Case <- c("B1F1","B1F2","B2F1","B2F2")
  
  write("",file="D:\\For_DebiB.csv")
  png("D:\\FigF2.png",width=600,height=800)
  par(mfrow=c(3,2),oma=c(2,2,2,2),mar=c(5,4,0,1)) 
  for (Iarea in 1:Narea)
   {
    ymax <- 0
    Ests <- matrix(0,nrow=length(Years),ncol=Nmodel)
    Sds <- matrix(0,nrow=length(Years),ncol=Nmodel)
    NQmodel <- 0
    for (Imodel in 1:Nmodel)
     if (sum(Resu[[Imodel]]$Breeds[Iarea,]) > 0)
      {
       NQmodel <- NQmodel + 1
       Use <- which(Resu[[Imodel]]$Years %in% Years)
       for (II in 1:2)
        if (Resu[[Imodel]]$Breeds[Iarea,II]>0)
         {
          Jarea <- Resu[[Imodel]]$Breeds[Iarea,II]
          Ests[,Imodel] <- Ests[,Imodel]+ Resu[[Imodel]]$Breed[Jarea,Use,1] 
          Sds[,Imodel] <- Sds[,Imodel]+ (Resu[[Imodel]]$Breed[Jarea,Use,1]*Resu[[Imodel]]$Breed[Jarea,Use,2])^2
         }
        Sds[,Imodel] <- sqrt(Sds[,Imodel])/Ests[,Imodel]
        if (NQmodel == 1) Lab <- cbind(Years,Ests[,Imodel],Ests[,Imodel]*exp(-1.96*Sds[,Imodel]),Ests[,Imodel]*exp(1.96*Sds[,Imodel]))
        if (NQmodel > 1) Lab <- cbind(Lab,Ests[,Imodel],Ests[,Imodel]*exp(-1.96*Sds[,Imodel]),Ests[,Imodel]*exp(1.96*Sds[,Imodel]))
        if (NQmodel == 1) { Tname <- c(" "); Lname <- "Year," }
        Tname <- c(Tname,rep(Case[Imodel],3))
        Lname <- c(Lname,c("Est,low95,upp95,"))
        ymax <- max(ymax,Ests)
      }
    ymax <- ymax * 1.2

    plot(Years,rep(-100,length(Years)),ylim=c(0,ymax),yaxs="i",type="n",xlab="",ylab="")
    for (Imodel in 1:Nmodel)
     if (sum(Resu[[Imodel]]$Breeds[Iarea,]) > 0)
      {
       lines(Years,Ests[,Imodel],lty=lty[Imodel],col=col[Imodel]) 
       if (Imodel==BaseM[Iarea])
        {
         lines(Years,Ests[,Imodel]*exp(-1.96*Sds[,Imodel]),lty=lty[Imodel],col=col[Imodel])
         lines(Years,Ests[,Imodel]*exp(1.96*Sds[,Imodel]),lty=lty[Imodel],col=col[Imodel])
        }
     }
    text(1980,ymax*0.9,Names[Iarea],adj=0,cex=1.4)
    write(Names[Iarea],file="D:\\For_DebiB.csv",append=T,sep=',')
    write(Tname,file="D:\\For_DebiB.csv",append=T,sep=',',ncol=3*NQmodel+1)
    write(Lname,file="D:\\For_DebiB.csv",append=T,ncol=3*NQmodel+1)
    write(t(Lab),file="D:\\For_DebiB.csv",append=T,sep=',',ncol=3*NQmodel+1)
    
    if (Iarea==2)
     legend("bottomright",leg=c("B1F1","B1F2","B2F1","B2F2"),lty=lty,col=col,cex=1.2,lwd=2)
    
   } # Iarea
    
  mtext("Year",side=1,line=0,outer=T)
  mtext("Abundance",side=2,line=0,outer=T)
  dev.off() 
 }

## ==================================================================================================================


Fig1C <- function(Resu,Years,BaseM=4,IsBreed,Names,Names2)
{
  Nmodel <-length(Resu)
  col <- c("red","red","blue","blue")
  lty <- c(1,3,1,3)
  Narea <- length(Names)
  Case <- c("B1F1","B1F2","B2F1","B2F2")
  
  write("",file="D:\\For_DebiC.csv")
  png("D:\\FigF3.png",width=600,height=800)
  par(mfrow=c(3,2),oma=c(2,2,2,2),mar=c(5,4,0,1)) 
  for (Iarea in 1:Narea)
  {
    ymax <- 0
    Ests <- matrix(0,nrow=length(Years),ncol=Nmodel)
    Sds <- matrix(0,nrow=length(Years),ncol=Nmodel)
    NQmodel <- 0
    for (Imodel in 1:Nmodel)
      if (sum(Resu[[Imodel]]$Feeds[Iarea,]) > 0)
       {
        NQmodel <- NQmodel + 1
        Use <- which(Resu[[Imodel]]$Years %in% Years)
        for (II in 1:2)
          if (Resu[[Imodel]]$Feeds[Iarea,II]>0)
          {
            Jarea <- Resu[[Imodel]]$Feeds[Iarea,II]
            Ests[,Imodel] <- Ests[,Imodel]+ Resu[[Imodel]]$Feed[Jarea,Use,1] 
            Sds[,Imodel] <- Sds[,Imodel]+ (Resu[[Imodel]]$Feed[Jarea,Use,1]*Resu[[Imodel]]$Feed[Jarea,Use,2])^2
          }
        Sds[,Imodel] <- sqrt(Sds[,Imodel])/Ests[,Imodel]
        if (NQmodel == 1) Lab <- cbind(Years,Ests[,Imodel],Ests[,Imodel]*exp(-1.96*Sds[,Imodel]),Ests[,Imodel]*exp(1.96*Sds[,Imodel]))
        if (NQmodel > 1) Lab <- cbind(Lab,Ests[,Imodel],Ests[,Imodel]*exp(-1.96*Sds[,Imodel]),Ests[,Imodel]*exp(1.96*Sds[,Imodel]))
        if (NQmodel == 1) { Tname <- c(" "); Lname <- "Year," }
        Tname <- c(Tname,rep(Case[Imodel],3))
        Lname <- c(Lname,c("Est,low95,upp95,"))
        ymax <- max(ymax,Ests)
      }
    ymax <- ymax * 1.2
    
    plot(Years,rep(-100,length(Years)),ylim=c(0,ymax),yaxs="i",type="n",xlab="",ylab="")
    for (Imodel in 1:Nmodel)
      if (sum(Resu[[Imodel]]$Feeds[Iarea,]) > 0)
      {
        lines(Years,Ests[,Imodel],lty=lty[Imodel],col=col[Imodel]) 
        if (Imodel==BaseM[Iarea])
        {
          lines(Years,Ests[,Imodel]*exp(-1.96*Sds[,Imodel]),lty=lty[Imodel],col=col[Imodel])
          lines(Years,Ests[,Imodel]*exp(1.96*Sds[,Imodel]),lty=lty[Imodel],col=col[Imodel])
        }
      }
    text(1980,ymax*0.9,Names[Iarea],adj=0,cex=1.4)
    write(Names2[Iarea],file="D:\\For_DebiC.csv",append=T,sep=',')
    write(Tname,file="D:\\For_DebiC.csv",append=T,sep=',',ncol=3*NQmodel+1)
    write(Lname,file="D:\\For_DebiC.csv",append=T,ncol=3*NQmodel+1)
    write(t(Lab),file="D:\\For_DebiC.csv",append=T,sep=',',ncol=3*NQmodel+1)
    
    if (Iarea==2)
      legend("bottomright",leg=c("B1F1","B1F2","B2F1","B2F2"),lty=lty,col=col,cex=1.2,lwd=2)
    
  } # Iarea
  
  mtext("Year",side=1,line=0,outer=T)
  mtext("Abundance",side=2,line=0,outer=T)
  dev.off() 
}


## ==================================================================================================================
Resu <- vector(mode="list",length=4)
Resu[[1]] <- ReadData("B1F1","BC",Yr1=1970,Yr2=2023,Breeds<-matrix(c(1,0, 2,0, 3,0, 0,0, 0,0, 4,0),ncol=2,byrow=T),Feeds<-matrix(c(1,0,2,3,4,0,5,0,6,0),ncol=2,byrow=T))
Resu[[2]] <- ReadData("B1F2","BC",Yr1=1970,Yr2=2023,Breeds<-matrix(c(1,0, 2,0, 3,0, 0,0, 0,0, 4,0),ncol=2,byrow=T),Feeds<-matrix(c(1,0,2,3,4,0,5,0,6,0),ncol=2,byrow=T))
Resu[[3]] <- ReadData("B2F1","BC",Yr1=1970,Yr2=2023,Breeds<-matrix(c(1,0, 2,0, 0,0, 3,0, 4,0, 5,0),ncol=2,byrow=T),Feeds<-matrix(c(1,0,2,3,4,0,5,0,6,0),ncol=2,byrow=T))
Resu[[4]] <- ReadData("B2F2","BC",Yr1=1970,Yr2=2023,Breeds<-matrix(c(1,0, 2,0, 0,0, 3,0, 4,0, 5,0),ncol=2,byrow=T),Feeds<-matrix(c(1,0,2,3,4,0,5,0,6,0),ncol=2,byrow=T))
Fig1A(Resu,Years=1980:2021)
Fig1B(Resu,Years=1980:2021,IsBreed=T,BaseM=c(4,4,2,4,4,4),Names=c("Asia","Hawaii","Mexico","Archipielago de Revillagigedo","Mainland Mexico","Central America"))
Fig1C(Resu,Years=1980:2021,IsBreed=F,BaseM=c(4,4,4,4,4),Names=c("RUSSIA+\nWESTERN ALEUTIAN ISLANDS",
                                                                "EASTERN ALEUTIAN+BERING SEA+\nGULF OF ALASKA",
                                                                "SOUTHEAST ALASKA+\nNORTHERN BRITISH COLUMBIA","SOUTHERN BRITISH COLUMBIA+\nWASHINGTON",
                                                                "OREGON+CALIFORNIA"),
      Names2=c("RUSSIA+WESTERN ALEUTIAN ISLANDS",
              "EASTERN ALEUTIAN+BERING SEA+GULF OF ALASKA",
              "SOUTHEAST ALASKA+NORTHERN BRITISH COLUMBIA","SOUTHERN BRITISH COLUMBIA+WASHINGTON",
              "OREGON+CALIFORNIA"))

