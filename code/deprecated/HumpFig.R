
Plots1 <- function(Dirn,FileName,FileName2,Yr1=1656,Yr2=2020,Nbreed=4,Nfeed=6,PDF=F)
{
 PDFFile <- paste("D:\\Results",Dirn,".PDF",sep="")
 print(PDFFile)

 if (PDF==T) pdf(PDFFile)
  
 Folder <- paste("C:\\Research\\iwc23\\IWC23B\\",Dirn,"\\",sep="")
 print(Folder)
 setwd(Folder)
 Nyear <- Yr2-Yr1+1
 print(Nyear)

 if (Case=="B1F1")
 {
    labels <- c("Asia","Hawaii","Mex","Cent Am")
    labelsF <- c("RUS+WAL","EAL+BER","GOA","SEA-NBC","SBC-WA","OR-CA")
 }
 if (Case=="B1F2")
 {
    labels <- c("Asia","Hawaii","Mex","Cent Am")
    labelsF <- c("RUS+WAL","EAL+BE+WGOA","NGOA","SEA-NBC","SBC-WA","OR-CA")
 }
 
 Nmat <- read.table(FileName,skip=4,nrow=Nyear)
 Years <- Nmat[,1]
 Nmat <- Nmat[,-1]

 par(mfrow=c(2,1),mar=c(5,4,2,2)+0.1,oma=c(3,3,3,3))
 NN <- apply(Nmat[,1:Nbreed],1,sum)
 ymax <- max(NN)
 plot(Years,NN,xlab="Year",ylab="Total population size",ylim=c(0,ymax),type="l")
 Y1 <- 1900-Yr1+1
 plot(Years[Y1:Nyear],NN[Y1:Nyear],xlab="Year",ylab="Total population size",ylim=c(0,ymax),type="l")
 
 #par(mfrow=c(2,1),mar=c(5,4,2,2)+0.1,oma=c(3,3,3,3))
 #NN <- apply(Nmat[,Nbreed+1:Nfeed],1,sum)
 #ymax <- max(NN)
 #plot(Years,NN,xlab="Year",ylab="Total population size",ylim=c(0,ymax),type="l")
 #Y1 <- 1900-Yr1+1
 #plot(Years[Y1:Nyear],NN[Y1:Nyear],xlab="Year",ylab="Total population size",ylim=c(0,ymax),type="l")
 
 
 par(mfrow=c(2,1),mar=c(5,4,2,2)+0.1,oma=c(3,3,3,3))
 ymax <- max(Nmat)*1.1
 plot(Years,Nmat[,1],xlab="Year",ylab="Population size",ylim=c(0,ymax),type="l")
 for (Ibreed in 2:Nbreed)
  lines(Years,Nmat[,Ibreed],lty=Ibreed) 
 
 plot(Years,Nmat[,1]/Nmat[1,1]*100,xlab="Year",ylab="Relative Population size",ylim=c(0,120),type="l")
 for (Ibreed in 2:Nbreed)
   lines(Years,Nmat[,Ibreed]/Nmat[1,Ibreed]*100,lty=Ibreed) 
 legend("bottomleft",legend=labels,lty=1:Nbreed)
 
 plot(Years,Nmat[,Nbreed+1],xlab="Year",ylab="Population size",ylim=c(0,ymax),type="l")
 for (Ifeed in 2:Nfeed)
   lines(Years,Nmat[,Nbreed+Ifeed],lty=Ifeed) 
 
 plot(Years,Nmat[,Nbreed+1]/Nmat[1,Nbreed+1]*100,xlab="Year",ylab="Relative Population size",ylim=c(0,105),type="l")
 for (Ifeed in 2:Nfeed)
   lines(Years,Nmat[,Nbreed+Ifeed]/Nmat[1,Nbreed+Ifeed]*100,lty=Ibreed) 
 legend("bottomleft",legend=labelsF,lty=1:Nfeed)

 
# =================================================================================== 
 
 print("Reading std file")
 StdStuff <- read.table(FileName2,skip=1)
 Index <- which(StdStuff[,2] %in% "LogNb" |  StdStuff[,2] %in% "LogNf" |  StdStuff[,2] %in% "LogNT")
 StdStuff <- StdStuff[Index,]
 NstdYr <- length(which(StdStuff[,2] %in% "LogNb"))/Nbreed

 skip <- 4 + Nyear+1 + 2+Nbreed
 Nindex <- scan(FileName,skip=skip,n=1,quiet=T)
 Index <- read.table(FileName,skip=skip+1,nrow=Nindex)
 print(Index)
 
 par(mfrow=c(4,3))

 for (Ibreed in 1:Nbreed)
  {
   Index2 <- Index[Index[,3]== Ibreed & Index[,4]==1,] 
   Low <- Index2[,7]*exp(-1.96*Index2[,8])
   Upp <- Index2[,7]*exp(1.96*Index2[,8])
   ymax2 <- max(Nmat[,Ibreed],max(Upp))*1.1
   plot(Years,Nmat[,Ibreed],xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
   title(labels[Ibreed])
  } 
 for (Ifeed in 1:Nfeed)
 {
   Index2 <- Index[Index[,3]== Ifeed+Nbreed & Index[,4]==1,] 
   Low <- Index2[,7]*exp(-1.96*Index2[,8])
   Upp <- Index2[,7]*exp(1.96*Index2[,8])
   ymax2 <- max(Nmat[,Nbreed+Ifeed],max(Upp))*1.1
   plot(Years,Nmat[,Nbreed+Ifeed],xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
   title(labelsF[Ifeed])
 } 

 par(mfrow=c(4,3))
 Index2 <- Index[Index[,3]== -1 & Index[,4]==1,] 
 Low <- Index2[,7]*exp(-1.96*Index2[,8])
 Upp <- Index2[,7]*exp(1.96*Index2[,8])
 ymax2 <- max(Nmat[,Ibreed],max(Upp))*1.1
 plot(Years[(Nyear-35):Nyear],Nmat[(Nyear-35):Nyear,Ibreed],xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
 for (II in 1:length(Index2[,1]))
 {
   Yr <- (Index2[II,1]+Index2[II,2])/2
   #cat(Yr,Low[II],Upp[II],"\n")
   if (Index2[II,5]==1) points(Yr,Index2[II,7],pch=16)  
   if (Index2[II,5]==0) points(Yr,Index2[II,7],pch=1)  
   lines(c(Yr,Yr),c(Low[II],Upp[II]),lty=1)
 }  
 title("All")
 rows <- (Ibreed-1)*NstdYr+c(1:NstdYr)
 rows <- NstdYr*Nbreed+Nfeed*(NstdYr-1)+c(1:(NstdYr))
 Tyears <- 1970:(Yr2+1)
 lines(Tyears,exp(StdStuff[rows,3]),lty=1)
 lines(Tyears,exp(StdStuff[rows,3]-1.645*StdStuff[rows,4]),lty=2)
 lines(Tyears,exp(StdStuff[rows,3]+1.645*StdStuff[rows,4]),lty=2)

 
  for (Ibreed in 1:Nbreed)
 {
   Index2 <- Index[Index[,3]== Ibreed & Index[,4]==1,] 
   Low <- Index2[,7]*exp(-1.96*Index2[,8])
   Upp <- Index2[,7]*exp(1.96*Index2[,8])
   ymax2 <- max(Nmat[,Ibreed],max(Upp))*1.1
   plot(Years[(Nyear-35):Nyear],Nmat[(Nyear-35):Nyear,Ibreed],xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
   for (II in 1:length(Index2[,1]))
   {
     Yr <- (Index2[II,1]+Index2[II,2])/2
     #cat(Yr,Low[II],Upp[II],"\n")
     if (Index2[II,5]==1) points(Yr,Index2[II,7],pch=16)  
     if (Index2[II,5]==0) points(Yr,Index2[II,7],pch=1)  
     lines(c(Yr,Yr),c(Low[II],Upp[II]),lty=1)
   }  
   title(labels[Ibreed])
   rows <- (Ibreed-1)*NstdYr+c(1:NstdYr)
   Tyears <- 1970:(Yr2+1)
   lines(Tyears,exp(StdStuff[rows,3]),lty=1)
   lines(Tyears,exp(StdStuff[rows,3]-1.645*StdStuff[rows,4]),lty=2)
   lines(Tyears,exp(StdStuff[rows,3]+1.645*StdStuff[rows,4]),lty=2)
 } 
 
 for (Ifeed in 1:Nfeed)
 {
   Index2 <- Index[Index[,3]== Ifeed+Nbreed & Index[,4]==1,] 
   Low <- Index2[,7]*exp(-1.96*Index2[,8])
   Upp <- Index2[,7]*exp(1.96*Index2[,8])
   ymax2 <- max(Nmat[,Nbreed+Ifeed],max(Upp))*1.1
   plot(Years[(Nyear-35):Nyear],Nmat[(Nyear-35):Nyear,Nbreed+Ifeed],xlab="Year",ylab="Population size",ylim=c(0,ymax2),type="l") 
   for (II in 1:length(Index2[,1]))
   {
     Yr <- (Index2[II,1]+Index2[II,2])/2
     #cat(Yr,Low[II],Upp[II],"\n")
     if (Index2[II,5]==1) points(Yr,Index2[II,7],pch=16)
     if (Index2[II,5]==0) points(Yr,Index2[II,7],pch=1)
     lines(c(Yr,Yr),c(Low[II],Upp[II]),lty=1)
   }  
   title(labelsF[Ifeed])
   rows <- NstdYr*Nbreed+(Ifeed-1)*(NstdYr-1)+c(1:(NstdYr-1))
   Tyears <- 1970:Yr2
   lines(Tyears,exp(StdStuff[rows,3]),lty=1)
   lines(Tyears,exp(StdStuff[rows,3]-1.645*StdStuff[rows,4]),lty=2)
   lines(Tyears,exp(StdStuff[rows,3]+1.645*StdStuff[rows,4]),lty=2)
 } 
 
 # Relative index
 par(mfrow=c(2,2),mar=c(8,3,2,2)+1)
 Nrel <- sum(Index[,4]==2)
 Index2 <- Index[Index[,4]==2,] 
 if (nrow(Index2)> 0)
  {
   Low <- Index2[,5]*exp(-1.96*Index2[,6])
   Upp <- Index2[,5]*exp(1.96*Index2[,6])
   ymax2 <- max(Index2[,7],max(Upp))*1.1
   plot(0,0,xlab="",ylab="Relative index",xlim=c(0.5,Nrel+0.5),ylim=c(0,ymax2),axes=F,type="n")
   box();  axis(2)
   XLAB <- NULL
   for (II in 1:Nrel)
    XLAB <- c(XLAB,paste(labelsF[Index2[II,3]-Nbreed],(Index2[II,1]+Index2[II,2])/2))  
   points(c(1:Nrel)+0.1,Index2[,7],pch=1)
   for (II in 1:Nrel)
    {
     points(II,Index2[II,5],pch=16) 
     lines(c(II,II),c(Low[II],Upp[II]))
    } 
   axis(1,1:Nrel,XLAB,las=2) 
  } 
 print("done index")

 # ================================================================================================
 
 skip <- skip +Nindex + 2
 Nprop <- scan(FileName,skip=skip,n=1,quiet=T)
 skip <- skip + 2

 par(mfrow=c(2,1),mar=c(10,4,2,2)+0.1,oma=c(1,1,1,1))
 Nprops <- c(15,13,16,20)
 
 Nprop2a <- Nprops[1]
 Propsa <- read.table(FileName,skip=skip,nrow=Nprop2a) 
 skip <- skip + Nprop2a
 Nprop2b <- Nprops[2]
 Propsb <- read.table(FileName,skip=skip,nrow=Nprop2b) 
 skip <- skip + Nprop2b
 VecOut <- matrix(0,ncol=9,nrow=100)
 Ipnt <- 0
 cols <- NULL
 for (Ibreed in 1:4)
  for (Ifeed in 1:6) 
   {
    IndexA <- which(Propsa[,2]==Ibreed & Propsa[,3]==Ifeed) 
    IndexB <- which(Propsb[,2]==Ibreed & Propsb[,3]==Ifeed) 
    #cat(Ibreed,Ifeed,IndexA,IndexB,"\n")
    if (length(IndexA) > 0 || length(IndexB) > 0)
     {
      Ipnt <- Ipnt + 1
      if (Ipnt %% 2 ==1) cols <- c(cols,"green")
      if (Ipnt %% 2 ==1) cols <- c(cols,"red")
      if (length(IndexA) > 0) Pred <- Propsa[IndexA,6] else Pred <- Propsb[IndexB,6]
      ObsA <- c(-1,-1,-1)
      if (length(IndexA) > 0) ObsA <- c(Propsa[IndexA,4],Propsa[IndexA,4]-1.96*Propsa[IndexA,5],Propsa[IndexA,4]+1.96*Propsa[IndexA,5])
      ObsB <- c(-1,-1,-1)
      if (length(IndexB) > 0) ObsB <- c(Propsb[IndexB,4],Propsb[IndexB,4]-1.96*Propsb[IndexB,5],Propsb[IndexB,4]+1.96*Propsb[IndexB,5])
      vec <- c(Ibreed,Ifeed,Pred,ObsA,ObsB)
      VecOut[Ipnt,] <- vec
     }
   }
 VecOut <- VecOut[1:Ipnt,]
 NpropV <- length(VecOut[,1])
 #print(VecOut)
 #print(cols)

 plot(0,0,xlab="",ylab="Breeding -> Feeding proportions",xlim=c(0.5,NpropV+0.5),ylim=c(0,1.05),axes=F,type="n")
 points(c(1:NpropV),VecOut[,3],pch=17,col=cols)
 box();  axis(2)
 XLAB <- NULL
 for (II in 1:NpropV)
   XLAB <- c(XLAB,paste(labels[VecOut[II,1]],labelsF[VecOut[II,2]]))  
 for (II in 1:NpropV) 
  { points(II-0.2,VecOut[II,4],pch=1,col=cols[II]); lines(c(II-0.2,II-0.2),c(VecOut[II,5],VecOut[II,6]),col=cols[II]); }
 for (II in 1:NpropV) 
  { points(II+0.2,VecOut[II,7],pch=16,col=cols[II]); lines(c(II+0.2,II+0.2),c(VecOut[II,8],VecOut[II,9]),col=cols[II]); }
 axis(1,1:NpropV,XLAB,las=2,cex.axis=1) 
 abline(h=0,col="black")
 print(skip)

 Nprop2a <- Nprops[3]
 Propsa <- read.table(FileName,skip=skip,nrow=Nprop2a) 
 skip <- skip + Nprop2a
 Nprop2b <- Nprops[4]
 Propsb <- read.table(FileName,skip=skip,nrow=Nprop2b) 
 skip <- skip + Nprop2b
 VecOut <- matrix(0,ncol=9,nrow=100)
 Ipnt <- 0
 cols <- NULL
 for (Ifeed in 1:6) 
  for (Ibreed in 1:4)
   {
    IndexA <- which(Propsa[,2]==Ifeed & Propsa[,3]==Ibreed) 
    IndexB <- which(Propsb[,2]==Ifeed & Propsb[,3]==Ibreed) 
    #cat(Ibreed,Ifeed,IndexA,IndexB,"\n")
    if (length(IndexA) > 0 || length(IndexB) > 0)
     {
      Ipnt <- Ipnt + 1
      if (Ipnt %% 2 ==1) cols <- c(cols,"green")
      if (Ipnt %% 2 ==1) cols <- c(cols,"red")
      if (length(IndexA) > 0) Pred <- Propsa[IndexA,6] else Pred <- Propsb[IndexB,6]
      ObsA <- c(-1,-1,-1)
      if (length(IndexA) > 0) ObsA <- c(Propsa[IndexA,4],Propsa[IndexA,4]-1.96*Propsa[IndexA,5],Propsa[IndexA,4]+1.96*Propsa[IndexA,5])
      ObsB <- c(-1,-1,-1)
      if (length(IndexB) > 0) ObsB <- c(Propsb[IndexB,4],Propsb[IndexB,4]-1.96*Propsb[IndexB,5],Propsb[IndexB,4]+1.96*Propsb[IndexB,5])
      vec <- c(Ibreed,Ifeed,Pred,ObsA,ObsB)
      VecOut[Ipnt,] <- vec
     }
   }
 VecOut <- VecOut[1:Ipnt,]
 NpropV <- length(VecOut[,1])

 plot(0,0,xlab="",ylab="Feeding -> Breeding proportions",xlim=c(0.5,NpropV+0.5),ylim=c(0,1.05),axes=F,type="n")
 points(c(1:NpropV),VecOut[,3],pch=17,col=cols)
 box();  axis(2)
 XLAB <- NULL
 for (II in 1:NpropV)
   XLAB <- c(XLAB,paste(labels[VecOut[II,1]],labelsF[VecOut[II,2]]))  
 for (II in 1:NpropV) 
 { points(II-0.2,VecOut[II,4],pch=1,col=cols[II]); lines(c(II-0.2,II-0.2),c(VecOut[II,5],VecOut[II,6]),col=cols[II]); }
 for (II in 1:NpropV) 
 { points(II+0.2,VecOut[II,7],pch=16,col=cols[II]); lines(c(II+0.2,II+0.2),c(VecOut[II,8],VecOut[II,9]),col=cols[II]); }
 axis(1,1:NpropV,XLAB,las=2,cex.axis=1) 
 abline(h=0,col="red")
 print(skip)
 

 
 
 if (PDF==T) dev.off()
 
} 

Case <- "B1F2"
#Plots1("","Hump.rep","Hump.std",PDF=F)
Plots1("","report.out","Hump.std",PDF=F)
#Plots1("2","Hump.rep",PDF=T)
#Plots1("3","Hump.rep",PDF=T)
