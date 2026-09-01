setwd("C:\\Research\\Iwc26\\NP_humpbacks\\")

ReadData <- function(Code,Exts,Yr1=1965,Yr2=2023,Yr1U)
{
 Nyear <- Yr2-Yr1+1 
 Years <- 1980:2021; Nyear <- length(Years)
 Use <- Years-Yr1+1
 Summary <- array(0,dim=c(length(Exts),Nyear,2))  
 for (II in 1:length(Exts))
  {
   FileName <- paste0("Diags/",Code,Exts[II],".out")
   print(FileName)
   DatFile <- read.table(FileName,fill=T,col.names=c(1:200),comment.char="?")
   Index <- which(DatFile[,1]=="Total" & DatFile[,2]=="abundance:")
   
   Skip <- Yr2+1-Yr1U[II]+1
   for (Iyr in 1:Nyear)  Summary[II,Iyr,1] <- as.numeric(DatFile[Index+Iyr+Skip+Years[1]-Yr1U[II],3]) 
   for (Iyr in 1:Nyear)  Summary[II,Iyr,2] <- as.numeric(DatFile[Index+Iyr+Skip+Years[1]-Yr1U[II],4])
   }
 ActT <- Summary[1,,1]
 SDsT <- Summary[1,,2]
 print(ActT)

 # Plot for total 
 ymax <-max(Summary[,,1],ActT*exp(1.645*SDsT))
 plot(Years,ActT,xlab="",ylab="",type="l",ylim=c(0,ymax),yaxs="i",lwd=2,col="black")
 lines(Years,ActT*exp(1.645*SDsT),lty=1,col="red",lwd=2)
 lines(Years,ActT*exp(-1.645*SDsT),lty=1,col="red",lwd=2)
 mtext(Code,side=3,line=1)
 for (II in 2:length(Exts))
  lines(Years,Summary[II,,1],lty=3) 
}

# ======================================================================================

Exts <- c("BC","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S15","S16","S17","S18","S19","S20","S21","S22")
Yr1U <- c(rep(1970,14),1965,1975,rep(1970,6))
print(Yr1U)

Codes <- c("B1F1","B1F2","B2F1","B2F2")
png(filename="E:/Fig5.png",width=700,height=800)
par(mfrow=c(2,2),oma=c(3,4,3,4),mar=c(4,4,2,2))

for (Icode in 1:4)
 {
  Code <- Codes[Icode]  
  ReadData(Code,Exts,Yr1U=Yr1U)
}

dev.off()