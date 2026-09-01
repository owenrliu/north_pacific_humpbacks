setwd("C:\\Research\\Iwc25\\NP_Humpbacks\\")

ReadData <- function(Code,Exts,Yr1=1965,Yr2=2023)
{
 Nyear <- Yr2-Yr1+1 
 Years <- 1980:2021
 Use <- Years-Yr1+1
 Summary <- array(0,dim=c(length(Exts),Nyear,2))  
 for (II in 1:length(Exts))
  {
   FileName <- paste0("Diags/",Code,Exts[II],".out")
   print(FileName)
   DatFile <- read.table(FileName,fill=T,col.names=c(1:200),comment.char="?")
   Index <- which(DatFile[,1]=="Total" & DatFile[,2]=="abundance:")
   for (Iyr in Yr1U[II]:Yr2)
    {
     Summary[II,Iyr-Yr1+1,1] <- as.numeric(DatFile[Index+Iyr-Yr1U[II]+1,2]) 
     Summary[II,Iyr-Yr1+1,2] <- as.numeric(DatFile[Index+Iyr-Yr1U[II]+1,3]) 
    }
  }
 ymax <- max(Summary[,,1])*1
 plot(Years,Summary[1,Use,1],ylim=c(0,ymax),lty=1,type="l",lwd=2,axes=F,ylab="",xlab="")
 box(); axis(1); axis(2)
 for (II in 2:length(Exts)) lines(Years,Summary[II,Use,1],lty=2) 
 lines(Years,Summary[1,Use,1],lwd=2)
 upp <- Summary[1,Use,1] *exp(-1.96*Summary[1,Use,2])
 low <- Summary[1,Use,1] *exp(1.96*Summary[1,Use,2])
 lines(Years,low,col="red")
 lines(Years,upp,col="red")
 mtext(Code,side=3,line=0.5,cex=1.4)
 mtext("Year",side=1,line=0.5,cex=1.4,outer=T)
 mtext("Total abundance",side=2,line=0.5,cex=1.4,outer=T)
 }

# ======================================================================================


Exts <- c("BC","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S15","S16","S17","S18","S19","S20","S21","S22")
Yr1U <- c(rep(1970,14),1965,1975,rep(1970,5))  
#Exts <- c("BC","S1","S2","S5","S6","S7","S8","S9","S12","S13","S16","S17","S19","S20","S21","S22")
#Yr1U <- c(rep(1970,10),1965,1975,rep(1970,4))  


png("D:/FigF2.png",width=700,height=800)
par(mfrow=c(2,2),oma=c(3,3,1,1),mar=c(4,3,2,2))
Codes <- c("B1F1","B1F2","B2F1","B2F2")
for (Icode in 1:4)
 {
  Code <- Codes[Icode]  
  ReadData(Code,Exts)
}
dev.off()
