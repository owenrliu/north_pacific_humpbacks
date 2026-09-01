setwd("C:\\Research\\Iwc25\\NP_humpbacks\\")

ReadData <- function(Code,Exts,Yr1=1965,Yr2=2023)
{
 Nyear <- Yr2-Yr1+1 
 Years <- 1980:2021
 Use <- Years-Yr1+1
 Summary <- matrix(0,ncol=2,nrow=length(Exts))  
 for (II in 1:length(Exts))
  {
   FileName <- paste0("Diags/",Code,Exts[II],".out")
   print(FileName)
   DatFile <- read.table(FileName,fill=T,col.names=c(1:200),comment.char="?")
   Index <- which(DatFile[,1]=="Total" & DatFile[,2]=="objective")
   #print(Index)
   Summary[II,1] <- as.numeric(DatFile[Index,4]) 
   Summary[II,2] <- as.numeric(DatFile[Index+1,4]) 
  }
 print(Summary)
 write(t(cbind(Exts,Summary)),file=FileNameO,ncol=3,append=T,sep=",")
 
}

# ======================================================================================

Exts <- c("BC","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S15","S16","S17","S18","S19","S20","S21","S22")
Yr1U <- c(rep(1970,14),1965,1975)  

Codes <- c("B1F1","B1F2","B2F1","B2F2")
FileNameO <- "PaperPlots/Table4RTMB.csv"
write("",FileNameO)
for (Icode in 1:4)
 {
  Code <- Codes[Icode]  
  ReadData(Code,Exts)
}
