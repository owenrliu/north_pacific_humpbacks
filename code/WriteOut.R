# ==========================================================================================

WriteOut <- function(Code,Abbrev,Yr1,Yr2,BreedNames,FeedNames,rept,rep,rep2,StockDef,data)
{
  
  # save rept as rds
  readr::write_rds(rept,paste0("Diags/",Code,Abbrev,".rds"))
  
  ###
  Nyear <- Yr2-Yr1+1
  
  Nbreed <- length(BreedNames)
  Nfeed <- length(FeedNames)
  
  FileName <- paste0("Diags/",Code,Abbrev,".Out")
  write(paste0("Stock structure: ",Code,"; Model type = ",Abbrev),FileName)
  if (StockDef$StochSopt==0) write("Mortality impacts breeding stocks",FileName,append=T)
  if (StockDef$StochSopt==1) write("Mortality impacts feeding grounds",FileName,append=T)
  if (StockDef$DensDepOpt==0) write("Recuitment related to unfished",FileName,append=T)
  if (StockDef$DensDepOpt==1) write("Recuitment related to unfished and current",FileName,append=T)
  write(paste0("Staying rate: ",StockDef$StrayBase),FileName,append=T)
  write(paste0("Age-at-matrity: ",StockDef$IAmat),FileName,append=T)
  write(paste0("Base adult survival: ",StockDef$TimeLag),FileName,append=T)
  if (StockDef$WithMirror) write("Mirrored survival",FileName,append=T)
  write(paste0("Sigma for survival devs: ",StockDef$SigmaDevS),FileName,append=T)
  if (StockDef$AddCV) write("Additional variance estimated",FileName,append=T)
  write(paste0("Weight for mark-recapture proportions: ",StockDef$MixWeights[1]),FileName,append=T)
  write(paste0("Weight for genetics proportions: ",StockDef$MixWeights[2]),FileName,append=T)
  
  # Log-likelihoods
  write("\n",FileName,append=T)
  write(paste0("Total objective fn: ",rept$neglogL),FileName,append=T)
  write(paste0("Negative log-likelihood fn: ",rept$datalike),FileName,append=T)
  
  write("\nNatural Survival by feeding ground",FileName,append=T)
  write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
  for (iyr in Yr1:Yr2)
  { xx <- c(iyr,rept$SurvOutF[,iyr-Yr1+1]); write(xx,FileName,append=T,ncol=Nfeed+1); }
  write("\nNaturak Survival by breeding stock",FileName,append=T)
  write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
  for (iyr in Yr1:Yr2)
  { xx <- c(iyr,rept$SurvOutB[,iyr-Yr1+1]); write(xx,FileName,append=T,ncol=Nbreed+1); }
  write("\nFecundity by feeding ground",FileName,append=T)
  write(t(unlist(FeedNames)),FileName,append=T,sep=',',ncol=Nfeed)
  for (iyr in Yr1:Yr2)
  { xx <- c(iyr,rept$Fec[,iyr-Yr1+1]); write(xx,FileName,append=T,ncol=Nfeed+1); }
  
  BreedK <- rept$BreedK
  FeedK <- rept$FeedK
  
  Index <- which(rownames(rep)=="LogNT") 
  LogNT<- rep[Index,]
  Index <- which(rownames(rep)=="LogNb")
  LogNb<- rep[Index,]
  Index <- which(rownames(rep)=="LogNf")
  LogNf<- rep[Index,]
  write("\nTotal abundance: Estimate SD_log",FileName,append=T)
  for (iyr in Yr1:(Yr2+1))
  { xx <- c(iyr,exp(LogNT[iyr-Yr1+1,1]),LogNT[iyr-Yr1+1,2]); write(xx,FileName,append=T,ncol=3); }
  
  write("\nBreeding stock abundance: Estimate SD_log Estimate/K",FileName,append=T)
  for (Ibreed in 1:Nbreed)
  {
    write(unlist(BreedNames)[Ibreed],FileName,append=T)  
    Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
    for (iyr in Yr1:Yr2)
    { II <- Ipnt[iyr-Yr1+1]; xx <- c(iyr,exp(LogNb[II,1]),LogNb[II,2],exp(LogNb[II,1])/BreedK[Ibreed]*100); 
    write(xx,FileName,append=T,ncol=4); }
  }
  
  write("\nFeeding stock abundance: Estimate SD_log Estimate/K",FileName,append=T)
  for (Ifeed in 1:Nfeed)
  {
    write(unlist(FeedNames)[Ifeed],FileName,append=T)  
    Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
    for (iyr in Yr1:Yr2)
    { II <- Ipnt[iyr-Yr1+1]; xx <- c(iyr,exp(LogNf[II,1]),LogNf[II,2],exp(LogNf[II,1])/FeedK[Ifeed]*100); 
    write(xx,FileName,append=T,ncol=4); }
  }
  
  write("\nParameters: Estimate  Std. Error      z value    Pr(>|z^2|)",FileName,append=T)
  names <- row.names(rep2)
  for (Irow in 1:length(rep2[,1]))
    write(c(names[Irow],round(rep2[Irow,],8)),FileName,ncol=5,append=T) 
  
  write("\nMixing matrix",FileName,append=T)
  for (Ibreed in 1:length(rept$Mix[,1]))
  {
    xx <- paste(c(round(rept$Mix[Ibreed,],8)))
    xx <- c(xx,",#",unlist(BreedNames[Ibreed]))
    write(xx,FileName,append=T,ncol=Nfeed+2,sep=" ")
  }
  
  write("\nSplit of breeding stock to feeding grounds over time",FileName,append=T)
  for (Ibreed in 1:(length(BreedNames)))
  {
    write(paste0("Split for: ",unlist(BreedNames)[Ibreed]),FileName,append=T)
    write(c("Year",unlist(FeedNames)),FileName,append=T,ncol=Nfeed+1)
    for (Iyr in Yr1:(Yr2+1))
    {
      xx <- c(Iyr,round(rept$NN[Ibreed,,Iyr-Yr1+1]/sum(rept$NN[Ibreed,,Iyr-Yr1+1]),5))
      write(xx,FileName,append=T,ncol=Nfeed+1)
    }
  }
  
  write("\nDifference in mortality over time",FileName,append=T)
  for (Ibreed in 1:(length(BreedNames)))
  {
    write(paste0("Difference for: ",unlist(BreedNames)[Ibreed]),FileName,append=T)
    write(c("Year",unlist(FeedNames)),FileName,append=T,ncol=Nfeed+1)
    for (Iyr in Yr1:Yr2)
    {
      xx <- c(Iyr,round(rept$MortDiff[Ibreed,,Iyr-Yr1+1],0))
      write(xx,FileName,append=T,ncol=Nfeed+1)
    }
  }
  
  write("\nTotal survival over time",FileName,append=T)
  for (Ibreed in 1:(length(BreedNames)))
  {
    write(paste0("Values for breeding stock: ",unlist(BreedNames)[Ibreed]),FileName,append=T)
    write(c("Year",unlist(FeedNames)),FileName,append=T,ncol=Nfeed+1)
    for (Iyr in Yr1:Yr2)
    {
      xx <- c(Iyr,round(rept$SurvTot[Ibreed,,Iyr-Yr1+1],3))
      write(xx,FileName,append=T,ncol=Nfeed+1)
    }
  }
  
  write("\nNumbers over time by herd",FileName,append=T)
  for (Ibreed in 1:(length(BreedNames)))
  {
    write(paste0("Values for breeding stock: ",unlist(BreedNames)[Ibreed]),FileName,append=T)
    write(c("Year",unlist(FeedNames)),FileName,append=T,ncol=Nfeed+1)
    for (Iyr in Yr1:Yr2)
    {
      xx <- c(Iyr,round(rept$NN[Ibreed,,Iyr-Yr1+1],0))
      write(xx,FileName,append=T,ncol=Nfeed+1)
    }
  }
  
  write("\nCatch numbers by breeding stock",FileName,append=T)
  write(c("Year",unlist(BreedNames)),FileName,append=T,ncol=Nbreed+1)
  write(t(cbind(Yr1:Yr2,data$CatchB)),FileName,append=T,ncol=Nbreed+1)
  write("\nCatch numbers by feeding ground",FileName,append=T)
  write(c("Year",unlist(FeedNames)),FileName,append=T,ncol=Nfeed+1)
  write(t(cbind(Yr1:Yr2,data$CatchF)),FileName,append=T,ncol=Nfeed+1)
}
