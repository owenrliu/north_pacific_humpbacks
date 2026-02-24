# ==========================================================================================
require(lubridate)
require(here)

WriteOut <- function(Code,Abbrev,Yr1,Yr2,BreedNames,FeedNames,rept,rep,rep2,StockDef,data,subdir=""){
  # create file directory if it does not exist yet
  fildir <- here("Diags",subdir)
  if(!dir.exists(fildir)) dir.create(fildir)
  # save data and rept as rds for easier access to plotting later
  outs <- list(Code=Code,SensCase=Abbrev,input=data,report=rept,sdreport=rep,sdfixed=rep2)
  readr::write_rds(outs,paste0(fildir,"/",Code,Abbrev,".rds"))
  
  ### Write out relevant information into a basic text file
  Nyear <- Yr2-Yr1+1
  Nage <- data$Nage
  
  Nbreed <- length(BreedNames)
  Nfeed <- length(FeedNames)
  
  FileName <- here("Diags",subdir,paste0(Code,Abbrev,".Out"))
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
  # Changed indexing since moving to full age structure
  # Total abundance is now only the first row of LogNT
  whichLogNT <- which(rownames(rep)=="LogNT")
  whichLogNb <- which(rownames(rep)=="LogNb")
  whichLogNf <- which(rownames(rep)=="LogNf")
  # Reform the matrices/arrays from the reported vectors, to make sure we are
  # pulling the right data
  LogNT <- rep[whichLogNT,1] |> matrix(nrow=3)
  sdLogNT <- rep[whichLogNT,2] |> matrix(nrow=3)
  
  LogNb <- rep[whichLogNb,1] |> array(dim=c(3,Nbreed,length(Yr1:(Yr2+1))))
  sdLogNb <- rep[whichLogNb,2] |> array(dim=c(3,Nbreed,length(Yr1:(Yr2+1))))
  
  LogNf <- rep[whichLogNf,1] |> array(dim=c(3,Nfeed,length(Yr1:(Yr2+1))))
  sdLogNf <- rep[whichLogNf,2] |> array(dim=c(3,Nfeed,length(Yr1:(Yr2+1))))
  
  write("\nTotal abundance: Estimate SD_log",FileName,append=T)
  for (iyr in Yr1:(Yr2+1))
  { xx <- c(iyr,exp(LogNT[1,iyr-Yr1+1]),sdLogNT[1,iyr-Yr1+1]); write(xx,FileName,append=T,ncol=3); }
  
  write("\nBreeding stock abundance: Estimate SD_log Estimate/K",FileName,append=T)
  for (Ibreed in 1:Nbreed)
  {
    write(unlist(BreedNames)[Ibreed],FileName,append=T)  
    # Ipnt <- seq(from=Ibreed,by=Nbreed,length=Nyear)
    for (iyr in seq_along(Yr1:Yr2))
    { xx <- c(iyr,exp(LogNb[1,Ibreed,iyr]),sdLogNb[1,Ibreed,iyr],exp(LogNb[1,Ibreed,iyr])/BreedK[Ibreed]*100) 
    write(xx,FileName,append=T,ncol=4); }
  }
  
  write("\nFeeding stock abundance: Estimate SD_log Estimate/K",FileName,append=T)
  for (Ifeed in 1:Nfeed)
  {
    write(unlist(FeedNames)[Ifeed],FileName,append=T)  
    Ipnt <- seq(from=Ifeed,by=Nfeed,length=Nyear)
    for (iyr in seq_along(Yr1:Yr2))
    { xx <- c(iyr,exp(LogNf[1,Ifeed,iyr]),sdLogNf[1,Ifeed,iyr],exp(LogNf[1,Ifeed,iyr])/FeedK[Ifeed]*100) 
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
      Vec1 <- rep(0,Nfeed)
      for (Ifeed in 1:Nfeed) Vec1[Ifeed] <- sum(rept$NNN[Ibreed,Ifeed,Iyr-Yr1+1,1:Nage])
      xx <- c(Iyr,round(Vec1/sum(Vec1),5))
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
    # print(paste0("Values for breeding stock: "))
    # print(c("Year",unlist(FeedNames)))
    
    for (Iyr in Yr1:Yr2)
    {
      Vec1 <- rep(0,Nfeed)
      for (Ifeed in 1:Nfeed) Vec1[Ifeed] <- sum(rept$NNN[Ibreed,Ifeed,Iyr-Yr1+1,1:Nage])
      xx <- c(Iyr,round(Vec1,0))
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
