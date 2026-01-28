setwd("C:/courses/FISH 559_26/RTMB Workshop/In Class assignments/Ex5")

rm(list=ls()) # clean

require(RTMB)
library(rstan)
library(StanHeaders)
library(tmbstan)
library(shinystan)

DoMCMCEst <- F
DoMCMCProj <- T

################################################################################

Nyear <- scan("Ex5.dat",skip=1,n=1,quiet=T)
Nclass <- scan("Ex5.dat",skip=3,n=1,quiet=T)
Length <- scan("Ex5.dat",skip=5,n=Nclass,quiet=T)
Weight <- scan("Ex5.dat",skip=7,n=Nclass,quiet=T)
X <- matrix(scan("Ex5.dat",skip=9,n=Nclass*Nclass,quiet=T),ncol=Nclass,byrow=T)
M <- scan("Ex5.dat",skip=19,n=1,quiet=T)
CWObs <- scan("Ex5.dat",skip=22,n=Nyear,quiet=T)
CALObs <- matrix(scan("Ex5.dat",skip=49,n=Nyear*Nclass,quiet=T),ncol=Nclass,byrow=T)
Neff <- scan("Ex5.dat",skip=75,n=1,quiet=T)
BioIndex <- scan("Ex5.dat",skip=78,n=Nyear,quiet=T)
BioSig <- scan("Ex5.dat",skip=104,n=1,quiet=T)

# Selectivity parameter (convert to selectivity)
S50 <- scan("Ex5.dat",skip=15,n=2,quiet=T)[1]
S95 <- scan("Ex5.dat",skip=15,n=2,quiet=T)[2]
SS50 <- scan("Ex5.dat",skip=17,n=2,quiet=T)[1]
SS95 <- scan("Ex5.dat",skip=17,n=2,quiet=T)[2]
S <- 1.0/(1+exp(-log(19.0)*(Length-S50)/(S95-S50)))
SurveyS <- 1.0/(1+exp(-log(19.0)*(Length-SS50)/(SS95-SS50)))

# Normalize the catch-at-length data
for (Iyear in 1:Nyear) CALObs[Iyear,] <- CALObs[Iyear,]/sum(CALObs[Iyear,])

#========================================================================================================================

f2 <- function(parms) {
  getAll(data, parms, warn=FALSE)
  
  # Square function
  square <- function(x) {x*x}
  
  N <- matrix(0,nrow=Nyear+Nproj+1,ncol=Nclass)
  FF <- matrix(0,nrow=Nyear+Nproj,ncol=Nclass)
  CAL <- matrix(0,nrow=Nyear+Nproj,ncol=Nclass)
  CW <- rep(0,Nyear+Nproj)
  BioPred <- rep(0,Nyear+Nproj)
  
  # First set F and Z by size-classs (note that Fproj applies after year Nyear)
  for (Iyear in 1:(Nyear+Nproj))
    for (Iclass in 1:Nclass)
     {
      if (Iyear <= Nyear)
       FF[Iyear,Iclass] <- exp(LogFullF[Iyear])*S[Iclass]
      else
       FF[Iyear,Iclass] <- Fproj*S[Iclass]
     }
  Z <- M + FF
  
  # Now set the N matrix
  for (Iclass in 1:Nclass) N[1,Iclass] <- exp(LogNinit[Iclass]);
  for (Iyear in 1:(Nyear+Nproj))
   {
    # Catch-at-length 
    CAL[Iyear,] <- FF[Iyear,]/Z[Iyear,]*N[Iyear,]*(1.0-exp(-Z[Iyear,]));
    # Total catch in weight
    CW[Iyear] <- sum(CAL[Iyear,]*Weight)
    # Catch proportion-at-length
    CAL[Iyear,] <- CAL[Iyear,]/sum(CAL[Iyear,])
    
    # Numbers-at-length
    for (Iclass in 1:Nclass)
     {
      N[Iyear+1,Iclass] <- 0;
      for (Jclass in 1:Nclass)
        N[Iyear+1,Iclass] <- N[Iyear+1,Iclass] + N[Iyear,Jclass]*exp(-Z[Iyear,Jclass])*X[Jclass,Iclass];
     }
    
    # Recruitment (watch for the index for Eps - and N)
    N[Iyear+1,1] <- N[Iyear+1,1] + exp(LogRbar)*exp(Eps[Iyear]);
   } # Iyear

  # Catch Likelihood 
  LikeCatch <- -sum(dnorm(log(CWObs[1:Nyear]),log(CW[1:Nyear]),0.05,log=T));
  
  # Biomass predictions
  for (Iyear in 1:(Nyear+Nproj)) BioPred[Iyear] <- sum(N[Iyear,]*SurveyS*Weight)     
  Top <- 0; Bot <- 0
  for (Iyear in 1:Nyear)
    { Top <- Top + log(BioIndex[Iyear]/BioPred[Iyear]); Bot <- Bot + 1.0; }
  q = exp(Top/Bot);
  # Likelihood for the index
  LikeBio <- -sum( dnorm(log(BioIndex[1:Nyear]),log(q*BioPred[1:Nyear]),BioSig,log=T));
  
  # CAL Likelihood
  LikeCAL <- 0;
  for (Iyear in 1:Nyear)
   for (Iclass in 1:Nclass)
    if (CALObs[Iyear,Iclass] > 0)
     LikeCAL <- LikeCAL - Neff*CALObs[Iyear,Iclass]*log(CAL[Iyear,Iclass]/CALObs[Iyear,Iclass]);
  
  # Recruitment penalty (include years after Nyear)
  Penal <- - sum(dnorm(Eps,0,0.6,log=T));
  
  # Total objective function
  obj_fun <- dummy*dummy + LikeCatch+LikeBio+LikeCAL+Penal;

  # Stuff to report
  REPORT(FF);
  REPORT(Z);
  REPORT(N);
  REPORT(q);
  REPORT(LikeCatch);
  REPORT(LikeBio);
  REPORT(LikeCAL);
  REPORT(Penal);
  REPORT(CAL);
  REPORT(CW);
  REPORT(BioPred);
  REPORT(obj_fun);
  #print(c(obj_fun,LikeCatch,LikeBio,LikeCAL,Penal))
  
  return(obj_fun)
 }

################################################################################
# Hint Nproj = 0 for basic estimation
# Basic estimation
################################################################################

print("Doing basic estimation")
Fproj <- 0

# Set the parameters to the initial values
LogRbar <- scan("Ex5.pin",skip=1,n=1,quiet=T)
LogNinit <- scan("Ex5.pin",skip=3,n=Nclass,quiet=T)
LogFullF <- scan("Ex5.pin",skip=5,n=Nyear,quiet=T)
Eps <- scan("Ex5.pin",skip=9,n=Nyear,quiet=T)

# Data vector
data <- list(Nyear=Nyear,Nclass=Nclass,Length=Length,Weight=Weight,X=X,S=S,SurveyS=SurveyS,M=M,
             CWObs=CWObs,CALObs=CALObs,Neff=Neff,BioIndex=BioIndex,BioSig=BioSig,Nproj=0,Fproj=0)
parameters <- list(dummy=0,LogRbar=LogRbar,LogNinit=LogNinit,LogFullF=LogFullF,Eps=Eps)
# When I was testing the code
map<-list(LogRbar=factor(NA),LogNinit=rep(factor(NA),Nclass),LogFullF=rep(factor(NA),Nyear),
          Eps=rep(factor(NA),Nyear))
# Estimate everything
map<-list(dummy=factor(NA))

print("Test application - no estimation")
model <- MakeADFun(f2, parameters,silent=T,map=map,hessian=T)
# test code - for checking for minimization
xx <- model$fn(model$env$last.par)
rpt <-model$report()
#print(rpt)
cat(rpt$obj_fun,rpt$LikeCatch,rpt$LikeCAL,rpt$LikeBio,rpt$Penal,"\n")
#AAA

print("Actual application - with estimation")
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
#print(rep)
rpt <-model$report()
cat(rpt$obj_fun,rpt$LikeCatch,rpt$LikeCAL,rpt$LikeBio,rpt$Penal,"\n")
print("Done basic estimation")

################################################################################
# Hint Nproj = 20 for projections
# Basic estimation
################################################################################

print("Doing basic estimation but with a projection period of 20 years")
Nproj = 20

# Set the parameters to initial values
LogRbar <- scan("Ex5.pin",skip=1,n=1,quiet=T)
LogNinit <- scan("Ex5.pin",skip=3,n=Nclass,quiet=T)
LogFullF <- scan("Ex5.pin",skip=5,n=Nyear,quiet=T)
Eps <- scan("Ex5.pin",skip=9,n=Nyear,quiet=T)
Eps <- c(Eps,rep(0,Nproj))

# Data vector and map
data <- list(Nyear=Nyear,Nclass=Nclass,Length=Length,Weight=Weight,X=X,S=S,SurveyS=SurveyS,M=M,
             CWObs=CWObs,CALObs=CALObs,Neff=Neff,BioIndex=BioIndex,BioSig=BioSig,Nproj=Nproj,Fproj=0)
parameters <- list(dummy=0,LogRbar=LogRbar,LogNinit=LogNinit,LogFullF=LogFullF,Eps=Eps)
map<-list(dummy=factor(NA))

# Note hessian=T so we get the Hessian matrix
model <- MakeADFun(f2, parameters, silent=T,map=map,hessian=T)

# test code - for checking for minimization
xx <- model$fn(model$env$last.par)
rpt <-model$report()
cat(rpt$obj_fun,rpt$LikeCatch,rpt$LikeCAL,rpt$LikeBio,rpt$Penal,"\n")

fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
best <- model$env$last.par.best
#print(model$fn(best))
rpt <-model$report()
cat(rpt$obj_fun,rpt$LikeCatch,rpt$LikeCAL,rpt$LikeBio,rpt$Penal,"\n")
rep <- sdreport(model)
#print(rep)
print("Done basic estimation but with a projection period of 20 years")

################################################################################
# Now for MCMC sampling and projections
################################################################################

# Run the MCMC
if (DoMCMCEst==T)
{
  print("Doing MCMC sampling followed by SHinySTAN")
  
  # Now consider mcmc sampling (stan)
  mcmcout <- tmbstan(obj=model,iter=50000,warmup = 50000/10, chains=3,init = list(best,best,best), 
                     seed = 1916, thin = 100)
  
  #mcmcout <- tmbstan(obj=model,iter=500,warmup = 500/10, chains=3,init = list(best,best,best), 
  #                   seed = 1916, thin = 100)
  
  ## Key information from run. Including the two recommended
  ## convergence diagnostics:
  print(summary(mcmcout))
  save(mcmcout,file="Ex5.RData")
  
  ## Interactive tools (must close out browser to regain console)
  launch_shinystan(mcmcout)
}

# Now do the projections
# ==========================
if (DoMCMCProj==T)
 {
  load(file="Ex5.RData")
  launch_shinystan(mcmcout)
  
  # Could save this as a data.frame
  post2 <- rstan::extract(mcmcout)
  Nsim <- length(post2$lp__)
  print(attributes(post2))

  # Note to change Fproj to 0.1 use the following code
  # data$Fproj <- 0.1
  # model <- MakeADFun(data, parameters, DLL="Ex4",silent=T,map=map,hessian=T)
  
  # Set up for graphs
  par(mfrow=c(2,2))
  
  Npar <- 1+ length(post2$LogNinit[1,])+length(post2$LogFullF[1,])+length(post2$Eps[1,])
  
  # Extract the results from tmbstan into the correct format
  post1 <- matrix(0,nrow=Nsim,ncol=Npar)
  post1[,1] <- post2$LogRbar
  Ioff <- 1
  for (Ilen in 1:length(post2$LogNinit[1,])) post1[,Ilen+Ioff] <- post2$LogNinit[,Ilen] 
  Ioff <- Ioff + length(post2$LogNinit[1,])
  for (Iyear in 1:length(post2$LogFullF[1,])) post1[,Iyear+Ioff] <- post2$LogFullF[,Iyear]
  Ioff <- Ioff + length(post2$LogFullF[1,])
  for (Iyear in 1:length(post2$Eps[1,])) post1[,Iyear+Ioff] <- post2$Eps[,Iyear]
  
  # Good to check (compare with the sdreport)
  print("# std report check")
  for (II in 1:Npar)  
    cat(II,best[II],mean(post1[,II]),sd(post1[,II]),"\n") 
  print("")
  
  # Extract  biomass
  Biomass <- matrix(0,nrow=Nsim,ncol=Nyear+1)
  for (Isim in 1:Nsim)
  { xx <- model$fn(post1[Isim,]); Biomass[Isim,] <- model$report()$BioPred[1:(Nyear+1)]; }  
  quant <- matrix(0,nrow=5,ncol=Nyear+1)
  for (Iyear in 1:(Nyear+1))
    quant[,Iyear] <- quantile(Biomass[,Iyear],probs=c(0.05,0.25,0.5,0.75,0.95))  
  Years <- 1:(Nyear+1)
  
  # Plot of biomass posterior
  ymax <- max(quant)
  plot(Years,quant[3,],xlab="Year",ylab="Biomass",ylim=c(0,ymax),type='l')
  xx <- c(Years,rev(Years))
  yy <- c(quant[1,],rev(quant[5,]))
  polygon(xx,yy,col="gray10")
  xx <- c(Years,rev(Years))
  yy <- c(quant[2,],rev(quant[4,]))
  polygon(xx,yy,col="gray90")
  lines(Years,quant[3,],lwd=3,lty=3)
  
  # Projection    
  # ==========
  # Now thin (otherwise this is super slow)
  Use <- seq(from=1,to=Nsim,by=10)
  #Use <- seq(from=1,to=Nsim,by=1)
  post1 <- post1[Use,]
  Nsim <- length(Use)
  
  Biomass <- matrix(0,nrow=Nsim,ncol=Nyear+20)
  
  # profile of probability of recovery with time  
  Probs <- rep(0,101)
  fs <- rep(0,101)
  for (IP in 1:101)
  {  
    print(IP)
    fs[IP] <- (IP-1)*0.015
    data$Fproj <- fs[IP]
    model <- MakeADFun(f2, parameters, silent=T,map=map,hessian=T)
    for (Isim in 1:Nsim) Biomass[Isim,] <- model$report(post1[Isim,])$BioPred[1:(Nyear+20)]; 
    Probs[IP] <- sum(Biomass[,Nyear+20]>1000)/Nsim
  }  
  plot(fs,Probs,xlab="Fishing mortality",ylab="Probability > 1000t",ylim=c(0,1),type='l',lwd=2,lty=1) 
  
  # Bisection to find desire recovery probability
  Fmin <- 0
  Fmax <- 1
  for (II in 1:20)
  { 
    Fbar <- (Fmin+Fmax)/2.0
    data$Fproj <- Fbar
    model <- MakeADFun(f2,parameters, silent=T,map=map,hessian=T)
    for (Isim in 1:Nsim) Biomass[Isim,] <- model$report(post1[Isim,])$BioPred[1:(Nyear+20)]; 
    Probs <- sum(Biomass[,Nyear+20]>1000)/Nsim
    cat(Fbar,Probs,"\n")
    if (Probs > 0.5) Fmin <- Fbar else Fmax <- Fbar;
  }  
  cat(Fbar,Probs,"\n")
  points(Fbar,Probs,pch=16)
  abline(h=0.5,lty=2)
  abline(v=Fbar,lty=2)
  
}   
