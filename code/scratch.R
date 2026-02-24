library(tidyverse)
library(here)

#######################
# Run all this to seed the DoRun function
Code="B2F1";SensCase="BC";StrayBase=0;
IAmat=8;SA=0.96;SC=0.8;TimeLag=0;DensDepOpt=0;
SF=c(0,1,1,1,1,1); WithMirror=0;
UseKPrior=1; Kmax=60000;
YrSDevs=2000;
rvars="Kdev";
envOpt="env-K";envVars=c("sst","mld");splineK=7;envlag=0;
SigmaDevS=6;SigmaDevF=0.01;Yr1=1970;Yr2=2023;
AddCV=T;MixWeights=c(1,1);CatchSer="B";AllPlots=F;DoBoot=F;
ByCatchFile="BycatchActual_2024_04_24.csv";
DataFileName= here('data',"Hump.dat");
FullDiag=T;
WghtTotal=1;Idirichlet=1;MaxN=100;seed=19101;
Init=NULL;Nage=11;DoBayes=F;

subdir = "B2F1 FAenvDirect"

# =================================================================================================================================
# Checking and Debugging
# =================================================================================================================================

# A fitted model obj
model <- read_rds(here('Diags','final','env-K test','B2F1BC_TMB.rds'))
fit <- nlminb(model$par, model$fn, model$gr,verbose=T)

# Check gradients
grad <- model$gr(fit$par)
names(fit$par)[abs(grad) > 0.01]  # Should be empty at convergence

#get a list of parameters
parms <- model$env$parList()
lastpar <- model$env$last.par
lastparbest <- model$env$last.par.best

# Check the Hessian with the most recent fit
hess <- model$he(fit$par)
# (if a model with random effects, you could do this over the fixed effects only)
marg_hess <- optimHess(fit$par, model$fn, model$gr)
hess <- marg_hess

eigen(hess)$values  # Look for values near zero or negative

cat("Any NaN in Hessian:", any(is.nan(hess)), "\n")
cat("Any Inf in Hessian:", any(is.infinite(hess)), "\n")
eigs <- eigen(hess)$values
cat("Min eigenvalue:", min(eigs), "\n")
cat("Near-zero eigenvalues:", sum(abs(eigs) < 1e-6), "\n")

# 2. Find which parameters correspond to near-zero eigenvalues
eig_decomp <- eigen(hess)
# The smallest eigenvalue indicates the "flattest" direction in the likelihood
# Find the eigenvector for the smallest eigenvalue
problematic_evec <- eig_decomp$vectors[, which.min(abs(eigs))]
# Find the parameters with large loadings onto the problematic eigenvector
names(fit$par)[abs(problematic_evec) > 0.1]
# Parameters with large loadings on the flat direction are your problem

# For all of the smallest eigenvalues
nearzero_eig <- which(abs(eigs) < 1e-6)
parnames_temp <-paste(names(fit$par),1:length(fit$par))
map(nearzero_eig,\(x){
  problematic_evec <- eig_decomp$vectors[, x]
  parnames_temp[abs(problematic_evec) > 0.1]
}) |> unlist() |> unique() |> sort()


# Run the objective function with the parameters from a fit
list2env(parameters,envir=environment())
parNew <- model$env$parList() |> map(unname)
parNew <- split(fit$par,names(fit$par)) |> map(unname)
test <- f(parNew,dat)

# Fix some parameter and re-fit
mapNew <- model$env$map
mapNew$logK <- factor(rep(NA,length(obj$input$parameters$logK)))
parmsNew <- obj$input$parameters
parmsNew$logK <- c(8.449,10.81,9.23,9.99,8.49)
modelNew <- MakeADFun(cmb(f,obj$input), parmsNew, map=mapNew,DLL="Hump",silent=F)
fitFix <- nlminb(modelNew$par, modelNew$fn, modelNew$gr,verbose=T)

# =================================================================================================================================
# Load Results/Outputs
# =================================================================================================================================

obj <- read_rds(here('Diags','final','rS no RUS_WAL','B2F1BC.rds'))
obj <- read_rds(here('Diags','final','env-survival no RUS_WAL','B2F1BC.rds'))
obj <- read_rds(here('Diags','final','env-K test','B2F1BC.rds'))
obj <- read_rds(here('Diags','final','ddOnly','B2F1BC.rds'))
obj <- read_rds(here('Diags','final','env-survival generic env test','B2F1BC.rds'))
###

# Reports and outputs
# all fixed effects
sdr <- obj$sdfixed
# all REPORTed things
rept <- obj$report
# all ADREPORTed things
adr <- obj$sdreport

# Amounts of data we're actually fitting to
Surveys <- read_csv(here('Diags','SurveyUse',paste0(obj$Code,obj$SensCase,".csv")),show_col_types = F)
dt <- Surveys |> count(Area,Component)

### Comparing results
obj1 <- read_rds(here("Diags","B2F1 direct/no RUS_WAL","B2F1BC.rds"))
obj2 <- read_rds(here('Diags','B2F1 FAenvDirect','B2F1BC.rds'))
obj3 <- read_rds(here('Diags','B2F1 FAvarK','B2F1BC.rds'))
obj4 <- read_rds(here('Diags','B2F1 FArS','B2F1BC.rds'))
obj5 <- read_rds(here('Diags','B2F1 FAdd','B2F1BC.rds'))

objlist <- list(obj4,obj5,obj2)
scen.names <- c("Random","Density-dependence","Environment")

plot_compare_abundance(objlist,scen.names)
plot_compare_abundance(objlist,scen.names,opt='feed')
plot_compare_abundance(objlist,scen.names,opt='breed')

plot_compare_FEs(objlist,scen.names,effect = c("logK"))
plot_compare_FEs(objlist,scen.names,effect = c("rval"))
plot_compare_S(objlist,scen.names)

#-------------------EQUATION TESTING-------------------------#

# trying to understand shape of survival functions
foo1 <- function(S) log(1/S-1) # eq. B.6b in Appendix X

foo1(0.96)
test <- tibble(S=seq(0.3,0.99,length.out=50)) %>% 
  mutate(phi=foo1(S))
glimpse(test)

test %>% ggplot(aes(S,phi))+geom_line()

foo2 <- function(eps,S,sig){ # eq. B.6a in Appendix X
  1/(1+exp(foo1(S)+eps*sig))
}

test2 <- crossing(eps=seq(0.01,0.99,length.out=50),
                  S=c(0.8,0.96),
                  sig=seq(0,6,by=2)) %>% 
  mutate(Sout=foo2(eps,S,sig))
glimpse(test2)

test2 %>% 
  ggplot(aes(eps,Sout,color=ordered(S)))+
  geom_line(linewidth=1.25)+
  facet_wrap(~sig,labeller = label_bquote(sigma==.(sig)))+
  theme_minimal()+
  labs(color="Baseline S",x=expression(epsilon),
       y="Survival")+
  theme(strip.text = element_text(size=12),
        panel.border = element_rect(color='black',fill=NA))


# Calculating kdevs
foor <- function(Kd,h=1.5,l=0.5,medK=10000){
  # Bounded logit
  # l+(h-l)*1/(1+exp(-Kd))
  
  # B-H
  # h/(1+exp(-Kd))
  
  # Stewart et al.
  medK*exp(Kd)
}
tibble(Kd=rnorm(100,0,1)) |> 
  mutate(out=foor(Kd)) |> 
  ggplot(aes(Kd,out))+
  geom_point()

# density dependent survival calc?
fooh <- function(base_s,beta,depl){
  # Linear
  # survival to mort0
  # m0 <- -log(base_s)
  # Sy <- exp(-m0*(1+beta*depl))
  # Sy
  # # annual survival to instantaneous hazard
  # base_haz=log(1/base_s-1)
  # haz <- exp(base_haz+beta*depl)
  # exp(-haz)
  
  # B-H
  Sy <- base_s/(1+beta*depl)
  
  # Ricker
  # Sy <- base_s*exp(-beta*depl)
  Sy
}
crossing(depl=seq(0,1.2,by=0.05),base_s=0.96,beta=seq(0,6,by=1)) |> 
  mutate(surv=fooh(base_s,beta,depl)) |> 
  ggplot(aes(depl,surv,color=ordered(beta)))+
  geom_line()+
  theme_minimal()+
  labs(color=expression(beta),x="Depletion",y="Survival")+
  theme(panel.border=element_rect(color='black',fill=NA))

# density-dependent fecundity with varying K?
foof <- function(rval=0.09,Amat=8,SA=0.96,SC=0.8,betaK=4,Depl){
  fmax <- 2*(exp(rval*(Amat+1.0))-SA*exp(rval*Amat))/(SC*SA^(Amat))
  f0 <- 2*(1.0-SA)/(SC*SA^(Amat))
  # print(c(f0,fmax))
  # Fecundity at carrying capacity
  # ParA <- (fmax-f0)/f0
  # # ParA
  # # ft <- -log(fmax)
  # # fout <- exp(-ft*(1+betaK*Depl))
  # # fout <- exp(fmax*(1+betaK*Depl))
  # 
  # Term1 <- f0*(1.0+ParA*exp(1.0-Depl));
  # print(Term1)
  # Term1 <- 0.0001+(Term1-0.0001)/(1+exp(-30.0*Term1))
  # Term1
  
  # B-H
  # fout <- f0+(fmax-f0)/(1+betaK*Depl)
  fout <- fmax/(1+betaK*Depl)
  fout
  
  # 
}
foof(SA=0.9,Depl=50)
crossing(SA=seq(0.98,0.7,by=-0.02),Depl=seq(0,1.5,by=0.1)) |> 
  mutate(fec=foof(SA=SA,Depl=Depl,betaK=1)) |> 
  ggplot(aes(Depl,fec,color=SA))+
  geom_point()+
  theme_minimal()+
  theme(panel.border=element_rect(color='black',fill=NA))

tibble(Depl=seq(0,1.4,by=0.02)) |> 
  mutate(out=foof(Depl=Depl)) |> 
  ggplot(aes(Depl,out))+
  geom_point()

foofec <- function(f0=0.129,ParA=5.891,Depl=0.01){
  
  Term1 <- f0*(1.0+ParA*(1.0-Depl));
  print(Term1)
  Term1 <- 0.0001+(Term1-0.0001)/(1+exp(-30.0*Term1))
  # print(Term1)
  # # Term1 <- fmax
  # # print(paste("Term 1:",Term1," Depletion:",Depl))
  # LogitFec <- log(1.0/Term1-1.0);
  # Term1 <- 1.0/(1.0+exp(LogitFec));
  Term1
}
tibble(depl=seq(0,2,by=0.1)) |> 
  mutate(fec=foofec(Depl=depl)) |> 
  ggplot(aes(depl,fec))+
  geom_point()

#--------------------------------------------------------------------------------
# Alluvial plot of mixing
library(ggalluvial)
library(viridis)
BreedNames <- pluck(obj,'input','BreedNames') |> as.character()
FeedNames <- pluck(obj,'input','FeedNames') |> as.character()
m <- pluck(obj,'report','Mix')
rownames(m) <- BreedNames
colnames(m) <- FeedNames
test <- as_tibble(m,rownames="Breed") |> pivot_longer(-Breed,names_to="Feed",values_to="Proportion") |> mutate(id=row_number())
ggplot(test,aes(y=Proportion,axis1=Breed,axis2=Feed))+
  geom_alluvium(aes(fill=Breed))+
  geom_stratum(width=1/12,fill='black',color='white')+
  geom_text(stat="stratum",aes(label=after_stat(stratum)),angle=90,color='white')+
  theme_classic()+
  scale_fill_manual(values=viridis_pal(option="H")(5),guide='none')+
  theme(panel.border = element_rect(fill=NA,color='black'),
        axis.text=element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_text(size=16))+
  labs(y="Breeding Ground")+
  scale_y_continuous(expand=c(0,0),sec.axis = sec_axis(~., name = "Feeding Ground"))+
  scale_x_continuous(expand=c(0,0))

#--------------------------------------------------------------------------------

# SPLINES!
# first you need to set up the data
library(mgcv)
yrdat <- data.frame(year=as.numeric(obj$input$YrSDevs:obj$input$Yr2))
spline_setup <- smoothCon(
  s(year, k=3),
  data = yrdat,
  absorb.cons = TRUE
)[[1]]
X_yr_obs <- spline_setup$X
S_yr <- spline_setup$S[[1]]
n_spline_coef <- ncol(X_yr_obs)

# Then in RTMB, you'd have something like 
yday_effect_bin <- sum(X_yr_obs[i, ] * beta_yday_bin)

# Where
for (j in 1:n_spline_coef) {
  beta_yday_bin[j] %~% dnorm(0, sigma_spline_bin, log = TRUE)
}

mcycle = MASS::mcycle
smoothCon(
  s(times, k=3),
  data = mcycle,
  absorb.cons = TRUE
)[[1]]
