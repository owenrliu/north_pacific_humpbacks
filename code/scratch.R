library(tidyverse)
library(here)

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

test2 <- crossing(eps=seq(0.01,0.99,length.out=50),S=c(0.8,0.9,0.96,0.99)) %>% 
  mutate(Sout=foo2(eps,S,sig=6))
glimpse(test2)

test2 %>% ggplot(aes(eps,Sout,color=ordered(S)))+geom_line(linewidth=1.25)+theme_classic()+labs(color="Baseline S",x=expression(epsilon))


#######################
# Run all this to seed the DoRun function
Code="B2F1";SensCase="BC";StochSopt=1;StrayBase=0;IAmat=8;SA=0.96;SC=0.8;TimeLag=0;DensDepOpt=0;
SF=c(0,1,0,1,1,1); 
UseKPrior=1; Kmax=60000;
YrSDevs=2000;
envOpt="randomS";
SigmaDevS=6;SigmaDevF=0.01;WithMirror=1;Yr1=1970;Yr2=2023;
AddCV=T;MixWeights=c(1,1);CatchSer="B";AllPlots=F;DoBoot=F;
ByCatchFile="BycatchActual_2024_04_24.csv";
DataFileName= here('data',"Hump.dat");
FullDiag=T;
WghtTotal=1;Idirichlet=1;MaxN=100;seed=19101;
SetNew=0;Init=NULL;Nage=11;

subdir = "B2F2 FAbase"

# =================================================================================================================================
# Checking and Debugging
# =================================================================================================================================

# Gradient checking
fit$par
model$gr(fit$par)

#get a list of parameters
parms <- model$env$parList()
lastpar <- model$env$last.par
lastparbest <- model$env$last.par.best

# Checking the Hessian
sqrt(diag(solve(model$he())))
# Errors prob indicate non-invertible Hessian
# Error in solve.default(a, b) : 
#   system is computationally singular: reciprocal condition number = 8.08535e-31

# Run the model function with the parameters from a fit
parNew <- model$env$parList() |> map(unname)
parNew <- split(fit$par,names(fit$par)) |> map(unname)
test <- f(parNew,dat)

### Pulling a result
obj <- read_rds(here('Diags','B2F2 FAbase','B2F2BC.rds'))
obj <- read_rds(here('Diags','B2F1 FAenvIndex','B2F1BC.rds'))

obj <- read_rds(here('Diags','B2F1 FArS','B2F1BC.rds'))
sdr <- obj$sdfixed
rept <- obj$report
# all ADREPORTed things
adr <- obj$sdreport

# survival
survOut <- obj$report$SurvOutF |> as.numeric()
SFdevOut <- adr[grepl("SFdevYr",row.names(adr)),1]
survOutdf <- tibble(SurvOutF=survOut,
                    SFdevYr=SFdevOut,
                    zoneID=rep(as.character(obj$input$FeedNames),length(obj$input$Yr1:obj$input$Yr2)),
                    year=rep(obj$input$Yr1:obj$input$Yr2,each=length(obj$input$FeedNames)))
survOutdf |> 
  filter(year>1999) |> 
  ggplot(aes(SFdevYr))+
  geom_histogram()
survOutdf |> 
  ggplot(aes(SFdevYr,SurvOutF,color=zoneID))+
  geom_point()
survOutdf |> 
  ggplot(aes(year,SurvOutF,color=zoneID))+
  geom_line()+
  facet_wrap(~zoneID)
survOutdf |> 
  ggplot(aes(year,SFdevYr,color=zoneID))+
  geom_line()+
  facet_wrap(~zoneID)

# omegas
omegasst <- adr[grepl("omega_sst",row.names(adr)),]
omegachl <- adr[grepl("omega_chl",row.names(adr)),]


# Developing output writing code:
# Abbrev = SensCase;
# rept <- model$report()
# best <- model$env$last.par.best
# stdreport <-sdreport(model)
# rep <- summary(stdreport)
# rep2<- summary(stdreport, "fixed", p.value = TRUE)
# WriteOut(Code,Abbrev=SensCase,Yr1=Yr1,Yr2=Yr2,BreedNames=BreedNames,FeedNames=FeedNames,
#          rept=rept,rep=rep,rep2=rep2,StockDef=StockDef,data=dat,subdir = subdir)

### Comparing results

obj1 <- read_rds(here('Diags','B2F2 base','B2F2BC.rds'))
obj2 <- read_rds(here('Diags','B2F2 index','B2F2BC.rds'))
obj3 <- read_rds(here('Diags','B2F2 varK','B2F2BC.rds'))
obj4 <- read_rds(here('Diags','B2F2 direct','B2F2BC.rds'))
objlist <- list(obj1,obj4,obj2,obj3)
scen.names <- c("Base","Direct","Index","VarK")
objlist <- list(obj1,obj4,obj2,obj3)
scen.names <- c("Base","Direct","Index","VarK")
plot_compare_abundance(objlist,scen.names)
