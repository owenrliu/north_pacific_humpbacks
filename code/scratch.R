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
Code="B2F2";SensCase="BC";StochSopt=1;StrayBase=0;IAmat=8;SA=0.96;SC=0.8;TimeLag=0;DensDepOpt=0;
SF=c(1,1,1,1,1,1); YrSDevs=1995;
SigmaDevS=6;SigmaDevF=0.01;WithMirror=1;Yr1=1970;Yr2=2023;
AddCV=F;MixWeights=c(1,1);CatchSer="B";AllPlots=F;DoBoot=F;
ByCatchFile="BycatchActual_2024_04_24.csv";
WghtTotal=1;Idirichlet=1;MaxN=100;seed=19101;
SetNew=0;Init=NULL; envOpt="varK"


#### 
# Gradient checking
fit$par
model$gr(fit$par)

### Pulling a result
obj <- read_rds(here('Diags','B2F2 index','B2F2BC.rds'))
