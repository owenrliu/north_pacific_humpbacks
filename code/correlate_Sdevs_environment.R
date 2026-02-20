# Correlate random deviations with environmental variables
library(tidyverse)
library(mgcv)
library(here)
library(rstan)
library(RTMB)
# Load "data" (devs from a Bayesian model)
bayes <- read_rds(here('Diags','B2F1 FArS','B2F1BC_Bayes.rds'))
mle <- read_rds(here('Diags','B2F1 FArS','B2F1BC.rds'))
source(here('code','full age structure','base_randomSDevs.R'))

bayes_df <- summary(bayes)$summary |> as_tibble(rownames="parameter")
bayes_thin <- bayes_df |> 
  mutate(upper=`97.5%`,lower=`2.5%`) |> 
  dplyr::select(parameter,mean,upper,lower) |> 
  mutate(type="bayes")
bayes_devs <- bayes_thin |> 
  filter(grepl("SFdev",parameter)) |> 
  dplyr::select(-type)|>
  mutate(zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$YrSDevs:mle$input$Yr2)),
         year=rep(mle$input$YrSDevs:mle$input$Yr2,each=length(mle$input$FeedNames))) |> 
  dplyr::select(dev=mean,zoneID,year)

# for the survival values we need the posterior
post <- as.matrix(bayes)
post <- post[,-ncol(post)]
obj <- read_rds(here('Diags','B2F1 FArS','B2F1BC_TMB.rds'))
# let's write a function to pull this for a specific reported quantity
# the full reporting with a huge set of posteriors takes forever
# what about a 1000 draws cutoff for now
which_rows <- sample(1:nrow(post),size = 1000,replace = F)
all_derived_posteriors <- map(which_rows,\(x) obj$report(post[x,]),.progress="Calculating posteriors")
# all_derived_posteriors2 <- map(which_rows,\(x) obj$report(post[x,]),.progress="Calculating posteriors")

derive_post <- function(param) map(all_derived_posteriors,\(x)pluck(x,param))
testSurvmean <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),mean)
testSurvSD <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),sd)
testSurvlower <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),quantile,probs=0.025)
testSurvupper <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),quantile,probs=0.975)
sdat <- tibble(surv=as.numeric(testSurvmean),
                           surv_se = as.numeric(testSurvSD)/1000,
                           upper=as.numeric(testSurvupper),
                           lower=as.numeric(testSurvlower),
                           zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$Yr1:mle$input$Yr2)),
                           year=rep(mle$input$Yr1:mle$input$Yr2,each=length(mle$input$FeedNames)))

## Load environmental variables
edat <- read_rds(here("data", "processed covariates", "glorys_covars_by_zone.rds"))
edat <- map2_df(edat,names(edat),\(x,y)mutate(x,variable=y)) |> 
  mutate(zoneID = str_remove_all(zoneID, "_\\d*")) |> 
  filter(zoneID %in% unique(sdat$zoneID)) |> 
  pivot_wider(names_from=variable,values_from=c(value,sd))

cdat <- sdat |> left_join(edat) |> 
  # filter to only years since we let devs vary (2000)
  filter(year>=2000) |> 
  left_join(bayes_devs)

cdat |> 
  ggplot(aes(dev,surv))+
  geom_point()

cdat |> 
  ggplot(aes(year,surv))+
  geom_line()+facet_wrap(~zoneID)

# corplot
ecor <- cdat |> dplyr::select(surv,dev,contains("value")) |> cor()
corrplot::corrplot(ecor)

## Test
m2 <- gam(dev~0+s(year)+s(value_sst)+s(value_chl)+s(value_mld)+s(value_no3)+s(value_nppv),data = cdat)
margs <- map(names(m2$var.summary),\(x)visreg::visreg(m2,x,gg=T))
plot_grid(plotlist=margs)

## normalized
cdatn <- cdat |> 
  group_by(zoneID) |> 
  mutate(across(contains("value_"),\(x) (x-mean(x))/sd(x))) |> 
  ungroup()
m3 <- gam(dev~0+s(year)+s(value_sst)+s(value_chl)+s(value_mld)+s(value_no3)+s(value_nppv),data = cdatn)
margs <- map(names(m3$var.summary),\(x)visreg::visreg(m3,x,gg=T))
plot_grid(plotlist=margs)

# Some evidence for MLD here
# all combinations
vars_to_try <- c("s(value_sst)","s(value_chl)","s(value_mld)","s(value_no3)","s(value_nppv)")
combs <- map(1:length(vars_to_try), \(x) combn(vars_to_try,x,simplify = F)) |> 
  unlist(recursive = F)
gms <- map(combs,\(v){
  form <- paste("dev~0+",paste(v,collapse="+")) |> as.formula()
  gam(formula=form,data = cdatn)
  })

aics <- map_dbl(gms,AIC)
rsq <- map_dbl(gms,\(x)x |> summary() |> pluck("r.sq"))
out <- tibble(modelnum=1:length(gms),aic=aics,r.sq=rsq,preds=combs)
# "Best" of all these crappy models have MLD, SST, and maybe no3 or nppv
# Seems to be more "support" for the zone-normalized environmental data

# Random forest
# library(randomForest)
# rf <- randomForest(dev~value_sst+value_chl+value_mld+value_no3+value_nppv,data = cdatn,importance=T,na.action = na.roughfix)
# print(rf)
# # Make predictions on the test set
# prf <- predict(rf)
# 
# # Evaluate the model performance using RMSE
# rmse <- sqrt(mean((cdatn$dev - prf)^2))
# print(paste("Root Mean Squared Error (RMSE):", rmse))
# partialPlot(rf, pred.data = na.omit(cdatn), x.var = "value_sst", 
#             main = "Partial Dependence on SST")
