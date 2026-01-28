# Compare MLE and Bayes estimation of the baseline humpback assessment model
library(tidyverse)
require(RTMB)
library(rstan)
library(StanHeaders)
library(tmbstan)
library(shinystan)
library(here)
theme_set(theme_classic())

# Bayesian estimates of the model, with survival as a random effect
bayes <- read_rds(here('Diags','B2F1 FArS','B2F1BC_Bayes.rds'))
launch_shinystan(bayes)

# Parameter estimates as a df
bayes_df <- summary(bayes)$summary |> as_tibble(rownames="parameter")
bayes_thin <- bayes_df |> 
  mutate(upper=`97.5%`,lower=`2.5%`) |> 
  dplyr::select(parameter,mean,upper,lower) |> 
  mutate(type="bayes")

# MLE Estimates
mle <- read_rds(here('Diags','B2F1 FArS','B2F1BC.rds'))
sdr <- mle$sdfixed |> as_tibble(rownames='parameter') |> 
  set_names(c("parameter","mean","se_mean","z","pr.Z")) |> 
  mutate(upper=mean+1.96*se_mean,lower=mean-1.96*se_mean) |>
  dplyr::select(mean,upper,lower) |> 
  mutate(type="ML")
# all ADREPORTed things
adr <- mle$sdreport|> as_tibble(rownames="parameter")

# older, penalized likelihood version
mle2 <- read_rds(here('Diags','B2F1 FAbase','B2F1BC.rds'))
sdr2 <- mle2$sdfixed |> as_tibble(rownames='parameter') |> 
  set_names(c("parameter","mean","se_mean","z","pr.Z")) |> 
  mutate(upper=mean+1.96*se_mean,lower=mean-1.96*se_mean) |>
  dplyr::select(parameter,mean,upper,lower) |> 
  mutate(type="ML_penal") |> 
  filter(parameter!="SFdev")
adr2 <- mle2$sdreport|> as_tibble(rownames="parameter")

# Other parameters comparison
compare1 <- bayes_thin |> 
  filter(!grepl("SFdev",parameter),!grepl("lp_",parameter))
sdr$parameter <- compare1$parameter
sdr2$parameter <- compare1$parameter[-length(compare1$parameter)]
compare1 <- compare1 |> 
  bind_rows(sdr) |> 
  bind_rows(sdr2)

# all fixed effects
FEs_compare <- compare1 |> 
  ggplot()+
  geom_pointrange(aes(parameter,mean,ymin=lower,ymax=upper,color=type),position=position_dodge(width=0.5))+
  labs(x="Parameter",y="Estimate (95% CI)",title="Parameter Estimates\nML vs. Bayesian estimation")+
  theme_minimal()+
  theme(panel.border = element_rect(color='black',fill=NA),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1))

# Survival Devs Comparison
bayes_devs <- bayes_thin |> filter(grepl("SFdev",parameter)) |> dplyr::select(-type)
mle_devs <- adr |> filter(grepl('SFdev',parameter),!grepl("SFdevYr",parameter)) |> 
  set_names(c("parameter","mle_est","mle_se")) |> 
  mutate(mle_upper=mle_est+1.96*mle_se,mle_lower=mle_est-mle_se) |> 
  dplyr::select(-parameter)
mlep_devs <- adr2 |> filter(grepl('SFdev',parameter),!grepl("SFdevYr",parameter)) |> 
  set_names(c("parameter","mle_est","mle_se")) |> 
  mutate(mle_upper=mle_est+1.96*mle_se,mle_lower=mle_est-mle_se) |> 
  dplyr::select(-parameter)

# survOutdf <- tibble(SurvOutF=survOut,
#                     SFdevYr=SFdevOut,
#                     zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$Yr1:mle$input$Yr2)),
#                     year=rep(mle$input$Yr1:mle$input$Yr2,each=length(mle$input$FeedNames)))
compare2 <- bind_cols(bayes_devs,mle_devs,mlep_devs) |>
  mutate(zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$YrSDevs:mle$input$Yr2)),
         year=rep(mle$input$YrSDevs:mle$input$Yr2,each=length(mle$input$FeedNames)))

# Scatter against each other
Sdevs_compare_scatter <- compare2 |> 
  ggplot(aes(mean,mle_est,xmin=lower,xmax=upper,ymin=mle_lower,ymax=mle_upper,
             color=zoneID,group=zoneID))+
  geom_pointrange()+
  geom_errorbarh()+
  geom_abline(slope=1,intercept=0,linetype=2)+
  geom_hline(yintercept=0,linetype=3,linewidth=1)+
  geom_vline(xintercept=0,linetype=3,linewidth=1)+
  geom_smooth(method='lm',se=F)+
  labs(x="Bayes (95% CI)",y="MLE (95% CI)",title="Survival Deviates\nML vs. Bayesian estimation")+
  theme_minimal()+
  coord_fixed()+
  theme(panel.border = element_rect(color='black',fill=NA))

# Time series
mle_devs2 <- mle_devs |> 
  dplyr::select(-mle_se) |> 
  mutate(parameter=bayes_devs$parameter) |> 
  set_names(c("mean","upper","lower","parameter")) |> 
  mutate(type='mle')
compare3 <- bayes_devs |> 
  mutate(type="bayes") |> 
  bind_rows(mle_devs2) |> 
  mutate(zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$YrSDevs:mle$input$Yr2)*2),
         year=rep(rep(mle$input$YrSDevs:mle$input$Yr2,each=length(mle$input$FeedNames)),2))

Sdevs_compare_ts <- compare3 |> 
  ggplot(aes(year,mean,ymax=upper,ymin=lower,fill=type,linetype=type))+
  geom_ribbon(alpha=0.5)+
  geom_line()+
  geom_hline(yintercept=0,linetype=2)+
  facet_wrap(~zoneID)+
  labs(x="Year",y="SFDev (95% CI)",title="Survival Deviates\nML vs. Bayesian estimation")

# translate to actual survival values
# SFfoo1 <- function(S) log(1/S-1) # eq. B.6b from the specs
# SFfoo2 <- function(eps){ # eq. B.6a from the model specs, 
#   1/(1+exp(SFfoo1(0.96)+eps)) 
# }
# compare3 <- compare3 |> 
#   mutate(surv_mean=SFfoo2(mean),surv_lower=SFfoo2(upper),surv_upper=SFfoo2(lower))
# 
# Surv_compare_ts <- compare3 |> 
#   ggplot(aes(year,surv_mean,ymax=surv_upper,ymin=surv_lower,fill=type,linetype=type))+
#   geom_ribbon(alpha=0.5)+
#   geom_line()+
#   facet_wrap(~zoneID)+
#   # ylim(0.25,1)+
#   labs(x="Year",y="Survival (+/- 1SE)",title="Survival\nML vs. Bayesian estimation")

Sdevs_compare_scatter
Sdevs_compare_ts

# Try to get derived quantities from the Stan model
# Logic- need to run tmbobj$report() with the params from each Stan iteration
# from https://github.com/kaskr/tmbstan
# posteriors
post <- as.matrix(bayes)
# need to use a saved version of a TMB model (in this case named "model")
# last column is log_posterior, so we drop it
post <- post[,-ncol(post)]
obj <- read_rds(here('Diags','B2F1 FArS','B2F1BC_TMB.rds'))
obj_p <- read_rds(here('Diags','B2F1 FAbase','B2F1BC_TMB.rds'))
# let's write a function to pull this for a specific reported quantity
# the full reporting with a huge set of posteriors takes forever
# what about a 1000 draws cutoff for now
which_rows <- sample(1:nrow(post),size = 1000,replace = F)
all_derived_posteriors <- map(which_rows,\(x) obj$report(post[x,]),.progress="Calculating posteriors")
# all_derived_posteriors2 <- map(which_rows,\(x) obj$report(post[x,]),.progress="Calculating posteriors")

derive_post <- function(param) map(all_derived_posteriors,\(x)pluck(x,param))

test2 <- derive_post(all_derived_posteriors,"Qest")
test2 <- do.call(rbind,test2) |> 
  as_tibble(.name_repair='unique') |> 
  set_names(paste0("Q",1:7)) |> 
  mutate(iteration=row_number()) |> 
  pivot_longer(-iteration,names_to="Survey",values_to='Qest')
test2 |> 
  ggplot(aes(Qest,fill=Survey))+geom_histogram(bins=10)+facet_wrap(~Survey,scales='free')

# check survival again
testSurvmean <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),mean)
testSurvSD <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),sd)
testSurvlower <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),quantile,probs=0.025)
testSurvupper <- derive_post("SurvOutF") |> simplify2array() |> apply(c(1,2),quantile,probs=0.975)
compare_survival <- tibble(surv=as.numeric(testSurvmean),
       surv_se = as.numeric(testSurvSD)/1000,
       upper=as.numeric(testSurvupper),
       lower=as.numeric(testSurvlower),
       zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$Yr1:mle$input$Yr2)),
       year=rep(mle$input$Yr1:mle$input$Yr2,each=length(mle$input$FeedNames)))

# from the MLE
mle_surv <- adr |> filter(grepl('SurvOutF',parameter)) |> 
  set_names(c("parameter","mle_est","mle_se")) |> 
  mutate(upper=mle_est+1.96*mle_se,lower=mle_est-mle_se) |> 
  # set max at 1 for surivival
  mutate(upper=pmin(1,upper)) |> 
  dplyr::select(surv=mle_est,upper,lower) |> 
  mutate(type="mle",
         zoneID=rep(as.character(mle$input$FeedNames),length(mle$input$Yr1:mle$input$Yr2)),
         year=rep(mle$input$Yr1:mle$input$Yr2,each=length(mle$input$FeedNames)))
# from the MLE-penalized version
mle_surv2 <- adr2 |> filter(grepl('SurvOutF',parameter)) |> 
  set_names(c("parameter","mle_est","mle_se")) |> 
  mutate(upper=mle_est+1.96*mle_se,lower=mle_est-mle_se) |> 
  # set max at 1 for surivival
  mutate(upper=pmin(1,upper)) |> 
  dplyr::select(surv=mle_est,upper,lower) |> 
  mutate(type="mle-penal",
         zoneID=rep(as.character(mle2$input$FeedNames),length(mle2$input$Yr1:mle2$input$Yr2)),
         year=rep(mle2$input$Yr1:mle2$input$Yr2,each=length(mle2$input$FeedNames)))
compare_survival <- compare_survival |> 
  mutate(type='bayes') |> 
  bind_rows(mle_surv) |> 
  bind_rows(mle_surv2)

Surv_compare_ts <- compare_survival |> 
  ggplot(aes(year,surv,ymax=upper,ymin=lower,fill=type,linetype=type))+
  geom_ribbon(alpha=0.5)+
  geom_line()+
  facet_grid(type~zoneID)+
  labs(x="Year",y="Survival (95% CI)",title="Survival\nML vs. Bayesian estimation")
Surv_compare_ts

# MortDiff
md <- derive_post("MortDiff") |> simplify2array()
mdm <- apply(md,c(1,2,3),mean)
mdm2 <- apply(mdm,c(2,3),sum)
mdmout <- tibble(md=as.numeric(mdm2),
                 zoneID=rep(as.character(mle2$input$FeedNames),length(mle2$input$Yr1:mle2$input$Yr2)),
                 year=rep(mle2$input$Yr1:mle2$input$Yr2,each=length(mle2$input$FeedNames)))
mdmout |> 
  ggplot(aes(year,md))+
  geom_line()+
  facet_wrap(~zoneID)
