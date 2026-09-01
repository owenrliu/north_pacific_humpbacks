library(RTMB)
library(tidyverse)
library(here)
# ========================================================================================

# all the functions we need
source(here('code','final','read_data.R'))
source(here('code','final',"write_outputs.R"))
source(here('code','final',"plotting_functions.R"))
source(here('code','final','bayes_helpers.R'))
# source(here('code','final',"plotting_functions_Bayes.R"))

# Fitting Function
DoRun <- function(Code,SensCase,StrayBase=0,Nage=11,IAmat=8,SA=0.96,SC=0.8,TimeLag=0,DensDepOpt=0,
                  SF=c(0,1,0,1,1,1),WithMirror=0,Yr1=1970,Yr2=2023,YrSDevs=2000, 
                  UseKPrior=1, Kmax=60000,
                  rvars=NULL, # which parameters are random
                  AddCV=T,MixWeights=c(1,1),CatchSer="B",
                  envOpt="rS", splineK=5, envlag=0, envVars=c("sst","mld"),
                  AllPlots=T,DoBoot=F,BootUse,
                  ByCatchFile="BycatchActual_2024_04_24.csv",
                  DataFileName= here('data',"Hump.dat"),
                  FullDiag=T,
                  # PlotsDir= paste0(here('plots'),"/"),
                  WghtTotal=1,Idirichlet=1,MaxN=100,seed=19101,
                  DoBayes=F,Init=NULL,subdir="",return_model_obj=F)
{
  
  cat("Doing ",paste0(Code,"-",SensCase)," using data file ",DataFileName,"\n")
  # =================================================================================================================================
  # Load Model Version
  # =================================================================================================================================
  # source the correct model code, depending on "envOpt"
  # i.e., choice about how to drive survival with environmental variables
  if(envOpt== 'rS') source(here('code','final','rS.R'))
  if(envOpt== 'env-survival') source(here('code','final','env-survival.R'))
  if(envOpt== 'ddOnly') source(here('code','final','ddOnly.R'))
  if(envOpt== 'env-K') source(here('code','final','env-K.R'))
  
  # =================================================================================================================================
  # Make Model Data
  # =================================================================================================================================
  # Pass the correct parameters to the data-making function
  run_args <- as.list(environment())
  
  # make the dataset
  dat <- do.call(MakeDataScenario,run_args)
  
  # make the data objects available to this function
  list2env(dat, envir = environment())

  # =================================================================================================================================
  # Build TMB Model
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling MakeADM")
  model <- MakeADFun(cmb(f,dat), parameters,random=rvars, map=map,DLL="Hump",silent=F)

  # =================================================================================================================================
  # Optimize
  # =================================================================================================================================
  
  if (FullDiag==T) print("Calling Fit")
  FitBest  <- 1.0e20
  fit <- nlminb(model$par, model$fn, model$gr,verbose=T)
  while(abs(FitBest-fit$objective)>0.01){
      FitBest <- fit$objective
      fit <- nlminb(fit$par, model$fn, model$gr,verbose=T)
      if (FullDiag==T) print(fit$objective)
  }

  # Save model
  fildir <- here("Diags",subdir)
  if(!dir.exists(fildir)) dir.create(fildir,recursive = T)
  readr::write_rds(model,paste0(fildir,"/",Code,SensCase,"_TMB.rds"))
  
  # =================================================================================================================================
  # Extract Results
  # =================================================================================================================================
  
  rept <- model$report()
  best <- model$env$last.par.best
  adr <-sdreport(model)
  sdr <- summary(adr)
  sdrf <- summary(adr, "fixed", p.value = TRUE)
  
  cat("final nll:",rept$neglogL,"\nTotal data likelihood:",rept$datalike,"\nTotal penalty:",rept$Penal,
      "\nSurvey likelihood:",rept$LogLike1,"\nMixing BtF likelihood:",rept$LogLike2a,"\nMixing FtB likelihood",rept$LogLike2b,"\n")
  #print(exp(LogNT))
  
  # =================================================================================================================================
  # Write and Plot
  # =================================================================================================================================
  
  WriteOut(Code,Abbrev=SensCase,rept=rept,sdr=sdr,sdrf=sdrf,data=dat,subdir = subdir)
  
  obj <- list(Code=Code,SensCase=SensCase,input=dat,report=rept,sdreport=sdr,sdfixed=sdrf,converge=fit$convergence)
  
  if(AllPlots){
    fildir <- here("plots",subdir)
    if(!dir.exists(fildir)) dir.create(fildir,recursive = T)
    
    suppressWarnings({
      plot_title <- paste0(fildir,"/",Code,"-",SensCase)
      # basic fixed effect parameters (minus environmental covariates)
      p0 <- plot_fixed_p(obj = obj)
      # abundance plots
      p1 <- plot_abundance(obj=obj,opt = "total")
      p2 <- plot_abundance(obj=obj,opt= "breed")
      p3 <- plot_abundance(obj=obj,opt="feed")
      # Mixing parameters plots
      p4 <- plot_proportions(obj,direction="B-F")
      p5 <- plot_proportions(obj,direction="F-B")
      p13 <- plot_mixing(obj)
      # Survival plots
      p6 <- plot_survival(obj,opt = "F")
      p7 <- plot_survival_curve(obj,opt="F")
      # Mortality plots
      p8 <- plot_mortdiff(obj,opt="feed",type="raw")
      p9 <- plot_compare_mort(obj,opt="feed",type="rate")
      p10 <- plot_compare_mort(obj,opt="total",type="rate")
      p11 <- plot_compare_mort(obj,opt="total",type="cumulative")

      ggsave(paste(plot_title,"fixed effects.png"),p0,height=6,width=9)
      ggsave(paste(plot_title,"total abundance.png"),p1,height=4,width=6)
      ggsave(paste(plot_title,"breeding ground abundance.png"),p2,height=6,width=9)
      ggsave(paste(plot_title,"feeding ground abundance.png"),p3,height=6,width=9)
      ggsave(paste(plot_title,"proportions breed to feed.png"),p4,height=5,width=10)
      ggsave(paste(plot_title,"proportions feed to breed.png"),p5,height=5,width=10)
      ggsave(paste(plot_title,"feeding ground survival.png"),p6,height=6,width=8)
      ggsave(paste(plot_title,"environment vs survival.png"),p7,height=6,width=8)
      ggsave(paste(plot_title,"mortality feeding grounds.png"),p8,height=6,width=8)
      ggsave(paste(plot_title,"mortality rate feeding grounds.png"),p9,height=6,width=8)
      ggsave(paste(plot_title,"mortality rate total.png"),p10,height=4,width=6)
      ggsave(paste(plot_title,"cumulative total mortality.png"),p11,height=4,width=6)
      ggsave(paste(plot_title,"mixing.png"),p13,width=10,height=8)
      if(envOpt=="env-K"){
        p12 <- plot_varK(obj)
        ggsave(paste(plot_title,"varying K.png"),p12,height=6,width=8)
      }
      })
    
  }
  
  # =================================================================================================================================
  # Full Bayesian sampling
  # =================================================================================================================================
  if (DoBayes==T){
    require(tmbstan)
    # Now consider mcmc sampling (stan)
    print("trying Bayes")
    options(mc.cores=4)
    mcmcout <- tmbstan(obj=model,iter=10000,refresh=100,warmup = 2000, chains=4,cores=1,
                       init = list(best,best,best,best),
                       seed = 1916, thin = 1)
    ## Key information from run. Including the two recommended
    ## convergence diagnostics:
    # print(summary(mcmcout))
    # create file directory if it does not exist yet
    fildir <- here("Diags",subdir)
    if(!dir.exists(fildir)) dir.create(fildir,recursive = T)
    # save Bayes out
    path_bayes <- paste0(fildir,"/",Code,SensCase,"_Bayes.rds")
    path_out <- paste0(fildir,"/",Code,SensCase,"_Out.rds")
    path_tmb <- paste0(fildir,"/",Code,SensCase,"_TMB.rds")
    
    # write model to file
    readr::write_rds(mcmcout,path_bayes)
    # write model summary file, including rnning derived quantities
    WriteBayesOut(path_bayes, path_tmb, path_out, force=TRUE)
  }

  if(return_model_obj){
    return(obj)
  }else {
    return(best)
  }
}

# ==========================================================================================
# THIS IS WHERE ALL THE ACTION HAPPENS- ACTUALLY RUN THE MODEL
###################################################################################################

# ###################################################################################################
# # Bayesian Versions
# ###################################################################################################
# # test/temp
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/env-survival test",
#             envOpt="env-survival",
#             rvars=c("epsEnv"),
#             envVars=c("chl","no3"),
#             SF=c(0,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
#             DoBayes = T,Init=NULL)
# 
# # Survival as random effect (with Bayesian sampling)
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/rS Bayes test",
#             envOpt="rS",
#             rvars=c("SFdev"),
#             SF=c(0,1,1,1,1,1), WithMirror=0,
#             DoBayes = T)
# # Density dependence only
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/ddOnly Bayes",
#             envOpt="ddOnly",
#             SF=c(1,1,1,1,1,1), WithMirror=0,
#             DoBayes = T)
# # environment survival (with Bayesian sampling)
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/env-survival Bayes",
#             envOpt="env-survival",
#             rvars=c("epsEnv"),
#             envVars=c("chl","no3"),
#             SF=c(0,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
#             DoBayes = T,Init=NULL)
# # environment K (with Bayesian sampling)
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/env-K Bayes",
#             rvars=c("Kdev"),
#             envOpt="env-K",
#             envVars=c("chl","no3"),
#             UseKPrior = 1, Kmax=60000,
#             SF=c(0,1,1,1,0,0),WithMirror = 0,AllPlots=T,DoBoot=F,
#             DoBayes = T,Init=NULL)
# 
# ###################################################################################################
# # Post-hoc Bayes plotting
# ###################################################################################################
# # Plot the outputs from the Bayesian models
# msd <- "env-survival test"
# make_all_Bayes_plots(msd,mtype="env-survival")
# 
# ###################################################################################################
# # Optimize Environmental Covariates
# ###################################################################################################
# # Try env-survival model with all possible combinations of environmental variables
# # However, after doing a correlation analysis, chl and nppv are always tightly coupled.
# # So, we don't allow them in the same models
# # For GLORYS variables first
# ev <- c("sst","chl","mld","no3","nppv")
# env_var_combs <- unlist(
#   lapply(seq_along(ev), function(k) combn(ev, k, simplify = FALSE)),
#   recursive = FALSE
# )
# removes <- map_lgl(env_var_combs,\(x)"chl"%in%x&"nppv"%in%x)
# env_var_combs <- env_var_combs[!removes]
# 
# # for now, try looping these combinations without varying other params
# lldata <- list()
# llS <- list()
# penals <- list()
# conv <- list()
# safe_run <- possibly(DoRun,otherwise="Error.")
# for(i in 1:length(env_var_combs)){
#   xx <- safe_run(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#               SensCase="BC",subdir="final/temp",
#               envOpt="env-survival",
#               rvars=c("epsEnv"),envVars = env_var_combs[[i]],
#               SF=c(0,1,1,1,1,1), WithMirror=0,
#               AllPlots=F,return_model_obj = T)
#   if(is.list(xx)){
#     lldata[[i]] <- xx$report$datalike
#     conv[[i]] <- xx$converge
#     llS[[i]] <- xx$report$LogLikeS
#     penals[[i]] <- xx$report$Penal
#   } else{
#     lldata[[i]] <- NA
#     conv[[i]] <- NA
#     llS[[i]] <- NA
#     penals[[i]] <- NA
#   }
# }
# mdl_tbl <- tibble(vars=env_var_combs,lldata=unlist(lldata),converged=unlist(conv),llS=unlist(llS),penalty=unlist(penals))
# mdl_tbl <- mdl_tbl |> 
#   mutate(mn=row_number(),
#          nev=map_dbl(vars,length)*5) |> 
#   mutate(converge= ifelse(converged==0,"Yes","No")) |> 
#   mutate(pAICd=2*nev+2*(lldata+penalty),
#          pAICs=2*nev+2*llS)
# 
# write_rds(mdl_tbl,here('Diags','final','env-survival variable testing.rds'))
# 
# p_dll <- mdl_tbl |> 
#   ggplot(aes(mn,lldata,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="Data Likelihood",color="Converged?",
#        title="NLL all data")
# p_sll <- mdl_tbl |> 
#   ggplot(aes(mn,llS,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="Survival Likelihood",color="Converged?",
#        title="NLL, survival data only")
# p_pen <- mdl_tbl |> 
#   ggplot(aes(mn,penalty,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="Total Data Penalty",color="Converged?",
#        title="Summed Data Penalty")
# # pseudo AIC- how does likelihood relate to total number of environmental parameters?
# # logic- for a given set of X environmental variables, there are X*5 parameters,
# # because there are 5 feeding grounds (6, but no covars for RUS+WAL) with estimated parameters for each
# p_AIC <- mdl_tbl |> 
#   ggplot(aes(mn,pAICd,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="pseudo AIC",color="Converged?",
#        title="2*(Nparams)-2*Summed LL")
# 
# allp <- cowplot::plot_grid(p_dll,p_sll,p_pen,p_AIC,nrow=2)
# ggsave(here('plots','final','environmental covariates','env-survival variable testing.png'),allp,w=10,h=6)
# 
# ## For MOM6 variables
# # Try env-survival model with all possible combinations of MOM6 environmental variables
# ev <- c("sst","chl","mld","no3","nppv","nlgz")
# env_var_combs <- unlist(
#   lapply(seq_along(ev), function(k) combn(ev, k, simplify = FALSE)),
#   recursive = FALSE
# )
# # for now, try looping these combinations without varying other params
# lldata <- list()
# llS <- list()
# penals <- list()
# conv <- list()
# safe_run <- possibly(DoRun,otherwise="Error.")
# for(i in 1:length(env_var_combs)){
#   xx <- safe_run(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#                  SensCase="BC",subdir="final/temp",
#                  envOpt="env-survival",
#                  rvars=c("epsEnv"),envVars = env_var_combs[[i]],
#                  SF=c(0,1,1,1,1,1), WithMirror=0,
#                  AllPlots=F,return_model_obj = T)
#   if(is.list(xx)){
#     lldata[[i]] <- xx$report$datalike
#     conv[[i]] <- xx$converge
#     llS[[i]] <- xx$report$LogLikeS
#     penals[[i]] <- xx$report$Penal
#   } else{
#     lldata[[i]] <- NA
#     conv[[i]] <- NA
#     llS[[i]] <- NA
#     penals[[i]] <- NA
#   }
# }
# mdl_tbl <- tibble(vars=env_var_combs,lldata=unlist(lldata),converged=unlist(conv),llS=unlist(llS),penalty=unlist(penals))
# mdl_tbl <- mdl_tbl |> 
#   mutate(mn=row_number(),
#          nev=map_dbl(vars,length)*5) |> 
#   mutate(converge= ifelse(converged==0,"Yes","No")) |> 
#   mutate(nlldata=lldata*-1) |> 
#   mutate(nllS=llS*-1) |> 
#   mutate(pAIC=2*nev-2*nlldata)
# write_rds(mdl_tbl,here('Diags','final','env-survival variable testing MOM6.rds'))
# 
# p_dll <- mdl_tbl |> 
#   ggplot(aes(mn,nlldata,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="Data Likelihood",color="Converged?",
#        title="Summed Log Likelihood")
# p_sll <- mdl_tbl |> 
#   ggplot(aes(mn,nllS,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="Survival Likelihood",color="Converged?",
#        title="Summed Survival Likelihood")
# p_pen <- mdl_tbl |>  
#   ggplot(aes(mn,penalty,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="Total Data Penalty",color="Converged?",
#        title="Summed Data Penalty")
# # pseudo AIC- how does likelihood relate to total number of environmental parameters?
# # logic- for a given set of X environmental variables, there are X*5 parameters,
# # because there are 5 feeding grounds with estimated parameters for each
# p_AIC <- mdl_tbl |>  
#   ggplot(aes(mn,pAIC,color=factor(converge)))+
#   geom_point()+
#   labs(x="Model Number",y="pseudo AIC",color="Converged?",
#        title="2*(Nparams)-2*Summed LL")
# 
# allp <- cowplot::plot_grid(p_dll,p_sll,p_pen,p_AIC,nrow=2)
# ggsave(here('plots','final','environmental covariates','env-survival variable testing MOM6.png'),allp,w=10,h=6)
# 
# # Try best model?
# # lowest AIC was every variables except sst
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/env-survival MOM6 test",
#             envOpt="env-survival",
#             rvars=c("epsEnv"),
#             envVars=env_var_combs[[49]],
#             SF=c(0,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
#             DoBayes = F,Init=NULL)
# obj <- read_rds(here('Diags','final','env-survival MOM6 test','B2F1BC.rds'))
# plot_omegas(obj)
# 
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/env-survival MOM6 Bayes",
#             envOpt="env-survival",
#             rvars=c("epsEnv"),
#             envVars=env_var_combs[[49]],
#             SF=c(0,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
#             DoBayes = T,Init=NULL)
# 
# xx <- DoRun(Code="B2F1",Yr1=1970, Yr2=2023,YrSDevs=2000,
#             SensCase="BC",subdir="final/env-survival MOM6 Bayes 2",
#             envOpt="env-survival",
#             rvars=c("epsEnv"),
#             envVars=c("sst", "mld", "nlgz"),
#             SF=c(0,1,1,1,1,1),WithMirror = 0,AllPlots=T,DoBoot=F,
#             DoBayes = T,Init=NULL)
