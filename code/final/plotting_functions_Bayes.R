# Version of plotting functions for Bayesian outputs
library(tidyverse)
library(here)
library(cowplot)
library(viridis)
library(colorspace)
library(RTMB)
library(rstan)
library(tmbstan)

theme_set(theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA)))

##-----Load----
# Load a saved model object
# bayes <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC_Bayes.rds'))
# # and its ML equivalent
# # model obj
# TMBobj <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC_TMB.rds'))
# # inputs and outputs
# TMBout <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC.rds'))
calcKPrior <- function(Kmax){
  
  Ks <- seq(from=0,to=3*Kmax,length=100)
  
  join1 <- 1/(1+exp(30*(Ks-Kmax)))
  join2 <- 1/(1+exp(-30*(Ks-Kmax)))
  
  K1 <- 2*Kmax; Y1 <- log(1/0.0001-1)
  K2 <-   Kmax; Y2 <- log(1/0.9999-1)
  Prior_slope <- (Y2-Y1)/(K2-K1)
  Prior_int <- Y2 - Prior_slope*K2
  #cat(Prior_slope,Prior_int,"\n")
  Priors <- 1.0/(1+exp(Prior_int+Prior_slope*Ks))
  KPriors <- list(Prior_int=Prior_int,Prior_slope=Prior_slope)
  
  KPriors
}

##----MCMC Diagnostics----
# quick function to check some diagnostics
# these work from rstan
# traceplot(bayes,pars=c("rval"))
# rKpairs <- pairs(bayes,pars=c("rval","logK[1]","log_alphaK","log_betaK"))

run_Bayes_diags <- function(bayesobj){
  check_hmc_diagnostics(bayesobj)
  bayes_df <- summary(bayesobj)$summary |> 
    as_tibble(rownames="parameter") |> 
    filter(parameter!="lp__")
  Rhat_count <- bayes_df |>
    filter(Rhat>1.01) |> 
    pull(parameter)
  if(length(Rhat_count)==0){
    message("No parameters had an R_hat greater than 1.01.")
  } else{
    message(cat("Parameters",Rhat_count, "had an R_hat greater than 1.01"))
  }
  min_neff <- bayes_df |> slice_min(n_eff)
  max_neff <- bayes_df |> slice_max(n_eff)
  message(cat("Minimum ESS was",min_neff$n_eff,"for parameter",min_neff$parameter))
  message(cat("Maximum ESS was",max_neff$n_eff,"for parameter",max_neff$parameter))
  message(cat("Mean ESS across all parameters was",mean(bayes_df$n_eff)))
}

##----Parameter Fits----
plot_FEs <- function(bayesobj){
  df <- summary(bayesobj)$summary |> as_tibble(rownames="parameter")
  # split into -1:1 scale and larger
  df2 <- df |> 
    filter(parameter!="lp__",
           !grepl("epsEnv",parameter),
           !grepl("Kdev",parameter)) |> 
    mutate(type=case_when(
      grepl("MixPars",parameter)~"Mixing",
      grepl("logBK",parameter) ~ "Initial Depletion",
      grepl("envParams",parameter) ~ 'Environmental',
      grepl("logK",parameter)~"Carrying Capacity",
      TRUE ~ "Other"
        ))
  p1 <- df2 |> 
    ggplot()+
    geom_linerange(aes(parameter,mean,ymin=`2.5%`,ymax=`97.5%`),linewidth=1,color='lightblue')+
    geom_pointrange(aes(parameter,mean,ymin=`25%`,ymax=`75%`),linewidth = 1.5,color='gray50')+
    labs(x="Parameter",y="Estimate (95% CI)",title="Parameter Estimates")+
    theme_minimal()+
    facet_wrap(~type,scales='free')+
    theme(panel.border = element_rect(color='black',fill=NA),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  p1
}

##----## Derived Quantities ##----
# Get derived quantities from the Stan model
# Logic- need to run tmbobj$report() with the params from each Stan iteration
# from https://github.com/kaskr/tmbstan
# Because these calculations can take a long time with large numbers of samples,
# we include an `nsamp` option that determines how many (random) samples are taken 
# from the posterior with which to calculate the derived quantities
calc_dq <- function(bayesobj,TMBobj,nsamp=min(5000,nrow(as.matrix(bayesobj)))){
  # posteriors, where rows are samples and columns are sampled parameter values
  post <- as.matrix(bayesobj)
  # last column is log_posterior, so we drop it
  post <- post[,-ncol(post)]
  which_samps <- sample(1:nrow(post),size = nsamp,replace = F)
  dqp <- map(which_samps,\(x) TMBobj$report(post[x,]),
             .progress=paste0("Calculating derived quantities for ",nsamp," samples."))
  dqp
}
# Calculate the derived quantities-- this is the time consuming step
# dq <- calc_dq(bayesobj=bayes,TMBobj=TMBobj,nsamp=1000)

# Pull one quantity of interest from the posteriors
pull_dq <- function(dqlist,qname){
  map(dqlist,\(x)pluck(x,qname))
}

# Abundance
plot_abundance <- function(dqlist,TMBout,opt="total"){
  
  yrs <- pluck(TMBout,"input","Years")
  numyr <- length(yrs)
  
  # get the right hypothesis
  breednames <- TMBout$input$BreedNames
  feednames <- TMBout$input$FeedNames
  fopt <- ifelse("EAL+BER+WGOA"%in%feednames,"F2","F1")
  bopt <- ifelse("MX_ML"%in%breednames,"B2","B1")
  hyp <- paste0(bopt,fopt)
  if(fopt=="F1") hyp<-c(hyp,"F1 only")
  if(fopt=="F2") hyp<-c(hyp,"F2 only")
  if(bopt=="B1") hyp<-c(hyp,"B1 only")
  if(bopt=="B2") hyp<-c(hyp,"B2 only")
  hyp <- c(hyp,"All")
  
  # get raw survey data
  Surveys <- read_csv(here('Diags','SurveyUse',paste0(TMBout$Code,TMBout$SensCase,".csv")),show_col_types = F)
  # survey catchability for scaling (use ML estimate for now)
  Qvec <- TMBout$report$Qest
  Qvec[1] <- 1
  survd <- Surveys |> 
    filter(Hypothesis%in%hyp,Use=="Yes") |> 
    rename(year=Year2) |> 
    dplyr::select(year,Estimate,CV,Area,Class,Rel,Hypothesis,Class,component=Component) |> 
    mutate(sd.estimate=Estimate*CV,
           q = 1/Qvec[Class],
           rescaled.est = q*Estimate,
           rescaled.upper=q*(Estimate+sd.estimate),
           rescaled.lower=q*(Estimate-sd.estimate)) |> 
    mutate(component=ifelse(component==0,"0+","1+"))
  
  NmatF <- pull_dq(dqlist,"NfitFeed") |> simplify2array()
  NmatB <- pull_dq(dqlist,"NfitBreed") |> simplify2array()
  
  if(opt=="total") {
    NmatTot <- apply(NmatB,c(1,3,4),sum)
    quants <- apply(NmatTot,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    d <- tibble(year=rep(as.integer(yrs),each=3),
                component=rep(c("0+","1+","Mature"),numyr)) |>
      mutate(low=as.numeric(quants[1,,]),
             lowmid=as.numeric(quants[2,,]),
             median=as.numeric(quants[3,,]),
             uppermid=as.numeric(quants[4,,]),
             upper=as.numeric(quants[5,,]))
    surv <- survd |> filter(Area=="Total")
    d <- d |> left_join(surv,by=join_by(year,component))
    p <-ggplot()+
      geom_ribbon(data=d,aes(x=year,ymin=low,ymax=upper,fill=component),alpha=0.5)+
      geom_ribbon(data=d,aes(x=year,ymin=lowmid,ymax=uppermid,fill=component),alpha=0.7)+
      geom_line(data=d,aes(year,median,color=component))+
      geom_pointrange(data=surv,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5)+
      scale_fill_manual(values=c("#756BB1","#238B8B","#E6781E"))+
      scale_color_manual(values=c("#756BB1","#238B8B","#E6781E"))+
      labs(x="Year",y="Total Abundance",fill="Component",color="Component")
  }
  if(opt=='breed'){
    zn <- pluck(TMBout,'input','BreedNames') |> as.character()
    numz <- length(zn)
    
    quants <- apply(NmatB,c(1,2,3),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    
    d <- tibble(year=rep(yrs,each=numz*3)) |> 
      mutate(zone=rep(rep(zn,each=3),numyr)) |>
      mutate(component=rep(c("0+","1+","Mature"),numyr*numz)) |> 
      mutate(year=as.integer(year)) |> 
      mutate(low=as.numeric(quants[1,,,]),
             lowmid=as.numeric(quants[2,,,]),
             median=as.numeric(quants[3,,,]),
             uppermid=as.numeric(quants[4,,,]),
             upper=as.numeric(quants[5,,,]))
    
    surv <- survd |> filter(Area%in%unique(d$zone)) |> 
      rename(zone=Area)
    survRel <- surv |> filter(Rel=="Rel")
    survAbs <- surv |> filter(Rel=="Abs")
    
    p <-ggplot()+
      geom_ribbon(data=d,aes(x=year,ymin=low,ymax=upper,fill=component),alpha=0.5)+
      geom_ribbon(data=d,aes(x=year,ymin=lowmid,ymax=uppermid,fill=component),alpha=0.7)+
      geom_line(data=d,aes(year,median,color=component))+
      scale_fill_manual(values=c("#756BB1","#238B8B","#E6781E"))+
      scale_color_manual(values=c("#756BB1","#238B8B","#E6781E"))+
      labs(x="Year",y="Abundance",fill="Component",color="Component")+
      geom_pointrange(data=survRel,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5,color="gray50")+
      geom_pointrange(data=survAbs,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5,color="gray20")+
      facet_wrap(~zone,scales='free_y')
  }
  if(opt=='feed'){
    zn <- pluck(TMBout,'input','FeedNames') |> as.character()
    numz <- length(zn)
    
    quants <- apply(NmatF,c(1,2,3),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    
    d <- tibble(year=rep(yrs,each=numz*3)) |> 
      mutate(zone=rep(rep(zn,each=3),numyr)) |>
      mutate(component=rep(c("0+","1+","Mature"),numyr*numz)) |> 
      mutate(year=as.integer(year)) |> 
      mutate(low=as.numeric(quants[1,,,]),
             lowmid=as.numeric(quants[2,,,]),
             median=as.numeric(quants[3,,,]),
             uppermid=as.numeric(quants[4,,,]),
             upper=as.numeric(quants[5,,,]))
    
    surv <- survd |> filter(Area%in%unique(d$zone)) |> 
      rename(zone=Area)
    survRel <- surv |> filter(Rel=="Rel")
    survAbs <- surv |> filter(Rel=="Abs")
    
    p <-ggplot()+
      geom_ribbon(data=d,aes(x=year,ymin=low,ymax=upper,fill=component),alpha=0.5)+
      geom_ribbon(data=d,aes(x=year,ymin=lowmid,ymax=uppermid,fill=component),alpha=0.7)+
      geom_line(data=d,aes(year,median,color=component))+
      scale_fill_manual(values=c("#756BB1","#238B8B","#E6781E"))+
      scale_color_manual(values=c("#756BB1","#238B8B","#E6781E"))+
      labs(x="Year",y="Abundance",fill="Component",color="Component")+
      geom_pointrange(data=survRel,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5,color="gray50")+
      geom_pointrange(data=survAbs,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5,color="gray20")+
      facet_wrap(~zone,scales='free_y')
  }
  p
}

# Survival
plot_survival <- function(dqlist,TMBout=TMBout){
  
  # posterior draws
  sdat <- pull_dq(dqlist,"SurvOutF") |> simplify2array()
  # remove non-varying survival zone
  SF <- pluck(TMBout,"input","SF")
  sdat <- sdat[which(SF==1),,]
  quants <- apply(sdat,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
  
  # fitted data
  yrs <- pluck(TMBout,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  numz <- length(zn)
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  zn <- zn[which(SF==1)]
  numz <- length(zn)
  
  sdf <- tibble(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |>
    mutate(low=as.numeric(quants[1,,]),
           lowmid=as.numeric(quants[2,,]),
           median=as.numeric(quants[3,,]),
           uppermid=as.numeric(quants[4,,]),
           upper=as.numeric(quants[5,,])) |> 
    filter(year>1994)
  
  p <- sdf |> 
    ggplot()+
    geom_ribbon(aes(year,median,ymin=low,ymax=upper),fill="#DADAEB")+
    geom_ribbon(aes(year,median,ymin=lowmid,ymax=uppermid),fill="#756BB1")+
    geom_line(aes(year,median))+
    geom_hline(yintercept=0.96,linetype=2,color="orange")+
    labs(x="Year",y="Survival",title="Survival")+
    facet_wrap(~zone)
  p
}

# Variable K
plot_varK <- function(dqlist,TMBout=TMBout){
  
  # posterior draws
  kdat <- pull_dq(dqlist,"KYr") |> 
    simplify2array() |> 
    log()
  quants <- apply(kdat,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))

  yrs <- pluck(TMBout,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  
  # feedk <- derive_post("FeedK") |> 
  #   simplify2array() |> 
  #   # median estimated K by feeding ground
  #   apply(1,median)
  feedk <- pull_dq(dqlist,"FeedK") |> 
    simplify2array() |> 
    log() |> 
    apply(1,quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
  
  # names
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  numz <- length(zn)
  feedk <- as_tibble(t(feedk)) |> 
    set_names(c("low","lowmid","median","uppermid","upper")) |> 
    mutate(zone=zn) |> 
    # pivot_longer(-zone,names_to='quant',values_to='logK') |> 
    # filter to only AK feeding grounds?
    filter(zone %in% c("EAL+BER","GOA","SEA+NBC"))
  
  # ydevs <- pluck(obj,"input","YrSDevs")
  kdat <- tibble(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |>
    mutate(low=as.numeric(quants[1,,]),
           lowmid=as.numeric(quants[2,,]),
           median=as.numeric(quants[3,,]),
           uppermid=as.numeric(quants[4,,]),
           upper=as.numeric(quants[5,,])) |> 
    filter(year>1994) |> 
    filter(zone %in% c("EAL+BER","GOA","SEA+NBC"))
  
  p <-  ggplot()+
    geom_ribbon(data=kdat,aes(year,median,ymin=low,ymax=upper),fill="#DADAEB",alpha=0.7)+
    geom_ribbon(data=kdat,aes(year,median,ymin=lowmid,ymax=uppermid),fill="#756BB1",alpha=0.7)+
    geom_line(data=kdat,aes(year,median))+
    geom_hline(data=feedk,aes(yintercept=median),linetype=2)+
    geom_hline(data=feedk,aes(yintercept=low),linetype=3)+
    geom_hline(data=feedk,aes(yintercept=upper),linetype=3)+
    labs(x="Year",y="Survival")+
    facet_wrap(~zone)
  p
  
}

# Mortality by origin herd
plot_mort <- function(dqlist,TMBout,type='raw'){
  # posterior draws
  mdat <- pull_dq(dqlist,"MortDiff") |> simplify2array()
  mdat <- mdat*-1
  quants <- apply(mdat,c(1,2,3),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
  
  # 
  yrs <- pluck(TMBout,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  bn <-pluck(TMBout,'input','BreedNames') |> as.character()
  fn <- pluck(TMBout,'input','FeedNames') |> as.character()
  
  mdf <- tibble(breed=rep(rep(bn,length(fn)),length(yrs)),
                feed=rep(rep(fn,each=length(bn)),length(yrs)),
                year=rep(yrs,each=length(fn)*length(bn))) |> 
    mutate(across(where(is.list),as.character)) |> 
    unite(herd,breed,feed,remove = F) |> 
    mutate(low=as.numeric(quants[1,,,]),
           lowmid=as.numeric(quants[2,,,]),
           median=as.numeric(quants[3,,,]),
           uppermid=as.numeric(quants[4,,,]),
           upper=as.numeric(quants[5,,,])) |> 
    # only count positive relative mortalities?
    mutate(across(low:upper,\(x)ifelse(x<0,0,x))) |> 
    # filter years?
    filter(year>1994)
  
  if(type=='raw'){
    
    p <- ggplot()+
      # geom_ribbon(data=mdf,aes(x=year,ymin=low,ymax=upper,fill=breed),alpha=0.5)+
      # geom_ribbon(data=mdf,aes(x=year,ymin=lowmid,ymax=uppermid,fill=breed),alpha=0.7)+
      # geom_line(data=mdf,aes(year,median,color=breed))+
      geom_col(data=mdf,aes(year,upper),color=NA,fill='gray80')+
      geom_col(data=mdf,aes(year,median,color=breed,fill=breed))+
      scale_fill_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
      scale_color_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
      facet_grid(breed~feed)+
      labs(y="Mortality (ind.) relative to baseline",x="Year",fill="Breeding\nStock",color="Breeding\nStock")
    # 
    # if(opt=="breed"){
    #   p <- ggplot()+
    #     # geom_ribbon(data=mdf,aes(x=year,ymin=low,ymax=upper,fill=breed),alpha=0.5)+
    #     # geom_ribbon(data=mdf,aes(x=year,ymin=lowmid,ymax=uppermid,fill=breed),alpha=0.7)+
    #     # geom_line(data=mdf,aes(year,median,color=breed))+
    #     geom_col(data=mdf,aes(year,upper),color=NA,fill='gray80')+
    #     geom_col(data=mdf,aes(year,median,color=breed,fill=breed))+
    #     scale_fill_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
    #     scale_color_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
    #     facet_grid(breed~feed)+
    #     labs(y="Mortality (ind.) relative to baseline",x="Year",fill="Breeding\nStock",color="Breeding\nStock")
    # }
    # if(opt=="feed"){
    #   p <- ggplot()+
    #     # geom_ribbon(data=mdf,aes(x=year,ymin=low,ymax=upper,fill=breed),alpha=0.5)+
    #     # geom_ribbon(data=mdf,aes(x=year,ymin=lowmid,ymax=uppermid,fill=breed),alpha=0.7)+
    #     # geom_line(data=mdf,aes(year,median,color=breed))+
    #     geom_col(data=mdf,aes(year,upper),color=NA,fill='gray80')+
    #     geom_col(data=mdf,aes(year,median,color=feed,fill=feed))+
    #     scale_fill_manual(values=viridis_pal(option="D",direction = -1)(length(fn)))+
    #     scale_color_manual(values=viridis_pal(option="D",direction=-1)(length(fn)))+
    #     facet_wrap(~breed)+
    #     labs(y="Mortality (ind.) relative to baseline",x="Year",fill="Feeding\nStock",color="Feeding\nStock")
    # }
  }
  if(type=="cumulative"){
    mdat[mdat<0] <- 0
    cmdat <- apply(mdat,c(1,2,4),cumsum) |> aperm(c(2,3,1,4))
    cmquants <- apply(cmdat,c(1,2,3),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    cmdf <- tibble(breed=rep(rep(bn,length(fn)),length(yrs)),
                  feed=rep(rep(fn,each=length(bn)),length(yrs)),
                  year=rep(yrs,each=length(fn)*length(bn))) |> 
      mutate(across(where(is.list),as.character)) |> 
      unite(herd,breed,feed,remove = F) |> 
      mutate(low=as.numeric(cmquants[1,,,]),
             lowmid=as.numeric(cmquants[2,,,]),
             median=as.numeric(cmquants[3,,,]),
             uppermid=as.numeric(cmquants[4,,,]),
             upper=as.numeric(cmquants[5,,,]))
    
    p <- ggplot()+
      geom_ribbon(data=cmdf,aes(x=year,ymin=low,ymax=upper,fill=breed),alpha=0.5)+
      geom_ribbon(data=cmdf,aes(x=year,ymin=lowmid,ymax=uppermid,fill=breed),alpha=0.7)+
      geom_line(data=cmdf,aes(year,median,color=breed))+
      # geom_col(data=mdf,aes(year,upper),color=NA,fill='gray80')+
      # geom_col(data=mdf,aes(year,median,color=breed,fill=breed))+
      scale_fill_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
      scale_color_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
      facet_grid(breed~feed)+
      labs(y="Cumulative Mortality Difference",x="Year",fill="Breeding\nStock",color="Breeding\nStock")
    # if(opt=="breed"){
    #   p <- ggplot()+
    #     geom_ribbon(data=cmdf,aes(x=year,ymin=low,ymax=upper,fill=breed),alpha=0.5)+
    #     geom_ribbon(data=cmdf,aes(x=year,ymin=lowmid,ymax=uppermid,fill=breed),alpha=0.7)+
    #     geom_line(data=cmdf,aes(year,median,color=breed))+
    #     # geom_col(data=mdf,aes(year,upper),color=NA,fill='gray80')+
    #     # geom_col(data=mdf,aes(year,median,color=breed,fill=breed))+
    #     scale_fill_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
    #     scale_color_manual(values=viridis_pal(option="D",direction=-1)(length(bn)))+
    #     facet_wrap(feed)+
    #     labs(y="Cumulative Mortality Difference",x="Year",fill="Breeding\nStock",color="Breeding\nStock")
    # }
    # if(opt=="feed"){
    #   p <- ggplot()+
    #     geom_ribbon(data=cmdf,aes(x=year,ymin=low,ymax=upper,fill=feed),alpha=0.5)+
    #     geom_ribbon(data=cmdf,aes(x=year,ymin=lowmid,ymax=uppermid,fill=feed),alpha=0.7)+
    #     geom_line(data=cmdf,aes(year,median,color=feed))+
    #     # geom_col(data=mdf,aes(year,upper),color=NA,fill='gray80')+
    #     # geom_col(data=mdf,aes(year,median,color=breed,fill=breed))+
    #     scale_fill_manual(values=viridis_pal(option="D",direction=-1)(length(fn)))+
    #     scale_color_manual(values=viridis_pal(option="D",direction=-1)(length(fn)))+
    #     facet_wrap(~breed)+
    #     labs(y="Cumulative Mortality Difference",x="Year",fill="Feeding\nStock",color="Feeding\nStock")
  # }
  }
  p
}

# Mortality in relation to bycatches
# options:
# opt= breed, feed, or total: organize plots and sums by breeding/feeding stock
# type= "raw", "cumulative", or "rate": plot results as an annual value vs. cumulative across the time series
plot_compare_mort <- function(dqlist,TMBout,opt='total',type="raw"){
  ### Basic inputs ###
  yrs <- pluck(TMBout,"input","Years")
  yrs <- yrs[-length(yrs)]
  numyr <- length(yrs)
  bn <-pluck(TMBout,'input','BreedNames') |> as.character()
  fn <- pluck(TMBout,'input','FeedNames') |> as.character()
  
  
  ### Annual Mortality by Herd ##
  mdat <- pull_dq(dqlist,"MortDiff") |> simplify2array()
  mdat <- mdat*-1
  mdat[mdat<0] <- 0
  
  ### Breeding Ground Catches ##
  catchb <- TMBout$input$CatchB |> t()
  catchbdf <- tibble(catch=as.numeric(TMBout$input$CatchB),
                     zone=rep(bn,each=length(yrs)),
                     zonetype='breed',
                     year=rep(yrs,length(bn))) |> 
    # calculate cumulative catch
    group_by(zone) |> 
    arrange(year) |> 
    mutate(cume_catch=cumsum(catch)) |> 
    ungroup()
  ### Feeding Ground Catches ##
  catchf <- TMBout$input$CatchF |> t()
  catchfdf <- tibble(catch=as.numeric(TMBout$input$CatchF),
                     zone=rep(fn,each=length(yrs)),
                     zonetype='feed',
                     year=rep(yrs,length(fn))) |> 
    # calculate cumulative catch
    group_by(zone) |> 
    arrange(year) |> 
    mutate(cume_catch=cumsum(catch)) |> 
    ungroup()
  # total catch
  totcatchbdf <- catchbdf |> 
    group_by(year) |> 
    summarise(totcatchb=sum(catch),
              totcumeb=sum(cume_catch),.groups='drop')
  totcatchfdf <- catchfdf |> 
    group_by(year) |> 
    summarise(totcatchf=sum(catch),
              totcumef=sum(cume_catch),.groups='drop')
  
  # now make plots
  if(type=="raw"){
    if(opt=="total"){
      ### Annual Mortality Entire Population ##
      tmdat <- apply(mdat,c(3,4),sum)
      tquants <- apply(tmdat,1,quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      
      tmdf <- tibble(year=yrs) |>
        mutate(low=as.numeric(tquants[1,]),
               lowmid=as.numeric(tquants[2,]),
               median=as.numeric(tquants[3,]),
               uppermid=as.numeric(tquants[4,]),
               upper=as.numeric(tquants[5,]))
      tmdf <- tmdf |> 
        left_join(totcatchbdf,by=join_by(year)) |> 
        left_join(totcatchfdf,by=join_by(year))
      
      p <- tmdf |> 
        ggplot(aes(x=year))+
        geom_line(aes(y=totcatchb,color="Catch: Breeding"))+
        geom_line(aes(y=totcatchf,color="Catch: Feeding"))+
        geom_ribbon(aes(ymin=low,ymax=upper,fill="Natural"),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid,fill="Natural"),alpha=0.7)+
        geom_line(aes(y=median,color="Natural"),color="#756BB1")+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality (ind.)",color="Source",fill="") 
    }
    if(opt=="breed"){
      # breeding grounds (sum over feeding grounds)
      mdatB <- apply(mdat,c(1,3,4),sum)
      mBquants <- apply(mdatB,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      mBdf <- tibble(zone=rep(bn,length(yrs)),
                     year=rep(yrs,each=length(bn))) |> 
        mutate(across(where(is.list),as.character)) |>
        mutate(low=as.numeric(mBquants[1,,]),
               lowmid=as.numeric(mBquants[2,,]),
               median=as.numeric(mBquants[3,,]),
               uppermid=as.numeric(mBquants[4,,]),
               upper=as.numeric(mBquants[5,,]))
      # annual catch
      mBdf <- mBdf |> 
        left_join(catchbdf,by=join_by(year,zone))
      p <- mBdf |> 
        ggplot(aes(x=year))+
        geom_line(aes(y=catch,color="Catch: Breeding"))+
        geom_ribbon(aes(ymin=low,ymax=upper,fill="Natural"),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid,fill="Natural"),alpha=0.7)+
        geom_line(aes(y=median,color="Natural"),color="#756BB1")+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality (ind.)",color="Source",fill="")+
        facet_wrap(~zone)
    }
    if(opt=="feed"){
      # feeding grounds (sum over breeding grounds)
      mdatF <- apply(mdat,c(2,3,4),sum)
      mFquants <- apply(mdatF,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      mFdf <- tibble(zone=rep(fn,length(yrs)),
                     year=rep(yrs,each=length(fn))) |> 
        mutate(across(where(is.list),as.character)) |>
        mutate(low=as.numeric(mFquants[1,,]),
               lowmid=as.numeric(mFquants[2,,]),
               median=as.numeric(mFquants[3,,]),
               uppermid=as.numeric(mFquants[4,,]),
               upper=as.numeric(mFquants[5,,]))
      # annual catch
      mFdf <- mFdf |> 
        left_join(catchfdf,by=join_by(year,zone))
      p <- mFdf |> 
        ggplot(aes(x=year))+
        geom_line(aes(y=catch,color="Catch: Feeding"))+
        geom_ribbon(aes(ymin=low,ymax=upper,fill="Natural"),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid,fill="Natural"),alpha=0.7)+
        geom_line(aes(y=median,color="Natural"),color="#756BB1")+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality (ind.)",color="Source",fill="")+
        facet_wrap(~zone)
    }
  }
  if(type=="cumulative"){
    ### Cumulative Mortality by Herd ##
    cmdat <- apply(mdat,c(1,2,4),cumsum) |> aperm(c(2,3,1,4))
    if(opt=="total"){
      ### Cumulative Mortality Entire Population ##
      tcmdat <- apply(tmdat,2,cumsum)
      tcmquants <- apply(tcmdat,1,quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      
      tcmdf <- tibble(year=yrs) |>
        mutate(low=as.numeric(tcmquants[1,]),
               lowmid=as.numeric(tcmquants[2,]),
               median=as.numeric(tcmquants[3,]),
               uppermid=as.numeric(tcmquants[4,]),
               upper=as.numeric(tcmquants[5,]))
      tcmdf <- tcmdf |> 
        left_join(totcatchbdf,by=join_by(year)) |> 
        left_join(totcatchfdf,by=join_by(year))
      
      p <- tcmdf |> 
        ggplot(aes(x=year))+
        geom_line(aes(y=totcumeb,color="Catch: Breeding"))+
        geom_line(aes(y=totcumef,color="Catch: Feeding"))+
        geom_ribbon(aes(ymin=low,ymax=upper,fill="Natural"),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid,fill="Natural"),alpha=0.7)+
        geom_line(aes(y=median,color="Natural"),color="#756BB1")+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Cumulative Mortality (ind.)",color="Source",fill="")
    }
    if(opt=="breed"){
      # cumulative mortality
      cmdatB <- apply(cmdat,c(1,3,4),sum)
      cmBquants <- apply(cmdatB,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      cmBdf <- tibble(zone=rep(bn,length(yrs)),
                      year=rep(yrs,each=length(bn))) |> 
        mutate(across(where(is.list),as.character)) |>
        mutate(low=as.numeric(cmBquants[1,,]),
               lowmid=as.numeric(cmBquants[2,,]),
               median=as.numeric(cmBquants[3,,]),
               uppermid=as.numeric(cmBquants[4,,]),
               upper=as.numeric(cmBquants[5,,]))
      cmBdf <- cmBdf |> 
        left_join(catchbdf,by=join_by(year,zone))
      
      p <- cmBdf |> 
        ggplot(aes(x=year))+
        geom_line(aes(y=cume_catch,color="Catch: Breeding"))+
        geom_ribbon(aes(ymin=low,ymax=upper,fill="Natural"),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid,fill="Natural"),alpha=0.7)+
        geom_line(aes(y=median,color="Natural"),color="#756BB1")+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality (ind.)",color="Source",fill="")+
        facet_wrap(~zone)
    }
    if(opt=="feed"){
      # cumulative mortality
      cmdatF <- apply(cmdat,c(2,3,4),sum)
      cmFquants <- apply(cmdatF,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      cmFdf <- tibble(zone=rep(fn,length(yrs)),
                      year=rep(yrs,each=length(fn))) |> 
        mutate(across(where(is.list),as.character)) |>
        mutate(low=as.numeric(cmFquants[1,,]),
               lowmid=as.numeric(cmFquants[2,,]),
               median=as.numeric(cmFquants[3,,]),
               uppermid=as.numeric(cmFquants[4,,]),
               upper=as.numeric(cmFquants[5,,]))
      cmFdf <- cmFdf |> 
        left_join(catchfdf,by=join_by(year,zone))
      p <- cmFdf |> 
        ggplot(aes(x=year))+
        geom_line(aes(y=cume_catch,color="Catch: Feeding"))+
        geom_ribbon(aes(ymin=low,ymax=upper,fill="Natural"),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid,fill="Natural"),alpha=0.7)+
        geom_line(aes(y=median,color="Natural"),color="#756BB1")+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality (ind.)",color="Source",fill="")+
        facet_wrap(~zone)
    }
  }
  if(type=='rate'){
    ### Breeding Ground Catch Rates ##
    catchrateb <- map(1:dim(NbS)[3],\(x)catchb/NbS[,,x]) |> simplify2array()
    catchratebquants <- apply(catchrateb,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    
    catchbrdf <- tibble(zone=rep(bn,each=length(yrs)),
                        year=rep(yrs,length(bn))) |> 
      mutate(low=as.numeric(catchratebquants[1,,]),
             lowmid=as.numeric(catchratebquants[2,,]),
             median=as.numeric(catchratebquants[3,,]),
             uppermid=as.numeric(catchratebquants[4,,]),
             upper=as.numeric(catchratebquants[5,,])) |> 
      mutate(source="Catch: Breeding")
    # catch rate across all breeding grounds
    catchBtot <- colSums(catchb)
    catchrateBtot <- map(1:dim(TotN)[2],\(x)catchBtot/TotN[,x]) |> simplify2array()
    totcatchratebquants <- apply(catchrateBtot,1,quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    totcatchbrdf <- tibble(year=yrs) |>
      mutate(low=as.numeric(totcatchratebquants[1,]),
             lowmid=as.numeric(totcatchratebquants[2,]),
             median=as.numeric(totcatchratebquants[3,]),
             uppermid=as.numeric(totcatchratebquants[4,]),
             upper=as.numeric(totcatchratebquants[5,])) |> 
      pivot_longer(low:upper,names_to="quant",values_to="rate") |> 
      mutate(source="CatchB")
    ### Feeding Ground Catch Rates ##
    catchratef <- map(1:dim(NfS)[3],\(x)catchf/NfS[,,x]) |> simplify2array()
    catchratefquants <- apply(catchratef,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    catchfrdf <- tibble(zone=rep(fn,each=length(yrs)),
                        year=rep(yrs,length(fn))) |> 
      mutate(low=as.numeric(catchratefquants[1,,]),
             lowmid=as.numeric(catchratefquants[2,,]),
             median=as.numeric(catchratefquants[3,,]),
             uppermid=as.numeric(catchratefquants[4,,]),
             upper=as.numeric(catchratefquants[5,,])) |> 
      mutate(source="Catch: Feeding")
    catchFtot <- colSums(catchf)
    catchrateFtot <- map(1:dim(TotN)[2],\(x)catchFtot/TotN[,x]) |> simplify2array()
    totcatchratefquants <- apply(catchrateFtot,1,quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
    totcatchfrdf <- tibble(year=yrs) |>
      mutate(low=as.numeric(totcatchratefquants[1,]),
             lowmid=as.numeric(totcatchratefquants[2,]),
             median=as.numeric(totcatchratefquants[3,]),
             uppermid=as.numeric(totcatchratefquants[4,]),
             upper=as.numeric(totcatchratefquants[5,])) |> 
      pivot_longer(low:upper,names_to="quant",values_to="rate") |> 
      mutate(source="CatchF")
    if(opt=="total"){
      # total N by herd, start of year
      NNS <- pull_dq(dqlist,"NNS") |> simplify2array()
      # Remove last year
      NNS <- NNS[,,-dim(NNS)[2],]
      # Sum 
      TotN <- apply(NNS,c(3,4),sum)
      Totrate <- tmdat/TotN
      Totrquants <- apply(Totrate,1,quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      trdf <- tibble(year=yrs) |>
        mutate(low=as.numeric(Totrquants[1,]),
               lowmid=as.numeric(Totrquants[2,]),
               median=as.numeric(Totrquants[3,]),
               uppermid=as.numeric(Totrquants[4,]),
               upper=as.numeric(Totrquants[5,])) |> 
        pivot_longer(low:upper,names_to="quant",values_to="rate") |> 
        mutate(source="Natural")
      # add catches
      trdf <- trdf |> 
        bind_rows(totcatchbrdf) |> 
        bind_rows(totcatchfrdf)|>
        pivot_wider(names_from='quant',values_from = rate)
      p <- trdf |> 
        mutate(source=case_when(
          source=="Natural" ~ "Natural",
          source=="CatchB" ~ "Catch: Breeding",
          source=="CatchF" ~ "Catch: Feeding"
        )) |> 
        ggplot(aes(x=year,fill=source))+
        geom_ribbon(aes(ymin=low,ymax=upper),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid),alpha=0.7)+
        geom_line(aes(y=median,color=source))+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality Rate",color="Source",fill="Source")
    }
    if(opt=="breed"){
      ### natural mortality by zone ##
      NbS <- pull_dq(dqlist,"NbS") |> simplify2array()
      # Remove last year
      NbS <- NbS[,-dim(NbS)[2],]
      # Mortality rate by breeding and feeding ground
      # to get mortality rate, divide mortality diff (from above) by total abundance
      rateBd <- mdatB/NbS
      rateBd[is.na(rateBd)] <-rateFd[is.na(rateFd)]<-0
      rBquants <- apply(rateBd,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      rBdf <- tibble(zone=rep(bn,length(yrs)),
                     year=rep(yrs,each=length(bn))) |> 
        mutate(across(where(is.list),as.character)) |>
        mutate(low=as.numeric(rBquants[1,,]),
               lowmid=as.numeric(rBquants[2,,]),
               median=as.numeric(rBquants[3,,]),
               uppermid=as.numeric(rBquants[4,,]),
               upper=as.numeric(rBquants[5,,])) |> 
        mutate(source="Natural")
      rBdf <- rBdf |> 
        bind_rows(catchbrdf)
      p <- rBdf |> 
        ggplot(aes(x=year,fill=source))+
        geom_ribbon(aes(ymin=low,ymax=upper),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid),alpha=0.7)+
        geom_line(aes(y=median,color=source))+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality Rate",color="Source",fill="Source")+
        facet_wrap(~zone)
    }
    if(opt=='feed'){
      ### natural mortality by zone ##
      NfS <- pull_dq(dqlist,"NfS") |> simplify2array()
      NfS <- NfS[,-dim(NfS)[2],]
      rateFd <- mdatF/NfS
      rFquants <- apply(rateFd,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
      rFdf <- tibble(zone=rep(fn,length(yrs)),
                     year=rep(yrs,each=length(fn))) |> 
        mutate(across(where(is.list),as.character)) |>
        mutate(low=as.numeric(rFquants[1,,]),
               lowmid=as.numeric(rFquants[2,,]),
               median=as.numeric(rFquants[3,,]),
               uppermid=as.numeric(rFquants[4,,]),
               upper=as.numeric(rFquants[5,,])) |> 
        mutate(source="Natural")
      rFdf <- rFdf |> 
        bind_rows(catchfrdf)
      p <- rFdf |> 
        ggplot(aes(x=year,fill=source))+
        geom_ribbon(aes(ymin=low,ymax=upper),alpha=0.5)+
        geom_ribbon(aes(ymin=lowmid,ymax=uppermid),alpha=0.7)+
        geom_line(aes(y=median,color=source))+
        scale_fill_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        scale_color_manual(values=c("Catch: Breeding"="#238B8B","Catch: Feeding"="#E6781E","Natural"="#756BB1"))+
        labs(x="Year",y="Mortality Rate",color="Source",fill="Source")+
        facet_wrap(~zone)
    } # close opt
  } # close type
  p
} # close fxn

##----## Environmental Index ##----
plot_omegas <- function(bayesobj){
  df <- summary(bayesobj)$summary |> 
    as_tibble(rownames="parameter")
  # split into -1:1 scale and larger
  # parameter names
  envVars <- pluck(TMBout,"input","envVars")
  SF <- pluck(TMBout,"input","SF")
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  zn <- zn[which(SF==1)]
  numz <- length(zn)
  
  df2 <- df |> 
    filter(grepl("envParams",parameter)) |> 
    mutate(zone=rep(zn,length(envVars)),
           variable=rep(envVars,each=numz))
  p <- df2 |> 
    ggplot()+
    geom_linerange(aes(variable,mean,ymin=`2.5%`,ymax=`97.5%`),linewidth=1,color='lightblue')+
    geom_pointrange(aes(variable,mean,ymin=`25%`,ymax=`75%`),linewidth = 1.5,color='gray50')+
    labs(x="Parameter",y="Estimate (95% CI)",title="Parameter Estimates")+
    theme_minimal()+
    geom_hline(yintercept=0,linetype=2)+
    facet_wrap(~zone)+
    theme(panel.border = element_rect(color='black',fill=NA),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  p
}
# Environmental Index
plot_envIndex <- function(dqlist,TMBout){
  # posterior draws
  edat <- pull_dq(dqlist,"env_index") |> simplify2array()
  quants <- apply(edat,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
  
  # fitted data
  yrEnd <- pluck(TMBout,"input","Years") |> last()
  yrSdevs <- pluck(TMBout,"input","YrSDevs")
  yrs <- yrSdevs:(yrEnd-1)
  numyr <- length(yrs)
  SF <- pluck(TMBout,"input","SF")
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  zn <- zn[which(SF==1)]
  numz <- length(zn)
  
  edf <- tibble(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |>
    mutate(low=as.numeric(quants[1,,]),
           lowmid=as.numeric(quants[2,,]),
           median=as.numeric(quants[3,,]),
           uppermid=as.numeric(quants[4,,]),
           upper=as.numeric(quants[5,,])) |> 
    filter(year>1994)
  
  p <- edf |> 
    ggplot()+
    geom_ribbon(aes(year,median,ymin=low,ymax=upper),fill="#2B7BBA",alpha=0.5)+
    geom_ribbon(aes(year,median,ymin=lowmid,ymax=uppermid),fill="#2B7BBA",alpha=0.7)+
    geom_line(aes(year,median))+
    geom_hline(yintercept=0,linetype=2)+
    labs(x="Year",y="Environmental Index",title="Index of Survival")+
    facet_wrap(~zone)
  p
}

# compare the environmental index against survival values
plot_envIndex_surv <- function(dqlist,TMBout){
  # fitted data
  yrEnd <- pluck(TMBout,"input","Years") |> last()
  yrSdevs <- pluck(TMBout,"input","YrSDevs")
  yrs <- yrSdevs:(yrEnd-1)
  numyr <- length(yrs)
  SF <- pluck(TMBout,"input","SF")
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  zn <- zn[which(SF==1)]
  numz <- length(zn)
  
  # posterior draws
  edat <- pull_dq(dqlist,"env_index") |> simplify2array()
  sdat <- pull_dq(dqlist,"SurvOutF") |> simplify2array()
  sdat <- sdat[which(SF==1),(dim(sdat)[2]-dim(edat)[2]+1):dim(sdat)[2],]
  
  df <- map(1:1000,\(x){
    tibble(year=rep(yrs,each=numz)) |> 
      mutate(zone=rep(zn,numyr)) |> 
      mutate(year=as.integer(year)) |>
      mutate(env_index=as.numeric(edat[,,x]),
             survival=as.numeric(sdat[,,x]))
  }) |> list_rbind()
  expected <- tibble(env_index=seq(min(df$env_index),max(df$env_index)),length.out=1000) |> 
    mutate(surv_expected = 1/(1+exp(log(1/0.96-1)+env_index)))
  p <- df |> 
    ggplot()+
    geom_point(data=df,aes(env_index,survival),size=0.2,color='gray60')+
    geom_line(data=expected,aes(env_index,surv_expected),color='#238B8B',linewidth = 1.2)+
    geom_hline(yintercept=0.96,linetype=2,color="orange")+
    labs(x="Environmental Index",y="Survival",title="Index vs. Survival")+
    facet_wrap(~zone)
  p
}

# plot the extra variance terms on survival
plot_epsEnv <- function(bayesobj,TMBout){
  
  df <- summary(bayesobj)$summary |> 
    as_tibble(rownames="parameter") |> 
    filter(grepl("epsEnv",parameter))
  
  # fitted data
  yrEnd <- pluck(TMBout,"input","Years") |> last()
  yrSdevs <- pluck(TMBout,"input","YrSDevs")
  yrs <- yrSdevs:(yrEnd-1)
  numyr <- length(yrs)
  SF <- pluck(TMBout,"input","SF")
  zn <- pluck(TMBout,'input','FeedNames') |> as.character()
  zn <- zn[which(SF==1)]
  numz <- length(zn)
  
  edf <- tibble(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |>
    mutate(low=df$`2.5%`,
           lowmid=df$`25%`,
           median=df$`50%`,
           uppermid=df$`75%`,
           upper=df$`97.5%`)
  p <- edf |> 
    ggplot()+
    geom_ribbon(aes(year,median,ymin=low,ymax=upper),fill="#238B8B",alpha=0.5)+
    geom_ribbon(aes(year,median,ymin=lowmid,ymax=uppermid),fill="#238B8B",alpha=0.7)+
    geom_line(aes(year,median))+
    geom_hline(yintercept=0,linetype=2)+
    labs(x="Year",y="Epsilon",title="Unexplained Variability")+
    facet_wrap(~zone)
  p
}
# plot a panel of 4 of these environmental index plots
plot_env_panel <- function(bayesobj,TMBout,dqlist){
  p1 <- plot_omegas(bayesobj)
  p2 <- plot_envIndex(dqlist,TMBout)
  p3 <- plot_epsEnv(bayesobj,TMBout)
  p4 <- plot_envIndex_surv(dqlist,TMBout)
  p5 <- plot_survival(dqlist,TMBout)
  cowplot::plot_grid(p1,p2,p3,p4,p5,nrow=2,rel_heights = c(1,1.2))
}

##----## Wrapper ##----
# take a file folder and make all of the above plots
# for now, it is hardcoded to pull 5000 posterior samples
library(tictoc)
make_all_Bayes_plots <- function(subdir){
  tic("Loading model objects")
  # Stanfit
  bayes <- read_rds(here('Diags','final',subdir,'B2F1BC_Bayes.rds'))
  # TMB model
  TMBobj <- read_rds(here('Diags','final',subdir,'B2F1BC_TMB.rds'))
  # inputs and outputs
  TMBout <- read_rds(here('Diags','final',subdir,'B2F1BC.rds'))
  toc()
  # calculate derived quantities
  tic("Calculating posterior samples")
  dq <- calc_dq(bayes,TMBobj,nsamp=5000)
  toc()
  # make plots:
  # fixed effects:
  tic("Making plots")
  pFEs <- plot_FEs(bayes)
  # abundance
  pabuntot <- plot_abundance(dq,TMBout,opt = 'total')
  pabunB <- plot_abundance(dq,TMBout,opt = 'breed')
  pabunF <- plot_abundance(dq,TMBout,opt = 'feed')
  # survival
  psurv <- plot_survival(dq,TMBout)
  if(grepl("env-K",subdir)){
    pvarK <- plot_varK(dq,TMBout)
  }
  # mortality
  pmortraw <- plot_mort(dq,TMBout,type='raw')
  pmortcume <- plot_mort(dq,TMBout,type='cumulative')
  # relative mortality
  p2mortrawT <- plot_compare_mort(dq,TMBout,opt="total",type = 'raw')
  p2mortrawB <- plot_compare_mort(dq,TMBout,opt="breed",type = 'raw')
  p2mortrawF <- plot_compare_mort(dq,TMBout,opt="feed",type = 'raw')
  p2mortrateT <- plot_compare_mort(dq,TMBout,opt="total",type = 'rate')
  p2mortrateB <- plot_compare_mort(dq,TMBout,opt="breed",type = 'rate')
  p2mortrateF <- plot_compare_mort(dq,TMBout,opt="feed",type = 'rate')
  p2mortcumeT <- plot_compare_mort(dq,TMBout,opt="total",type = 'cumulative')
  p2mortcumeB <- plot_compare_mort(dq,TMBout,opt="breed",type = 'cumulative')
  p2mortcumeF <- plot_compare_mort(dq,TMBout,opt="feed",type = 'cumulative')
  #environmental inex
  penv <- plot_env_panel(bayes,TMBout,dq)
  toc()
  tic("Saving plots")
  fdir <- here('plots','final',subdir)
  ggsave(paste0(fdir,"/Bayes FEs.png"),pFEs,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes total abundance.png"),pabuntot,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes breed abundance.png"),pabunB,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes feed abundance.png"),pabunF,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes feed survival.png"),psurv,w=10,h=8)
  if(grepl("env-K",subdir)){
    ggsave(paste0(fdir,"/Bayes variable K.png"),pvarK,w=10,h=8)
  }
  ggsave(paste0(fdir,"/Bayes mortality feed x breed.png"),pmortraw,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes cume mortality feed x breed.png"),pmortcume,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort total raw.png"),p2mortrawT,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort breed raw.png"),p2mortrawB,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort feed raw.png"),p2mortrawF,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort total cume.png"),p2mortcumeT,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort breed cume.png"),p2mortcumeB,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort feed cume.png"),p2mortcumeF,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort total rate.png"),p2mortrateT,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort breed rate.png"),p2mortrateB,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes relmort feed rate.png"),p2mortrateF,w=10,h=8)
  ggsave(paste0(fdir,"/Bayes envir index.png"),penv,w=12,h=8)
  toc()
  message(paste("Plots finished and saved in ",fdir))
}