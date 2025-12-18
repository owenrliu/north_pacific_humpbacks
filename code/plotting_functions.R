library(tidyverse)
library(here)
library(circlize)
library(cowplot)
library(viridis)
library(colorspace)
theme_set(theme_classic())
# ------------------------------------------------------------------------------------

# As of October 7, 2025, models are saved as .rds files as lists
# First list element is the input data from MakeDataScenario()
# Second list element is the rept from TMB
# Third and fourth are the summary sdreport for all effects and fixed effects only, respectively

# makes a chord diagram of best-estimate mixing proportions between feeding and breeding grounds
plot_mixing <- function(obj){
  BreedNames <- pluck(obj,'input','BreedNames') |> as.character()
  FeedNames <- pluck(obj,'input','FeedNames') |> as.character()
  m <- pluck(obj,'report','Mix')
  rownames(m) <- BreedNames
  colnames(m) <- FeedNames
  chordDiagram(m,col=viridis_pal()(length(m)),grid.col="grey")
  circos.clear()
}

# ------------------------------------------------------------------------------------
# Compare the observed to model-predicted breeding to feeding
# and feeding to breeding proportions.
# Option- direction must be either "B-F" or "F-B"
plot_proportions <- function(obj,direction="B-F"){
  
  if(!(direction%in% c("B-F","F-B"))) stop("direction must be B-F or F-B")
  
  BreedNames <- pluck(obj,'input','BreedNames') |> as.character()
  FeedNames <- pluck(obj,'input','FeedNames') |> as.character()
  # extract breeding to feeding ground proportions
  d <- pluck(obj,"report","ObsMixProp")
  preds <- pluck(obj,"report","PredMix")
  
  # identifiers:
  # column 1: 1= breeding to feeding; 2= feeding to breeding
  # column 2: mixing dataset (1 or 2, mark-recapture vs. genetics)
  # column 3: which breeding ground
  # column 4: which feeding ground
  dwhich <- pluck(obj,"report","ObsMixPropI") |> 
    as_tibble(.name_repair="minimal") |> 
    set_names(c("direction","dataset","breed","feed"))
  
  dp <- dwhich |> 
    mutate(direction=ifelse(direction==1,"B-F","F-B"),
           dataset=ifelse(dataset==1,"mark-recapture","genetics"),
           breed=BreedNames[breed],
           feed=FeedNames[feed]) |>  
    #smash the names together
    unite(labBF,breed,feed,remove=F) |> 
    unite(labFB,feed,breed,remove=F)
  
  dpred <- dp |> mutate(est=preds) |> 
    mutate(dataset="model predictions")
  
  dout <- dp |> 
    mutate(est=d[,1],sd=d[,2]) |> 
    mutate(upper=est+1.96*sd,lower=est-1.96*sd) |>
    # add the preds
    bind_rows(dpred) |>
    distinct()
  
  if(direction=="B-F"){
    p <- dout |> 
      filter(direction=="B-F") |> 
      ggplot(aes(labBF,est,ymin=lower,ymax=upper,shape=dataset,color=dataset))+
      geom_pointrange(position = position_dodge(width=0.3))+
      scale_y_continuous(labels=seq(0,1,by=0.2),breaks=seq(0,1,by=0.2))+
      scale_color_manual(values=qualitative_hcl(3,palette="Dark2"))+
      scale_shape_manual(values=c(16,1,17))+
      geom_hline(yintercept=0,linetype=2)+
      labs(x="Breeding to Feeding",y="Proportion",shape="",color="")+
      theme(axis.text.x=element_text(angle=45,vjust=1,hjust=0.9))
  }
  
  if(direction=="F-B"){
    p <- dout |> 
      filter(direction=="F-B") |> 
      ggplot(aes(labFB,est,ymin=lower,ymax=upper,shape=dataset,color=dataset))+
      geom_pointrange(position = position_dodge(width=0.3))+
      scale_y_continuous(labels=seq(0,1,by=0.2),breaks=seq(0,1,by=0.2))+
      scale_color_manual(values=qualitative_hcl(3,palette="Dark2"))+
      scale_shape_manual(values=c(16,1,17))+
      geom_hline(yintercept=0,linetype=2)+
      labs(x="Feeding to Breeding",y="Proportion",shape="",color="")+
      theme(axis.text.x=element_text(angle=45,vjust=1,hjust=0.9))
  }
  p
}

# ------------------------------------------------------------------------------------
# Plot estimated survival timeseries
# Option- either "B" or "F" for survival on breeding or feeding grounds
plot_survival <- function(obj,opt="F"){
  
  # fitted data
  if(opt=="B") sdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="SurvOutB")
  if(opt=="F") sdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="SurvOutF")
  yrs <- pluck(obj,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  
  # names
  if(opt=="B") zn <- pluck(obj,'input','BreedNames') |> as.character()
  if(opt=="F") zn <- pluck(obj,'input','FeedNames') |> as.character()
  numz <- length(zn)
  
  sdat <- sdat |> 
    mutate(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |> 
    set_names(c("param","survival","sd.survival","year","zone")) |> 
    mutate(upper=survival+1.96*sd.survival,lower=survival-1.96*sd.survival) |> 
    mutate(upper=pmin(upper,1),lower=pmax(0,lower))
  
  p <- sdat |> 
    ggplot(aes(year,survival,ymin=lower,ymax=upper))+
    geom_ribbon(fill="lightblue")+
    geom_line()+
    labs(x="Year",y="Survival")+
    facet_wrap(~zone)
  p
}

# ------------------------------------------------------------------------------------
# Plot estimated survival values relative to environmental data
# Option- either "B" or "F" for survival on breeding or feeding grounds
plot_survival_curve <- function(obj,opt="F"){
  
  # fitted data
  if(opt=="B") sdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="SurvOutB")
  if(opt=="F") sdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="SurvOutF")
  yrs <- pluck(obj,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  
  # names
  if(opt=="B") zn <- pluck(obj,'input','BreedNames') |> as.character()
  if(opt=="F") zn <- pluck(obj,'input','FeedNames') |> as.character()
  numz <- length(zn)
  
  sst <- pluck(obj,"input","sst")
  chl <- pluck(obj,"input","chl")
  ydevs <- pluck(obj,"input","YrSDevs")
  envdf <- tibble(year=rep(ydevs:max(yrs),each=numz),
                  sst=as.vector(sst),
                  chl=as.vector(chl)) |> 
    mutate(zone=rep(zn,length(ydevs:max(yrs))))
  
  sdat <- sdat |> 
    mutate(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |> 
    set_names(c("param","survival","sd.survival","year","zone")) |> 
    mutate(upper=survival+1.96*sd.survival,lower=survival-1.96*sd.survival) |> 
    mutate(upper=pmin(upper,1),lower=pmax(0,lower)) |>
    left_join(envdf,by=join_by(year,zone))
  
  psst <- sdat |> 
    ggplot(aes(sst,survival,color=year,ymin=lower,ymax=upper))+
    geom_point()+
    labs(x="SST",y="Survival")+
    scale_color_viridis()+
    facet_wrap(~zone)
  pchl <- sdat |> 
    ggplot(aes(chl,survival,color=year,ymin=lower,ymax=upper))+
    geom_point()+
    labs(x="Chlorophyll",y="Survival")+
    scale_color_viridis()+
    facet_wrap(~zone)
  cowplot::plot_grid(psst,pchl,nrow=2)
}

# ------------------------------------------------------------------------------------
# Plot estimated K timeseries, for models with time-varying K
plot_survival <- function(obj,opt="F"){
  
  # fitted data
  kdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(grepl('Kdev',param))
  feedk <- obj$report$FeedK
  yrs <- pluck(obj,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  
  # names
  zn <- pluck(obj,'input','FeedNames') |> as.character()
  numz <- length(zn)
  feedk <- tibble(zone=zn,feedk=feedk)
  
  ydevs <- pluck(obj,"input","YrSDevs")
  kdat <- kdat |> 
    mutate(year=rep(ydevs:max(yrs),each=numz))|> 
    mutate(zone=rep(zn,length(ydevs:max(yrs)))) |> 
    set_names(c("param","multK","sd.multK","year","zone")) |>
    left_join(feedk,by=join_by(zone)) |> 
    mutate(multK=multK+1) |> 
    mutate(upper.multK=multK+1.96*sd.multK,lower.multK=multK-1.96*sd.multK) |> 
    mutate(varK=multK*feedk,upper=upper.multK*feedk,lower=lower.multK*feedk)
  
  p <- kdat |> 
    ggplot(aes(year,varK,ymin=lower,ymax=upper))+
    geom_ribbon(fill="lightblue")+
    geom_line()+
    labs(x="Year",y="K")+
    facet_wrap(~zone)
  p
}

# ------------------------------------------------------------------------------------
# Plot index of abundance
# Right now this just plots output from one scenario. We'll write another to compare scenarios
# Options are: 
# "total" = total abundance across all populations
# "breed" = abundance in the breeding grounds
# "feed" = abundance in the feeding grounds

plot_abundance <- function(obj,opt="total"){

  yrs <- pluck(obj,"input","Years")
  numyr <- length(yrs)
  
  # get the right hypothesis
  breednames <- obj$input$BreedNames
  feednames <- obj$input$FeedNames
  fopt <- ifelse("EAL+BER+WGOA"%in%feednames,"F2","F1")
  bopt <- ifelse("MX_ML"%in%breednames,"B2","B1")
  hyp <- paste0(bopt,fopt)
  if(fopt=="F1") hyp<-c(hyp,"F1 only")
  if(fopt=="F2") hyp<-c(hyp,"F2 only")
  if(bopt=="B1") hyp<-c(hyp,"B1 only")
  if(bopt=="B2") hyp<-c(hyp,"B2 only")
  hyp <- c(hyp,"All")
  
  # get raw survey data
  Surveys <- read.csv(here('data',"SurveyAll.duringWorkshop1.csv"),fill=T,comment.char="?",header=T,row.names=NULL)[,1:12]
  colnames <- c("Year1","Year2","Estimate","CV","Area","Rel","Use","Add.cv","Hypothesis","Class","SensUse","Reference")
  colnames(Surveys) <- colnames
  # survey bias for scaling
  Qvec <- obj$report$Qest
  Qvec[1] <- 1
  survd <- Surveys |> 
    filter(Hypothesis%in%hyp,Use=="Yes") |> 
    rename(year=Year2) |> 
    dplyr::select(year,Estimate,CV,Area,Class,Rel,Hypothesis,Class) |> 
    mutate(sd.estimate=Estimate*CV,
           q = 1/Qvec[Class],
           rescaled.est = q*Estimate,
           rescaled.upper=q*(Estimate+sd.estimate),
           rescaled.lower=q*(Estimate-sd.estimate))
  
  if(opt=="total") {
    d <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNT")
    d <- d |> 
      mutate(year=as.integer(yrs)) |> 
      set_names(c("param","total.ln","sd.abun","year")) |> 
      mutate(total=exp(total.ln)) |> 
      mutate(upper=exp(total.ln+1.96*sd.abun),lower=exp(total.ln-1.96*sd.abun))
    surv <- survd |> filter(Area=="Total")
    d <- d |> left_join(surv,by=join_by(year))
    p <- d |> 
      ggplot()+
      geom_ribbon(aes(year,total,ymax=upper,ymin=lower),fill='lightblue')+
      geom_line(aes(year,total))+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5)+
      labs(x="Year",y="Total Abundance")
  }
  if(opt=='breed'){
    d <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNb")
    zn <- pluck(obj,'input','BreedNames') |> as.character()
    numz <- length(zn)
    
    d <- d |> 
      mutate(year=rep(yrs,each=numz)) |> 
      mutate(zone=rep(zn,numyr)) |> 
      mutate(year=as.integer(year)) |> 
      set_names(c("param","ln.abun","sd.abun","year","zone")) |> 
      mutate(abun=exp(ln.abun)) |> 
      mutate(upper=exp(ln.abun+1.96*sd.abun),lower=exp(ln.abun-1.96*sd.abun))
    
    surv <- survd |> filter(Area%in%unique(d$zone))
    d <- d |> left_join(surv,by=join_by(year,zone==Area))
    
    p <- d |> 
      ggplot()+
      geom_ribbon(aes(year,abun,ymax=upper,ymin=lower),fill='lightblue')+
      geom_line(aes(year,abun))+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=Rel),
                      linetype=2,position=position_jitter(),size=0.5)+
      facet_wrap(~zone,scales='free_y')+
      scale_color_manual(na.translate=FALSE,values=c("black","orange"))+
      labs(x="Year",y="Abundance",color="Observation\nType")
  }
  if(opt=='feed'){
    d <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNf")
    zn <- pluck(obj,'input','FeedNames') |> as.character()
    numz <- length(zn)
    
    d <- d |> 
      mutate(year=rep(yrs,each=numz)) |> 
      mutate(zone=rep(zn,numyr)) |> 
      mutate(year=as.integer(year)) |> 
      set_names(c("param","ln.abun","sd.abun","year","zone")) |> 
      mutate(abun=exp(ln.abun)) |> 
      mutate(upper=exp(ln.abun+1.96*sd.abun),lower=exp(ln.abun-1.96*sd.abun))
    surv <- survd |> filter(Area%in%unique(d$zone))
    d <- d |> left_join(surv,by=join_by(year,zone==Area))
    
    p <- d |> 
      ggplot()+
      geom_ribbon(aes(year,abun,ymax=upper,ymin=lower),fill='lightblue')+
      geom_line(aes(year,abun))+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=Rel),
                      linetype=2,position = position_jitter(),size=0.5)+
      facet_wrap(~zone,scales='free_y')+
      scale_color_manual(na.translate=FALSE,values=c("black","orange"))+
      labs(x="Year",y="Abundance",color="Observation\nType")
  }
  p
}

#------------------------------------------------------------------------------------
# Compare total abundance index across scenarios
# Instead of one output object, this function takes a list of multiple objects
# and a character vector of the names of the scenarios
# MAKE SURE THE NAMES MATCH THE ORDER OF THE LIST
# 
# obj1 <- read_rds(here('Diags','B1F1BC.rds'))
# obj2 <- read_rds(here('Diags','B2F1BC.rds'))
# obj3 <- read_rds(here('Diags','B1F2BC.rds'))
# obj4 <- read_rds(here('Diags','B2F2BC.rds'))
# objlist <- list(obj1,obj2,obj3,obj4)
# scen.names <- c("B1F1","B2F1","B1F2","B2F2")

plot_compare_abundance <- function(objlist,scen.names){
  yrs <- pluck(objlist,1,"input","Years")
  numyr <- length(yrs)
  
  dall <- map2(objlist,scen.names,\(x,y){
    d <- pluck(x,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNT")
    d <- d |> 
      mutate(year=as.integer(yrs),scenario=y) |> 
      set_names(c("param","total.ln","sd.abun","year","scenario")) |> 
      mutate(total=exp(total.ln)) |> 
      mutate(upper=exp(total.ln+1.96*sd.abun),lower=exp(total.ln-1.96*sd.abun))
    d
  }) |> list_rbind()

  p <- dall |> 
    ggplot(aes(year,total,ymax=upper,ymin=lower,
               fill=scenario,linetype = scenario))+
    geom_ribbon(alpha=0.2,color='gray20')+
    geom_line(linewidth=0.8)+
    labs(x="Year",y="Total Abundance")
}
