library(tidyverse)
library(here)
# library(circlize)
library(ggalluvial)
library(cowplot)
library(viridis)
library(ggarrow)
library(colorspace)
theme_set(theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA)))

# As of October 7, 2025, models are saved as .rds files as lists
# First list element is the input data from MakeDataScenario()
# Second list element is the rept from TMB
# Third and fourth are the summary sdreport for all effects and fixed effects only, respectively

# Plot the fixed effect estimates, except the environmental coefficients
plot_fixed_p <- function(obj){
  # zn <- pluck(obj,'input','FeedNames') |> as.character()
  # numz <- length(zn)
  
  dat <- pluck(obj,"sdfixed") |> 
    as_tibble(rownames="param") |> 
    # don't include env. coefficients (we'll look at those separately)
    filter(!(grepl("envParam",param))) |> 
    filter(!(grepl("SFdev",param))) |> 
    dplyr::select(1:3) |> 
    set_names(c("param","mean","se")) |> 
    mutate(pn=row_number()) |> 
    unite(name,pn,param,sep=".",remove = F) |>
    mutate(lower=mean-1.96*se,upper=mean+1.96*se)
  # dat1 <- dat |> filter(mean>2|mean< -2)
  # dat2 <- dat |> filter(mean<2,mean>-2)
  # p1 <- dat1 |> 
  #   ggplot(aes(fct_reorder(name,pn),mean,ymax=upper,ymin=lower))+
  #   geom_pointrange()+
  #   labs(y="Estimate",x="Parameter")+
  #   geom_hline(yintercept=0,linetype=2)+
  #   theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  # p2 <- dat2 |> 
  #   ggplot(aes(fct_reorder(name,pn),mean,ymax=upper,ymin=lower))+
  #   geom_pointrange()+
  #   labs(y="Estimate",x="Parameter")+
  #   geom_hline(yintercept=0,linetype=2)+
  #   theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  # 
  # cowplot::plot_grid(p1,p2,nrow=1)
  p <- dat |> 
    ggplot(aes(fct_reorder(name,pn),mean,ymax=upper,ymin=lower))+
    geom_pointrange()+
    labs(y="Estimate",x="Parameter")+
    geom_hline(yintercept=0,linetype=2)+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  p
}

# makes a chord diagram of best-estimate mixing proportions between feeding and breeding grounds
plot_mixing <- function(obj){
  # BreedNames <- pluck(obj,'input','BreedNames') |> as.character()
  # FeedNames <- pluck(obj,'input','FeedNames') |> as.character()
  # m <- pluck(obj,'report','Mix')
  # rownames(m) <- BreedNames
  # colnames(m) <- FeedNames
  # chordDiagram(m,col=viridis_pal()(length(m)),grid.col="grey")
  # circos.clear()
  BreedNames <- pluck(obj,'input','BreedNames') |> as.character()
  FeedNames <- pluck(obj,'input','FeedNames') |> as.character()
  m <- pluck(obj,'report','Mix')
  rownames(m) <- BreedNames
  colnames(m) <- FeedNames
  test <- as_tibble(m,rownames="Breed") |> pivot_longer(-Breed,names_to="Feed",values_to="Proportion") |> mutate(id=row_number())
  p <- ggplot(test,aes(y=Proportion,axis1=Breed,axis2=Feed))+
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
  p
}

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
  
  # extract environmental data
  envData <- pluck(obj,"input","envData")
  envdf <- map(1:length(envData),\(x) as.numeric(envData[[x]])) |> 
    set_names(names(envData)) |> as_tibble()
  ydevs <- pluck(obj,"input","YrSDevs")
  envdf <- envdf |> 
    mutate(year=rep(ydevs:max(yrs),each=numz)) |> 
    mutate(zone=rep(zn,length(ydevs:max(yrs))))
  
  sdat <- sdat |> 
    mutate(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |> 
    set_names(c("param","survival","sd.survival","year","zone")) |> 
    mutate(upper=survival+1.96*sd.survival,lower=survival-1.96*sd.survival) |> 
    mutate(upper=pmin(upper,1),lower=pmax(0,lower)) |>
    left_join(envdf,by=join_by(year,zone))
  nzone=length(unique(sdat$zone))
  envnames <- pluck(obj,"input","envVars")
  pl <- map(envnames,\(x){
    sdat |> 
      ggplot(aes(.data[[x]],survival))+
      geom_point(aes(color=zone))+
      # geom_text(aes(label=year,color=year),size=2)+
      # geom_arrow_chain(color='gray80',linewidth=0.25,arrow_head="head_wings")+
      labs(x=toupper(x),y="Survival")+
      scale_color_manual(values=viridis_pal()(nzone))
  })
  cowplot::plot_grid(plotlist = pl,ncol=2)
}

# Plot relationship between Sdevs and environmental data
# This is the linear relationship estimated in the model
# Option- either "B" or "F" for survival on breeding or feeding grounds
plot_Sdevs_curve <- function(obj){
  
  sdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="SFdevYr")
  
  # fitted data sdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="SFdevYr")
  yrs <- pluck(obj,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  # names
  zn <- pluck(obj,'input','FeedNames') |> as.character()
  numz <- length(zn)
  
  # extract environmental data
  envData <- pluck(obj,"input","envData")
  envdf <- map(1:length(envData),\(x) as.numeric(envData[[x]])) |> 
    set_names(names(envData)) |> as_tibble()
  ydevs <- pluck(obj,"input","YrSDevs")
  envdf <- envdf |> 
    mutate(year=rep(ydevs:max(yrs),each=numz)) |> 
    mutate(zone=rep(zn,length(ydevs:max(yrs))))
  sdat <- sdat |> 
    mutate(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |> 
    set_names(c("param","sdev","sd.sdev","year","zone")) |> 
    mutate(upper=sdev+1.96*sd.sdev,lower=sdev-1.96*sd.sdev) |>
    left_join(envdf,by=join_by(year,zone)) |> 
    filter(year>=ydevs)
  nzone=length(unique(sdat$zone))
  pl <- map(envnames,\(x){
    sdat |> 
      ggplot(aes(.data[[x]],sdev))+
      geom_point(aes(color=zone))+
      # geom_text(aes(label=year,color=year),size=2)+
      # geom_arrow_chain(color='gray80',linewidth=0.25,arrow_head="head_wings")+
      labs(x=toupper(x),y="Survival Deviate")+
      scale_color_manual(values=viridis_pal()(nzone))
  })
  cowplot::plot_grid(plotlist = pl,ncol=2)
}

# Plot coefficient estimates for environmental covariates
plot_omegas <- function(obj){
  
  zn <- pluck(obj,'input','FeedNames') |> 
    as.character()
  whichz <- which(obj$input$SF==1)
  numz <- length(whichz)
  zn <- zn[whichz]
  envnames <- pluck(obj,"input","envVars")
  
  dat <- pluck(obj,"sdfixed") |> 
    as_tibble(rownames="param") |> 
    filter(grepl('envParam',param)) |>
    dplyr::select(1:3) |> 
    set_names(c("param","mean","se")) |> 
    mutate(lower=mean-1.96*se,upper=mean+1.96*se) |> 
    mutate(zone=rep(zn,length(envnames)),
           evar=rep(envnames,each=numz))

  p <- dat |> 
    ggplot(aes(zone,mean,ymax=upper,ymin=lower,color=param))+
    geom_pointrange()+
    facet_wrap(~evar)+
    scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2),guide='none')+
    labs(x="Feeding Ground",y="Parameter Estimate")+
    geom_hline(yintercept=0,linetype=2)
  p
}

# Plot estimated K timeseries, for models with time-varying K
plot_varK <- function(obj){
  
  # fitted data
  kdat <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> 
    filter(grepl('KYr',param)) 
    # filter(!grepl('KdevYr',param))
  feedk <- obj$report$FeedK
  yrs <- pluck(obj,"input","Years")
  yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
  numyr <- length(yrs)
  
  # names
  zn <- pluck(obj,'input','FeedNames') |> as.character()
  numz <- length(zn)
  feedk <- tibble(zone=zn,feedk=feedk)
  
  # ydevs <- pluck(obj,"input","YrSDevs")
  kdat <- kdat |> 
    mutate(year=rep(yrs,each=numz)) |> 
    mutate(zone=rep(zn,numyr)) |> 
    mutate(year=as.integer(year)) |> 
    # mutate(year=rep(ydevs:max(yrs),each=numz))|> 
    # mutate(zone=rep(zn,length(ydevs:max(yrs)))) |> 
    set_names(c("param","Kyr","sd.Kyr","year","zone")) |>
    left_join(feedk,by=join_by(zone)) |> 
    mutate(CV = sd.Kyr/Kyr) |> 
    mutate(vlog = log(1+CV^2),sdlog=sqrt(vlog),Klog= log(Kyr)-vlog/2) |> 
    mutate(upper.K=Klog+1.96*sdlog,lower.K=Klog-1.96*sdlog)
  
  p <- kdat |> 
    ggplot(aes(year,Klog,ymin=lower.K,ymax=upper.K))+
    geom_ribbon(fill="lightblue")+
    geom_line()+
    labs(x="Year",y="K")+
    facet_wrap(~zone)
  p
}

# Plot index of abundance
# Right now this just plots output from one scenario. We'll write another to compare scenarios
# Options are: 
# "total" = total abundance across all populations
# "breed" = abundance in the breeding grounds
# "feed" = abundance in the feeding grounds

plot_abundance <- function(obj,opt="total",include_age0=T){

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
  Surveys <- read_csv(here('Diags','SurveyUse',paste0(obj$Code,obj$SensCase,".csv")),show_col_types = F)
  # colnames <- c("Year1","Year2","Estimate","CV","Area","Rel","Use","Add.cv","Hypothesis","Class","SensUse","Reference")
  # colnames(Surveys) <- colnames
  # survey bias for scaling
  Qvec <- obj$report$Qest
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
  
  if(opt=="total") {
    # As of Dec 2025, LogNT should be re-formed into a 3x(Nyr+1)
    d <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNT")
    d <- d |> 
      mutate(year=rep(as.integer(yrs),each=3),component=rep(c("0+","1+","Mature"),numyr)) |> 
      set_names(c("param","total.ln","sd.abun","year","component")) |> 
      mutate(total=exp(total.ln)) |> 
      mutate(upper=exp(total.ln+1.96*sd.abun),lower=exp(total.ln-1.96*sd.abun))
    surv <- survd |> filter(Area=="Total")
    d <- d |> left_join(surv,by=join_by(year,component))
    p <- d |> 
      filter(year<max(year)) |> 
      ggplot()+
      geom_ribbon(aes(year,total,ymax=upper,ymin=lower,group=component,fill=component),alpha=0.5)+
      geom_line(aes(year,total,color=component))+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=component),
                      linetype=2,size=0.5)+
      labs(x="Year",y="Total Abundance",color="Component",fill="Component")
  }
  if(opt=='breed'){
    # As of Dec 2025, LogNb is an array (3+NbreedxNyr+1), so we have to 
    # invert the way that TMB made LogNb into a long vector
    d <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNb")
    zn <- pluck(obj,'input','BreedNames') |> as.character()
    numz <- length(zn)
    
    d <- d |> 
      # make sure to check this dimensioning
      mutate(year=rep(yrs,each=numz*3)) |> 
      mutate(zone=rep(rep(zn,each=3),numyr)) |>
      mutate(component=rep(c("0+","1+","Mature"),numyr*numz)) |> 
      mutate(year=as.integer(year)) |> 
      set_names(c("param","ln.abun","sd.abun","year","zone","component")) |> 
      mutate(abun=exp(ln.abun)) |> 
      mutate(upper=exp(ln.abun+1.96*sd.abun),lower=exp(ln.abun-1.96*sd.abun))
    
    surv <- survd |> filter(Area%in%unique(d$zone))
    d <- d |> left_join(surv,by=join_by(year,component,zone==Area))
    if(include_age0){
      d <- d |> 
        filter(year<max(year)) |> 
        # for now, filter to 1+
        filter(component=="0+")
      ytitle <- "0+ Abundance"
    } else{
      d <- d |> 
        filter(year<max(year)) |> 
        # for now, filter to 1+
        filter(component=="1+")
      ytitle <- "1+ Abundance"
    }
    p <- d |>  
      ggplot()+
      geom_ribbon(aes(year,abun,ymax=upper,ymin=lower),fill='lightblue')+
      geom_line(aes(year,abun))+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=Rel),
                      linetype=2,position=position_jitter(),size=0.5)+
      facet_wrap(~zone,scales='free_y')+
      scale_color_manual(na.translate=FALSE,values=c("black","orange"))+
      labs(x="Year",y=ytitle,color="Observation\nType")
  }
  if(opt=='feed'){
    d <- pluck(obj,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNf")
    zn <- pluck(obj,'input','FeedNames') |> as.character()
    numz <- length(zn)
    
    d <- d |> 
      # make sure to check this dimensioning
      mutate(year=rep(yrs,each=numz*3)) |> 
      mutate(zone=rep(rep(zn,each=3),numyr)) |>
      mutate(component=rep(c("0+","1+","Mature"),numyr*numz)) |> 
      mutate(year=as.integer(year)) |> 
      set_names(c("param","ln.abun","sd.abun","year","zone","component")) |> 
      mutate(abun=exp(ln.abun)) |> 
      mutate(upper=exp(ln.abun+1.96*sd.abun),lower=exp(ln.abun-1.96*sd.abun))
    surv <- survd |> filter(Area%in%unique(d$zone))
    d <- d |> left_join(surv,by=join_by(year,component,zone==Area))
    if(include_age0){
      d <- d |> 
        filter(year<max(year)) |> 
        # for now, filter to 1+
        filter(component=="0+")
      ytitle <- "0+ Abundance"
    } else{
      d <- d |> 
        filter(year<max(year)) |> 
        # for now, filter to 1+
        filter(component=="1+")
      ytitle <- "1+ Abundance"
    }
    
    p <- d |>
      ggplot()+
      geom_ribbon(aes(year,abun,ymax=upper,ymin=lower),fill='lightblue')+
      geom_line(aes(year,abun))+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=Rel),
                      linetype=2,position = position_jitter(),size=0.5)+
      facet_wrap(~zone,scales='free_y')+
      scale_color_manual(na.translate=FALSE,values=c("black","orange"))+
      labs(x="Year",y=ytitle,color="Observation\nType")
  }
  p
}

## Mortality difference between realized and baseline survival
# options:
# opt= breed or feed: organize plots and sums by breeding/feeding stock
# type= "raw" or "cumulative": plot results as an annual value vs. cumulative across the time series
plot_mortdiff <- function(obj,opt="breed",type='raw'){
  # Get Mortdiff- the difference in mortality if survival was equal to baseline survival
  md <- obj$report$MortDiff*-1 # switch sign
  bn <-obj$input$BreedNames |> as.character()
  fn <- obj$input$FeedNames |> as.character()
  yrs <- obj$input$Yr1:obj$input$Yr2
  mddf <- tibble(md=as.numeric(md),
                 breed=rep(rep(bn,length(fn)),length(yrs)),
                 feed=rep(rep(fn,each=length(bn)),length(yrs)),
                 year=rep(yrs,each=length(fn)*length(bn))) |> 
    mutate(across(where(is.list),as.character)) |> 
    unite(herd,breed,feed,remove = F)
  # mddf |> 
  #   ggplot(aes(year,md,color=herd))+
  #   geom_line()+
  #   theme_minimal()
  
  # mddf |> 
  #   ggplot(aes(year,md))+
  #   geom_line()+
  #   facet_grid(feed~breed)+
  #   theme_minimal()
  if(type=='raw'){
    if(opt=="breed"){
      p <- mddf |> 
        filter(year>1999) |> 
        ggplot(aes(year,md,fill=breed))+
        geom_col()+
        scale_fill_manual(values=viridis_pal(option="G")(length(bn)))+
        facet_wrap(~feed)+
        labs(y="Mortality (ind.) relative to baseline",x="Year",fill="Breeding\nStock")
    }
    if(opt=="feed"){
      p <- mddf |> 
        filter(year>1999) |> 
        ggplot(aes(year,md,fill=feed))+
        geom_col()+
        scale_fill_manual(values=viridis_pal(option="G")(length(fn)))+
        facet_wrap(~breed)+
        labs(y="Mortality (ind.) relative to baseline",x="Year",fill="Feeding\nStock")
    }
  }
  
  if(type=="cumulative"){
    if(opt=="breed"){
      # cumulative by breeding
      cume_mort_df <- mddf |> 
        group_by(herd) |> 
        arrange(year) |> 
        mutate(cume_md=cumsum(md)) |> 
        filter(last(cume_md)!=0) |> 
        group_by(breed,year) |> 
        summarise(cume_md=sum(cume_md),.groups = 'drop')
      p <- cume_mort_df |> 
        filter(year>1999) |> 
        ggplot(aes(year,cume_md,color=breed,fill=breed))+
        # geom_line()+
        # geom_area()+
        geom_col()+
        scale_fill_manual(values=viridis_pal(option="G")(length(bn)))+
        scale_color_manual(values=viridis_pal(option="G")(length(bn)))+
        labs(x="Year",y="Cumulative Mortality Difference",
             fill="Breeding\nStock",color="Breeding\nStock")
    }
    if(opt=="feed"){
      # cumulative by breeding
      cume_mort_df <- mddf |> 
        group_by(herd) |> 
        arrange(year) |> 
        mutate(cume_md=cumsum(md)) |> 
        filter(last(cume_md)!=0) |> 
        group_by(feed,year) |> 
        summarise(cume_md=sum(cume_md),.groups = 'drop')
      p <- cume_mort_df |> 
        filter(year>1999) |> 
        ggplot(aes(year,cume_md,color=feed,fill=feed))+
        # geom_line()+
        # geom_area()+
        geom_col()+
        scale_fill_manual(values=viridis_pal(option="G")(length(fn)))+
        scale_color_manual(values=viridis_pal(option="G")(length(fn)))+
        labs(x="Year",y="Cumulative Mortality Difference",
             fill="Feeding\nStock",color="Feeding\nStock") 
    }
  }
  p
}

## Human (bycatch and ship strikes) vs. natural mortality
# options:
# opt= breed, feed, or total: organize plots and sums by breeding/feeding stock
# type= "raw", "cumulative", or "rate": plot results as an annual value vs. cumulative across the time series
plot_compare_mort <- function(obj,opt="breed",type='raw'){
  # Natural mortality from the model
  md <- obj$report$MortDiff*-1 # switch sign
  bn <-obj$input$BreedNames |> as.character()
  fn <- obj$input$FeedNames |> as.character()
  yrs <- obj$input$Yr1:obj$input$Yr2
  mddf <- tibble(md=as.numeric(md),
                 breed=rep(rep(bn,length(fn)),length(yrs)),
                 feed=rep(rep(fn,each=length(bn)),length(yrs)),
                 year=rep(yrs,each=length(fn)*length(bn))) |> 
    mutate(across(where(is.list),as.character)) |> 
    unite(herd,breed,feed,remove = F)
  # Catch data- breeding grounds
  catchb <- tibble(catchb=as.numeric(obj$input$CatchB),
                   breed=rep(bn,each=length(yrs)),
                   year=rep(yrs,length(bn))) |> 
    group_by(breed) |> 
    arrange(year) |> 
    mutate(cume_catchb=cumsum(catchb))
  # Catch data- feeding grounds
  catchf <- tibble(catchf=as.numeric(obj$input$CatchF),
                   feed=rep(fn,each=length(yrs)),
                   year=rep(yrs,length(fn))) |> 
    group_by(feed) |> 
    arrange(year) |> 
    mutate(cume_catchf=cumsum(catchf))
  # Total abundance (for rate calculations)
  NNS <- pluck(obj,"report","NNS")
  nyrs <- length(yrs)+1
  abun <- tibble(abun=as.numeric(NNS),
                 breed=rep(rep(bn,length(fn)),nyrs),
                 feed=rep(rep(fn,each=length(bn)),nyrs),
                 year=rep(c(yrs,max(yrs)+1),each=length(fn)*length(bn)))
  babun <- abun |> 
    group_by(year,breed) |> 
    summarise(tot=sum(abun),.groups='drop')
  fabun <- abun |> 
    group_by(year,feed) |> 
    summarise(tot=sum(abun),.groups='drop')
  tabun <- abun |> 
    group_by(year) |> 
    summarise(tot=sum(abun),.groups='drop')
  
  if(opt=="total"){
    # cumulative
    cume_mort_df <- mddf |> 
      group_by(herd) |> 
      arrange(year) |> 
      mutate(cume_md=cumsum(md)) |> 
      filter(last(cume_md)!=0) |> 
      group_by(year) |> 
      summarise(md=sum(md),cume_md=sum(cume_md),.groups = 'drop')
    cb <- catchb |> ungroup() |> 
      group_by(year) |> 
      summarise(cume_catchb=sum(cume_catchb),
                catchb=sum(catchb),
                .groups='drop')
    cf <- catchf |> ungroup() |> 
      group_by(year) |> 
      summarise(cume_catchf=sum(cume_catchf),
                catchf=sum(catchf),
                .groups='drop')
    mdc <- cume_mort_df |> 
      left_join(cb,by="year") |> 
      left_join(cf,by="year") |> 
      pivot_longer(contains('cume'),names_to='cumetype',values_to='cumemort') |> 
      pivot_longer(all_of(c("md","catchb","catchf")),names_to='type',values_to='mort')
    
    if(type=="raw"){
      p <- mdc |>
        mutate(Source=case_when(
          type=="md"~"Natural",
          type=="catchb"~"Catch: Breeding",
          type=="catchf" ~ "Catch: Feeding")) |> 
        ggplot(aes(year,mort,color=Source))+
        geom_line()+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(3))+
        labs(x="Year",y="Mortality (ind.)")
    }
    if(type=="cumulative"){
      p <- mdc |>
        mutate(Source=case_when(
          cumetype=="cume_md"~"Natural",
          cumetype=="cume_catchb"~"Catch: Breeding",
          cumetype=="cume_catchf" ~ "Catch: Feeding")) |> 
        dplyr::select(year,contains('cume'),Source) |> 
        distinct() |> 
        ggplot(aes(year,cumemort,color=Source,fill=Source))+
        geom_col()+
        scale_fill_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(3))+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(3))+
        labs(x="Year",y="Cumulative Mortality (ind.)")
    }
    if(type=="rate"){
      p <- mdc |>
        mutate(Source=case_when(
          type=="md"~"Natural",
          type=="catchb"~"Catch: Breeding",
          type=="catchf" ~ "Catch: Feeding")) |> 
        left_join(tabun,by=join_by(year)) |> 
        mutate(rate=mort/tot) |> 
        ggplot(aes(year,rate,color=Source))+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(3))+
        geom_line()+
        labs(x="Year",y="Mortality Rate (per ind.)")
    }
  }
  if(opt=="breed"){
    # cumulative by breeding
    cume_mort_df <- mddf |> 
      group_by(herd) |> 
      arrange(year) |> 
      mutate(cume_md=cumsum(md)) |> 
      filter(last(cume_md)!=0) |> 
      group_by(breed,year) |> 
      summarise(md=sum(md),cume_md=sum(cume_md),.groups = 'drop')
    d <- cume_mort_df |> 
      left_join(catchb,by=join_by(breed,year))
    if(type=="raw"){
      p <- d |> 
        pivot_longer(all_of(c("md","catchb")),names_to="Source",values_to="mort") |> 
        mutate(Source=ifelse(Source=="md","Natural","Catch")) |> 
        ggplot(aes(year,mort,color=Source))+
        geom_line()+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        facet_wrap(~breed)+
        labs(x="Year",y="Mortality (ind.)")
    }
    if(type=="cumulative"){
      p <- d |> 
        pivot_longer(all_of(c("cume_md","cume_catchb")),names_to="Source",values_to="mort") |> 
        mutate(Source=ifelse(Source=="cume_md","Natural","Catch")) |> 
        ggplot(aes(year,mort,color=Source,fill=Source))+
        geom_col()+
        scale_fill_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        facet_wrap(~breed)+
        labs(x="Year",y="Cumulative Mortality (ind.)")
    }
    if(type=="rate"){
      p <- d |> 
        pivot_longer(all_of(c("md","catchb")),names_to="Source",values_to="mort") |> 
        mutate(Source=ifelse(Source=="md","Natural","Catch")) |> 
        left_join(babun,by=join_by(year,breed)) |> 
        mutate(rate=mort/tot) |> 
        ggplot(aes(year,rate,color=Source))+
        geom_line()+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        facet_wrap(~breed)+
        labs(x="Year",y="Mortality Rate (per ind.)")
    }
  }
  if(opt=="feed"){
    # cumulative by breeding
    cume_mort_df <- mddf |> 
      group_by(herd) |> 
      arrange(year) |> 
      mutate(cume_md=cumsum(md)) |> 
      filter(last(cume_md)!=0) |> 
      group_by(feed,year) |> 
      summarise(md=sum(md),cume_md=sum(cume_md),.groups = 'drop')
    d <- cume_mort_df |> 
      left_join(catchf,by=join_by(feed,year))
    if(type=="raw"){
      p <- d |> 
        pivot_longer(all_of(c("md","catchf")),names_to="Source",values_to="mort") |> 
        mutate(Source=ifelse(Source=="md","Natural","Catch")) |> 
        ggplot(aes(year,mort,color=Source))+
        geom_line()+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        facet_wrap(~feed)+
        labs(x="Year",y="Mortality (ind.)")
    }
    if(type=="cumulative"){
      p <- d |> 
        pivot_longer(all_of(c("cume_md","cume_catchf")),names_to="Source",values_to="mort") |> 
        mutate(Source=ifelse(Source=="cume_md","Natural","Catch")) |> 
        ggplot(aes(year,mort,color=Source,fill=Source))+
        geom_col()+
        scale_fill_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        facet_wrap(~feed)+
        labs(x="Year",y="Cumulative Mortality (ind.)")
    }
    if(type=="rate"){
      p <- d |> 
        pivot_longer(all_of(c("md","catchf")),names_to="Source",values_to="mort") |> 
        mutate(Source=ifelse(Source=="md","Natural","Catch")) |> 
        left_join(fabun,by=join_by(year,feed)) |> 
        mutate(rate=mort/tot) |> 
        ggplot(aes(year,rate,color=Source))+
        geom_line()+
        scale_color_manual(values=viridis_pal(begin=0.3,end=0.8,option="G")(2))+
        facet_wrap(~feed)+
        labs(x="Year",y="Mortality Rate (per ind.)")
    }
  }
  p
}

# =================================Compare Models=========================================

# Compare total abundance index across scenarios
# Instead of one output object, this function takes a list of multiple objects
# and a character vector of the names of the scenarios
# Options are: 
# "total" = total abundance across all populations
# "breed" = abundance in the breeding grounds
# "feed" = abundance in the feeding grounds

# MAKE SURE THE NAMES MATCH THE ORDER OF THE LIST
# 
# obj1 <- read_rds(here('Diags','B1F1BC.rds'))
# obj2 <- read_rds(here('Diags','B2F1BC.rds'))
# obj3 <- read_rds(here('Diags','B1F2BC.rds'))
# obj4 <- read_rds(here('Diags','B2F2BC.rds'))
# objlist <- list(obj1,obj2,obj3,obj4)
# scen.names <- c("B1F1","B2F1","B1F2","B2F2")

# obj1 <- read_rds(here('Diags','B2F2 base','B2F2BC.rds'))
# obj2 <- read_rds(here('Diags','B2F2 index','B2F2BC.rds'))
# obj3 <- read_rds(here('Diags','B2F2 varK','B2F2BC.rds'))
# obj4 <- read_rds(here('Diags','B2F2 direct','B2F2BC.rds'))
# objlist <- list(obj1,obj4,obj2,obj3)
# scen.names <- c("Base","Direct","Index","VarK")

plot_compare_abundance <- function(objlist,scen.names,opt="total"){
  yrs <- pluck(objlist,1,"input","Years")
  numyr <- length(yrs)
  
  # get the right hypothesis
  obj <- objlist[[1]]
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
  # get raw survey data
  Surveys <- read_csv(here('Diags','SurveyUse',paste0(obj$Code,obj$SensCase,".csv")),show_col_types = F)
  # colnames <- c("Year1","Year2","Estimate","CV","Area","Rel","Use","Add.cv","Hypothesis","Class","SensUse","Reference")
  # colnames(Surveys) <- colnames
  # survey bias for scaling
  Qvec <- obj$report$Qest
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
  
  if(opt=="total"){
    dall <- map2(objlist,scen.names,\(x,y){
      d <- pluck(x,"sdreport") |> as_tibble(rownames="param") |> filter(param=="LogNT")
      d <- d |> 
        mutate(year=rep(as.integer(yrs),each=3),component=rep(c("0+","1+","Mature"),length(yrs)),scenario=y) |> 
        set_names(c("param","total.ln","sd.abun","year","component","scenario")) |> 
        mutate(total=exp(total.ln)) |> 
        mutate(upper=exp(total.ln+1.96*sd.abun),lower=exp(total.ln-1.96*sd.abun))
      surv <- survd |> filter(Area=="Total",)
      d <- d |> left_join(surv,by=join_by(year,component))
    }) |> list_rbind()

    p <- dall |> 
      ggplot(aes(year,total,ymax=upper,ymin=lower,
                 fill=scenario,linetype = scenario))+
      geom_ribbon(alpha=0.2,color='gray20')+
      geom_line(linewidth=0.8)+
      geom_pointrange(aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                      linetype=2,size=0.5)+
      facet_wrap(~component,nrow=3,scales='free_y')+
      labs(x="Year",y="Total Abundance",fill="Model",linetype="Model")
  }
  
  if(opt=='breed'){ # just look at 0+ abundance for now
    zn <- pluck(obj,'input','BreedNames') |> as.character()
    numz <- length(zn)
    dall <- map2(objlist,scen.names,\(x,y){
      d <- pluck(x,"sdreport") |> as_tibble(rownames="param") |> 
        filter(param=="LogNb") |> 
        # only 1st component
        filter(row_number()%%3==1)
      d <- d |> 
        mutate(year=rep(yrs,each=numz),scenario=y) |> 
        mutate(zone=rep(zn,numyr)) |> 
        mutate(year=as.integer(year)) |> 
        set_names(c("param","ln.abun","sd.abun","year","scenario","zone")) |> 
        mutate(abun=exp(ln.abun)) |> 
        mutate(upper=exp(ln.abun+1.96*sd.abun),lower=exp(ln.abun-1.96*sd.abun))
      surv <- survd |> filter(Area%in%unique(d$zone))
      d <- d |> left_join(surv,by=join_by(year,zone==Area))
    }) |> list_rbind()

    survd_filt <- survd |> filter(Area %in%unique(dall$zone)) |> rename(zone=Area)
    p <- ggplot()+
      geom_ribbon(data=dall,alpha=0.2,aes(year,abun,ymax=upper,ymin=lower,fill=scenario))+
      geom_line(data=dall,aes(year,abun,linetype=scenario))+
      geom_pointrange(data=survd_filt,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=Rel),
                      linetype=2,position=position_jitter(),size=0.5)+
      facet_wrap(~zone,scales='free_y')+
      scale_color_manual(na.translate=FALSE,values=c("black","orange"))+
      labs(x="Year",y="Abundance",color="Observation\nType")
  }
  if(opt=='feed'){
    zn <- pluck(obj,'input','FeedNames') |> as.character()
    numz <- length(zn)
    dall <- map2(objlist,scen.names,\(x,y){
      d <- pluck(x,"sdreport") |> as_tibble(rownames="param") |> 
        filter(param=="LogNf") |> 
        # only 1st component
        filter(row_number()%%3==1)
      d <- d |> 
        mutate(year=rep(yrs,each=numz),scenario=y) |> 
        mutate(zone=rep(zn,numyr)) |> 
        mutate(year=as.integer(year)) |> 
        set_names(c("param","ln.abun","sd.abun","year","scenario","zone")) |> 
        mutate(abun=exp(ln.abun)) |> 
        mutate(upper=exp(ln.abun+1.96*sd.abun),lower=exp(ln.abun-1.96*sd.abun))
    }) |> list_rbind()
    
    survd_filt <- survd |> 
      filter(Area %in%unique(dall$zone),component=="0+") |> 
      rename(zone=Area)
    p <- ggplot()+
      geom_pointrange(data=survd_filt,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower,color=Rel),
                      linetype=2,position=position_jitter(),size=0.5)+
      geom_ribbon(data=dall,alpha=0.2,aes(year,abun,ymax=upper,ymin=lower,fill=scenario))+
      geom_line(data=dall,aes(year,abun,linetype=scenario))+
      facet_wrap(~zone,scales='free_y')+
      scale_color_manual(na.translate=FALSE,values=c("black","orange"))+
      labs(x="Year",y="Abundance",color="Observation\nType",fill="Model",linetype="Model")
  }
  p
}

# Compare fixed effects among a suite of models
plot_compare_FEs <- function(objlist,scen.names,effect="all"){
  dall <- map2(objlist,scen.names,\(x,y){
    dat <- pluck(x,"sdfixed") |> 
      as_tibble(rownames="param") |> 
      filter(!(grepl("envParam",param))) |> 
      filter(!(grepl("SFdev",param))) |> 
      dplyr::select(1:3) |> 
      set_names(c("param","mean","se")) |> 
      mutate(pn=row_number(),scenario=y) |> 
      unite(name,pn,param,sep=".",remove = F) |>
      mutate(lower=mean-1.96*se,upper=mean+1.96*se)
    if(all(grepl("all",effect))==F){
      dat <- dat |> filter(param%in%effect)
    }
    dat
  }) |> list_rbind()
  p <- dall |> 
    ggplot(aes(fct_reorder(name,pn),mean,ymax=upper,ymin=lower,color=scenario))+
    geom_pointrange(position = position_dodge(width=0.3))+
    labs(y="Estimate",x="Parameter",color="Model")+
    geom_hline(yintercept=0,linetype=2)+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  p
}

# Compare survival deviates among models
plot_compare_S <- function(objlist,scen.names){
  
  dall <- map2(objlist,scen.names,\(x,y){
    sdat <- pluck(x,"sdreport") |> as_tibble(rownames="param") |> 
      filter(param=="SurvOutF")
    yrs <- pluck(x,"input","Years")
    yrs <- yrs[-length(yrs)] #remove the last year because there's no survival calculated
    numyr <- length(yrs)
    
    # names
    zn <- pluck(x,'input','FeedNames') |> as.character()
    numz <- length(zn)
    
    sdat <- sdat |> 
      mutate(year=rep(yrs,each=numz)) |> 
      mutate(zone=rep(zn,numyr)) |> 
      mutate(year=as.integer(year)) |> 
      set_names(c("param","survival","sd.survival","year","zone")) |>
      mutate(scenario=y) |> 
      mutate(upper=survival+1.96*sd.survival,lower=survival-1.96*sd.survival) |> 
      mutate(upper=pmin(upper,1),lower=pmax(0,lower))
    sdat
  }) |> list_rbind()
  # timeseries
  p <- dall |> 
    # filter(scenario=="Index") |> 
    ggplot(aes(year,survival,ymin=lower,ymax=upper))+
    geom_ribbon(aes(fill=scenario),alpha=0.7)+
    geom_line(aes(linetype=scenario))+
    labs(x="Year",y="Survival")+
    facet_wrap(~zone)
  
  p
  
}


