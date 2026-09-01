# Make manuscript figures
library(here)
library(sf)
library(rnaturalearth)
source(here('code','final','plotting_functions_Bayes.R'))

bayes <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC_Bayes.rds'))
# # and its ML equivalent
TMBobj <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC_TMB.rds'))
# # inputs and outputs
TMBout <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC.rds'))
# calculate derived quantities
dq <- calc_dq(bayes,TMBobj,nsamp = 1000)

#### Figure 1- map, mixing estimates, and feeding ground abundances ####
# shapefile of zones
zones <- read_sf(here('data','spatial','NPhump_zones_scenarios.shp')) |> 
  filter(B2|F1)
zones_merc <- zones |> 
  # PDC Mercator
  st_transform(3832)|> unite("zoneID",Name,id,remove=F) |> 
  mutate(label=c("Asia","Central\nAmerica","Hawai'i","Asia","Mainland\nMexico",
                 "Mexican\nIslands","Asia","Asia","Southern B.C. and\nWashington",
                 "Southeast Alaska\nand Northern B.C.","Oregon and\nCalifornia",
                 "Russia and\nWestern Aleutians","Eastern Aleutians\nand Bering Sea","Gulf of Alaska")) |> 
  mutate(znum=c(1,5,2,1,3,4,1,1,10,9,11,6,7,8)) |> 
  mutate(type=ifelse(Feeding,"Feeding","Breeding"))

bb <- st_bbox(zones_merc) # bounding box
bbp <- bb |> st_as_sfc() |> st_as_sf() # bounding box as a polygon
# coastline (just for background/plotting)
coast <- ne_countries(scale=50,continent = c("North America","Asia","Europe"))
coast <- coast |> 
  filter(geounit!="Greenland") |> # greenland really messing up the transformation
  st_transform(3832) |> 
  st_filter(bbp) # filter to only countries that touch the zones

# plot
pz <- ggplot()+
  geom_sf(data=coast,fill="#BFB89E",color=NA)+
  geom_sf(data=zones_merc,fill=NA,aes(linetype=type))+
  geom_sf_text(data=zones_merc,aes(label=znum))+ 
  # scale_color_manual(values=c("#238B8B","#8B7355"))+
  xlim(bb[1],bb[3]-1e6)+ylim(bb[2],bb[4]-2e6)+
  labs(linetype="",x="",y="")+
  scale_linetype_manual(values=c(3,1))+
  theme(panel.background = element_rect(fill="lightblue"),legend.position=c(0.9,0.6),
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",color='gray80'),
        panel.grid.major = element_blank(),panel.grid.minor=element_blank())

pabun <- plot_abundance(dq,TMBout,opt = 'feed')+
  facet_wrap(~zonelab,nrow=2)+
  coord_cartesian(ylim=c(0,15000))+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
pabun
pmix <- plot_mixing(dq)+theme(margins = unit(c(0,0.5,0.2,1),units='in'),
                              text=element_text(size=6),
                              axis.title = element_text(size=10))
p1 <- plot_grid(pz,pmix,nrow=2,rel_heights=c(1.1,1),labels=c("a","b"))

fig1 <- plot_grid(p1,pabun,ncol=2,labels=c("","c"),rel_widths=c(1.1,1))
fig1
ggsave(here('writing','fig1.png'),fig1,w=15,h=7,bg="white")

# 
# fabun <- plot_abundance(dq,TMBout,opt="feed",return_df = T)
# fabun$zone <- factor(fabun$zone,levels=unique(fabun$zone))
# 
# # survey data
# Surveys <- read_csv(here('Diags','SurveyUse',paste0(TMBout$Code,TMBout$SensCase,".csv")),show_col_types = F)
# # rescaling relative abundance
# Qvec <- TMBout$report$Qest
# Qvec[1] <- 1
# survd <- Surveys |> 
#   filter(Hypothesis%in%c("All","B2F1","B2 only","F1 only"),Use=="Yes") |> 
#   rename(year=Year2) |> 
#   dplyr::select(year,Estimate,CV,Area,Class,Rel,Hypothesis,Class,component=Component) |> 
#   mutate(sd.estimate=Estimate*CV,
#          q = 1/Qvec[Class],
#          rescaled.est = q*Estimate,
#          rescaled.upper=q*(Estimate+sd.estimate),
#          rescaled.lower=q*(Estimate-sd.estimate)) |> 
#   mutate(component=ifelse(component==0,"0+","1+"))
# 
# surv <- survd |> filter(Area%in%unique(fabun$zone)) |> 
#   rename(zone=Area)
# survRel <- surv |> filter(Rel=="Rel")
# survAbs <- surv |> filter(Rel=="Abs")
# # make individual plots for the feeding grounds
# pabun <- map(unique(fabun$zone),\(x){
#   d <- fabun |> filter(zone==x)
#   survRel <- survRel |> filter(zone==x)
#   survAbs <- survAbs |> filter(zone==x)
#   p <-ggplot()+
#     geom_ribbon(data=d,aes(x=year,ymin=low,ymax=upper,fill=component),alpha=0.5)+
#     geom_ribbon(data=d,aes(x=year,ymin=lowmid,ymax=uppermid,fill=component),alpha=0.7)+
#     geom_line(data=d,aes(year,median,color=component))+
#     scale_fill_manual(values=c("#756BB1","#238B8B","#E6781E"))+
#     scale_color_manual(values=c("#756BB1","#238B8B","#E6781E"))+
#     labs(x="",y="",fill="",color="",title=unique(d$zone))+
#     theme(legend.position = c(0.2,0.7),plot.background = element_rect(color='black'))+
#     geom_pointrange(data=survRel,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
#                     linetype=2,size=0.5,color="tan")+
#     geom_pointrange(data=survAbs,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
#                     linetype=2,size=0.5,color="gray20")
#   p
# })
# library(patchwork)
# coord_sf(xlim=c(NA,250),datum=NA)

# png(here('writing','fig1.png'), width = 4000, height = 4000, res = 200) 
# pz+inset_element(pabun[[1]],0.05,0.6,0.2,0.8)+
#   inset_element(pabun[[2]],0.35,0.7,0.5,0.9)+
#   inset_element(pabun[[3]],0.55,0.7,0.7,0.9)+
#   inset_element(pabun[[4]],0.75,0.65,0.9,0.85)+
#   inset_element(pabun[[5]],0.84,0.4,0.99,0.6)+
#   inset_element(pabun[[6]],0.84,0.1,0.99,0.3)
# dev.off()

#### Figure 2- Model comparison ####
# Gather total abundance timeseries from each model
# Estimates of r, K, and sigma
rSobj <- read_rds(here('Diags','final','rS Bayes test','B2F1BC_Bayes.rds'))
envSobj <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC_Bayes.rds'))
envKobj <- read_rds(here('Diags','final','env-K Bayes','B2F1BC_Bayes.rds'))
ddobj <- read_rds(here('Diags','final','ddOnly Bayes','B2F1BC_Bayes.rds'))

rSTMBobj <- read_rds(here('Diags','final','rS Bayes test','B2F1BC_TMB.rds'))
envSTMBobj <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC_TMB.rds'))
envKTMBobj <- read_rds(here('Diags','final','env-K Bayes','B2F1BC_TMB.rds'))
ddTMBobj <- read_rds(here('Diags','final','ddOnly Bayes','B2F1BC_TMB.rds'))

rSTMBout <- read_rds(here('Diags','final','rS Bayes test','B2F1BC.rds'))
envSTMBout <- read_rds(here('Diags','final','env-survival Bayes','B2F1BC.rds'))
envKTMBout <- read_rds(here('Diags','final','env-K Bayes','B2F1BC.rds'))
ddTMBout <- read_rds(here('Diags','final','ddOnly Bayes','B2F1BC.rds'))
# 
# MOM6envSobj <- read_rds(here('Diags','final','env-survival test','B2F1BC_Bayes.rds'))

# compare survival "uncertainty" (variance terms)
SFs <- extract(rSobj,pars=c("log_SFsigma")) |> unlist() |> exp()
Ks <- extract(envKobj,pars=c("log_Ksigma")) |> unlist() |> exp()
eps <- extract(envSobj,pars=c("log_sigmaEnv")) |> unlist() |> exp()
df <- tibble(random=SFs,envK=Ks,envS=eps) |> 
  pivot_longer(everything(),names_to="type",values_to="value")
psigma <- df |> 
  ggplot(aes(value,fill=type,color=type))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("#756BB1","#238B8B","#E6781E"),
                    labels=c("envK"=expression(paste("env-K: ", sigma[K]^2)),
                    "envS"=expression(paste("env-S: ", sigma[Epsilon]^2)),
                    "random"=expression(paste("random: ",sigma[epsilon]^2))))+
  scale_color_manual(values=c("#756BB1","#238B8B","#E6781E"),
                    labels=c("envK"=expression(paste("env-K: ", sigma[K]^2)),
                             "envS"=expression(paste("env-S: ", sigma[Epsilon]^2)),
                             "random"=expression(paste("random: ",sigma[epsilon]^2))))+
  labs(fill="Term",color="Term",x="Uncertainty",y="Density")+
  theme(legend.position=c(0.8,0.7),
        legend.key.size = unit(0.5,'cm'),
        legend.box.background = element_rect(color='black'),
        legend.text = element_text(size=10))
psigma
ggsave(here('plots','final','survival variability sigma comparison.png'),psigma,w=6,h=5)
# compare rmax (max growth rate)
r1 <-extract(ddobj,pars=c("rval")) |> unlist()
r2 <-extract(rSobj,pars=c("rval")) |> unlist()
r3 <-extract(envSobj,pars=c("rval")) |> unlist()
r4 <-extract(envKobj,pars=c("rval")) |> unlist()

df <- tibble(dd=r1,random=r2,envK=r4,envS=r3) |> 
  pivot_longer(everything(),names_to="type",values_to="value")

pr.1 <- df |>
  mutate(typelab=case_when(
    type=="dd" ~"Density\ndependence",
    type=="envK" ~"Environment\nK",
    type=="envS" ~"Environment\nsurvival",
    type=="random" ~"Random",
  )) |> 
  ggplot(aes(typelab,value,fill=typelab,color=typelab))+
  geom_violin(alpha=0.5)+
  geom_hline(yintercept=0.0905,linetype=2)+
  # annotate("text",x=-0.2,y=0.0905,label="0.0905")+
  guides(fill="none",color="none")+
  scale_fill_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  scale_color_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  labs(x="Model",y=expression(paste("Maximum Growth Rate,    ", r[max])),title="")
pr.1

# compare K
ls <- list(ddobj,rSobj,envSobj,envKobj)

zn <- pluck(rSTMBout,'input','BreedNames') |> as.character()
zlabs <- zones_merc |> st_set_geometry(NULL) |> distinct(Name,label)
numz <- length(zn)

df <- imap(ls,\(x,i){
  kdf <- extract(x,pars=c("logK"))[[1]] |> 
    as_tibble(.name_repair = "unique_quiet") |> 
    set_names(c("Asia","Hawai'i","Mexican Islands","Mainland Mexico","Central America")) |> 
    pivot_longer(everything(),names_to="breed",values_to="K") |> 
    mutate(type=i)
  kdf
}) |> list_rbind()|>
  mutate(type=as.factor(c("dd","random","envS","envK")[type]))
pK <- df |> 
  mutate(typelab=case_when(
    type=="dd" ~"Density\ndependence",
    type=="envK" ~"Environment\nK",
    type=="envS" ~"Environment\nsurvival",
    type=="random" ~"Random",
  )) |> 
  ggplot(aes(breed,K,fill=typelab,color=typelab))+
  geom_violin(alpha=0.5,quantile.linetype = 1,quantiles=c(0.5),quantile.color = 'black')+
  scale_fill_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  scale_color_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  theme(legend.position = 'top')+
  labs(fill="Model",color="Model",x="Breeding Ground",y="Carrying Capacity, ln (K)",title="")
pK

# Calc range of total K?

dfK <- df |> 
  mutate(K=exp(K),iter=rep(rep(1:32000,each=5),4)) |> 
  group_by(type,iter) |> 
  summarise(totK=sum(K)) |> 
  group_by(type) |> 
  summarise(Kmean=mean(totK),K2.5=quantile(totK,0.025),K97.5=quantile(totK,0.975))

# compare BK
df <- imap(ls,\(x,i){
  bkdf <- extract(x,pars=c("logBK"))[[1]] |> 
    as_tibble(.name_repair = "unique_quiet") |> 
    set_names(c("Asia","Hawai'i","Mexican Islands","Mainland Mexico","Central America")) |> 
    pivot_longer(everything(),names_to="breed",values_to="logBK") |> 
    mutate(type=i)
  bkdf
}) |> list_rbind()|>
  mutate(type=as.factor(c("dd","random","envS","envK")[type]))
pBK <- df |> 
  mutate(typelab=case_when(
    type=="dd" ~"Density\ndependence",
    type=="envK" ~"Environment\nK",
    type=="envS" ~"Environment\nsurvival",
    type=="random" ~"Random",
  )) |> 
  mutate(BK=1/exp(logBK)) |> 
  ggplot(aes(breed,logBK,fill=typelab,color=typelab))+
  geom_violin(alpha=0.5,quantile.linetype = 1,quantiles=c(0.5),quantile.color = 'black')+
  scale_fill_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  scale_color_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  # theme(legend.position = 'top')+
  guides(color='none',fill='none')+
  coord_cartesian(ylim=c(2,8))+
  labs(fill="Model",color="Model",x="Breeding Ground",y="Initial Depletion, ln(K/B)",title="")
pBK

# df <- imap(ls,\(x,i){
#   kdf <- summary(x) |>
#     pluck('summary') |> 
#     as_tibble(rownames="parameter") |> 
#     set_names(c("parameter","mean","se_mean","sd","low","lowmid","median","uppermid","upper","n_eff","Rhat")) |> 
#     filter(grepl("logK",parameter)) |> 
#     mutate(type=i)
#   kdf
# }) |> list_rbind()
#   
# df2 <- df |> 
#   mutate(breed=rep(zn,4)) |> 
#   left_join(zlabs,by=join_by(breed==Name)) |> 
#   mutate(type=as.factor(c("dd","random","envS","envK")[type]))
  
# pK <- df2 |> 
#   ggplot(aes(label,fill=type,color=type))+
#   geom_linerange(aes(ymin=low,ymax=upper),linewidth=0.5,position=position_dodge(width=0.4),linetype=2)+
#   geom_pointrange(aes(y=median,ymin=lowmid,ymax=uppermid),linewidth = 1.5,position=position_dodge(width=0.4))+
#   scale_fill_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
#   scale_color_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
#   labs(fill="Model",color="Model",x="Breeding Ground",y="ln (K)",title="Estimates of K")
# pK

# compare total abundance
rS_dq <- calc_dq(rSobj,rSTMBobj,nsamp = 1000)
envS_dq <- calc_dq(envSobj,envSTMBobj,nsamp = 1000)
envK_dq <- calc_dq(envKobj,envKTMBobj,nsamp = 1000)
dd_dq <- calc_dq(ddobj,ddTMBobj,nsamp = 1000)

rSabun <- plot_abundance(rS_dq,TMBout = rSTMBout,opt = 'total',return_df = T) |> 
  filter(component=="0+") |> mutate(type="rS")
envSabun <- plot_abundance(envS_dq,TMBout = envSTMBout,opt = 'total',return_df = T) |> 
  filter(component=="0+") |> mutate(type="envS")
envKabun <- plot_abundance(envK_dq,TMBout = envKTMBout,opt = 'total',return_df = T) |> 
  filter(component=="0+") |> mutate(type="envK")
ddabun <- plot_abundance(dd_dq,TMBout = ddTMBout,opt = 'total',return_df = T) |> 
  filter(component=="0+") |> mutate(type="dd")

# survey data
all_abun <- list(rSabun,envSabun,envKabun,ddabun) |> list_rbind() |> 
  mutate(typelab=case_when(
    type=="dd" ~"Density dependence",
    type=="envK" ~"Environment K",
    type=="envS" ~"Environment survival",
    type=="rS" ~"Random"))
all_abun_trunc <- all_abun |> 
  filter(year>1994)
surv <- all_abun |> filter(!is.na(rescaled.est)) |> distinct()

pabun <-ggplot()+
  geom_ribbon(data=all_abun_trunc,aes(x=year,ymin=low,ymax=upper,fill=typelab),alpha=0.5)+
  geom_ribbon(data=all_abun_trunc,aes(x=year,ymin=lowmid,ymax=uppermid,fill=typelab),alpha=0.7)+
  geom_line(data=all_abun_trunc,aes(year,median,color=typelab))+
  geom_pointrange(data=surv,aes(year,rescaled.est,ymax=rescaled.upper,ymin=rescaled.lower),
                  linetype=2,size=0.1)+
  scale_y_continuous(labels = scales::label_number(scale = 1/1000))+
  scale_fill_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  scale_color_manual(values=c("#E63946","#756BB1","#238B8B","#E6781E"))+
  guides(fill="none",color="none")+
  facet_wrap(~typelab)+
  labs(x="Year",y="Total Abundance (1000s)",fill="",color="")
pabun

fig2 <- cowplot::plot_grid(pabun,pK,pr.1,pBK,nrow=2,rel_heights = c(1.2,1),labels='auto')
fig2
ggsave(here('writing','fig2.png'),fig2,w=12,h=7,bg="white")


#### Figure 3: Time-varying Survival ####
# We want to have environmental coefficients, feeding ground abundance, survival, and mortality

pomegas <- plot_omegas(envSobj,TMBout)+facet_wrap(~zone,nrow=1)

# probability coefs<0?
# from plot_omegas dataframe:
# test <- df2 |> group_by(parameter,zone) |> summarise(negs=sum(value<0),tot=n()) |> mutate(pr0=negs/tot)

penvInd <- plot_envIndex(envS_dq,TMBout)+facet_wrap(~zone,nrow=1)+labs(title="")
psurv <- plot_survival(envS_dq,TMBout)+facet_wrap(~zone,nrow=1)+labs(title="")

plot_feedmort <- function(dqlist,TMBout){
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
  
  p <- mFdf |> 
      filter(zone!="RUS+WAL") |> 
      filter(year>1994) |> 
      ggplot(aes(x=year))+
      # geom_line(aes(y=catch,color="Catch: Feeding"))+
      geom_ribbon(aes(ymin=low,ymax=upper),alpha=0.5,fill="#756BB1")+
      geom_ribbon(aes(ymin=lowmid,ymax=uppermid),alpha=0.7,fill="#756BB1")+
      geom_line(aes(y=median,color="Natural"),color="black")+
      labs(x="Year",y="Mortality (ind.)")+
      facet_wrap(~zone)+
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  }
  
pmort <- plot_feedmort(envS_dq,envSTMBout)+
  facet_wrap(~zone,nrow=1)+
  scale_y_continuous(labels = scales::label_number(scale = 1/1000))+
  labs(y="Mortality (thousands)")

fig3 <- cowplot::plot_grid(pomegas,penvInd,psurv,pmort,ncol=2,labels='auto')
fig3
ggsave(here('writing','fig3.png'),fig3,w=12,h=7,bg="white")

#### Figure 4: Demographic outlook ####
# Realized growth rate and time to recovery
pr <- plot_r(envSobj,envS_dq,TMBout)
prec <- plot_recovery_t(envSobj,envS_dq,TMBout)+
  xlim(1995,NA)+
  labs(y="Minimum Recovery Year")
prec
fig4 <- cowplot::plot_grid(pr,prec,nrow=2,labels='auto')
fig4
ggsave(here('writing','fig4.png'),fig4,w=4,h=8,bg="white")

#### Supplementary ####
# Compare natural and human direct mortality
# dq <- calc_dq(envSobj,TMBobj = envSTMBobj,nsamp = 3000)
# 
# plot_mort_scatter <- function(dqlist,TMBout){
#   yrs <- pluck(TMBout,"input","Years")
#   yrs <- yrs[-length(yrs)]
#   numyr <- length(yrs)
#   bn <-pluck(TMBout,'input','BreedNames') |> as.character()
#   fn <- pluck(TMBout,'input','FeedNames') |> as.character()
#   
#   
#   ### Annual Mortality by Herd ##
#   mdat <- pull_dq(dqlist,"MortDiff") |> simplify2array()
#   mdat <- mdat*-1
#   mdat[mdat<0] <- 0
#   
#   ### natural mortality by zone ##
#   # breeding
#   # NbS <- pull_dq(dqlist,"NbS") |> simplify2array()
#   # NbS <- NbS[,-dim(NbS)[2],]
#   # feeding
#   NfS <- pull_dq(dqlist,"NfS") |> simplify2array()
#   NfS <- NfS[,-dim(NfS)[2],]
#   catchf <- TMBout$input$CatchF |> t()
#   catchfdf <- tibble(catch=as.numeric(TMBout$input$CatchF),
#                      zone=rep(fn,each=length(yrs)),
#                      zonetype='feed',
#                      year=rep(yrs,length(fn))) |> 
#     # calculate cumulative catch
#     group_by(zone) |> 
#     arrange(year) |> 
#     mutate(cume_catch=cumsum(catch)) |> 
#     ungroup()
#   # feeding grounds (sum over breeding grounds)
#   mdatF <- apply(mdat,c(2,3,4),sum)
#   mFquants <- apply(mdatF,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
#   mFdf <- tibble(zone=rep(fn,length(yrs)),
#                  year=rep(yrs,each=length(fn))) |> 
#     mutate(across(where(is.list),as.character)) |>
#     mutate(low=as.numeric(mFquants[1,,]),
#            lowmid=as.numeric(mFquants[2,,]),
#            median=as.numeric(mFquants[3,,]),
#            uppermid=as.numeric(mFquants[4,,]),
#            upper=as.numeric(mFquants[5,,]))
#   # annual catch
#   mFdf <- mFdf |> 
#     left_join(catchfdf,by=join_by(year,zone))
#   p <- mFdf |> 
#     filter(zone!="RUS+WAL") |> 
#     filter(year>1994) |> 
#     mutate(mhw=ifelse(year %in% 2013:2017,"2013-2017", "other\nyears")) |> 
#     ggplot(aes(x=catch,y=median,ymin=low,ymax=upper,color=mhw))+
#     geom_point()+
#     facet_wrap(~zone,scales='free')+
#     labs(x="Direct Mortality",y="Excess Natural Mortality",color=NA)+
#     scale_color_manual(values=c("#E6781E","#238B8B"))
#   p
# }


# cumulative mortality for 2013-2023
# cmdat <- apply(mdat,c(1,2,4),cumsum) |> aperm(c(2,3,1,4))
# cmdatF <- apply(cmdat,c(2,3,4),sum)
# cmFquants <- apply(cmdatF,c(1,2),quantile,probs=c(0.025,0.25,0.50,0.75,0.975))
# cmFdf <- tibble(zone=rep(fn,length(yrs)),
#                 year=rep(yrs,each=length(fn))) |> 
#   mutate(across(where(is.list),as.character)) |>
#   mutate(low=as.numeric(cmFquants[1,,]),
#          lowmid=as.numeric(cmFquants[2,,]),
#          median=as.numeric(cmFquants[3,,]),
#          uppermid=as.numeric(cmFquants[4,,]),
#          upper=as.numeric(cmFquants[5,,]))
# cmFdf <- cmFdf |> 
#   left_join(catchfdf,by=join_by(year,zone))
# # Table
# cmFdf |> filter(year %in% c(2012,2023)) |> select(year,zone,median) |> group_by(zone) |> mutate(diff=median[year==2023]-median[year==2012])


# environmental coefficients
# df <- summary(envSobj)$summary |> as_tibble(rownames="parameter")
# # split into -1:1 scale and larger
# df2 <- df |> 
#   filter(parameter!="lp__",
#          !grepl("epsEnv",parameter),
#          !grepl("Kdev",parameter),
#          !grepl("SFdev",parameter)) |> 
#   mutate(type=case_when(
#     grepl("MixPars",parameter)~"Mixing",
#     grepl("logBK",parameter) ~ "Initial Depletion",
#     grepl("envParams",parameter) ~ 'Environmental',
#     grepl("logK",parameter)~"Carrying Capacity",
#     TRUE ~ "Other"
#   ))
