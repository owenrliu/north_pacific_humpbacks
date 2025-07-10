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

# try reading some outputs?
totabun_B2F2 <- read_table(here('Diags','B2F2BC.Out'),skip=188,n_max=55,col_names=c("Year","est","SD_log"))
totabun_B2F2 %>% ggplot(aes(Year,est))+geom_line()+
  theme_classic()+
  xlim(1980,NA)

# survival
b2f2_surv_asia <- read_table(here('Diags','B2F2BC.Out'),skip=1560,na='NaN',n_max=54,col_names=c("Year","RUS+WAL", "EAL+BER+WGOA", "NGOA", "SEA+NBC", "SBC+WA", "OR+CA"),col_types='idddddd') %>% 
  pivot_longer(`RUS+WAL`:`OR+CA`,names_to="feed",values_to="surv") %>% mutate(breed="Asia")
b2f2_surv_hawaii <- read_table(here('Diags','B2F2BC.Out'),skip=1616,na='NaN',n_max=54,col_names=c("Year","RUS+WAL", "EAL+BER+WGOA", "NGOA", "SEA+NBC", "SBC+WA", "OR+CA"),col_types='idddddd') %>% 
  pivot_longer(`RUS+WAL`:`OR+CA`,names_to="feed",values_to="surv") %>% mutate(breed="Hawaii")
b2f2_surv_mxar <- read_table(here('Diags','B2F2BC.Out'),skip=1672,na='NaN',n_max=54,col_names=c("Year","RUS+WAL", "EAL+BER+WGOA", "NGOA", "SEA+NBC", "SBC+WA", "OR+CA"),col_types='idddddd') %>% 
  pivot_longer(`RUS+WAL`:`OR+CA`,names_to="feed",values_to="surv") %>% mutate(breed="MX_AR")
b2f2_surv_mxml <- read_table(here('Diags','B2F2BC.Out'),skip=1728,na='NaN',n_max=54,col_names=c("Year","RUS+WAL", "EAL+BER+WGOA", "NGOA", "SEA+NBC", "SBC+WA", "OR+CA"),col_types='idddddd') %>% 
  pivot_longer(`RUS+WAL`:`OR+CA`,names_to="feed",values_to="surv") %>% mutate(breed="MX_ML")
b2f2_surv_cenam <- read_table(here('Diags','B2F2BC.Out'),skip=1784,na='NaN',n_max=54,col_names=c("Year","RUS+WAL", "EAL+BER+WGOA", "NGOA", "SEA+NBC", "SBC+WA", "OR+CA"),col_types='idddddd') %>% 
  pivot_longer(`RUS+WAL`:`OR+CA`,names_to="feed",values_to="surv") %>% mutate(breed="CenAm")

b2f2_surv_all <- list(b2f2_surv_asia,b2f2_surv_hawaii,b2f2_surv_cenam,b2f2_surv_mxar,b2f2_surv_mxml) %>% bind_rows()

b2f2_surv_all %>% 
  # group_by(Year,feed) %>% 
  # summarise(mean_surv=mean(surv,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(Year,surv,color=breed),alpha=0.7)+
  geom_line(linewidth=1.25)+
  geom_hline(yintercept=0.96,linetype=2)+
  geom_vline(xintercept=2015,linetype=2)+
  theme_classic()+
  facet_wrap(~feed)

# Check out MOM6 SST
library(tidync)
library(sf)
library(ncdf4)
library(viridis)

ncfn <- here('data','mom6','tos.nep.full.hcast.monthly.regrid.r20250509.199301-202412.nc')
nc_open(ncfn)
x<-tidync(ncfn)

times <- x %>% activate("D0") %>% hyper_tibble() %>%mutate(ind=row_number()) %>% mutate(date=as_date(time)) %>% mutate(year=year(date),month=month(date))

# find a random month
july2014ind <- times %>% filter(year==2014,month==7) %>% pull(ind)

july2014 <- tidync(ncfn) %>% hyper_filter(time=index==july2014ind) %>% hyper_tibble() %>% 
  mutate(across(everything(),as.numeric))
july2014 %>% 
  st_as_sf(coords=c("lon","lat"),crs=4326) %>% 
  ggplot(aes(color=tos))+
  geom_sf()+
  scale_color_viridis(option='turbo')+
  theme_classic()+
  labs(color="Surface T\nJuly 2014")

# Breeding/feeding circular diagram
mixing <- read_table(here('Diags','B2F2BC.Out'),skip=982,n_max=5,col_names=c("RUS+WAL", "EAL+BER+WGOA", "NGOA", "SEA+NBC", "SBC+WA", "OR+CA"))
mixing <- mixing %>% mutate(from=c("Asia","Hawaii","MX_AR","MX_ML","CenAm")) %>% pivot_longer(-from,names_to="to",values_to="prop")

library(circlize)
chordDiagram(mixing)
