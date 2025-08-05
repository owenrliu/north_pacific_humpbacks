# Goal: extract data from already-downloaded GLORYS ncdf files (01_download_GLORYS.R)

library(ncdf4) 
library(tidync)
library(data.table)
library(tidyverse)
library(here)
library(geosphere)
library(sf)
library(viridis)
library(cowplot)
library(rnaturalearth)
library(tictoc) # for timing code
sf_use_s2(FALSE)

# NOTE: TRYING TO DECIDE ON A MAP PROJECTION TO USE; NEED TO BE CAREFUL ABOUT AREA OR DISTANCE CALCS
# USING PDC MERCATOR (EPGS 3832) FOR NOW

# shapefile of humpback zones
zones <- read_sf(here('data','spatial','NPhuwhRegions.1Apr2024.shp')) 

zones_merc <- zones %>% 
  # PDC Mercator
  st_transform(3832)%>% unite("zoneID",Name,id,remove=F)
bb <- st_bbox(zones_merc) # bounding box
bbp <- bb %>% st_as_sfc() %>% st_as_sf() # bounding box as a polygon

# coastline (just for background/plotting)
coast <- ne_countries(scale=50,continent = c("North America","Asia","Europe")) %>% st_transform(3832)
coast <- coast %>% filter(geounit!="Greenland") %>% # greenland really messing up the transformation
  st_filter(bbp) # filter to only countries that touch the zones

# plot
ggplot()+geom_sf(data=coast,fill="#dacf80")+
  geom_sf(data=zones_merc,fill='lightblue')+
  geom_sf_text(data=zones_merc,aes(label=zoneID))+ 
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+theme_minimal()

## Index of GLORYS grid cells <=300m (continental shelf), produced in `define spatial zones.qmd`
shelf <- read_rds(here('data','glorys','phys','glorys_bathy_crop_whalezones_300m.rds'))
shelf_merc <- shelf %>% st_transform(3832)

# plot again?

zones_shelf_p <- ggplot()+
  geom_sf(data=coast,fill="#dacf80")+
  geom_sf(data=zones_merc,fill=NA)+
  geom_sf(data=shelf_merc,aes(color=deptho),size=0.1)+
  geom_sf_text(data=zones_merc,aes(label=zoneID))+
  scale_color_viridis()+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+theme_minimal()

ggsave(here('plots','humpback zones with GLORYS bathymetry.png'),zones_shelf_p,w=10,h=7)

# Let's try extracting data from the large GLORYS netCDFs
# Write a function that uses the zoneID, year, and month and loads the correct GLORYS data
# filtered by the shelf/300m dataset

extr_glorys_yr_mth_zone <- function(y,m,z,plot=F){
  # first, find the right file
  if(y%in%1993:2020) fn <- paste0("data/glorys/phys/raw/glorys_physics_",z,"_1993-2021.nc")
  if(y %in% 2022:2024) fn <- paste0("data/glorys/phys/raw/glorys_physics_",z,"_2021-2024.nc")
  if(y==2021 & m %in% 1:6) fn <- paste0("data/glorys/phys/raw/glorys_physics_",z,"_1993-2021.nc")
  if(y==2021 & m %in% 7:12) fn <- paste0("data/glorys/phys/raw/glorys_physics_",z,"_2021-2024.nc")
  
  nc <- tidync(fn)
  
  # shelf dataset, filtered
  shelffilt <- shelf_merc %>% filter(zoneID==z)
  
  # datetimes
  dt <- nc %>% activate('time') %>% hyper_tibble() %>% 
    mutate(date=as_date(time)) %>% 
    mutate(year=year(date),month=month(date)) %>% 
    mutate(ind=row_number()) # the time indices, for filtering the nc
  
  # filter to the correct year/month and pull the time index
  tind <- dt %>% filter(year==y,month==m) %>% pull(ind)
  
  # now pull the data
  # temperature
  xsst <- nc %>% hyper_filter(time=index==tind) %>% hyper_tibble()
  # MLD
  xmld <- nc %>% activate('mlotst') %>% hyper_filter(time=index==tind) %>% hyper_tibble()
  
  # join and then crop
  x <- xsst %>% left_join(xmld,by=join_by(latitude,longitude)) %>% mutate(across(latitude:longitude,as.numeric))
  xsf <- x %>% st_as_sf(coords=c("longitude","latitude"),crs=4326) %>% st_transform(3832)
  xsf_crop <- xsf %>% 
    # join the shelf dataset, and allow points to join that are within 100m
    st_join(shelffilt,join=st_is_within_distance,dist=100) %>% 
    mutate(include=replace_na(include,0L)) %>% 
    # spatial filter- is the GLORYS 
    st_filter(zones %>% filter(zoneID==z)) %>% 
    # then, we have the "include" variable from the shelf dataset to filter later
    mutate(year=y,month=m,glorys_time_ind=tind)
  
  if(plot){
    mldp <- ggplot()+
      geom_sf(data=xsf_crop,aes(color=mlotst))+scale_color_viridis(option="turbo")+
      theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA))
    sstp <- ggplot()+
      geom_sf(data=xsf_crop,aes(color=thetao))+scale_color_viridis(option="turbo")+
      theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA))
    p <- cowplot::plot_grid(mldp,sstp,nrow = 1)
    print(p)
  }
  xsf_crop
}

# try this
# WGOA, July 2020
test <- extr_glorys_yr_mth_zone(2023,12,z="RUS+WAL_15",plot=T)
# it works!

# calculate an anomaly value, from a defined baseline period (start and end years)
calc_anom_ts <- function(z,m,base_start,base_end){
  dat <- map(1993:2024,\(yr) extr_glorys_yr_mth_zone(yr,m,z=z)) %>% 
    bind_rows()
  dat_out <- dat %>% 
    # drop spatial info for now
    st_set_geometry(NULL) %>% 
    # filter for shelf only
    filter(include==1) %>% 
    rename(sst=thetao,mld=mlotst) %>% 
    mutate(period=ifelse(year %in% base_start:base_end,"baseline","non-baseline")) %>% 
    group_by(year) %>% 
    # make a unique zone/grid cell id number
    mutate(grid_idx=row_number()) %>% 
    mutate(grid_idx=paste(grid_idx,z)) %>% 
    ungroup() %>% 
    group_by(grid_idx,month) %>% 
    # make derived calculations
    # baseline mean sst and mld
    mutate(sst_baseline=mean(sst[period=="baseline"]),mld_baseline=mean(mld[period=="baseline"])) %>%
    ungroup() %>% 
    # and the anomalies off of the baselines
    mutate(sst_anom=sst-sst_baseline,mld_anom=mld-mld_baseline)
  
  dat_out
    
}

# try to extract July mean field, 1993:2013
tic("Extracting July EAL+BER_16: ")
july_EAL_BER <- calc_anom_ts("EAL+BER_16",m=7,base_start=1993,base_end=2013)
toc()
#this took 109s and produced 1.95M data rows

tic("Extracting July NGOA_8: ")
july_NGOA <- calc_anom_ts("NGOA_8",m=7,base_start=1993,base_end=2013)
toc()
#this took 75s and produced 380K data rows


# and then, an average anomaly (across the spatial domain) in each year
july_eal_ber_ts <- july_EAL_BER %>% 
  group_by(year) %>% 
  summarise(sst_anom=mean(sst_anom),mld_anom=mean(mld_anom))

july_eal_ber_ts %>% ggplot(aes(year,sst_anom))+geom_line()

july_NGOA_ts <- july_NGOA %>% 
  group_by(year) %>% 
  summarise(sst_anom=mean(sst_anom),mld_anom=mean(mld_anom),
            sst_q90=quantile(sst,0.9))

july_NGOA_ts %>% ggplot(aes(year,sst_anom))+geom_line()+geom_smooth(method='lm')