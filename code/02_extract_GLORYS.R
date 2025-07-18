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

# # get the GLORYS grids from one of the nc files
# # this is for physics, which is higher resolution than the biogeochemistry
# gr <- tidync(here('data','glorys','phys','raw','glorys_physics_1995.nc'))
# # this is for the bgc
# gr_bgc <- tidync(here('data','glorys','bgc','raw','glorys_bgc_1995.nc'))
# 
# # Footprint of the FEAT survey data, as a polygon (created in 00_make_FEAT_grid.R)
# feat_foot <- read_rds(here('data','grids','FEAT_footprint.rds'))
# # in lat/lon coords
# feat_foot_ll <- feat_foot %>% st_transform(4326)
# # bounding box for the FEAT footprint
# bb <- st_bbox(feat_foot_ll)
# 
# ## GLORYS PHYSICS
# # grab first time/depth slice, just to get grid
# # note, this first part is only for the 4D variables (i.e. indexed by time, lat, lon, depth)
# # other variables are stored in different grids because they are depth-invariant:
# # i.e., mixed layer depth, ssh, and bottom temperature are 3D only (lat/lon/time)
# slice1 <- gr %>% hyper_filter(time=index==2) %>% hyper_tibble()
# # converted from dataframe to sf object
# slice1_sf <- slice1 %>% st_as_sf(coords=c("longitude","latitude"),crs=4326)
# 
# # plot, relative to FEAT
# ggplot()+
#   geom_sf(data=coast,fill='gray80')+
#   geom_sf(data=slice1_sf %>% sample_n(10000))+
#   geom_sf(data=feat_foot,color='red',fill=NA)+
#   xlim(bb[1],bb[3])+ylim(bb[2],bb[4])
# 
# # see, we only need a small subset of these values. Let's do a simple overlay and 
# # keep track of which GLORYS points overlap with the FEAT footprint
# glorys_in_feat <- slice1_sf %>% 
#   # make a grid index. if we're careful, this should be consistent for every year of data
#   mutate(glorID=row_number()) %>% 
#   st_filter(feat_foot_ll)
# 
# # these are the grid indices for the rows of GLORYS data we will extract for each year/month
# glorIDs <- glorys_in_feat$glorID
# 
# # lets try plotting again
# ggplot()+
#   geom_sf(data=coast,fill='gray80')+
#   geom_sf(data=glorys_in_feat,size=0.25)+
#   geom_sf(data=feat_foot,color='red',fill=NA)+
#   xlim(bb[1],bb[3])+ylim(bb[2],bb[4])
# 
# # the spatial filter/overlay is slow for large data, so BECAUSE GLORYS is consistent in its grid, 
# # we should be able to use non-spatial indexing to filter instead
# 
# # Do this for all year/months
# # 1. open a year's ncdf
# # 2. for each month, extract all data
# # 3. filter by row number using the glorIDs from above
# # then 
# # 4. recombine across months
# 
# extract_glorys_feat_footprint <- function(nc,mth){
#   # extract the data
#   temp <- tidync(nc) %>% hyper_filter(time=index==mth) %>% hyper_tibble()
#   # filter using the grid IDs above
#   tempsub <- temp %>% slice(glorIDs)
#   tempsub %>% 
#     #convert to meaningful dates (instead of hours since 1950-01-01)
#     mutate(date=as_datetime("1950-01-01 00:00:00")+hours(time))
# }
# 
# # glorys physics files lists
# fls <- list.files(here('data','glorys','phys','raw'),full.names = T)
# 
# tic("trying first year, first month")
# test <- extract_glorys_feat_footprint(fls[1],mth=1)
# toc()
# 
# tic("trying first year, all months")
# test <- purrr::map(1:12,extract_glorys_feat_footprint,nc=fls[1]) %>% 
#   list_rbind()
# toc()
# # ~10s for one year
# 
# # finally, apply to all years (except 2021 is a special case)
# for(i in c(1:27,30,31)){
#   tic(paste("processing",fls[i]))
#   extr <- purrr::map(1:12,extract_glorys_feat_footprint,nc=fls[i]) %>% 
#     list_rbind()
#   fn_out <- fls[i] %>% str_replace(".nc","_spatial_filtered.rds") %>% str_replace("raw","filtered") 
#   write_rds(extr,fn_out)
#   toc()
# }
# 
# # 2021 is weird because it splits two datasets
# # the first has january to july, the second, august to december
# tic(paste("processing",fls[28],"and",fls[29]))
# extr2021.1 <- purrr::map(1:6,extract_glorys_feat_footprint,nc=fls[28]) %>% 
#     list_rbind()
# extr2021.2 <- purrr::map(1:6,extract_glorys_feat_footprint,nc=fls[29]) %>% 
#   list_rbind()
# extr <- bind_rows(extr2021.1,extr2021.2)
# fn_out <- fls[27] %>% str_replace(".nc","_spatial_filtered.rds") %>% str_replace("raw","filtered")
# write_rds(extr,fn_out)
# toc()
# 
# ### 3D VARIABLES #
# # for the 3D variables, the process can be the same
# # find gridIDs
# slice3d <- gr %>% 
#   activate("D2,D1,D3") %>% 
#   hyper_filter(time=index==1) %>% 
#   hyper_tibble()
# 
# # these are the row numbers (i.e., the indices) of the data we will extract for each year/month
# slice3d_sf <- slice1 %>% st_as_sf(coords=c("longitude","latitude"),crs=4326)
# glorys3d_in_feat <- slice3d_sf %>% 
#   mutate(glor3dID=row_number()) %>% 
#   st_filter(feat_foot_ll)
# glor3dIDs <- glorys3d_in_feat$glor3dID
# 
# # same idea, extract filtered data for each year/month
# extract_glorys3D_feat_footprint <- function(nc,mth){
#   # extract the data
#   temp <- tidync(nc) %>% activate("D2,D1,D3") %>%
#     hyper_filter(time=index==mth) %>% hyper_tibble()
#   # filter using the grid IDs above
#   tempsub <- temp %>% slice(glor3dIDs)
#   tempsub %>% 
#     #real times (hours since 1950-01-01)
#     mutate(date=as_datetime("1950-01-01 00:00:00")+hours(time))
# }
# 
# 
# for(i in c(1:27,30,31)){
#   tic(paste("processing",fls[i]))
#   extr <- purrr::map(1:12,extract_glorys3D_feat_footprint,nc=fls[i]) %>% 
#     list_rbind()
#   fn_out <- fls[i] %>% str_replace(".nc","_3d_filtered.rds") %>% str_replace("raw","filtered")
#   write_rds(extr,fn_out)
#   toc()
# }
# 
# # for 2021
# tic(paste("processing",fls[28],"and",fls[29]))
# extr2021.1 <- purrr::map(1:6,extract_glorys3D_feat_footprint,nc=fls[27]) %>% 
#   list_rbind()
# extr2021.2 <- purrr::map(1:6,extract_glorys3D_feat_footprint,nc=fls[28]) %>% 
#   list_rbind()
# extr <- bind_rows(extr2021.1,extr2021.2)
# fn_out <- fls[27] %>% str_replace(".nc","_3d_filtered.rds") %>% str_replace("raw","filtered")
# write_rds(extr,fn_out)
# toc()
# 
# ## GLORYS BGC, with a different resolution than the physics
# # grab first time/depth slice, just to get grid
# slice_bgc <- gr_bgc %>% hyper_filter(time=index==1) %>% hyper_tibble()
# # converted from dataframe to sf object
# slice_bgc_sf <- slice_bgc %>% st_as_sf(coords=c("longitude","latitude"),crs=4326)
# 
# # plot, relative to FEAT
# ggplot()+
#   geom_sf(data=coast,fill='gray80')+
#   geom_sf(data=slice_bgc_sf %>% sample_n(10000))+
#   geom_sf(data=feat_foot,color='red',fill=NA)+
#   xlim(bb[1],bb[3])+ylim(bb[2],bb[4])
# 
# glorys_bgc_in_feat <- slice_bgc_sf %>% 
#   # make a grid index. if we're careful, this should be consistent for every year of data
#   mutate(glorID_bgc=row_number()) %>% 
#   st_filter(feat_foot_ll)
# 
# # these are the grid indices for the rows of GLORYS data we will extract for each year/month
# glorIDs_bgc <- glorys_bgc_in_feat$glorID_bgc
# 
# # lets try plotting again
# ggplot()+
#   geom_sf(data=coast,fill='gray80')+
#   geom_sf(data=glorys_bgc_in_feat,size=0.25)+
#   geom_sf(data=feat_foot,color='red',fill=NA)+
#   xlim(bb[1],bb[3])+ylim(bb[2],bb[4])
# 
# # apply to all bgc year/months
# extract_glorys_bgc_feat_footprint <- function(nc,mth){
#   # extract the data
#   temp <- tidync(nc) %>% hyper_filter(time=index==mth) %>% hyper_tibble()
#   # filter using the grid IDs above
#   tempsub <- temp %>% slice(glorIDs_bgc)
#   tempsub %>% 
#     #convert to meaningful dates (instead of hours since 1950-01-01)
#     mutate(date=as_datetime("1950-01-01 00:00:00")+hours(time))
# }
# 
# # apply to all years
# fls_bgc <- list.files(here('data','glorys','bgc','raw'),full.names = T)
# 
# for(i in 1:length(fls_bgc)){
#   tic(paste("processing",fls_bgc[i]))
#   extr <- purrr::map(1:12,extract_glorys_bgc_feat_footprint,nc=fls_bgc[i]) %>% 
#     list_rbind()
#   fn_out <- fls_bgc[i] %>% str_replace(".nc","_spatial_filtered.rds") %>% str_replace("raw","filtered") 
#   write_rds(extr,fn_out)
#   toc()
# }
# 
# ## At the end of this script, we should have extracted and coarsely spatially filtered both the physical and biogeochemical data
# ## In the next script, we can actually join these data to both our observations and to the projection grid