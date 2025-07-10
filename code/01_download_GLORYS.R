# Go Get GLORYS data and download it, using python workaround
# from https://github.com/pepijn-devries/CopernicusMarine/issues/42#issuecomment-2080179599
library(tidyverse)
library(sf)
library(here)
library(reticulate)
library(tidync)
library(viridis)
library(tictoc)

# install_python() use this instead of the next line if this is the first time using reticulate/python
# install_python()
# use_python("C:/Users/owenrliu/AppData/Local/Programs/Python/Python313/python.exe")

virtualenv_create(envname = "thiscopernicus")
# virtualenv_install("thiscopernicus", packages = c("copernicusmarine"))
use_virtualenv("thiscopernicus", required = TRUE)
# py_install("copernicusmarine")

cm <- (import("copernicusmarine"))
# cm$login() # I was having trouble with this, so I did it in the R Terminal with cmd 'copernicusmarine login'
# then i could put in my credentials, which saved to my registry

# copernicusmarine subset -i cmems_mod_glo_phy_my_0.083deg_P1M-m -x -126 -X -115 -y 32 -Y 50 -z 120. -Z 150. -v thetao -t 1994-01-01 -T 2024-08-20 -o ./copernicus-data -f glorys_data_for_hake_temperature.nc
# copernicusmarine subset -i cmems_mod_glo_phy_myint_0.083deg_P1M-m -x -126 -X -115 -y 32 -Y 50 -z 120. -Z 150. -v thetao -t 1994-01-01 -T 2024-08-20 -o ./copernicus-data -f glorys_data_for_hake_temperature2.nc

# bounding boxes for humpback areas
zones <- read_sf(here('data','spatial','NPhuwhRegions.1Apr2024.shp'))
zones_bbs <- map(zones$geometry,st_bbox)

# need to shift the longitude so it doesn't plot backwards across the intl dateline (i.e. move to the 0-360 longitude instead of -180-180)
zones_360 <- zones %>% st_shift_longitude() %>% unite("zoneID",Name,id)
zones_360_bbs <- map(zones_360$geometry,st_bbox) %>% set_names(zones_360$zoneID)

# quick plot
zones_360 %>% ggplot()+geom_sf(fill='lightblue')+theme_minimal()+geom_sf_text(aes(label=Name))

## DOWNLOAD PHYSICS DATA FROM GLORYS HINDCAST
# 1993 to 06/2021: cmems_mod_glo_phy_my_0.083deg_P1M-M
# 07/2021-present: cmems_mod_glo_phy_myint_0.083deg_P1M-m
# variables: (not exhaustive, check back if something important is missing later)
# mixed layer depth mlotst
# salinity so
# temperature thetao
# velocities, uo and vo
# ssh zos
# bottom temperature, bottomT


if(!dir.exists(here('data','glorys','phys','raw'))) {dir.create(here('data','glorys','phys','raw'),recursive = T)}

#quick function to download mixed layer depth and temperature for the surface layer
# for a given bounding box and year
dl_glorys_yr <- function(year,zoneID){
  
  bbox <- zones_360_bbs %>% pluck(zoneID)
  
  minlon <- bbox$xmin
  maxlon <- bbox$xmax
  minlat <- bbox$ymin
  maxlat <- bbox$ymax
  
  tic(paste("Downloading GLORYS: ",year," ",zoneID))
  dat <- cm$subset(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
    variables=list("mlotst","thetao"), #mld and sst
    minimum_longitude=minlon,
    maximum_longitude=maxlon,
    minimum_latitude=minlat,
    maximum_latitude=maxlat,
    start_datetime=paste0(year,"-01-01T00:00:00"),
    end_datetime=paste0(year,"-12-31T23:59:00"),
    minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
    maximum_depth=0.49402499198913574,
    output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_",year,".nc")
  )
  toc()
}

# a version to download all years
# for the version up to 2021:
dl_glorys_1993_2021_allyrs <- function(zoneID){
  
  bbox <- zones_360_bbs %>% pluck(zoneID)
  
  minlon <- bbox$xmin
  maxlon <- bbox$xmax
  minlat <- bbox$ymin
  maxlat <- bbox$ymax
  
  tic(paste("Downloading GLORYS: ",zoneID))
  dat <- cm$subset(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
    variables=list("mlotst","thetao"),
    minimum_longitude=minlon,
    maximum_longitude=maxlon,
    minimum_latitude=minlat,
    maximum_latitude=maxlat,
    start_datetime=paste0(1993,"-01-01T00:00:00"),
    end_datetime=paste0(2021,"-12-31T23:59:00"),
    minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
    maximum_depth=0.49402499198913574,
    output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_1993-2021.nc")
  )
  toc()
}

# for 2021-2024
dl_glorys_2021_2024_allyrs <- function(zoneID){
  
  bbox <- zones_360_bbs %>% pluck(zoneID)
  
  minlon <- bbox$xmin
  maxlon <- bbox$xmax
  minlat <- bbox$ymin
  maxlat <- bbox$ymax
  
  tic(paste("Downloading GLORYS: ",zoneID))
  dat <- cm$subset(
    dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1M-m",
    variables=list("mlotst","thetao"),
    minimum_longitude=minlon,
    maximum_longitude=maxlon,
    minimum_latitude=minlat,
    maximum_latitude=maxlat,
    start_datetime=paste0(2021,"-07-01T00:00:00"),
    end_datetime=paste0(2024,"-12-31T23:59:00"),
    minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
    maximum_depth=0.49402499198913574,
    output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_2021-2024.nc")
  )
  toc()
}

# test
dl_glorys_yr(year=1993,zoneID="EAL+BER_16")
xsst <- tidync(here('data','glorys','phys','raw','glorys_physics_EAL+BER_16_1993.nc')) %>% hyper_tibble() %>% 
  mutate(across(longitude:latitude,as.numeric)) %>% 
  mutate(time=as_date(time))
xmld <- tidync(here('data','glorys','phys','raw','glorys_physics_EAL+BER_16_1993.nc')) %>% 
  activate("mlotst") %>% 
  hyper_tibble() %>% 
  mutate(across(longitude:latitude,as.numeric)) %>% 
  mutate(time=as_date(time))
xsst %>% filter(time=="1993-01-01") %>% ggplot(aes(longitude,latitude,color=thetao))+geom_point()+scale_color_viridis(option='turbo')
xmld %>% filter(time=="1993-01-01") %>% ggplot(aes(longitude,latitude,color=mlotst))+geom_point()+scale_color_viridis(option='turbo')

# try downloading everything by zone
# for 1993:2021
walk(zones_360$zoneID,\(z) dl_glorys_1993_2021_allyrs(z))
# this took about ~30 minutes

# for 2021:2024
walk(zones_360$zoneID,\(z) dl_glorys_2021_2024_allyrs(z))

### END

### OLDER 
# for(i in startyr:2021){
#   tic(paste("Downloading GLORYS: ",i))
#   dat <- cm$subset(
#     dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
#     variables=list("mlotst","thetao"),
#     minimum_longitude=-140,
#     maximum_longitude=-117,
#     minimum_latitude=32,
#     maximum_latitude=59,
#     start_datetime=paste0(i,"-01-01T00:00:00"),
#     end_datetime=paste0(i,"-12-31T23:59:00"),
#     minimum_depth=0.49402499198913574,
#     maximum_depth=2000,
#     output_filename = paste0("data/glorys/phys/raw/glorys_physics_",i,".nc"),
#     force_download = TRUE
#   )
#   toc()
# }
# 
# # dataset switch since 07-01-2021
# # cmems_mod_glo_phy_myint_0.083deg_P1M
# for(i in 2021:endyr){
#   tic(paste("Downloading GLORYS: ",i))
#   dat <- cm$subset(
#     dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1M-m",
#     variables=list("mlotst", "so", "thetao", "uo", "vo", "zos", "bottomT"),
#     minimum_longitude=-140,
#     maximum_longitude=-117,
#     minimum_latitude=32,
#     maximum_latitude=59,
#     start_datetime=paste0(i,"-01-01T00:00:00"),
#     end_datetime=paste0(i,"-12-31T23:59:00"),
#     minimum_depth=0.49402499198913574,
#     maximum_depth=2000,
#     output_filename = paste0("data/glorys/phys/raw/glorys_physics_",i,".nc"),
#     force_download = TRUE
#   )
#   toc()
# }
# 
# ## DOWNLOAD BIO-GEO-CHEMICAL DATA FROM GLORYS HINDCAST
# # for BGC data, different resolution and different datasets
# # 1993 to end of 2022: cmems_mod_glo_bgc_my_0.25deg_P1M-m
# # 2023-present: cmems_mod_glo_bgc_myint_0.25deg_P1M-m
# # variables: (not exhaustive, check back if something important is missing later)
# # chlorophyll
# # no3
# # nppv
# # o2
# # ph
# # phytoplankton
# 
# for(i in startyr:endyr){
#   tic(paste("Downloading GLORYS BGC: ",i))
#   dat <- cm$subset(
#     dataset_id=ifelse(i %in% 1993:2022,"cmems_mod_glo_bgc_my_0.25deg_P1M-m","cmems_mod_glo_bgc_myint_0.25deg_P1M-m"),
#     variables=list("chl", "no3", "nppv", "o2", "ph", "phyc"),
#     minimum_longitude=-140,
#     maximum_longitude=-117,
#     minimum_latitude=32,
#     maximum_latitude=59,
#     start_datetime=paste0(i,"-01-01T00:00:00"),
#     end_datetime=paste0(i,"-12-31T23:59:00"),
#     minimum_depth=0.49402499198913574,
#     maximum_depth=2000,
#     output_filename = paste0("data/glorys/bgc/raw/glorys_bgc_",i,".nc"),
#     force_download = TRUE
#   )
#   toc()
# }
## END