# Go Get GLORYS data and download it, using python workaround
# from https://github.com/pepijn-devries/CopernicusMarine/issues/42#issuecomment-2080179599
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(viridis)
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

# bounding boxes for humpback areas, updated with assessment scenarios
zones <- read_sf(here('data','spatial','NPhump_zones_scenarios.shp'))
zones_bbs <- map(zones$geometry,st_bbox) %>% 
  set_names(zones$zoneID)

# quick plot (B2F2 scenario)
zones %>%
  filter(B2|F2) %>% 
  ggplot()+
  geom_sf(aes(fill=factor(Feeding)))+
  labs(fill="Feeding\nGround")+
  theme_classic()+
  geom_sf_text(aes(label=Name))

# feeding zones only, because we are only going to download data for feeding grounds for now
feeding_zones <- zones %>% 
  filter(Feeding)

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
if(!dir.exists(here('data','glorys','bgc','raw'))) {dir.create(here('data','glorys','bgc','raw'),recursive = T)}

# static dataset with bathymetry
dat_static <- cm$subset(
  dataset_id="cmems_mod_glo_phy_my_0.083deg_static",
  minimum_longitude=118,
  maximum_longitude=281,
  minimum_latitude=5,
  maximum_latitude=72,
  output_filename = paste0("data/glorys/phys/raw/glorys_physics_statics.nc")
)

# for the 0.25deg biogeochemistry
dat_static_bgc <- cm$subset(
  dataset_id="cmems_mod_glo_bgc_my_0.25deg_static",
  dataset_part="mask",
  minimum_longitude=118,
  maximum_longitude=281,
  minimum_latitude=5,
  maximum_latitude=72,
  output_filename = paste0("data/glorys/bgc/raw/glorys_bgc_statics.nc")
)

# a version to download all years for a chosen whale zone
# can choose physics or bgc, which will point to different GLORYS/CMEMS datasets

dl_glorys_zone_allyrs <- function(zoneID,
                                  which_model= "physics", # physics or bgc are the options
                                  vars){
  # find the spatial bounding box
  bbox <- zones_bbs %>% pluck(zoneID)
  
  minlon <- bbox$xmin
  maxlon <- bbox$xmax
  minlat <- bbox$ymin
  maxlat <- bbox$ymax
  
  # find and download the right dataset
  if(which_model=="physics"){
    
    tic(paste("Downloading GLORYS Physics:",zoneID))
    # for the physics hindcast, it is split into two different datasets,
    # with the breakpoint between June and July, 2021
      dat1 <- cm$subset(
        dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
        variables=as.list(vars),
        minimum_longitude=minlon,
        maximum_longitude=maxlon,
        minimum_latitude=minlat,
        maximum_latitude=maxlat,
        start_datetime=paste0(1993,"-01-01T00:00:00"),
        end_datetime=paste0(2021,"-06-01T00:00:00"),
        minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
        maximum_depth=0.49402499198913574,
        output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_1993-2021.nc")
      )
      dat2 <- cm$subset(
        dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1M-m",
        variables=as.list(vars),
        minimum_longitude=minlon,
        maximum_longitude=maxlon,
        minimum_latitude=minlat,
        maximum_latitude=maxlat,
        start_datetime=paste0(2021,"-07-01T00:00:00"),
        end_datetime=paste0(2024,"-12-01T00:00:00"),
        minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
        maximum_depth=0.49402499198913574,
        output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_2021-2024.nc")
      )
      toc()
    }
  
  # the biogeochemistry is also in two different datasets,
  # but the split is between December 2022 and Jan 2023
  if(which_model=="bgc"){
    
    tic(paste("Downloading GLORYS BGC:",zoneID))
    dat1 <- cm$subset(
      dataset_id="cmems_mod_glo_bgc_my_0.25deg_P1M-m",
      variables=as.list(vars),
      minimum_longitude=minlon,
      maximum_longitude=maxlon,
      minimum_latitude=minlat,
      maximum_latitude=maxlat,
      start_datetime=paste0(1993,"-01-01T00:00:00"),
      end_datetime=paste0(2022,"-12-01T00:00:00"),
      minimum_depth=0.51, # this is the shallowest GLORYS level (i.e., the surface)
      maximum_depth=0.51,
      output_filename = paste0("data/glorys/bgc/raw/glorys_bgc_",zoneID,"_1993-2022.nc")
    )
    dat2 <- cm$subset(
      dataset_id="cmems_mod_glo_bgc_myint_0.25deg_P1M-m",
      variables=as.list(vars),
      minimum_longitude=minlon,
      maximum_longitude=maxlon,
      minimum_latitude=minlat,
      maximum_latitude=maxlat,
      start_datetime=paste0(2023,"-01-01T00:00:00"),
      end_datetime=paste0(2024,"-12-01T00:00:00"),
      minimum_depth=0.51, # this is the shallowest GLORYS level (i.e., the surface)
      maximum_depth=0.51,
      output_filename = paste0("data/glorys/bgc/raw/glorys_bgc_",zoneID,"_2023-2024.nc")
    )
    toc()
    
  }
}

# test
dl_glorys_zone_allyrs(zoneID="Hawaii_3",which_model='bgc',vars=c("chl","no3","nppv"))

xbgcfn <- here('data','glorys','bgc','raw','glorys_bgc_Hawaii_3_1993-2022.nc')
# random year/month to plot
xym <- as.Date("2016-08-01")

# chl
xchl <- rast(xbgcfn,subds="chl")
xchl <- xchl[[time(xchl)==xym]]
chl_p <- ggplot()+geom_spatraster(data=xchl)+scale_fill_viridis()

# no3
xno3 <- rast(xbgcfn,subds="no3")
xno3 <- xno3[[time(xno3)==xym]]
no3_p <- ggplot()+geom_spatraster(data=xno3)+scale_fill_viridis()

# nppv
xnppv <- rast(xbgcfn,subds="nppv")
xnppv <- xnppv[[time(xnppv)==xym]]
nppv_p <- ggplot()+geom_spatraster(data=xnppv)+scale_fill_viridis()

cowplot::plot_grid(chl_p,no3_p,nppv_p,nrow=1)

# try downloading everything by zone
# for biogeochemisty
walk(zones$zoneID,\(z) dl_glorys_zone_allyrs(zoneID=z,which_model="bgc",vars=c("chl","no3","nppv")))
# this took about ~30 minutes

# for physics
walk(feeding_zones$zoneID,\(z) dl_glorys_zone_allyrs(zoneID=z,which_model="physics",vars=c("thetao","mlotst")))

### END

### OLDER 
#quick function to download mixed layer depth and temperature for the surface layer
# for a given bounding box and year
# dl_glorys_yr <- function(year,zoneID,
#                          datID="cmems_mod_glo_phy_my_0.083deg_P1M-m",
#                          # dataset id. monthly physics: cmems_mod_glo_phy_my_0.083deg_P1M-m
#                          vars){
#   
#   bbox <- zones_bbs %>% pluck(zoneID)
#   
#   minlon <- bbox$xmin
#   maxlon <- bbox$xmax
#   minlat <- bbox$ymin
#   maxlat <- bbox$ymax
#   
#   # determine output file name/path
#   cmems_mod_glo_phy_my_0.083deg_P1M-m
#   
#   tic(paste("Downloading GLORYS: ",year," ",zoneID))
#   dat <- cm$subset(
#     dataset_id=datID,
#     variables=as.list(vars), #mld and sst
#     minimum_longitude=minlon,
#     maximum_longitude=maxlon,
#     minimum_latitude=minlat,
#     maximum_latitude=maxlat,
#     start_datetime=paste0(year,"-01-01T00:00:00"),
#     end_datetime=paste0(year,"-12-31T23:59:00"),
#     minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
#     maximum_depth=0.49402499198913574,
#     output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_",year,".nc")
#   )
#   toc()
# }

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

# for 2021-2024
# dl_glorys_2021_2024_allyrs <- function(zoneID){
#   
#   bbox <- zones_360_bbs %>% pluck(zoneID)
#   
#   minlon <- bbox$xmin
#   maxlon <- bbox$xmax
#   minlat <- bbox$ymin
#   maxlat <- bbox$ymax
#   
#   tic(paste("Downloading GLORYS: ",zoneID))
#   dat <- cm$subset(
#     dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1M-m",
#     variables=list("mlotst","thetao"),
#     minimum_longitude=minlon,
#     maximum_longitude=maxlon,
#     minimum_latitude=minlat,
#     maximum_latitude=maxlat,
#     start_datetime=paste0(2021,"-07-01T00:00:00"),
#     end_datetime=paste0(2024,"-12-31T23:59:00"),
#     minimum_depth=0.49402499198913574, # this is the shallowest GLORYS level (i.e., the surface)
#     maximum_depth=0.49402499198913574,
#     output_filename = paste0("data/glorys/phys/raw/glorys_physics_",zoneID,"_2021-2024.nc")
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