library(tidyverse)
library(sf)
library(tidync)
library(terra)
library(here)
library(tictoc)

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

## DOWNLOAD DATA FROM MOM6

if(!dir.exists(here('data','mom6','raw'))) {dir.create(here('data','mom6','raw'),recursive = T)}

# Table of variables
library(jsonlite)
cefi_dims <- read_json("https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northeast_pacific.full_domain.hindcast.json") |> 
  bind_rows(.id="data")
# monthly grid variables
cefi_monthly_vars <- cefi_dims |> 
  filter(cefi_output_frequency=="monthly")
# "search"
cefi_monthly_vars$cefi_long_name |> str_subset("Sea Surface Temperature")
# this looks right; pull the correct url
sst_url <- cefi_monthly_vars |> filter(cefi_long_name=="Sea Surface Temperature") |> pull(cefi_opendap)
# seems like the first one is the one we want
sst_url <- sst_url[4]
sst_nc <- tidync(sst_url)
sst_nc

#     4 dimensions:
# time  Size:384   *** is unlimited *** 
#   units: days since 1993-01-01 00:00:00
# long_name: time
# axis: T
# calendar_type: GREGORIAN
# calendar: gregorian
# bounds: time_bnds
# _ChunkSizes: 100
# nv  Size:2 
# long_name: vertex number
# _ChunkSizes: 2
# xh  Size:342 
# units: degrees_east
# long_name: h point nominal longitude
# axis: X
# _ChunkSizes: 200
# yh  Size:816 
# units: degrees_north
# long_name: h point nominal latitude
# axis: Y
# _ChunkSizes: 200

# dimensions tables with indices (i.e., dimension indices mapped to dimension values)
sst_dimtbls <- hyper_transforms(sst_nc)
xh_id_df <- sst_dimtbls$xh
yh_id_df <- sst_dimtbls$yh
xh_yh_sf <- crossing(x=xh_id_df$xh,y=yh_id_df$yh) |> 
  st_as_sf(coords=c('x','y'),crs=4326)

mom6_zones <- feeding_zones |>
  st_filter(st_union(xh_yh_sf))

pull_sst <- function(y,zoneID,saverast=T){
  
  bb <- zones_bbs[[zoneID]]
  #find the right indices to pull
  xh_id <- which(between(xh_id_df$xh,bb[1],bb[3])) #"latitudes" to pull
  yh_id <- which(between(yh_id_df$yh,bb[2],bb[4])) #'longitudes' to pull

  tic(paste("Pulling sst for",y))
  time_inds= which(year(sst_dimtbls$time$timestamp)==y) #time slices to pull (there should be 12 in a year for monthly vars)
  
  sst_pull <- sst_nc %>% 
    hyper_tibble(xh= index %in% xh_id,
                 yh= index %in% yh_id,
                 time = index %in% time_inds)
  
  sst_sf <- sst_pull %>% 
    st_as_sf(coords=c("xh","yh"),crs=4326,remove=F)
  
  sst_sf
}

sst93 <- pull_sst(1993,zoneID=zoneID)
glimpse(sst93)
sst93 |> 
  filter(time=="1993-01-16T12:00:00") |> 
  ggplot()+geom_sf(aes(col=tos))+
  geom_sf(data=zones |> filter(Name=="SEA+NBC"),fill=NA)+
  scale_color_gradient2()
