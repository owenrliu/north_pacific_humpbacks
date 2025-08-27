library(rstac)
library(earthdatalogin)
library(tidyverse)
library(sf)
library(terra)
library(here)
library(tictoc)

# edl_netrc(netrc_path = "C:/Users/owenrliu/_netrc",
#           username = "XXXXX",
#           password="XXXXX") # fill credentials

edl_netrc() # authenticate

### EXAMPLE cmd line
## AFTER edl_netrc is run
## Dataset: https://podaac.jpl.nasa.gov/dataset/MUR25-JPL-L4-GLOB-v04.2#

# spatial domain
zones <- read_sf(here('data','spatial','NPhump_zones_scenarios.shp'))
npext <- ext(zones)

y <- 2002
m <- 9
startdate <- ym(paste(y,m,sep="-")) |> paste0("T00:00:00Z")
enddate <- as_date(startdate)+months(1)
enddate <- enddate |> paste0("T00:00:00Z")

s <- paste("podaac-data-downloader -c MUR25-JPL-L4-GLOB-v04.2 -d ./data/podaac/temp --start-date", startdate,"--end-date",enddate)
# do the download
system(s)

# crop extent and combine
fls <- list.files(here('data','podaac','temp'),full.names = T)
ncfls <- fls |> str_subset(".nc")
outr <- map(fls,\(m){
  rast(m,subds="analysed_sst") |> 
    rotate() |> 
    crop(npext)
}) |> rast() |> 
  mean()
plot(outr) # yay this works!
unlink(fls)

# write a function and apply to all year-months since 2002 September (start of the dataset)
download_crop_podaac_sst <- function(y,m){
  tic(paste("Downloading/ Year:",y,"Month:",m))
  startdate <- ym(paste(y,m,sep="-")) |> paste0("T00:00:00Z")
  enddate <- as_date(startdate)+months(1)
  enddate <- enddate |> paste0("T00:00:00Z")
  
  s <- paste("podaac-data-downloader -c MUR25-JPL-L4-GLOB-v04.2 -d ./data/podaac/temp --start-date", startdate,"--end-date",enddate)
  
  # do the download
  system(s)
  
  # crop extent and combine
  fls <- list.files(here('data','podaac','temp'),full.names = T)
  ncfls <- fls |> str_subset(".nc")
  outr <- map(ncfls,\(m){
    rast(m,subds="analysed_sst") |> 
      rotate() |> 
      crop(npext)
  }) |> rast() |> 
    mean()
  writeRaster(outr,here('data','podaac',paste0("GHRSST_0.25deg_",y,"_",m,".tif")),overwrite=T)
  unlink(fls)
  toc()
}

# download_crop_podaac_sst(y=2015,m=3)
months_to_download <- crossing(year=2002:2024,month=1:12) |> 
  filter(!(year==2002&month<9))

walk2(months_to_download$year,months_to_download$month,download_crop_podaac_sst)

# test-2015
fls15 <- list.files(here('data','podaac'),full.names=T) |> str_subset("2005")
fls15 <- c(fls15[1],fls15[5:12],fls15[2:4])
r15 <- rast(fls15)
plot(r15)
