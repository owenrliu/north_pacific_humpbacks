# podaac-data-subscriber

# library(tidyverse)
# library(sf)
# library(terra)
# library(tidyterra)
# library(viridis)
# library(here)
# library(reticulate)
# library(tidync)
# library(viridis)
# library(tictoc)

# install_python() use this instead of the next line if this is the first time using reticulate/python
# install_python()
# use_python("C:/Users/owenrliu/AppData/Local/Programs/Python/Python313/python.exe")

# virtualenv_create(envname = "thispodaac")
# # virtualenv_install("thispodaac", packages = c("earthaccess"))
# use_virtualenv("thispodaac", required = TRUE)
# # 
# earthaccess <- reticulate::import("earthaccess") 
# auth = earthaccess$login(strategy="netrc")
# 
# granules <- earthaccess$search_data(
#   doi = "10.5067/SLREF-CDRV3",
#   temporal = reticulate::tuple("2017-01", "2017-02") # with an earthaccess update, this can be simply c() or list()
# )

# system("podaac-data-downloader -c MUR-JPL-L4-GLOB-v4.1 -d ./data --start-date 2002-05-31T21:00:00Z --end-date 2002-06-07T21:00:00Z -e '' ")

# PO.DAAC bulk-data downloader [-h] -c COLLECTION -d OUTPUTDIRECTORY [--cycle SEARCH_CYCLES] [-sd STARTDATE] [-ed ENDDATE] [-f] [-b BBOX] [-dc] [-dydoy] [-dymd] [-dy]
# [--offset OFFSET] [-e EXTENSIONS] [-gr GRANULENAME] [--process PROCESS_CMD] [--version] [--verbose] [-p PROVIDER] [--limit LIMIT] [--dry-run]

# Using R
# https://nasa-openscapes.github.io/earthdata-cloud-cookbook/how-tos/find-data/find-r.html
# install.packages("pak")
# pak::pak("boettiger-lab/earthdatalogin")
# pak::pak("rstac")

### AUTHENTICATION HAS BEEN ANNOYING ###
## try using environmental variables: https://earthaccess.readthedocs.io/en/latest/user_guide/authenticate/#login-using-environment-variables


library(rstac)
library(earthdatalogin)
library(tidyverse)
library(here)

edl_netrc(netrc_path = "C:/Users/owenrliu/_netrc",
          username = "XXXXX",
          password="XXXXX") # fill credentials

stac_source <- rstac::stac(
  "https://cmr.earthdata.nasa.gov/stac/POCLOUD/"
)
stac_source
get_request(stac_source)

collections_query <- stac_source %>% 
  rstac::collections()
collections_query
get_request(collections_query)

# this line takes a while to run.
allc <- collections_query %>%
  rstac::get_request() %>%
  earthdatalogin::collections_fetch()

allcc <- allc$collections

# available_collections
ids <- map_chr(allcc,"id")
descriptions <- map_chr(allcc,"description")
sst_ids <- which(grepl("MUR",ids)) #407 is GHRSST 0.1 degree, 408 is 0.25deg, 409 is 0.01deg
ids[sst_ids]
descriptions[sst_ids]

t <- rstac::stac_search(
  q = stac_source,
  collections = ids[408],
  datetime = "2018-02-12T00:00:00Z/2018-03-18T12:31:12Z",
  # collections="MUR-JPL-L4-GLOB-v4.1",
  # collections="MUR25-JPL-L4-GLOB-v04.2_4.2",
  limit=10
) %>% 
  get_request()

find_urls <- function(item_collection){
  feats <- item_collection %>% pluck("features")
  assets <- map(feats,"assets")
  
  # Direct Download?
  # is_dl <- map_lgl(assets,\(x) grepl("Download",x[[1]]$title))
  urls <- map_chr(assets,\(x)x[[1]]$href)
  
  # via s3
  
  # is_dl <- map_lgl(assets,\(x) grepl("Download",x[[1]]$title))
  # urls <- map_chr(assets,\(x)x[["s3_2"]]$href)
  
  urls
}

tu <- find_urls(t)
ncdf4::nc_open(tu[1])

tr <- terra::rast(tu[1],vsi=T) #doesn't work
download.file(tu[1],destfile=here::here('data','podaac','test1.nc'))
# tr <- terra::rast(here::here('data','podaac','test.nc'))
terra::plot(tr,1)

# ##
# 
# allc <- allc$collections
# echo "machine urs.earthdata.nasa.gov login <username> password <password>" >> %HOME%\_netrc
# # Get shortnames from the 'id' field and search for a match:
# ids <- map_chr(allc,"id")
# descriptions <- map_chr(allc,"description")
# sst_ids <- which(grepl("MUR",ids)) #407 is GHRSST 0.1 degree, 408 is 0.25deg, 409 is 0.01deg
# ids[sst_ids]
# descriptions[sst_ids]
# 
# start <- "2015-01-01T00:00:00Z"
# end   <- "2015-03-31T00:00:00Z" 
# trange <- paste(start,end, sep = "/")
# 

# testing with tutorial https://stacspec.org/en/tutorials/1-download-data-using-r/
stac_source <- rstac::stac(
  "https://planetarycomputer.microsoft.com/api/stac/v1"
)
stac_source
collections_query <- stac_source |>
  rstac::collections()

collections_query
available_collections <- rstac::get_request(collections_query)
available_collections

rstac::stac_search(
  q = stac_source,
  collections = "usgs-lcmap-conus-v13",
  datetime = "2021-01-01/2021-12-31",
  limit = 999
)

stac_query <- rstac::stac_search(
  q = stac_source,
  collections = "usgs-lcmap-conus-v13",
  bbox = c(-81.74091,36.23448, -81.23970,36.58977),
  datetime = "2021-01-01/2021-12-31"
)

executed_stac_query <- rstac::get_request(stac_query)
executed_stac_query

rstac::assets_download(executed_stac_query, "lcpri", output_dir = tempdir())

output_file <-file.path(
  tempdir(),
  "lcmap",
  "CU",
  "V13",
  "025011",
  "2021",
  "LCMAP_CU_025011_2021_20220721_V13_CCDC",
  "LCMAP_CU_025011_2021_20220629_V13_LCPRI.tif"
) |>
  terra::rast()

### try again
library(terra)
s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
it_obj <- s_obj %>%
  stac_search(collections = "landsat-c2-l2",
              bbox = c(-47.02148, -17.35063, -42.53906, -12.98314)) %>%
  get_request()

print(it_obj)
it_obj <- s_obj %>%
  stac_search(collections = "landsat-c2-l2",
              bbox = c(-47.02148, -17.35063, -42.53906, -12.98314)) %>%
  get_request() %>%
  items_sign(sign_fn = sign_planetary_computer())
url <- paste0("/vsicurl/", it_obj$features[[1]]$assets$blue$href)

data <- rast(url)
plot(data)

# try to get SST from planetary computer
# https://planetarycomputer.microsoft.com/dataset/noaa-cdr-sea-surface-temperature-whoi#Example-Notebook

s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
it_obj <- s_obj %>%
  stac_search(collections = "noaa-cdr-sea-surface-temperature-whoi",
              datetime="1988-01-01T22:00:00Z") %>%
  get_request()%>%
  items_sign(sign_fn = sign_planetary_computer())

print(it_obj)

url <- paste0("/vsicurl/", it_obj$features[[1]]$assets$sea_surface_temperature$href)

data <- rast(url) # wow this worked
plot(data)


### EXAMPLE cmd line
## AFTER edl_netrc is run
y <- 2002
m <- 9
startdate <- ym(paste(y,m,sep="-")) |> paste0("T00:00:00Z")
enddate <- as_date(startdate)+months(1)
enddate <- enddate |> paste0("T00:00:00Z")

s <- paste("podaac-data-downloader -c MUR25-JPL-L4-GLOB-v04.2 -d ./data/podaac/temp --start-date", startdate,"--end-date",enddate)
system(s)



minimum_longitude=118
maximum_longitude=281
minimum_latitude=5
maximum_latitude=72

npext <- ext(c(minimum_longitude,maximum_longitude,minimum_latitude,maximum_latitude))

# crop extent and combine
fls <- list.files(here('data','podaac','temp'),full.names = T) |> str_subset(".nc")
outr <- map(fls,\(m){
  rast(m,subds="analysed_sst") |> 
    crop(npext)
}) |> rast() |> 
  mean()

