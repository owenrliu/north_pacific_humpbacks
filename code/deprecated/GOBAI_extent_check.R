# Check GOBAI extent

library(here)
library(terra)

dat <- sds(here('data','gobaiO2','GOBAI-O2-clim-v2.3.nc'))
nlyr(dat)
# 4 variables of 696 layers each
# oxygen, temperature, salinity

plot(dat$temp[[1]])
plot(dat$oxy[[1]])
plot(dat$uncer[[1]])

# zoom in on the NEP
nep_crop <- ext(225,250,25,60)

nep_test <- dat$oxy[[1]] |> crop(nep_crop)
plot(nep_test)

library(rnaturalearth)
bg <- ne_countries(continent = "North America",scale = 'medium') |> 
  vect() |> crop(nep_test) |> rotate(left=F)
  
plot(nep_test);polys(bg)
