library(tidyverse)
library(rnaturalearth)
library(sf)

bb <- c(-175,23,-105,65)
names(bb) <- c("xmin","ymin","xmax","ymax")
out <- ne_countries(continent = "North America",scale=50,returnclass = 'sf') |> 
  st_crop(bb)
pout <- ggplot(out)+geom_sf(fill='tan')+theme_void()
ggsave(here::here('plots','coastline_background.png'),pout,h=10,w=10,bg='transparent')
