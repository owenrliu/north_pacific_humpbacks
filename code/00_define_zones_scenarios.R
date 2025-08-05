## Define all of the unique zones and combinations of zones in the 2024 NP Humpback Assessment scenarios (B1,B2,B3,B4)

library(tidyverse)
library(sf)
library(here)

# bounding boxes for humpback areas
zones <- read_sf(here('data','spatial','NPhuwhRegions.1Apr2024.shp'))

zones %>% st_set_geometry(NULL)
# manually add the zones that are combinations of other zones (in some scenarios)
mexico <- zones %>% filter(Name %in% c("MX_ML","MX_AR","MX_BJ")) %>%
  summarise() %>% 
  mutate(Name="Mexico",id="17")
goa <- zones %>% filter(Name %in% c("WGOA","NGOA")) %>% 
  summarise() %>% 
  mutate(Name="GOA",id="18")
ealberwgoa <- zones %>% filter(Name %in% c("WGOA","EAL+BER")) %>% 
  st_shift_longitude() %>%  # crossing the dateline issue
  summarise() %>% 
  mutate(Name="EAL+BER+WGOA",id="19")
zones_scen <- zones %>% 
  bind_rows(mexico) %>% 
  bind_rows(goa) %>% 
  bind_rows(ealberwgoa)
zones_scen %>% st_set_geometry(NULL)

# Scenarios mapped to names and polygon numbers
B1 <- c(1,4,7,9,3,17,2,14) %>% as.character()
B2 <- c(1,4,7,9,3,5,6,2) %>% as.character()
F1 <- c(15,16,18,11,12,13) %>% as.character()
F2 <- c(15,19,8,11,12,13) %>% as.character()
feeding_grounds <- c(8,10:16,18,19) %>% as.character()

zones_scen <- zones_scen %>% 
  mutate(B1=ifelse(id %in% B1,T,F)) %>% 
  mutate(B2=ifelse(id %in% B2,T,F)) %>% 
  mutate(F1=ifelse(id %in% F1,T,F)) %>% 
  mutate(F2=ifelse(id %in% F2,T,F)) %>% 
  mutate(Feeding=ifelse(id %in% feeding_grounds,T,F))
zones_scen %>% st_set_geometry(NULL)

# need to shift the longitude so it doesn't plot backwards across the intl dateline
# (i.e. move to the 0-360 longitude instead of -180-180)
zones_360 <- zones_scen %>% st_shift_longitude() %>% unite("zoneID",Name,id,remove = F)

# quick plot, B1F1
zones_360 %>% 
  filter(B1|F1) %>% 
  ggplot()+
  geom_sf(aes(fill=factor(Feeding)))+
  labs(fill="Feeding\nGround")+
  theme_classic()+
  geom_sf_text(aes(label=Name))

# save this updated version
write_sf(zones_360,here('data','spatial','NPhump_zones_scenarios.shp'))
