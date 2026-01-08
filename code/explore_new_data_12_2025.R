# Explore updated abundance data
# and compare to old data
# Dec 2025
library(tidyverse)
library(here)
library(viridis)
fls <- list.files(here('data','reprogressonhumpbacks'),full.names = T)
fls <- fls[grepl(".csv",fls)]

# MIXING DATA
# Feeding to wintering (mixing)
mix_f2b <- read_csv(fls[1])|> 
  mutate(Direction="FeedingtoBreeding")
# wintering to feeding (mixing)
mix_b2f <- read_csv(fls[4])|> 
  mutate(Direction="BreedingtoFeeding")

# OLD mixing data
oldmix <- read_csv(here('data',"Mark-Recapture_mixing_data_allscenarios_table_Long.csv")) |> 
  filter(Hypothesis=="B2F1") |> 
  dplyr::select(Feeding,Breeding,Direction,Estimate)
# change names and combine
# NOT SURE HOW TO COMBINE feed2breed from Kamchatka+WBer
mixnew <- bind_rows(mix_f2b,mix_b2f)

names_ref <- tibble(cheeseman_name=unique(c(mixnew$Wintering_Region,mixnew$Feeding_Region))) |> 
  mutate(assess_name=c("Asia","Hawaii","MX_AR","MX_ML","Central_Am","RUS+WAL","RUS+WAL","EAL+BER","GOA","SEA+NBC","SBC+WA","OR+CA"))

mixnew <- mixnew |> 
  left_join(names_ref,by=join_by(Wintering_Region==cheeseman_name)) |> 
  rename(assessment_breed_name=assess_name) |> 
  left_join(names_ref,by=join_by(Feeding_Region==cheeseman_name)) |> 
  rename(assessment_feed_name=assess_name) |> 
  # combine Kamchatka and WBer
  group_by(assessment_breed_name,assessment_feed_name,Direction) |> 
  summarise(new_est=sum(Overlap_Percentage),.groups = 'drop')

  # join the old data
mixall <- mixnew |> 
  left_join(oldmix,by=join_by(Direction,assessment_breed_name==Breeding,assessment_feed_name==Feeding)) |> 
  rename(old_est=Estimate) |> 
  mutate(old_est=old_est*100)

old_mat <- mixall |> 
  filter(Direction=="BreedingtoFeeding") |> 
  ggplot(aes(assessment_breed_name,assessment_feed_name,fill=old_est))+
  geom_tile()+
  geom_text(aes(label=old_est))+
  scale_fill_gradient(low = "white", high = "red",limits=c(0,100))+
  # scale_fill_viridis(option="D")+
  labs(x="Breed",y="Feed",fill="")
old_mat
new_mat <- mixall |> 
  filter(Direction=="BreedingtoFeeding") |> 
  ggplot(aes(assessment_breed_name,assessment_feed_name,fill=new_est))+
  geom_tile()+
  geom_text(aes(label=round(new_est,1)))+
  scale_fill_gradient(low = "white", high = "red",limits=c(0,100))+
  # scale_fill_viridis(option='D')+
  labs(x="Breed",y="Feed",fill="")
new_mat
cowplot::plot_grid(old_mat,new_mat,nrow=2,labels=c("Old","New"))

# ABUNDANCE DATA
# NPac total abundance, 
# year by year for 2002-2024 with the modified Chapman Peterson estimate
tot <- read_csv(fls[2])
tot |> 
  mutate(ymax=BiasCorrAbund+2*SE,ymin=BiasCorrAbund-2*SE) |> 
  ggplot(aes(MidYear,BiasCorrAbund,ymax=ymax,ymin=ymin))+
  geom_line()+geom_pointrange(linetype=2)

# Hawaii relative abundance estimate
rel <- read_csv(fls[3])
rel |> 
  ggplot()+
  geom_line(aes(MidYear,Abundance,ymax=U_CI,ymin=L_CI))+
  geom_pointrange(aes(MidYear,Abundance,ymax=U_CI,ymin=L_CI),linetype=2)

# OLD abundance data
# B2F1BC Survey
old <- read_csv(here('Diags','SurveyUse','B2F1BC.csv'))
old_HI <- old |>
  filter(Area=="Hawaii",Year2 %in% unique(rel$MidYear)) |> 
  dplyr::select(Year2,Estimate,CV)
HI_combined <- rel |>
  left_join(old_HI,by=join_by("MidYear"=="Year2"))

# compare
HI_combined |> 
  ggplot(aes(x=MidYear))+
  geom_line(aes(MidYear,Abundance))+
  geom_pointrange(aes(MidYear,Abundance,ymax=U_CI,ymin=L_CI),linetype=2)+
  geom_point(aes(y=Estimate),color='red')
