library(reshape2)
library(dplyr)

setwd("~/EPFL/MA3/MS/project")

data <- read.csv('fishes+nearest_env.csv')

env <- data %>%
  group_by(site_code) %>%
  summarise(CHL=mean(CHL_me_002), DHW= mean(DHW_me_002), FE=mean(FE_me_002),
          O2=mean(O2_me_002), PH=mean(PH_me_002), NO3=mean(NO3_me_002),
          PO4 = mean(PO4_me_002))

spe <- data %>%
  group_by(site_code, site_name, species_na) %>%
  summarise(abundance = sum(total, na.rm = TRUE)) %>%
  dcast(site_code ~ species_na, value.var="abundance")
spe[is.na(spe)] = 0
