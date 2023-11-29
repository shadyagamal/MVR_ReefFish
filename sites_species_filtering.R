library(reshape2)
library(dplyr)
library(vegan)

setwd(dir = "C:/Users/clem7/EPFL/MA3/MS_R")

# read csv file
data <- read.csv('fishes+nearest_env.csv')

# subsetting environmental variables

env <- data %>%
  group_by(site_code) %>%
  summarise(CHL=mean(CHL_me_002), DHW= mean(DHW_me_002), FE=mean(FE_me_002),
            O2=mean(O2_me_002), PH=mean(PH_me_002), NO3=mean(NO3_me_002),
            PO4 = mean(PO4_me_002), SCV=mean(SCV_me_002), SPM=mean(SPM_me_002), 
            SSS=mean(SSS_me_002), SST=mean(SST_me_002), BATHY=mean(BATHY_X_01),
            CROP=mean(CROP_me_01), LAND=mean(LAND_X_010), PDEN=mean(PDEN_me_01),
            URBA=mean(URBA_me_01), VBD=mean(VBD_me_010))
env[is.na(env)] = 0
rownames(env) <- env$site_code


# ============= Filter the sites : 
# remove sites with same environmental characteristics

env <- env %>%
  distinct(CHL, DHW, FE, O2, PH, NO3, PO4, SCV, SPM, SSS, SST, BATHY, CROP, LAND, PDEN, URBA, VBD, .keep_all = TRUE)

# ============= Filter the species

# subsetting fish abundances
spe <- data %>%
  group_by(site_code, site_name, species_na) %>%
  summarise(abundance = sum(total, na.rm = TRUE)) %>%
  dcast(site_code ~ species_na, value.var="abundance")
spe[is.na(spe)] = 0
rownames(spe) <- spe$site_code

# remove useless sites from species dataframe
spe <- spe[spe$site_code %in% env$site_code, ]

# selecting species with abundance>15%
# Calculate the occurrence percentage for each species
occurrence_percentage <- colMeans(spe > 0) * 100
# Select species with occurrence > 15%
selected_species <- colnames(spe)[occurrence_percentage > 15]
# Subset the original data with selected species
spe <- spe[ ,selected_species ]
spe <- spe[, -1]


env <- env[, -1]
rownames(env) <- spe$site_code

# ============= Principal coordinate analysis (PCoA) ============================

spe.bray <- vegdist(spe)
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe) - 1), eig = TRUE)
# Plot of the sites
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "n",
         main = "PCoA with species weighted averages")

# Ajouter des points aux coordonnées des sites
points(scores(spe.b.pcoa, choices = c(1, 2)),
       col = "blue",  # Couleur des points
       pch = 16)  # Type de point (utilisez pch = 16 pour des cercles pleins, mais vous pouvez ajuster selon vos préférences)



