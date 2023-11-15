library(reshape2)
library(dplyr)
library(vegan)

setwd("~/EPFL/MA3/MS/project")

# read csv file
data <- read.csv('fishes+nearest_env.csv')

# subsetting environmental variables 
env <- data %>%
  group_by(site_code, site_name) %>%
  summarise(CHL=mean(CHL_me_002), DHW= mean(DHW_me_002), FE=mean(FE_me_002),
          O2=mean(O2_me_002), PH=mean(PH_me_002), NO3=mean(NO3_me_002),
          PO4 = mean(PO4_me_002))
env[is.na(env)] = 0
rownames(env) <- env$site_code
env <- env[,-1]

# subsetting fish abundances
spe <- data %>%
  group_by(site_code, site_name, species_na) %>%
  summarise(abundance = sum(total, na.rm = TRUE)) %>%
  dcast(site_code ~ species_na, value.var="abundance")
spe[is.na(spe)] = 0
rownames(spe) <- spe$site_code
spe <- spe[,-1]

# subsetting fish biomass
spe.mass <- data %>%
  group_by(site_code, site_name, species_na) %>%
  summarise(biomass = sum(biomass, na.rm = TRUE)) %>%
  dcast(site_code ~ species_na, value.var="biomass")
spe.mass[is.na(spe.mass)] = 0
rownames(spe.mass) <- spe.mass$site_code
spe.mass <- spe.mass[,-1]

# selecting species with abundance>15%
# Calculate the occurrence percentage for each species
occurrence_percentage <- colMeans(spe > 0) * 100
# Select species with occurrence > 15%
selected_species <- colnames(spe)[occurrence_percentage > 15]
# Subset the original data with selected species
spe_15 <- spe[ ,selected_species ]

# Number of absences
sum(spe == 0)
# Proportion of zeros in the community data set
sum(spe == 0) / (nrow(spe) * ncol(spe))

# subsetting fish families
spe.family <- data %>%
  group_by(site_code, site_name, family) %>%
  summarise(abundance = sum(total, na.rm = TRUE)) %>%
  dcast(site_code ~ family, value.var="abundance")
spe.family[is.na(spe.family)] = 0
spe.family <- spe.family[,-1]  
# matrix visualization
heatmap(as.matrix(spe.family), dendrogram = "none", trace="none")
  

# Barplot of the distribution, all species families confounded
apply(spe_15, 2, range)
ab <- table(unlist(spe_15))
barplot(ab, 
        las = 1,
        xlab = "Abundance class",
        ylab = "Frequency",
        col = gray(5 : 0 / 5),
        horiz=F,
)

# test - classification
spe <- spe[,-1]
spe.norm <- decostand(spe, "normalize") # normalize date
spe.ch <- vegdist(spe.norm, "euc") # eucledian distance
heatmap(as.matrix(spe.ch))

spe.ch.ward <- hclust(spe.ch, method = "ward.D2")
plot(spe.ch.ward,  main = "Chord - Ward")

# Diversity index
div <- diversity(spe, index = "shannon")

# Rarefaction
S <- specnumber(spe) # observed number of species
(raremax <- min(rowSums(spe)))
Srare <- rarefy(spe, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(spe, step = 20, sample = raremax, col = "blue", cex = 0.6)


# Correspondence analysis (CA) ====================================

# Compute CA
env.standard <- scale(env)
env.ca <- cca(env.standard)
summary(env.ca)		# default scaling 2
summary(env.ca, scaling = 1)

# Scree plot and broken stick model using vegan's screeplot.cca()
screeplot(env.ca, bstick = TRUE, npcs = length(env.ca$CA$eig))

# CA biplots
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of species
plot(env.ca, 
     scaling = 1, 
     main = "CA fish abundances - biplot scaling 1"
)
# Scaling 2 (default): species are centroids of sites
plot(env.ca, main = "CA fish abundances - biplot scaling 2")


# ============= Principal coordinate analysis (PCoA) ============================

spe.bray <- vegdist(spe)
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe) - 1), eig = TRUE)
# Plot of the sites
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")

# project weighted average projection of species
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe)
text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")

# A posteriori projection of environmental variables
spe.b.pcoa.env <- envfit(spe.b.pcoa, env)
spe.b.pcoa.env

# Plot significant variables with a user-selected colour
plot(spe.b.pcoa.env, p.max = 0.05, col = 3)

