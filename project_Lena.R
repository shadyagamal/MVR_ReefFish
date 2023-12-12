library(reshape2)
library(dplyr)
library(vegan)

setwd("~/EPFL/MA3/MS/project")

# read csv file
data <- read.csv('fishes+nearest_env.csv')

# subsetting environmental variables 
env <- data %>%
  group_by(site_code) %>%
  summarise(CHL=mean(CHL_me_002), DHW= mean(DHW_me_002), FE=mean(FE_me_002),
          O2=mean(O2_me_002), PH=mean(PH_me_002), NO3=mean(NO3_me_002),
          PO4 = mean(PO4_me_002))
env[is.na(env)] = 0
rownames(env) <- env$site_code
env <- env[,-1]

env <- data %>%
  group_by(site_code) %>%
  summarise(CHL=mean(CHL_me_002), DHW= mean(DHW_me_002), FE=mean(FE_me_002),
            O2=mean(O2_me_002), PH=mean(PH_me_002), NO3=mean(NO3_me_002),
            PO4 = mean(PO4_me_002), SCV=mean(SCV_me_002), SPM=mean(SPM_me_002), 
            SSS=mean(SSS_me_002), SST=mean(SST_me_002), BATHY=mean(BATHY_X_01),
            CROP=mean(CROP_me_01), LAND=mean(LAND_X_010), PDEN=mean(PDEN_me_01),
            URBA=mean(URBA_me_01), VBD=mean(VBD_me_010)
            )
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

#subsetting spatial coordinates
spatial <- data %>%
  group_by(site_code) %>%
  summarise(X=mean(longitude), Y=mean(latitude), depth=mean(depth))
rownames(spatial) <- spatial$site_code
spatial <- spatial[,-1]

# subsetting species with abundance > 15%
occurrence_percentage <- colMeans(spe > 0) * 100
selected_species <- colnames(spe)[occurrence_percentage > 15]
spe_15 <- spe[ ,selected_species ]

# subsetting fish families
spe.family <- data %>%
  group_by(site_code, site_name, family) %>%
  summarise(abundance = sum(total, na.rm = TRUE)) %>%
  dcast(site_code ~ family, value.var="abundance")
spe.family[is.na(spe.family)] = 0
spe.family <- spe.family[,-1]  
# matrix visualization
heatmap(as.matrix(spe.family), dendrogram = "none", trace="none")

# subsetting species families with abundance > 15%
occurrence_percentage <- colMeans(spe.family > 0) * 100
selected_species <- colnames(spe.family)[occurrence_percentage > 15]
spe.family.15 <- spe.family[ ,selected_species ]

# subsetting fish biomass
spe.mass <- data %>%
  group_by(site_code, site_name, species_na) %>%
  summarise(biomass = sum(biomass, na.rm = TRUE)) %>%
  dcast(site_code ~ species_na, value.var="biomass")
spe.mass[is.na(spe.mass)] = 0
rownames(spe.mass) <- spe.mass$site_code
spe.mass <- spe.mass[,-1]  

# Number of absences
sum(spe == 0)
# Proportion of zeros in the community data set
sum(spe == 0) / (nrow(spe) * ncol(spe))

# Barplot of the distribution, all species families confounded
apply(spe_50, 2, range)
ab <- table(unlist(spe.family.15))
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

# Changing points transparancy
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "n",
         main = "PCoA with species weighted averages")
point_color1 <- rgb(0, 0, 1, alpha = 0.5)
point_color2 <- rgb(1, 0, 1, alpha = 0.4)
points(scores(spe.b.pcoa), display = "sites", pch = 20, col = point_color1, cex = 2)
points(spe.wa, display = "sites",pch = 17, cex = 1.5, col = point_color2)
legend("topright", legend = c("Species","Sites"), pch = c(20,17), col = c(point_color1, point_color2), cex = 1, bg = "white", box.lwd = 0)

# Changing symbol for depth variable
spatial$depth_cat <- cut(spatial$depth, breaks = c(0, 5, 10, 20), labels = c(0, 1, 2), include.lowest = TRUE)
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "n",
         main = "PCoA with species weighted averages")
points(scores(spe.b.pcoa), display = "sites", pch = as.numeric(spatial$depth_cat)+14, col = as.numeric(spatial$depth_cat)+4, alpha=0.5, cex = 1)
legend("topright", legend = c("depth < 5", "5 < depth < 10", "depth > 10"), pch = c(15,16,17), col = c(5,6,7), cex = 1, bg = "white", box.lwd = 0)
plot(spe.b.pcoa.env, p.max = 0.05, col = 3)



# ======================= Correspondence analysis (CA) ========================

# Compute CA
spe.ca <- cca(spe)
summary(spe.ca)		# default scaling 2
summary(spe.ca, scaling = 1)

# Scree plot and broken stick model using vegan's screeplot.cca()
screeplot(spe.ca, bstick = TRUE, npcs = length(spe.ca$CA$eig))

# CA biplots
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of species
plot(spe.ca, 
     scaling = 1, 
     main = "CA fish abundances - biplot scaling 1"
)
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main = "CA fish abundances - biplot scaling 2")


# ======================= DCA on species =======================

?decorana
spe.dca<-decorana(spe, iweigh=0, iresc=4, ira=0, mk=26, short=0, before=NULL, after=NULL)
summary(spe.dca)
spe.dca.env <- envfit(spe.dca, env)

par(mfrow = c(1, 2))
plot(spe.ca, 
     scaling = 1, 
     main = "CA fish abundances"
)

plot(spe.dca, 
     scaling = 1, 
     main = "DCA fish abundances"
)
plot(spe.dca.env, p.max = 0.001, col = 3)

# ======================= Mantel correlogram =======================
# on fish species with abundance > 15%

spe_15.hel <- decostand(spe_15, "hel")
spe.correlog <-mantel.correlog(spe.hel, XY=spatial, n.class=10, r.type="pearson", nperm=999, cutoff=FALSE)
spe.correlog
plot(spe.correlog)
