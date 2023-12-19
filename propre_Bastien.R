# working directory -------------------------------------------------------
setwd("C:/Users/Bastien Amez-Droz/Documents/EPFL/Ma3/Multivariate stats in R/project/github/")


# Install, load packages and data -----------------------------------------
library(vegan)
library(ade4)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(forcats)
library(tidyr)

# read data from file -----------------------------------------------------
data <- read.csv('fishes+nearest_env.csv')

# subsetting environmental variables
env <- data %>%
  group_by(site_code) %>%
  summarise(DEPTH=mean(depth), CHL=mean(CHL_me_002), DHW= mean(DHW_me_002), FE=mean(FE_me_002),
            O2=mean(O2_me_002), PH=mean(PH_me_002), NO3=mean(NO3_me_002),
            PO4 = mean(PO4_me_002), SCV=mean(SCV_me_002), SPM=mean(SPM_me_002), 
            SSS=mean(SSS_me_002), SST=mean(SST_me_002), BATHY=mean(BATHY_X_01),
            CROP=mean(CROP_me_01), LAND=mean(LAND_X_010), PDEN=mean(PDEN_me_01),
            URBA=mean(URBA_me_01), VBD=mean(VBD_me_010))
#env[is.na(env)] = 0 
rownames(env) <- env$site_code


# ============= Filter the sites :=====================
# remove sites with same environmental characteristics

env <- env %>%
  distinct(CHL, DHW, FE, O2, PH, NO3, PO4, SCV, SPM, SSS, SST, BATHY, CROP, LAND, PDEN, URBA, VBD, .keep_all = TRUE)

# ============= Filter the species=====================

# subsetting fish abundances
spe <- data %>%
  group_by(site_code, site_name, species_na) %>%
  summarise(abundance = sum(total, na.rm = TRUE)) %>%
  dcast(site_code ~ species_na, value.var="abundance")
spe[is.na(spe)] = 0
rownames(spe) <- spe$site_code

# remove useless sites from species dataframe
spe <- spe[spe$site_code %in% env$site_code, ]

rownames(env) <- spe$site_code
spe <- spe[, -1]
env <- env[, -1]

# selecting species with abundance>15%
# Calculate the occurrence percentage for each species
occurrence_percentage <- colMeans(spe > 0) * 100
# Select species with occurrence > 15%
selected_species <- colnames(spe)[occurrence_percentage > 15]
# Subset the original data with selected species
spe <- spe[ ,selected_species ]

# distribution of abundance (reefFish) ---------------------------------------
# Minimum and maximum of abundance values in the whole data set
range(spe)

# Minimum and maximum value for each species
apply(spe, 2, range)

# Count the cases for each abundance class
ab <- table(unlist(spe))
ab

# Barplot of the distribution, all species confounded
barplot(ab, 
        las = 1,
        xlab = "Abundance class",
        ylab = "Frequency",
        col = gray(5 : 0 / 5),
        horiz=F,
)

# ================= BASTIEN =======================

#Pour certaines analyses, je dois enlever les lignes avec des NaN (y'en a pas beaucoup)
rownames(env) <- rownames(spe)
row_names <- rownames(env)
rows_with_nans <- which(rowSums(is.na(env)) > 0)
env.nonan <- env[-rows_with_nans, ]
spe.nonan<- spe[-rows_with_nans, ]
rownames(env.nonan) = rownames(spe.nonan)



## Explanatory variables selection

# Hellinger-transform the species dataset
spe.hel <- decostand(spe.nonan, "hellinger")

## RDA of the Hellinger-transformed fish species data, constrained
## by all (~.) environmental variables contained in env3
spe.rda <- rda(spe.hel ~ ., env.nonan)
summary(spe.rda)


#Forward selection of variables
mod0 <- rda(spe.hel ~ 1, data = env.nonan) # starting with unconstrained RDA
summary(mod0)
step.forward <- ordistep(mod0, scope = formula(spe.rda), direction = "forward", permutations = how(nperm = 999))
RsquareAdj(step.forward)


#RDA 
# Hellinger-transform the species dataset
spe.hel <- decostand(spe.nonan, "hellinger")

## RDA of the Hellinger-transformed fish species data, constrained
## Même pas toutes les variables sélectionnées par la forward selection
RDA.result <- rda(spe.hel ~ DEPTH + CHL + DHW + FE + O2 + PH + NO3 + PO4 + SCV + SPM + SSS + SST, data = env.nonan)
result.anova <- anova(RDA.result)

result.anova$Variance[1]/sum(result.anova$Variance)

sc_bp <- scores(RDA.result, display="bp", choices=c(1, 2), scaling=1)


#data about which species are coral dependent
coral_dep <- read.csv('coral_dependent.csv', sep=";")
species_names <- colnames(spe.hel)
colors <- ifelse(coral_dep$dep[match(species_names, coral_dep$species)] == 1, "deeppink", "blue")
#Pink means coral dependent

# Plot RDA with colored fishes by dependency
plot(RDA.result,type="n")
#text(RDA.result, "sites", col="red", cex=0.5)
points(RDA.result, pch=21, col="black", cex=1)
points(RDA.result, "species", pch=19, col=colors, cex=1)

#Code ci-dessous pour ajouter les flèches des variables environnementales

# arrows(0,0, # start them from (0,0)
#        sc_bp[,1], sc_bp[,2], # end them at the score value
#        col = "red", 
#        lwd = 1)
# text(x = sc_bp[,1] -0.05, # adjust text coordinate to avoid overlap with arrow tip
#      y = sc_bp[,2] - 0.03, 
#      labels = rownames(sc_bp), 
#      col = "red", 
#      cex = 0.5, 
#      font = 0.8)




# GRAPHES AVEC FACTOEXTRA ========

library(InPosition)
library(factoextra)

#Compute PCA
res.pca<-princomp(env.nonan,cor=TRUE)
biplot(res.pca,xlabs=rownames(env.nonan))

#Plot des valeurs propres. Show the percentage of variances explained by each principal component.
fviz_eig(res.pca)

#PCA sur les sites, colors by quality of representation
fviz_pca_ind(res.pca,col.ind = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE,geom="point")

#Axes des variables environmentales sur les sites en PCA
fviz_pca_var(res.pca,col.var = "contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

#Combinaison des 2
fviz_pca_biplot(res.pca, repel = TRUE,col.var = "#2E9FDF",col.ind = "#696969", geom="point")


#k-means pour faire des clusters avec les variables 
env.nonan.std <- scale(env.nonan) 
kmeans.result<-kmeans(env.nonan.std,2)
groups <- as.factor(kmeans.result$cluster)
fviz_pca_ind(res.pca,col.ind = groups,palette = c("#00AFBB",  "#FC4E07","#E7B800"),addEllipses = TRUE, ellipse.type = "convex", legend.title = "Groups",repel = TRUE, geom="point")

#Analyse supp par clusters ? Pas sûr...
