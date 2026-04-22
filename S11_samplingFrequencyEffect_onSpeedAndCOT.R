
#____________________________________________________________________
# TEST EFFECT OF SAMPLING FREQUENCY ON SPEED ESTIMATES AND COT MODEL
#____________________________________________________________________

library(phylolm)
library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(fields)
library(ggridges)
library(rr2) #R2 for phylo models

#______________________________________________________________
# Explain MR and COT ~ body mass, separate for soarers/flappers
#______________________________________________________________


# folder where the data files from the repository were downloaded and where the results from script 8 were saved
setwd("...")

# Import model dataset with COT calculations and new soar/flap categories
allSegmDfs <- readRDS("finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")
# Import the majority-rule consensus (MRC) of the 1000 Ericson trees
tree <- read.nexus("MRCtree_DendroPy_from1000_Ericson_Feb2025_pruned.tre")

#_________________________________________________
# Explore sampling frequency and speed per segment and effects on species' average

mod_seg <- lm(avg_grSpeed_ms~min_timeLag_min, allSegmDfs)
summary(mod_seg)

hist(allSegmDfs$avg_timeLag_min)
hist(allSegmDfs$min_timeLag_min)

allSegmDfs_HR <- allSegmDfs[allSegmDfs$min_timeLag_min <= 20,]

length(unique(allSegmDfs$species)) - length(unique(allSegmDfs_HR$species))
unique(allSegmDfs$species)[!unique(allSegmDfs$species)%in%unique(allSegmDfs_HR$species)]
length(unique(allSegmDfs$individualID)) - length(unique(allSegmDfs_HR$individualID))

allSegmDfs_HHR <- allSegmDfs[allSegmDfs$min_timeLag_min <= 15 & allSegmDfs$avg_timeLag_min <= 15,]
excludedSegms <- allSegmDfs[!(allSegmDfs$min_timeLag_min <= 15 & allSegmDfs$avg_timeLag_min <= 15),]

nrow(allSegmDfs) - nrow(allSegmDfs_HHR)
nrow(excludedSegms)

(nrow(excludedSegms)/nrow(allSegmDfs))*100 # as percentage of total n segments removed

length(unique(allSegmDfs$species)) - length(unique(allSegmDfs_HHR$species)) # n lost species and individuals
unique(allSegmDfs$species)[!unique(allSegmDfs$species)%in%unique(allSegmDfs_HHR$species)]
length(unique(excludedSegms$individualID))
table(unique(excludedSegms$individualID) %in% unique(allSegmDfs_HHR$individualID))
# 4 species and 76 individuals get dropped entirely if we only consider high-res data.

# compare mean speed between species that had both low and high sampling freq
specSpeed_LR <- group_by(excludedSegms, species) %>% summarise(meanCCspeed = mean(avg_grSpeed_ms),
                                               minCCspeed = min(avg_grSpeed_ms),
                                               maxCCspeed = max(avg_grSpeed_ms),
                                               samplFreq = unique(min(min_timeLag_min)),
                                               n_segm = n())
specSpeed_HR <- group_by(allSegmDfs_HHR, species) %>% summarise(meanCCspeed = mean(avg_grSpeed_ms),
                                                            minCCspeed = min(avg_grSpeed_ms),
                                                            maxCCspeed = max(avg_grSpeed_ms),
                                                            samplFreq = unique(min(min_timeLag_min)),
                                                            n_segm = n())
summary(specSpeed_LR$n_segm)
summary(specSpeed_HR$n_segm)
specSpeed_LR[specSpeed_LR$n_segm<5,]
# drop species with less than 5 segments
specSpeed_LR <- specSpeed_LR[specSpeed_LR$n_segm>5,]
specSpeed_HR <- specSpeed_HR[specSpeed_HR$n_segm>5,]
# Ensure the same set of species
common_species <- intersect(specSpeed_HR$species, specSpeed_LR$species)
specSpeed_HR <- specSpeed_HR[specSpeed_HR$species %in% common_species, ]
specSpeed_LR <- specSpeed_LR[specSpeed_LR$species %in% common_species, ]
# Sort both datasets by species name (alphabetically)
specSpeed_HR <- specSpeed_HR[order(specSpeed_HR$species), ]
specSpeed_LR <- specSpeed_LR[order(specSpeed_LR$species), ]
all(specSpeed_HR$species == specSpeed_LR$species)

cor(specSpeed_LR$maxCCspeed,specSpeed_HR$maxCCspeed)
cor(specSpeed_LR$meanCCspeed,specSpeed_HR$meanCCspeed)
plot(specSpeed_LR$meanCCspeed,specSpeed_HR$meanCCspeed)
hist(specSpeed_LR$meanCCspeed - specSpeed_HR$meanCCspeed)


#______________________________________________
# Calculate species' mean of individuals' mean 
# for the complete and filtered dataset

data <- allSegmDfs

data_HR <- allSegmDfs_HHR
summary(data_HR$min_timeLag_min)
summary(data_HR$avg_timeLag_min)

## Full dataset "data"
# Transform dependent and predictor to log scale
data$Y_sp <- log(data$avg_grSpeed_ms)
data$Y <- log(data$obs_tot_COT_J_KgM) # per unit of mass and distance, COT in Joules/(kg m)
data$X1 <- log(data$Body_mass_kg)
data$WL <- log(data$wingArea_ellipse_cm2)
data$X2 <- factor(data$soarFlap_pgls, levels=c("flap","soar"))
# individuals' means
data <- group_by(data, individualID) %>%
  summarise(Y = mean(Y, na.rm=T),
            Y_sp = mean(Y_sp, na.rm=T),
            X1 = unique(X1),
            X2 = unique(X2),
            WL = unique(WL),
            Z = unique(species_phy),
            n_segm = n())
# species' means
data <- group_by(data, Z) %>%
  summarise(Y = mean(Y, na.rm=T),
            Y_sp = mean(Y_sp, na.rm=T),
            X1 = unique(X1),
            X2 = unique(X2),
            WL = unique(WL),
            n_segm = sum(n_segm)) %>%
  as.data.frame()

# check if any species has < 5 segments
summary(data$n_segm)

# Match tree to data and data to tree
attr(data,"na.action") <- NULL
anyNA(data)
rownames(data) <- gsub(" ", "_", data$Z) # replace space with _ for matching
setdiff(rownames(data), tree$tip.label)  # Species in data but not in tree
setdiff(tree$tip.label, rownames(data))  # Species in tree but not in data
data <- data[rownames(data) %in% tree$tip.label,] # filter out bats
tree <- keep.tip(tree, rownames(data))
data <- data[tree$tip.label, ]
identical(rownames(data), tree$tip.label)

## Filtered dataset with sampl freq <= 15 min "data_HR"
# Transform dependent and predictor to log scale
data_HR$Y_sp <- log(data_HR$avg_grSpeed_ms)
data_HR$Y <- log(data_HR$obs_tot_COT_J_KgM) # per unit of mass and distance, COT in Joules/(kg m)
data_HR$X1 <- log(data_HR$Body_mass_kg)
data_HR$WL <- log(data_HR$wingLoading_kgm2)
data_HR$WA <- log(data_HR$wingArea_ellipse_cm2)
data_HR$X2 <- factor(data_HR$soarFlap_pgls, levels=c("flap","soar"))
# individuals' means
data_HR <- group_by(data_HR, individualID) %>%
  summarise(Y = mean(Y, na.rm=T),
            Y_sp = mean(Y_sp, na.rm=T),
            X1 = unique(X1),
            X2 = unique(X2),
            WL = unique(WL),
            WA = unique(WA),
            Z = unique(species_phy),
            n_segm = n())
# species' means
data_HR <- group_by(data_HR, Z) %>%
  summarise(Y = mean(Y, na.rm=T),
            Y_sp = mean(Y_sp, na.rm=T),
            X1 = unique(X1),
            X2 = unique(X2),
            WL = unique(WL),
            WA = unique(WA),
            n_segm = sum(n_segm)) %>%
  as.data.frame()

# check if any species has < 5 segments
summary(data_HR$n_segm)

# Match tree to data_HR and data_HR to tree
attr(data_HR,"na.action") <- NULL
anyNA(data_HR)
rownames(data_HR) <- gsub(" ", "_", data_HR$Z) # replace space with _ for matching
setdiff(rownames(data_HR), tree$tip.label)  # Species in data_HR but not in tree
setdiff(tree$tip.label, rownames(data_HR))  # Species in tree but not in data_HR
data_HR <- data_HR[rownames(data_HR) %in% tree$tip.label,] # filter out bats
tree_HR <- keep.tip(tree, rownames(data_HR))
data_HR <- data_HR[tree_HR$tip.label, ]
identical(rownames(data_HR), tree_HR$tip.label)

## Correlation with wing morphology
# par(mfrow=c(2,1))
# plot(data_HR$X1, data_HR$WL)
# plot(data_HR$X1, data_HR$WA)
# cor(data_HR$X1, data_HR$WL)
# cor(data_HR$X1, data_HR$WA)

#_______________________________
# Run model cc_speed ~ body mass
# Note that the interaction effect X1:X2soar represents the shift in slope between soarers and flappers. The main effect X2soar represents the shift in intercept.

# First with a simple linear model:
# a. with complete dataset
summary(lm(Y_sp ~ X1 * X2, data))
# vflap = 7.52 Mass^(0.10); vsoar = 4.98 Mass^(0.22); Adj R2 0.432

# b. with high res dataset
summary(lm(Y_sp ~ X1 * X2, data_HR))
# vflap = 7.56 Mass^(0.12); vsoar = 5.05 Mass^(0.23); Adj R2 0.386

# Then with a phylolm:
# a. with complete dataset
summary(phylolm(Y_sp ~ X1 * X2, data, phy = tree, model="lambda"))
# vflap = 7.52 Mass^(0.10); vsoar = 4.98 Mass^(0.22); Adj R2 0.432

# b. with high res dataset
summary(phylolm(Y_sp ~ X1 * X2, data_HR, phy = tree_HR, model="lambda"))
# vflap = 7.56 Mass^(0.12); vsoar = 5.05 Mass^(0.23); Adj R2 0.386

# Plot of phylo model predicting the allometry of cc-speed on full dataset

plot_full <- ggplot(data, aes(x=X1, y=Y_sp, col=X2)) +
  xlab("Body mass as log(x in kg)") + ylab("Cross-country speed (m/s)") +
  coord_cartesian(xlim = c(-1, 2.8), ylim = c(1.09, 2.7)) +
  scale_y_continuous(
    labels = function(x) round(exp(x), 1),
    breaks = log(c(3, 6, 9, 12, 15))) +  # choose breaks on real scale, log-transform them
  # Add regression lines from your model
  geom_abline(intercept = 2.017353, slope = 0.102671, color = "gold2") +  # flappers
  geom_abline(intercept = 2.017353 - 0.411328,
              slope = 0.102671 + 0.118794,
              color = "dodgerblue2") +  # soarers
  # Add species labels
  geom_text(aes(label = Z),
            hjust = -0.1, vjust = 0.5,
            size = 3, show.legend = FALSE) +
  geom_point() + scale_color_manual(values=c("gold2","dodgerblue2")) +
  theme_classic(base_size = 9) +
  theme(
    plot.title       = element_text(size = 9, face = "plain", hjust = 0, margin = margin(b = 4)),
    axis.title       = element_text(size = 8.5),
    axis.text        = element_text(size = 7.5, colour = "black"),
    axis.line        = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.35),
    legend.position  = "none",   # collected by patchwork below
    plot.margin      = margin(6, 10, 4, 4)
  )

plot_HR <- ggplot(data_HR, aes(x=X1, y=Y_sp, col=X2)) +
  xlab("Body mass as log(x in kg)") + ylab("Cross-country speed (m/s)") + # log scale
  coord_cartesian(xlim = c(-1, 2.8), ylim = c(1.09, 2.7)) +
  scale_y_continuous(
    labels = function(x) round(exp(x), 1),
    breaks = log(c(3, 6, 9, 12, 15))) +
    #breaks = log(c(3, 5, 9, 15.5))) +  # choose breaks on real scale, log-transform them
  # Add regression lines from your model
  geom_abline(intercept = 2.022854, slope = 0.123005, color = "gold2") +  # flappers
  geom_abline(intercept = 2.102601 - 0.402502,
              slope = 0.123005 + 0.108576,
              color = "dodgerblue2") +  # soarers
  # Add species labels
  geom_text(aes(label = Z),
            hjust = -0.1, vjust = 0.5,
            size = 3, show.legend = FALSE) +
  geom_point() + scale_color_manual(values=c("gold2","dodgerblue2")) +
  theme_classic(base_size = 9) +
  theme(
    plot.title       = element_text(size = 9, face = "plain", hjust = 0, margin = margin(b = 4)),
    axis.title       = element_text(size = 8.5),
    axis.text        = element_text(size = 7.5, colour = "black"),
    axis.line        = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.35),
    legend.position  = "none",   # collected by patchwork below
    plot.margin      = margin(6, 10, 4, 4)
  )

library(patchwork)
pdf("Revision/NewSupplFigures/phylolm_ccspeed_bodymass_highRes.pdf", width = 18,height = 6)
(plot_full | plot_HR) + plot_layout(guides = "collect") & theme(legend.position = "none")
dev.off()


#___________________________
# Run model COT ~ body mass

# a. on the full dataset
model_cot <- phylolm(Y ~ X1 * X2, 
                 data, phy=tree, model="lambda") # Pagel's lambda model
(summ <- summary(model_cot))
hist(residuals(model_cot))
R2_lik(model_cot)
# eCOT flap = 4.89 M^(−0.33); eCOT soar = 2.30 M^(−0.70)
# R2 Adj 0.809; R2 lik 0.905
# only flappers have a slightly steeper relationship in the full dataset

# For an 11 kg species (e.g. cignus buccinator 11.07 and vultur gryphus 11.24 kg)
cot_swan = 4.89 * 11^(-0.33) # 2.22 J / kg m for flapper HR model
cot_vult = 2.30 * 11^(-0.70) # 0.43 J / kg m for soarer HR model
(cot_swan - cot_vult)/cot_swan # 80% saving of vulture relative to swan

summ[["coefficients"]]
intFlap <- summ[["coefficients"]]["(Intercept)",1] # model intercept for flappers
intSoar <- intFlap + summ[["coefficients"]]["X2soar",1] # model intercept for swimmers
slopeFlap <- summ[["coefficients"]]["X1",1] # effect of body mass on runners
slopeSoar <- slopeFlap + summ[["coefficients"]]["X1:X2soar",1] # effect of body mass on swimmers

modelCoefficients <- cbind(intFlap,intSoar,slopeFlap,slopeSoar)

# b. on the high resolution data subset
model_cot_hr <- phylolm(Y ~ X1 * X2, 
                     data_HR, phy=tree_HR, model="lambda") # Pagel's lambda model
(summ_hr <- summary(model_cot_hr))
hist(residuals(model_cot_hr))
R2_lik(model_cot_hr)
# eCOT flap = 5.45 M^(−0.47); eCOT soar = 2.29 M^(−0.74)
# R2 Adj 0.697; R2 lik 0.857

# For an 11 kg species (e.g. cignus buccinator 11.07 and vultur gryphus 11.24 kg)
cot_swan = 5.45 * 11^(-0.47) # 1.77 J / kg m for flapper HR model
cot_vult = 2.29 * 11^(-0.74) # 0.39 J / kg m for soarer HR model
(cot_swan - cot_vult)/cot_swan # 78% saving of vulture relative to swan

summ_hr[["coefficients"]]
intFlap <- summ_hr[["coefficients"]]["(Intercept)",1] # model intercept for flappers
intSoar <- intFlap + summ_hr[["coefficients"]]["X2soar",1] # model intercept for swimmers
slopeFlap <- summ_hr[["coefficients"]]["X1",1] # effect of body mass on runners
slopeSoar <- slopeFlap + summ_hr[["coefficients"]]["X1:X2soar",1] # effect of body mass on swimmers

modelCoefficients <- cbind(intFlap,intSoar,slopeFlap,slopeSoar)

# Compute body mass ranges per group from data
flap_xrange <- range(data$X1[data$X2 == "flap"])
soar_xrange <- range(data$X1[data$X2 == "soar"])


plot_full_COT <- ggplot(data, aes(x=X1, y=Y, col=X2)) +
  xlab("Body mass as log(x in kg)") + ylab("Effective COT \n log(geometric mean of X in J kg⁻¹ m⁻¹)") + 
  xlim(c(-1, 2.8)) + coord_cartesian(ylim = c(-1.3, 3)) +
  scale_y_continuous(
    labels = function(x) round(exp(x), 1),
    breaks = log(c(4, 5, 7, 10, 15))  # choose breaks on real scale, log-transform them
  ) + # tick labels on real unit but plotted on log scale
  scale_x_continuous(
    labels = function(x) round(exp(x), 1),
    breaks = log(c(4, 5, 7, 10, 15))  # choose breaks on real scale, log-transform them
  ) +
  # Add regression lines manually from model
  geom_abline(intercept = 1.587056, slope = -0.328447, color = "goldenrod2") +  # flappers
  geom_abline(intercept = 1.587056 - 0.328447,
              slope = -0.328447 -0.376582,
              color = "dodgerblue2") +  # soarers
  
  # Draw segments instead of ablines to respect the x-y range of the observations
  geom_segment(aes(
    x    = flap_xrange[1],
    xend = flap_xrange[2],
    y    = int_flap + slp_flap * flap_xrange[1],
    yend = int_flap + slp_flap * flap_xrange[2]),
    colour = "dodgerblue2", linewidth = 0.7, inherit.aes = FALSE) +
  geom_segment(aes(
    x    = soar_xrange[1],
    xend = soar_xrange[2],
    y    = int_soar + slp_soar * soar_xrange[1],
    yend = int_soar + slp_soar * soar_xrange[2]), 
    colour = "goldenrod2", linewidth = 0.7, inherit.aes = FALSE) +
  # Add species labels
  geom_text(aes(label = Z),
            hjust = -0.1, vjust = 0.5,
            size = 3, show.legend = FALSE) +
  geom_point() + scale_color_manual(values=c("goldenrod2","dodgerblue2")) +
  theme_classic(base_size = 9) +
  theme(
    plot.title       = element_text(size = 9, face = "plain", hjust = 0, margin = margin(b = 4)),
    axis.title       = element_text(size = 8.5),
    axis.text        = element_text(size = 7.5, colour = "black"),
    axis.line        = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.35),
    legend.position  = "none",   # collected by patchwork below
    plot.margin      = margin(6, 10, 4, 4)
  )

plot_HR_COT <- ggplot(data_HR, aes(x=X1, y=Y, col=X2)) +
  xlab("Body mass as log(x in kg)") + ylab("Effective COT \n log(geometric mean of X in J kg⁻¹ m⁻¹)") + 
  xlim(c(-1, 2.8)) +
  # Add regression lines from your model
  geom_abline(intercept = 1.694600, slope = -0.469935, color = "gold2") +  # flappers
  geom_abline(intercept = 1.694600 - 0.865445,
              slope = -0.469935 -0.265913,
              color = "dodgerblue2") +  # soarers
  # Add species labels
  geom_text(aes(label = Z),
            hjust = -0.1, vjust = 0.5,
            size = 3, show.legend = FALSE) +
  geom_point() + scale_color_manual(values=c("gold2","dodgerblue2")) +
  theme_classic(base_size = 9) +
  theme(
    plot.title       = element_text(size = 9, face = "plain", hjust = 0, margin = margin(b = 4)),
    axis.title       = element_text(size = 8.5),
    axis.text        = element_text(size = 7.5, colour = "black"),
    axis.line        = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.35),
    legend.position  = "none",   # collected by patchwork below
    plot.margin      = margin(6, 10, 4, 4)
  )

library(patchwork)
pdf("Revision/NewSupplFigures/phylolm_ccspeedCOT_bodymass_highRes_panelsABCD.pdf", width = 18,height = 12)
combined_plot <- 
  (plot_full + plot_HR) /
  (plot_full_COT + plot_HR_COT)
combined_plot +
  plot_annotation(tag_levels = "A")
dev.off()

