
#_____________________
# PHYLOGENETIC MODELS
#_____________________

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


# setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")

# Import model dataset with COT calculations and new soar/flap categories
allSegmDfs <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")
# Import the majority-rule consensus (MRC) of the 1000 Ericson trees
tree <- read.nexus("Phylogeny/Trees_species_2025/MRCtree_DendroPy_from1000_Ericson_Feb2025_pruned.tre")


#______________________________________________
# Calculate species' mean of individuals' mean

data <- allSegmDfs

# Transform dependent and predictor to log scale
data$Y_w <- log(data$tot_MR_W) # total flight power in Watts
data$Y_wkg <- log(data$tot_MR_W_kg) # flight power per unit of mass
data$Y <- log(data$obs_tot_COT_J_KgM) # per unit of mass and distance, COT in Joules/(kg m)
data$X1 <- log(data$Body_mass_kg)

length(unique(data$individualID))

# individuals' means
data <- group_by(data, individualID) %>%
  summarise(Y = mean(Y, na.rm=T),
            Y_w = mean(Y_w, na.rm=T),
            Y_wkg = mean(Y_wkg, na.rm=T),
            X1 = unique(X1),
            X2 = unique(soarFlap_pgls),
            Z = unique(species_phy))

# species' means
#se <- function(x) {sd(x) / sqrt(length(x))} # s.e. function
data <- group_by(data, Z) %>%
  summarise(Y = mean(Y, na.rm=T),
            Y_w = mean(Y_w, na.rm=T),
            Y_wkg = mean(Y_wkg, na.rm=T),
            X1 = unique(X1),
            X2 = unique(X2)) %>%
  as.data.frame()

# Plot difference between the power per kilo and the COT (the difference is driven by speed)
data$diff_W_COT <- data$Y_wkg - data$Y
plot(diff_W_COT~X1, data)
data[order(data$diff_W_COT, decreasing = T),c("diff_W_COT","Y_wkg","Y")]

# Match tree to data and data to tree
attr(data,"na.action") <- NULL
anyNA(data)
rownames(data) <- gsub(" ", "_", data$Z) # replace space with _ for matching
setdiff(rownames(data), tree$tip.label)  # Species in data but not in tree
setdiff(tree$tip.label, rownames(data))  # Species in tree but not in data

data_bats <- data[!rownames(data) %in% tree$tip.label,] # filter out bats
data <- data[rownames(data) %in% tree$tip.label,] # filter out bats
tree <- keep.tip(tree, rownames(data))

# Sort species in data by tip.label in tree
data <- data[tree$tip.label, ]
identical(rownames(data), tree$tip.label)

# flap soar as factor
data$X2 <- factor(data$X2, levels=c("flap","soar"))
table(data$X2) # only 15 flappers because the 2 bat species get excluded


#___________________________
# Run model tot_MR_W ~ body mass
# Note that the interaction effect X1:X2soar represents the shift in slope between soarers and flappers. The main effect X2soar represents the shift in intercept.
model_w <- phylolm(Y_w ~ X1 * X2, 
                 data, phy=tree, model="lambda") # Pagel's lambda model
(summ_w <- summary(model_w))
hist(residuals(model_w))
R2_lik(model)

# # only for soarers
# model_w_soar <- phylolm(Y_w ~ X1, 
#                    data[data$X2=="soar",], phy=tree, model="lambda") # Pagel's lambda model
# summary(model_w_soar)

# Extract coefficients for plotting
summ_w[["coefficients"]]
intFlap <- summ_w[["coefficients"]]["(Intercept)",1] # model intercept for flappers
intSoar <- intFlap + summ_w[["coefficients"]]["X2soar",1] # model intercept for swimmers
slopeFlap <- summ_w[["coefficients"]]["X1",1] # effect of body mass on runners
slopeSoar <- slopeFlap + summ_w[["coefficients"]]["X1:X2soar",1] # effect of body mass on swimmers

modelCoefficients_w <- cbind(intFlap,intSoar,slopeFlap,slopeSoar)

#___________________________
# Run model COT ~ body mass
model <- phylolm(Y ~ X1 * X2, 
                  data, phy=tree, model="lambda") # Pagel's lambda model
(summ <- summary(model))
hist(residuals(model))
R2_lik(model)

# Extract coefficients for plotting
summ[["coefficients"]]
intFlap <- summ[["coefficients"]]["(Intercept)",1] # model intercept for flappers
intSoar <- intFlap + summ[["coefficients"]]["X2soar",1] # model intercept for swimmers
slopeFlap <- summ[["coefficients"]]["X1",1] # effect of body mass on runners
slopeSoar <- slopeFlap + summ[["coefficients"]]["X1:X2soar",1] # effect of body mass on swimmers

modelCoefficients <- cbind(intFlap,intSoar,slopeFlap,slopeSoar)

# Save models
save(model_w, modelCoefficients_w, file="DataFinalSummary/phyloModels_FPW-bodymass_Feb2025.rdata")
save(model, modelCoefficients, file="DataFinalSummary/phyloModels_COT-bodymass_Feb2025.rdata")

#_______
# Plots 
 
birds_bats_abbNames <- rbind(data[,c("Y","Y_w","Y_wkg","X1","X2","Z")], data_bats[,c("Y","Y_w","Y_wkg","X1","X2","Z")])
birds_bats_abbNames$Z <- paste0(substr(birds_bats_abbNames$Z, 1,1),
                                ". ",
                                sapply(strsplit(birds_bats_abbNames$Z, " "),"[",2))

# Plots flight power MR ~ body mass
pdf("Plots/finalPlots/soarersFlappers_FPW-Wmodel_Feb2025_onlyGPS.pdf", 8,5)
plot(data$X1, data$Y_w, pch=19, type="n", xlim=c(-1.9,3.5), ylim=c(2,7),
     xlab="log(Body mass in Kg)", ylab="log(Flight power in W)")
abline(modelCoefficients_w[,"intFlap"], modelCoefficients_w[,"slopeFlap"], lty=1, lwd=6, col="gold2")
abline(modelCoefficients_w[,"intSoar"], modelCoefficients_w[,"slopeSoar"], lty=1, lwd=6, col="dodgerblue3")
points(data$X1, data$Y_w, col=alpha(c("gold2","dodgerblue3")[data$X2],0.6), pch=19)
points(data$X1, data$Y_w, col=c("gold2","dodgerblue3")[data$X2], pch=1)
points(data_bats$X1, data_bats$Y_w, pch=19, col=alpha("black",0.6))
points(data_bats$X1, data_bats$Y_w, pch=1, col="black")
text(Y_w~X1, birds_bats_abbNames, labels=birds_bats_abbNames$Z, pos=4, cex=.7, offset=0.5, srt=30, font=3)
dev.off()
# Clearly the osprey lays among the runners when it comes to flight power. See script 11 for an exploration of Osprey's flight over water.

# Plots COT ~ body mass, adding theoretical lines from Alexander
load("DataFinalSummary/thereticalCOTs_Nielsen-Alexander.RData") # Object theoreticalCOTalex

pdf("Plots/finalPlots/soarersFlappers_COTmodel_Feb2025_onlyGPS_theorLines.pdf", 8,5)
plot(data$X1, data$Y, pch=19, type="n", xlim=c(-1.9,3.1), ylim=c(-1,2.7),
     xlab="log(Body mass in Kg)", ylab="log(COT in J/kg m)")
lines(log(COTrun_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex, 
      lty=2, lwd=2, col=alpha("brown", 0.7))
lines(log(COTflight_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex,
      lty=2, lwd=2, col=alpha("forestgreen", 0.7))
lines(log(COTswim_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex,
      lty=2, lwd=2, col=alpha("dodgerblue3", 0.7))
abline(modelCoefficients[,"intFlap"], modelCoefficients[,"slopeFlap"], lty=1, lwd=6, col="gold2")
abline(modelCoefficients[,"intSoar"], modelCoefficients[,"slopeSoar"], lty=1, lwd=6, col="dodgerblue3")
points(data$X1, data$Y, col=alpha(c("gold2","dodgerblue3")[data$X2],0.6), pch=19)
points(data$X1, data$Y, col=c("gold2","dodgerblue3")[data$X2], pch=1)
points(data_bats$X1, data_bats$Y, pch=19, col=alpha("black",0.6))
points(data_bats$X1, data_bats$Y, pch=1, col="black")
text(Y~X1, birds_bats_abbNames, labels=birds_bats_abbNames$Z, pos=4, cex=.7, offset=0.5, srt=30, font=3)
dev.off()

#_________________________________________
# Calculate residual variation per species
# As the difference between the fitted COT value per species and the actual value of each segment

# Re-add the bats to the dataset and predict the expected COT
data_all <- rbind(data,data_bats)
data_all$Yfitted <- predict(model, newdata = data_all[,c("X1","X2")], type="response")
# check that bats are included!

# Merge the fitted species-specific COT value to the original segment data to calculate differences
allSegmDfs_mode <- merge(allSegmDfs, data_all[,c("Z","Yfitted")], by.x="species_phy", by.y = "Z")

allSegmDfs_mode$COT_residuals <- log(allSegmDfs_mode$obs_tot_COT_J_KgM) - allSegmDfs_mode$Yfitted

hist(allSegmDfs_mode$COT_residuals[allSegmDfs_mode$soarFlap_pgls=="flap"], col="gold2")
hist(allSegmDfs_mode$COT_residuals[allSegmDfs_mode$soarFlap_pgls=="soar"], col="dodgerblue3")
pdf("Plots/finalPlots/soarersFlappers_boxplotCOTresiduals_Feb2025.pdf", 6,5)
boxplot(COT_residuals~soarFlap_pgls, data=allSegmDfs_mode, col=alpha(c("gold2","dodgerblue3"), 0.6))
dev.off()

# Re-save the dataset with the additional columns Yfitted and COTresiduals
saveRDS(allSegmDfs_mode, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")

#____________________________________
## Density plots and residual models
#____________________________________

library(ggplot2)
library(ggridges)
library(scales)
library(gridExtra)
library(dplyr)

setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")
allSegmDfs_mode <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")

#_______________________________________
# Plot distribution of COT per species

# Add COT of running and swimming for reference in plot
allSegmDfs_mode$swimCOT_alex_J_KgM <- 1.1*(allSegmDfs_mode$Body_mass_kg)^(-0.38)
allSegmDfs_mode$runCOT_alex_J_KgM <- 10.7*(allSegmDfs_mode$Body_mass_kg)^(-0.32)

# Sort species by body mass and extract predicted COT per species and median observed COT
species_order_soar <- allSegmDfs_mode %>%
  filter(soarFlap_pgls == "soar") %>%
  group_by(species_phy) %>%
  summarize(median_cot = median(obs_tot_COT_J_KgM), 
            bodymass=unique(Body_mass_kg)) %>%
  #arrange(desc(median_cot)) %>%
  arrange(desc(bodymass)) %>%
  mutate(species_phy=as.character(species_phy)) %>%
  pull(species_phy)
species_order_flap <- allSegmDfs_mode %>%
  filter(soarFlap_pgls == "flap") %>%
  group_by(species_phy) %>%
  summarize(median_cot = median(obs_tot_COT_J_KgM), 
            bodymass=unique(Body_mass_kg)) %>%
  #arrange(desc(median_cot)) %>%
  arrange(desc(bodymass)) %>%
  mutate(species_phy=as.character(species_phy)) %>%
  pull(species_phy)

soarers <- allSegmDfs_mode %>% filter(soarFlap_pgls == "soar")
flappers <- allSegmDfs_mode %>% filter(soarFlap_pgls == "flap")
soarers$species_phy <- factor(soarers$species_phy, levels = species_order_soar) 
flappers$species_phy <- factor(flappers$species_phy, levels = species_order_flap)

summary(soarers$obs_tot_COT_J_KgM)
summary(flappers$obs_tot_COT_J_KgM)

soarYfitted <- allSegmDfs_mode %>%
  filter(soarFlap_pgls == "soar") %>%
  group_by(species_phy) %>%
  summarize(predictedCOT=unique(Yfitted),
            alexFlight=unique(theor_COT_alex_J_KgM),alexSwim=unique(swimCOT_alex_J_KgM),alexRun=unique(runCOT_alex_J_KgM)) %>%
  mutate(species_phy = factor(species_phy, levels = levels(soarers$species_phy)))
#mutate(y = as.numeric(factor(species_phy, levels = levels(soarers$species_phy))))
#soarYfitted$species_phy <- factor(soarYfitted$species_phy, levels = species_order_soar)
flapYfitted <- allSegmDfs_mode %>%
  filter(soarFlap_pgls == "flap") %>%
  group_by(species_phy) %>%
  summarize(predictedCOT=unique(Yfitted),
            alexFlight=unique(theor_COT_alex_J_KgM),alexSwim=unique(swimCOT_alex_J_KgM),alexRun=unique(runCOT_alex_J_KgM)) %>%
  mutate(y = as.numeric(factor(species_phy, levels = levels(soarers$species_phy))))
#flapYfitted$species_phy <- factor(flapYfitted$species_phy, levels = species_order_flap)
 

# separately for flappers and soarers
cot1 <- ggplot(soarers, aes(x = obs_tot_COT_J_KgM, y = species_phy, group = species_phy)) +
  xlim(-0.5, 30) +
  labs(x = "Effective COT (J/kg m)", y = "Species") +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, color = "black", fill = NA, linewidth=0.9) +
  #geom_vline(xintercept = mean(soarers$obs_tot_COT_J_KgM), linetype = "dashed", color = "black") +
  # Add three vertical segments per species indicating the run,swim.flight reference values from Alexander
  geom_segment(
    data = soarYfitted,
    aes(x = alexRun, xend = alexRun, y = as.numeric(factor(species_phy)) - 0.01, yend = as.numeric(factor(species_phy)) + 0.9),
    linetype = "solid", color = "gold2", inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(
    data = soarYfitted,
    aes(x = alexFlight, xend = alexFlight, y = as.numeric(factor(species_phy)) - 0.01, yend = as.numeric(factor(species_phy)) + 0.9),
    linetype = "solid", color = "magenta", inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(
    data = soarYfitted,
    aes(x = alexSwim, xend = alexSwim, y = as.numeric(factor(species_phy)) - 0.01, yend = as.numeric(factor(species_phy)) + 0.9),
    linetype = "solid", color = "dodgerblue3", inherit.aes = FALSE, linewidth = 0.6) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = alpha("white",0.7), color = NA),  # background color
    plot.background = element_rect(fill = alpha("dodgerblue3",0.2), color = NA))   # around the panel


cot2 <- ggplot(flappers, aes(x = obs_tot_COT_J_KgM, y = species_phy, group = species_phy)) +
  xlim(-0.5, 30) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, color = "black", fill = NA, linewidth=0.9) +
  labs(x = "Effective COT (J/kg m)", y = "Species") +
  #geom_vline(xintercept = mean(flappers$obs_tot_COT_J_KgM), linetype = "dashed", color = "black") +
  geom_segment(
    data = flapYfitted,
    aes(x = alexRun, xend = alexRun, y = as.numeric(factor(species_phy)) - 0.01, yend = as.numeric(factor(species_phy)) + 0.9),
    linetype = "solid", color = "gold2", inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(
    data = flapYfitted,
    aes(x = alexFlight, xend = alexFlight, y = as.numeric(factor(species_phy)) - 0.01, yend = as.numeric(factor(species_phy)) + 0.9),
    linetype = "solid", color = "magenta", inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(
    data = flapYfitted,
    aes(x = alexSwim, xend = alexSwim, y = as.numeric(factor(species_phy)) - 0.01, yend = as.numeric(factor(species_phy)) + 0.9),
    linetype = "solid", color = "dodgerblue3", inherit.aes = FALSE, linewidth = 0.6) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = alpha("white",0.7), color = NA),  # background color
        plot.background = element_rect(fill = alpha("gold2",0.2), color = NA))   # around the panel

pdf("Plots/finalPlots/soarersFlappers_COTdensitiesPerSpecies_Nov2025.pdf", width=8, height=5)
grid.arrange(cot1, cot2, ncol = 2)
dev.off()



#_______________________________
# Plot distribution of residuals

# Other colors options
# "#D35400""#3F51B5""#FFC107""#673AB7"
# "#E69F00", "#56B4E9", "#009E73", "#CC79A7"
# magma	"#000004", "#51127C", "#B63779", "#FB8861", "#FCFDBA" #	High-contrast data
# plasma	"#0D0887", "#7E03A8", "#CC4678", "#F89441", "#F0F921"	#Sequential data
# inferno	"#000004", "#56106E", "#BB3754", "#F98C09", "#F0F921"	Highlight extremes
# mako "#0B0405", "#1F2C5C", "#3C4F8A", "#5FC3A6","#FAFCCB"	Cool-toned gradients
# rocket	"#03051A", "#CB1B4F", "#F48952", "#FCCE2B", "#E3F6FC"	Warm/cool contrast
# library(colorspace)
# muted_magenta <- lighten("magenta", amount = 0.2)
# my_gradient <- c("#009688", "white", muted_magenta)

# separately for soarers and flappers
species_order_soar <- allSegmDfs_mode %>%
  filter(soarFlap_pgls == "soar") %>%
  group_by(species_phy) %>%
  summarize(median_residual = median(COT_residuals)) %>%
  arrange(desc(median_residual)) %>%
  pull(species_phy)
species_order_flap <- allSegmDfs_mode %>%
  filter(soarFlap_pgls == "flap") %>%
  group_by(species_phy) %>%
  summarize(median_residual = median(COT_residuals)) %>%
  arrange(desc(median_residual)) %>%
  pull(species_phy)

soarers <- allSegmDfs_mode %>% filter(soarFlap_pgls == "soar")
flappers <- allSegmDfs_mode %>% filter(soarFlap_pgls == "flap")
soarers$species_phy <- factor(soarers$species_phy, levels = species_order_soar) 
flappers$species_phy <- factor(flappers$species_phy, levels = species_order_flap) 

my_gradient <- c("#0D0887","#5FC3A6","white","#B63779","#51127C")


ggplot(allSegmDfs_mode, aes(x = COT_residuals, y = species_phy, fill = ..x.., group = species_phy)) +
  xlim(-3,3) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = my_gradient,
                       values = scales::rescale(c(-3, 0, 3)),limits = c(-3, 3),
                       name = "COT residuals") +
  labs(x = "COT Residuals", y = "Species", fill = "COT Residuals") +
  facet_wrap(~soarFlap_pgls) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70")
ggsave("Plots/finalPlots/all_residualsDensities_Feb2025_greenMag_facetWrap.pdf", width=5, height=6)

# separately for flappers and soarers
r1 <- ggplot(soarers, aes(x = COT_residuals, y = species_phy, fill = ..x.., group = species_phy)) +
  xlim(-3,3) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = my_gradient,
                       values = scales::rescale(c(-3, 0, 3)),limits = c(-3, 3),
                       name = "COT residuals") +
  labs(x = "COT Residuals", y = "Species", fill = "COT Residuals") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70")
  
r2 <- ggplot(flappers, aes(x = COT_residuals, y = species_phy, fill = ..x.., group = species_phy)) +
  xlim(-3,3) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = my_gradient,
    values = scales::rescale(c(-3, 0, 3)),limits = c(-3, 3),
    name = "COT residuals") +
  labs(x = "COT Residuals", y = "Species", fill = "COT Residuals") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70")

pdf("Plots/finalPlots/soarersFlappers_residualsDensities_Feb2025_greenMag.pdf", width=5, height=8)
grid.arrange(r1, r2, ncol = 1)
dev.off()

# Plot all species together
species_order <- allSegmDfs_mode %>%
  group_by(species_phy) %>%
  summarize(median_residual = median(COT_residuals)) %>%
  arrange(desc(median_residual)) %>%
  pull(species_phy)
allSegmDfs_mode$species_phy <- factor(allSegmDfs_mode$species_phy, levels = species_order) 

ggplot(allSegmDfs_mode, aes(x = COT_residuals, y = species_phy, fill = ..x.., group = species_phy)) +
  xlim(-3,3) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = my_gradient,
                       values = scales::rescale(c(-3, 0, 3)),limits = c(-3, 3),
                       name = "COT residuals") +
  labs(x = "COT Residuals", y = "Species", fill = "COT Residuals") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70")
ggsave("Plots/finalPlots/all_residualsDensities_Feb2025_greenMag.pdf", width=5, height=6)

# Plot distribution of residuals per group
ggplot(allSegmDfs_mode, aes(x = COT_residuals, y = soarFlap_pgls, 
                            fill = after_stat(x))) +
  geom_density_ridges_gradient(alpha = 0.7, scale=3) +
  scale_fill_gradientn(
    colours = my_gradient,
    limits = c(-3, 3),
    name = "COT residuals"
  ) +
  xlim(-3, 3) +
  theme_minimal() +
  # theme(panel.grid.major.x = element_blank(),
  #       panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  labs(y = "Density", x = "COT Residuals")
ggsave("Plots/finalPlots/gps_residualsDensities_soarFlap.pdf", width=3, height=3)


#______________________________________________________
# Plot distribution of atmospheric energy


ws <- ggplot(allSegmDfs_mode, aes(x = avg_windSupp_1000m, y = soarFlap_pgls, 
                            fill = after_stat(x))) +
  geom_density_ridges_gradient(alpha = 0.7, scale=3) +
  scale_fill_gradientn(
    colours = alpha("cyan",0.6),
    name = "Average Wind Support (m/s)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  labs(y = "Density", x = "Average Wind Support")
tu <- ggplot(allSegmDfs_mode, aes(x = avg_thermUplift, y = soarFlap_pgls, 
                                  fill = after_stat(x))) +
  geom_density_ridges_gradient(alpha = 0.7, scale=3) +
  scale_fill_gradientn(
    colours = alpha("orange",0.6),
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  labs(y = "Density", x = "Average Thermal Uplift")

pdf("Plots/finalPlots/atmosphericEnergy_soarFlap.pdf", width=8, height=5)
grid.arrange(ws, tu, ncol = 2)
dev.off()

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
group_by(allSegmDfs_mode, soarFlap_pgls) %>%
  summarise(meanWS=mean(avg_windSupp_1000m, na.rm=T), seWS=se(avg_windSupp_1000m), 
            meanTU=mean(avg_thermUplift, na.rm=T), seTU=se(avg_thermUplift))


#______________________________________________________
# Model residuals as a function of the energy landscape
# (that is the atmospheric variables encountered en route)

library(lme4)
library(MuMIn)

cor(allSegmDfs_mode[,c("avg_thermUplift","avg_heatFlux","avg_windSupp_1000m")], 
    use = "complete.obs")
# thermal uplift and sensible heat flux are highly correlated, we keep the thermal uplift

allSegmDfs_mode$segm_timeDuration_min <- allSegmDfs_mode$segm_timeDuration_h*60

# Linear mixed models separately for flappers and soarers
modRes_lm_flap <- lmer(COT_residuals ~ 
                         avg_windSupp_1000m +
                         avg_thermUplift +
                         segm_timeDuration_min +
                         avg_timeLag_min +
                         # avg_windSupp_1000m * segm_timeDuration_h +
                         # avg_thermUplift * segm_timeDuration_h +
                         (1 | individualID), # species accounted for by previous model
                       data = allSegmDfs_mode[allSegmDfs_mode$soarFlap_pgls == "flap",])
summary(modRes_lm_flap)
r.squaredGLMM(modRes_lm_flap)
hist(residuals(modRes_lm_flap))

modRes_lm_soar <- lmer(COT_residuals ~ 
                         avg_windSupp_1000m +
                         avg_thermUplift +
                         segm_timeDuration_min +
                         avg_timeLag_min +
                         # avg_windSupp_1000m * segm_timeDuration_h +
                         # avg_thermUplift * segm_timeDuration_h +
                         (1 | individualID), # species accounted for by previous model
                       data = allSegmDfs_mode[allSegmDfs_mode$soarFlap_pgls == "soar",])
summary(modRes_lm_soar)
r.squaredGLMM(modRes_lm_soar)
hist(residuals(modRes_lm_soar))

save(modRes_lm_flap, modRes_lm_soar, file="DataFinalSummary/models_COTresiduals-energyLandscape_Feb2025.rdata")


# Plot the effects
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the range of the covariates
windSupp_range <- seq(min(allSegmDfs_mode$avg_windSupp_1000m, na.rm = TRUE),
                      max(allSegmDfs_mode$avg_windSupp_1000m, na.rm = TRUE),
                      length.out = 20)
thermUplift_range <- seq(min(allSegmDfs_mode$avg_thermUplift, na.rm = TRUE),
                         max(allSegmDfs_mode$avg_thermUplift, na.rm = TRUE),
                         length.out = 20)

# Create a grid of predictor values
grid_flap <- expand.grid(avg_windSupp_1000m = windSupp_range,
                    avg_thermUplift = thermUplift_range,
                    segm_timeDuration_h = median(allSegmDfs_mode$segm_timeDuration_h[allSegmDfs_mode$soarFlap_pgls == "flap"]),
                    individualID = unique(allSegmDfs_mode$individualID[allSegmDfs_mode$soarFlap_pgls == "flap"]))
grid_soar <- expand.grid(avg_windSupp_1000m = windSupp_range,
                         avg_thermUplift = thermUplift_range,
                         segm_timeDuration_h = median(allSegmDfs_mode$segm_timeDuration_h[allSegmDfs_mode$soarFlap_pgls == "soar"]),
                         individualID = unique(allSegmDfs_mode$individualID[allSegmDfs_mode$soarFlap_pgls == "soar"]))

# Predict COT_residuals for the "flap" model
grid_flap$COT_residuals <- predict(modRes_lm_flap, newdata = grid_flap)

# Predict COT_residuals for the "soar" model
grid_soar$COT_residuals <- predict(modRes_lm_soar, newdata = grid_soar)


# Plot heatmaps
summary(allSegmDfs_mode$COT_residuals)

p1 <- ggplot(grid_flap, aes(x = avg_windSupp_1000m, y = avg_thermUplift, fill = COT_residuals)) +
  geom_tile() +
  # geom_hline(yintercept = seq(min(grid_flap$avg_thermUplift), max(grid_flap$avg_thermUplift), 
  #                             length.out = 20), 
  #            color = "white", linewidth = 0.1) +
  # geom_vline(xintercept = seq(min(grid_soar$avg_windSupp_1000m), max(grid_soar$avg_windSupp_1000m), 
  #                             length.out = 20), 
  #            color = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = my_gradient,
    # force to the same range for both plots
    values = scales::rescale(c(-3, 0, 3)),
    limits = c(-3, 3),
    #values = scales::rescale(c(min(allSegmDfs_mode$COT_residuals), 0, max(allSegmDfs_mode$COT_residuals))),
    name = "COT residuals"
  ) +
  # facet_wrap(~ species_phy, ncol = 4) +  # Adjust ncol for layout
  labs(
    x = "Average Wind Support (at 1000 m)",
    y = "Average Thermal Uplift",
    title = "Flappers"
  ) +
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

p2 <- ggplot(grid_soar, aes(x = avg_windSupp_1000m, y = avg_thermUplift, fill = COT_residuals)) +
  geom_tile() +
  # geom_hline(yintercept = seq(min(grid_soar$avg_thermUplift), max(grid_soar$avg_thermUplift), 
  #                             length.out = 20), 
  #            color = "white", linewidth = 0.1) +
  # geom_vline(xintercept = seq(min(grid_soar$avg_windSupp_1000m), max(grid_soar$avg_windSupp_1000m), 
  #                             length.out = 20), 
  #            color = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = my_gradient,
    # force to the same range for both plots
    values = scales::rescale(c(-3, 0, 3)),
    limits = c(-3, 3),
    #values = scales::rescale(c(min(allSegmDfs_mode$COT_residuals), 0, max(allSegmDfs_mode$COT_residuals))),
    name = "COT residuals"
  ) +
  # facet_wrap(~ species_phy, ncol = 4) +  # Adjust ncol for layout
  labs(
    x = "Average Wind Support (at 1000 m)",
    y = "Average Thermal Uplift",
    title = "Soarers"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#pdf("Plots/finalPlots/soarersFlappers_residualsVSenergyLandscape_Feb2025_greenMag.pdf", width=10, height=5)
png("Plots/finalPlots/soarersFlappers_residualsVSenergyLandscape_Feb2025_greenMag.png", width=10, height=5, units = "in", res=400)
grid.arrange(p1, p2, ncol = 2)
dev.off()












