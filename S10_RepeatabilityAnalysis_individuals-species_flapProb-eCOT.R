
#______________________________________________________________________________________
# Quantify repeatability of propFlap and eCOT within individual and within species ####
#______________________________________________________________________________________

library(rptR)
library(dplyr)

setwd("...")

# Load the dataset created in script 8
# This includes both the flapping probability and eCOT calculations averaged per commuting segment
allSegmDfs_withOutliers <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025.rds")
# same but excludes the outlier individuals (the dataset that was used in all analyses after script 8)
allSegmDfs_noOutliers <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")

#_______________________
# Proportion of flapping
#_________________________________________________
# Segment-level repeatability WITHIN individuals 
# Does an individual show consistent flapping proportion across its segments?
# Response: log-transformed flapping probability (as in script 8)

summary(allSegmDfs_noOutliers$avg_probFlap)
summary(log(allSegmDfs_noOutliers$avg_probFlap))

allSegmDfs_noOutliers$log_pFlap <- log(allSegmDfs_noOutliers$avg_probFlap)
rpt_seg <- rpt(
  log_pFlap ~ 1 + (1 | individualID), # random effect = individual
  grname= 'individualID',
  data= allSegmDfs_noOutliers,
  datatype = 'Gaussian',
  nboot= 1000,# bootstrapped CIs
  npermut= 0# permutation test optional
)
summary(rpt_seg)
# R (Intraclass Correlation Coefficient ICC = 0.392), CI (2.5% 0.366 - 97.5% 0.418), p-value < 2.2 × 10-16 (smallest number representable)

# for comparison, when using the dataset with outlier individuals
allSegmDfs_withOutliers$log_pFlap <- log(allSegmDfs_withOutliers$avg_probFlap)
rpt_seg2 <- rpt(
  log_pFlap ~ 1 + (1 | individualID), # random effect = individual
  grname= 'individualID',
  data= allSegmDfs_withOutliers,
  datatype = 'Gaussian',
  nboot= 1000,# bootstrapped CIs
  npermut= 0# permutation test optional
)
summary(rpt_seg2)
#ICC = 0.398, CI 2.5% 0.372 - 97.5% 0.424, p-value < 2.2 × 10-16 (smallest number representable)


#______________________________________________
# Individual-level repeatability WITHIN species
# Does knowing the species predict an individual's mean flapping proportion?
# First compute individual means (log scale)
perInd <- allSegmDfs_noOutliers %>%
  group_by(individualID, species_phy) %>%
  summarise(log_pFlap_mean = mean(log(avg_probFlap), na.rm = TRUE),
            .groups = 'drop')
rpt_ind <- rpt(
  log_pFlap_mean ~ 1 + (1 | species_phy),   # random effect = species
  grname= 'species_phy',
  data= perInd,
  datatype = 'Gaussian',
  nboot= 1000,
  npermut= 0
)
summary(rpt_ind)
# ICC = 0.895, CI 2.5% 0.839 - 97.5% 0.927, p-value < 2.2 × 10-16

# therefore regarding proportion of flapping, individuals are moderately consistent across segments (justifying individual means), and species are highly consistent across individuals (justifying species means as the unit of the analysis in script 8 and 9)
icc_table <- data.frame(
  Level= c('Segments within individuals', 'Individuals within species'),
  Response= 'log(proportion flapping)',
  ICC= c(rpt_seg$R$individualID, rpt_ind$R$species_phy),
  CI_low= c(rpt_seg$CI_emp$individualID[1], rpt_ind$CI_emp$species_phy[1]),
  CI_high= c(rpt_seg$CI_emp$individualID[2], rpt_ind$CI_emp$species_phy[2]),
  p_value= c(rpt_seg$P$individualID[1],
             rpt_ind$P$species_phy[1]))
print(round(icc_table[, c('Level','ICC','CI_low','CI_high','p_value')], 3))

# for comparison, when using the dataset with outlier individuals
perInd2 <- allSegmDfs_withOutliers %>%
  group_by(individualID, species_phy) %>%
  summarise(log_pFlap_mean = mean(log(avg_probFlap), na.rm = TRUE),
            .groups = 'drop')
rpt_ind2 <- rpt(
  log_pFlap_mean ~ 1 + (1 | species_phy),   # random effect = species
  grname= 'species_phy',
  data= perInd2,
  datatype = 'Gaussian',
  nboot= 1000,
  npermut= 0
)
summary(rpt_ind2)
# ICC = 0.753, CI 2.5% 0.634 - 97.5% 0.82, p-value 1.3e-210

#_______
# eCOT
#__________________________________________________________
# Segment-level repeatability of eCOT WITHIN individuals

summary(allSegmDfs_noOutliers$obs_tot_COT_J_KgM)
summary(log(allSegmDfs_noOutliers$obs_tot_COT_J_KgM))

allSegmDfs_noOutliers$log_eCOT <- log(allSegmDfs_noOutliers$obs_tot_COT_J_KgM)
rpt_eCOT_seg <- rpt(
  log_eCOT ~ 1 + (1 | individualID),
  grname= 'individualID',
  data= allSegmDfs_noOutliers,
  datatype = 'Gaussian',
  nboot= 1000,
  npermut= 0
)
summary(rpt_eCOT_seg)
# ICC 0.807, CI 0.79-0.821, p < 2.2 × 10-16

#______________________________________________________
# Individual-level repeatability of eCOT WITHIN species

perInd_eCOT <- allSegmDfs_noOutliers %>%
  group_by(individualID, species_phy) %>%
  summarise(log_eCOT_mean = mean(log_eCOT, na.rm = TRUE),
            .groups = 'drop')
rpt_eCOT_ind <- rpt(
  log_eCOT_mean ~ 1 + (1 | species_phy),
  grname= 'species_phy',
  data= perInd_eCOT,
  datatype = 'Gaussian',
  nboot= 1000,
  npermut= 0
)
summary(rpt_eCOT_ind)
# ICC 0.933, CI 0.89-0.955, p < 2.2 × 10-16
