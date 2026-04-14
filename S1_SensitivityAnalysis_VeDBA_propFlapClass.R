
#_________________________________________________________________
# Robustness of the propFlap classification based on VeDBA ####
#_________________________________________________________________

library(mixR)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("...")
dir.create("Revision/NewSupplFigures", recursive = T)

# Load dataset saved in step 3A
load("DataFinalSummary/allStudies_allTags_allFlightSegments_binded_birdsBats_thresholdClass_transfGs_March2024_noDupl.RData") #object allStudies_noDups

# Load the mixfit model produced in step 3B, line 48
mod <- readRDS("DataFinalSummary/mixR_model_VeDBALognormBimodal_final.rds")

# Load the dataset associated with flapProp saved in step 3B, line 143
finalDf <- readRDS("DataFinalSummary/FinalDf_perPoint_VedbaGs_flappingProbs.rds")

# Extract means of the two distributions and antimode
mu1 <- mod$mu[1]
mu2 <- mod$mu[2]
equalProb <- which.min(abs(mod$comp.prob[,1] - mod$comp.prob[,2])) #point of almost equal probability of the two distributions
antimode <- mod$data[equalProb]

#____________________________________
## Explore classification uncertainty

# Maximum posterior probability for each observation
# (i.e. how confidently is each point assigned to its most likely component)
str(mod)
head(mod$comp.prob)
max_prob <- apply(mod$comp.prob, 1, max)
quantile(max_prob)

# Proportion of observations that are 'ambiguous'
# (i.e. where neither component reaches the confidence threshold)
thr <- 0.8
# adjust as desired; 0.90 is more conservative
prop_ambiguous <- mean(max_prob < thr)
cat(sprintf('Proportion of ambiguous observations (max prob < %.2f): %.1f%%\n',
            thr, prop_ambiguous * 100))

# Proportion of observation crossing the antimode boundary
# (i.e. difference between observation and antimode value being within 0.05 G)
near_antimode <- abs(finalDf$meanVedba_Gs - antimode) < 0.05
cat(sprintf('Observations within 0.05 g of antimode (%.3f g): %d (%.1f%%)\n',
            antimode, sum(near_antimode), mean(near_antimode)*100))

# Plot histogram of max posterior probability for Supplementary material
df_prob <- data.frame(max_prob = max_prob)
p_uncertainty <- ggplot(df_prob, aes(x = max_prob)) +
  geom_histogram(binwidth = 0.01, fill = 'darkgrey', colour = 'darkgrey', alpha = 0.85) +
  geom_vline(xintercept = thr, linetype = 'dashed',
             colour = '#C00000', linewidth = 0.6) +
  annotate('text', x = thr - 0.01, y = Inf,
           label = sprintf('%% Ambiguous \nObservations: %.1f%%', prop_ambiguous*100),
           hjust = 1, vjust = 1.3, colour = '#C00000', size = 3.5) +
  labs(#title = 'Classification confidence per observation',
    x = 'Maximum posterior probability from mixture model',
    y = 'Number of observations') +
  theme_classic(base_size = 11)
print(p_uncertainty)
# ggsave("Revision/NewSupplFigures/Fig_Sensitivity_propFlap_posteriorProb.pdf",
#        p_uncertainty, width = 5, height = 4)

#______________________________________
## Test sensitivity to threshold choice

print(antimode)
names(finalDf)

# try different thresholds
thresholds <- c(0.35, 0.50, 0.65)
thresh_labels <- paste0('Threshold ', thresholds, ' g')

# Binary flapping classification under each hard threshold
for (thr in thresholds) {
  col_name <- paste0('flap_thr_', gsub('\\.', '_', thr))
  finalDf[[col_name]] <- as.integer(finalDf$meanVedba_Gs >= thr)
}

# species mean flapping proportion for each threshold
speciesSummary <- finalDf %>%
  group_by(species) %>%
  summarise(
    prop_flap_model = mean(flapping_prob),
    prop_flap_thr035 = mean(flap_thr_0_35),
    prop_flap_thr050 = mean(flap_thr_0_5),
    prop_flap_thr065 = mean(flap_thr_0_65),
    .groups = 'drop'
  )

# Test correlation
cor_table <- speciesSummary %>%
  summarise(
    r_thr035 = cor(prop_flap_model, prop_flap_thr035, use = 'complete'),
    r_thr050 = cor(prop_flap_model, prop_flap_thr050, use = 'complete'),
    r_thr065 = cor(prop_flap_model, prop_flap_thr065, use = 'complete')
  )
print(cor_table)

# Scatter plots of species' mean flapProb of model vs each threshold
sp_long <- speciesSummary %>%
  pivot_longer(cols = c(prop_flap_thr035, prop_flap_thr050, prop_flap_thr065),
               names_to = 'method', values_to = 'prop_threshold') %>%
  mutate(method = recode(method,
                         prop_flap_thr035 = 'Threshold 0.35 g',
                         prop_flap_thr050 = 'Threshold 0.50 g',
                         prop_flap_thr065 = 'Threshold 0.65 g'))

p_sensitivity <- ggplot(sp_long,
                        aes(x = prop_flap_model, y = prop_threshold)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed',
              colour = 'grey50', linewidth = 0.6) +
  geom_point(alpha = 0.7, colour = '#4472C4', size = 2) +
  facet_wrap(~ method, nrow = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = 'Species mean flapping probability - mixture model)',
       y = 'Species mean flapping probability - pre-defined thresholds') +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'))
print(p_sensitivity)
# ggsave('Revision/NewSupplFigures/Fig_Sensitivity_propFlap_differentThresholds.pdf',
#        p_sensitivity, width = 9, height = 3.5)


#__________________________________________________
# Combine both figures together for the supplementary material

library(patchwork)

p_combined <- p_uncertainty + p_sensitivity +
  plot_layout(ncol = 1, heights = c(1.5, 2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 13))

ggsave('Revision/NewSupplFigures/Fig_Sensitivity_propFlap_VeDBAclassification.pdf',
       p_combined, width = 7, height = 10)
ggsave('Revision/NewSupplFigures/Fig_Sensitivity_propFlap_VeDBAclassification.png',
       p_combined, width = 7, height = 10, dpi = 300, units = "in")

