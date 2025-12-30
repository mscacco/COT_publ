
setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
#setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")

#___________________________________________
# EFFECT of BMR estimation on COT estimation
#___________________________________________

#______________________________________________________________
## 1. Plot different calculations of BMR relative to body mass

# Sequence of possible body masses from our dataset
mass_gr <- round(seq(12, 11000, length.out=100)) # 12 g is the minimum mass in our radar dataset; 11 kg the maximum mass in the our dataset

# Different BMR calculations, all in W (i.e. J/s)
BMR_Londono_tropical <- 0.044 * (mass_gr)^0.589
BMR_Londono_temperate <- 0.023 * (mass_gr)^0.729
BMR_Londono_pass <- 0.045 * (mass_gr)^0.627
BMR_Londono_nonPass <- 0.021 * (mass_gr)^0.724
BMR_McKechnie_birds <- 0.0346 * (mass_gr)^0.669 # same as: 10^(-1.461 + 0.669 * log10(mass_gr))
BMR_Speakman_bats <- exp(1.0895 + 0.744 * log(mass_gr)) * 0.005583

## Plot the different BMR lines against body mass values
## We add an inset to show how the position of the lines change depending on body masses
dev.off() # restore plotting device
pdf("Plots/finalPlots/BMR_lines.pdf", 7,7)
# Main plot
plot(mass_gr/1000, BMR_McKechnie_birds, type = "l", lwd=2, xlab="Body mass (kg)", ylab="BMR (W)" )
lines(mass_gr/1000, BMR_Londono_pass, col = "darkred", lty = 2, lwd=2)
lines(mass_gr/1000, BMR_Londono_nonPass, col = "forestgreen", lty = 2, lwd=2)
lines(mass_gr/1000, BMR_Londono_tropical, col = "darkblue", lty = 2, lwd=2)
lines(mass_gr/1000, BMR_Londono_temperate, col = "darkgrey", lty = 2, lwd=2)
lines(mass_gr/1000, BMR_Speakman_bats, col = "orange", lty = 3, lwd=3)
legend("bottomright",c("McKechnie 2004 (birds)","Londono 2015 (passerines)","Londono 2015 (non-passerines)",
                       "Londono 2015 (tropical)","Londono 2015 (temperate)","Speakman 2003 (bats)"),
       lty=c(1,2,2,2,2,3),col=c("black","darkred","forestgreen","darkblue","darkgrey","orange"), lwd=2,
       bty="n", cex=0.7)
# Inset plot, smaller mass range
# Zoomed-in inset plot (focusing between 12 grams and 3 kg of body mass)
# par(fig = c(0.05, 0.5, 0.5, 1), new = TRUE)  # Define inset plot area (left, right, bottom, top)
# plot(mass_gr/1000, BMR_McKechnie_birds, type = "l", lwd=2, 
#      xlab="", ylab="", xaxt='n', yaxt='n', bty='n', xlim=c(0,1), ylim=c(0,3.5))
# lines(mass_gr/1000, BMR_Londono_pass, col = "darkred", lty = 2, lwd=2)
# lines(mass_gr/1000, BMR_Londono_nonPass, col = "forestgreen", lty = 2, lwd=2)
# lines(mass_gr/1000, BMR_Londono_tropical, col = "darkblue", lty = 2, lwd=2)
# lines(mass_gr/1000, BMR_Londono_temperate, col = "darkgrey", lty = 2, lwd=2)
# lines(mass_gr/1000, BMR_Speakman_bats, col = "orange", lty = 3, lwd=2)
# axis(1, cex.axis=0.7, cex.lab=0.8)  # x-axis with smaller font
# axis(2, cex.axis=0.7, cex.lab=0.8)  # y-axis with smaller font
# box() 
dev.off()

#______________________________________________________________
## 2. Plot effect of BMR calculation on COT calculation
## for different body masses and proportions of flapping

dummyDf <- as.data.frame(cbind(mass_gr, BMR_Londono_nonPass, BMR_Londono_pass, BMR_Londono_temperate, BMR_Londono_tropical,
                               BMR_McKechnie_birds, BMR_Speakman_bats))

propFlap <- c(seq(0,1,length.out=5), 0.05) # add these additional two extreme soaring/flapping values

dummyDf_prop <- expand.grid(propFlap=propFlap, mass_gr=dummyDf$mass_gr)
dummyDf_prop <- merge(dummyDf_prop, dummyDf, by="mass_gr")
dummyDf_prop$mass_kg <- dummyDf_prop$mass_gr/1000
head(dummyDf_prop)

# calculate MR of flapping based on body mass
load("DataFinalSummary/flappingModel_KylesData.RData") #object modBirds, from script 0C
dummyDf_prop$MR_flap_guigueno2019 <- exp(predict(modBirds, newdata = data.frame(log_bodyMass_Kg=log(dummyDf_prop$mass_kg))))
dummyDf_prop$MR_flap_alexander2003 <- 3.6*(dummyDf_prop$mass_kg)^(-0.31)

# MR in Guigueno is in W, while in Alexander is in J/(kg m)
# To compare them we multiply the Guigueno values by the body mass and the median speed in our dataset (8 m/s)
dummyDf_prop$MR_flap_alexander2003 <- dummyDf_prop$MR_flap_alexander2003 * dummyDf_prop$mass_kg * 8

pdf("Plots/finalPlots/flapMR_lines.pdf", 7,7)
plot(dummyDf_prop$mass_kg, dummyDf_prop$MR_flap_guigueno2019, type = "l", lwd=2, 
     xlab="Body mass (kg)", ylab="Flapping MR (W)", ylim=c(1,380))
lines(dummyDf_prop$mass_kg, dummyDf_prop$MR_flap_alexander2003, col = "darkgrey", lty = 2, lwd=2)
legend("bottomright",c("Guigueno 2019","Alexander 2003"),
       lty=c(1,2),col=c("black","darkgrey"), lwd=2,
       bty="n", cex=0.7)
dev.off()


## Calculate passive costs based different BMRs
cost_pass_McKechnie_birds <- (2 * dummyDf_prop$BMR_McKechnie_birds) * (1 - dummyDf_prop$propFlap)
cost_pass_Londono_tropical <- (2 * dummyDf_prop$BMR_Londono_tropical) * (1 - dummyDf_prop$propFlap)
cost_pass_Londono_temperate <- (2 * dummyDf_prop$BMR_Londono_temperate) * (1 - dummyDf_prop$propFlap)
cost_pass_Londono_pass <- (2 * dummyDf_prop$BMR_Londono_pass) * (1 - dummyDf_prop$propFlap)
cost_pass_Londono_nonPass <- (2 * dummyDf_prop$BMR_Londono_nonPass) * (1 - dummyDf_prop$propFlap)
cost_pass_Speakman_bats <- (2 * dummyDf_prop$BMR_Speakman_bats) * (1 - dummyDf_prop$propFlap)

## Calculate active cost based on guigueno and alexander
cost_flap_guigueno <- dummyDf_prop$MR_flap_guigueno2019 * dummyDf_prop$propFlap
cost_flap_alexander <- dummyDf_prop$MR_flap_alexander2003 * dummyDf_prop$propFlap

# Sum all combinations of the two, and divide by body mass and speed to obtain COT in J/kg*m
# We consider a range of speeds between 3 and 15 m/s (median of the dataset of this study)
speeds <- c(seq(3, 15, length.out=5), 8) # add 8 m/s as it is the median speed in both datasets
dummyDf_prop <- do.call(rbind, lapply(speeds, function(s){
  dummyDf_prop$speed <- s
  return(dummyDf_prop)}))

# for Guigueno
dummyGuigueno <- dummyDf_prop
dummyGuigueno$MRsource <- "Guigueno 2019"
dummyGuigueno$overall_COT_McKechnie_birds_Jkgm <- (cost_flap_guigueno + cost_pass_McKechnie_birds) / dummyGuigueno$mass_kg / dummyGuigueno$speed
dummyGuigueno$overall_COT_Londono_tropical_Jkgm <- (cost_flap_guigueno + cost_pass_Londono_tropical) / dummyGuigueno$mass_kg / dummyGuigueno$speed
dummyGuigueno$overall_COT_Londono_temperate_Jkgm <- (cost_flap_guigueno + cost_pass_Londono_temperate) / dummyGuigueno$mass_kg / dummyGuigueno$speed
dummyGuigueno$overall_COT_Londono_pass_Jkgm <- (cost_flap_guigueno + cost_pass_Londono_pass) / dummyGuigueno$mass_kg / dummyGuigueno$speed
dummyGuigueno$overall_COT_Londono_nonPass_Jkgm <- (cost_flap_guigueno + cost_pass_Londono_nonPass) / dummyGuigueno$mass_kg / dummyGuigueno$speed
dummyGuigueno$overall_COT_Speakman_bats_Jkgm <- (cost_flap_guigueno + cost_pass_Speakman_bats) / dummyGuigueno$mass_kg / dummyGuigueno$speed

# for Alexander
dummyAlexander <- dummyDf_prop
dummyAlexander$MRsource <- "Alexander 2003"
dummyAlexander$overall_COT_McKechnie_birds_Jkgm <- (cost_flap_alexander + cost_pass_McKechnie_birds) / dummyAlexander$mass_kg / dummyAlexander$speed
dummyAlexander$overall_COT_Londono_tropical_Jkgm <- (cost_flap_alexander + cost_pass_Londono_tropical) / dummyAlexander$mass_kg / dummyAlexander$speed
dummyAlexander$overall_COT_Londono_temperate_Jkgm <- (cost_flap_alexander + cost_pass_Londono_temperate) / dummyAlexander$mass_kg / dummyAlexander$speed
dummyAlexander$overall_COT_Londono_pass_Jkgm <- (cost_flap_alexander + cost_pass_Londono_pass) / dummyAlexander$mass_kg / dummyAlexander$speed
dummyAlexander$overall_COT_Londono_nonPass_Jkgm <- (cost_flap_alexander + cost_pass_Londono_nonPass) / dummyAlexander$mass_kg / dummyAlexander$speed
dummyAlexander$overall_COT_Speakman_bats_Jkgm <- (cost_flap_alexander + cost_pass_Speakman_bats) / dummyGuigueno$mass_kg / dummyAlexander$speed

# Bind the two
dummyDf_prop_final <- rbind(dummyGuigueno, dummyAlexander)

# Simple plot with two extreme cases (flapping = 0.05 and 0.95) fixing the speed to its median (8 m/s)
sub <- dummyDf_prop_final[dummyDf_prop_final$propFlap %in% c(0.05, 0.75) & dummyDf_prop_final$speed == 8,]

library(scales)
pdf("Plots/finalPlots/EffectOnCOT_differentBMR_differentMR_speed8.pdf", 13, 7)
par(mfrow=c(1,2))
plot(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], type="l", 
     xlab="Log of body mass (kg)", ylab="Log of overall COT (J/(kg*m))", ylim=c(-1.5,3.5), col="black", lwd=3, main = "Flapping MR from Alexander 2003")
lines(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Alexander 2003",], col = alpha("black", 0.8), lty = 2, lwd=3)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("darkred", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Alexander 2003",], col = alpha("darkred", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("forestgreen", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Alexander 2003",], col = alpha("forestgreen", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("darkblue", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Alexander 2003",], col = alpha("darkblue", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("darkgrey", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Alexander 2003",], col = alpha("darkgrey", 0.8), lty = 2, lwd=2)
legend("bottomleft",c("McKechnie 2004 (birds)","Londono 2015 (passerines)","Londono 2015 (non-passerines)",
                       "Londono 2015 (tropical)","Londono 2015 (temperate)",
                       "soarers (pFlap = 0.05)","flappers (pFlap = 0.95)"),
       lty=c(1,1,1,1,1,1,2),col=c("black","darkred","forestgreen","darkblue","darkgrey","black","black"), lwd=2,
       bty="n", cex=0.7)
plot(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], type="l", 
     xlab="Log of body mass (kg)", ylab="Log of overall COT (J/(kg*m))", ylim=c(-1.5,3.5), col="black", lwd=3, main="Flapping MR from Guigueno 2019")
lines(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Guigueno 2019",], col = alpha("black", 0.8), lty = 2, lwd=3)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("darkred", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Guigueno 2019",], col = alpha("darkred", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("forestgreen", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Guigueno 2019",], col = alpha("forestgreen", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("darkblue", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Guigueno 2019",], col = alpha("darkblue", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("darkgrey", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.75 & sub$MRsource == "Guigueno 2019",], col = alpha("darkgrey", 0.8), lty = 2, lwd=2)
dev.off()

# Reshape df to fit ggplot (different COT estimations in rows instead of columns)
library(tidyr)
library(dplyr)
library(ggplot2)
df_long <- dummyDf_prop_final %>%
  pivot_longer(cols = starts_with("overall_COT"), 
               names_to = "BMR_source", 
               values_to = "overall_COT_Jkgm")

# summarise COT values for each group
df_summarized <- df_long %>%
  group_by(BMR_source, MRsource, log(mass_kg)) %>%
  summarize(
    mean_COT = mean(log(overall_COT_Jkgm)),   # Average COT per bodymass per group
    sd_COT = sd(log(overall_COT_Jkgm)),       # Standard deviation
    se_COT = sd_COT / sqrt(n()),              # Standard error
    ci_lower = mean_COT - 1.96 * se_COT,      # Lower 95% CI
    ci_upper = mean_COT + 1.96 * se_COT       # Upper 95% CI
  )
table(df_summarized$BMR_source)
table(df_summarized$MRsource)

df_summarized$BMR_source <- gsub("overall_COT_|_Jkgm","",df_summarized$BMR_source)

# Plot
ggplot(df_summarized, aes(x = `log(mass_kg)`, y = mean_COT, color = BMR_source, group = BMR_source)) +
  geom_line(size = 2) +                         # Line plot for each COT group
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = BMR_source), 
              alpha = 0.05, linetype = "blank") + # Confidence intervals
  labs(x = "log of Body Mass", y = "log of overall COT (J/(kg * m))", 
       #title = "COT vs Body Mass with Confidence Intervals",
       color = "COT Group", fill = "COT Group") +
  theme_minimal() +
  facet_wrap(~MRsource) +           # Create facets for each metabolic rate source
  theme(legend.position = "top")    # Optional: position legend at the top
ggsave("Plots/finalPlots/EffectOnCOT_differentBMR_differentMR_allPropsAllSpeeds_gg.pdf", width =13,height=7) 


#_______________________________________________
# EXPLORE EFFECT OF CALCULATION IN OUR DATASET!
#_______________________________________________

#__________________________________
## Summary of our species, both gps and radar
## in terms of body mass, proportion of flapping and speed

# Import our data, both gps and radar
gps <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Aug2024_max10hours.rds") #dataset with one row per segment
radar <- readRDS("DataFinalSummary/RADARdata_finalSummaryDataset_perEcho_COTvariables.rds")

# Import kile's MR model
load("DataFinalSummary/flappingModel_KylesData.RData") #object modBirds, from script 0C

# summarise speed prop of flapping and body mass
speciesSummary <- rbind(group_by(gps, species=species) %>% summarise(mass_kg=unique(Body_mass_kg), 
                                                                     avgSpeed=mean(avg_grSpeed_ms),
                                                                     avgPflap=mean(avg_probFlap),
                                                                     group="Non-passerines"),
                        group_by(radar, species=WFFmonth_species) %>% summarise(mass_kg=unique(meanMass_kg), 
                                                                                avgSpeed=mean(Speed.ms),
                                                                                avgPflap=mean(propPulse),
                                                                                group="Passerines"))
speciesSummary$logMass_kg <- log(speciesSummary$mass_kg)
summary(speciesSummary)
table(speciesSummary$group)

speciesSummary$BMR <- NA
speciesSummary_cot <- do.call(rbind, lapply(c("McKechnie","Londono_tropical","Londono_temperate","Londono_pass"), function(source){
  if(source=="McKechnie"){speciesSummary$BMR <- 0.0346 * (speciesSummary$mass_kg*1000)^0.669}
  if(source=="Londono_tropical"){speciesSummary$BMR <- 0.044 * (speciesSummary$mass_kg*1000)^0.589}
  if(source=="Londono_temperate"){speciesSummary$BMR <- 0.023 * (speciesSummary$mass_kg*1000)^0.729}
  if(source=="Londono_pass"){
    speciesSummary$BMR[speciesSummary$group=="Passerines"] <- 0.045 * (speciesSummary$mass_kg*1000)[speciesSummary$group=="Passerines"]^0.627
    speciesSummary$BMR[speciesSummary$group=="Non-passerines"] <- 0.045 * (speciesSummary$mass_kg*1000)[speciesSummary$group=="Non-passerines"]^0.627
  }
  speciesSummary$source <- source
  
  costPass <- 2 * speciesSummary$BMR * (1 - speciesSummary$avgPflap)
  costFlap_guig <- (exp(predict(modBirds, newdata = data.frame(log_bodyMass_Kg=speciesSummary$logMass_kg)))) * speciesSummary$avgPflap
  costFlap_alex <- (3.6*(speciesSummary$mass_kg)^(-0.31)) * speciesSummary$avgPflap
  
  speciesSummary$overallCOT_guig <- (costPass + costFlap_guig) / speciesSummary$mass_kg / speciesSummary$avgSpeed
  speciesSummary$overallCOT_alex <- (costPass + costFlap_alex) / speciesSummary$mass_kg / speciesSummary$avgSpeed
  return(speciesSummary)
}))
head(speciesSummary_cot)
speciesSummary_cot$source <- factor(speciesSummary_cot$source, levels = c("McKechnie","Londono_tropical","Londono_temperate","Londono_pass"))
saveRDS(speciesSummary_cot, file="DataFinalSummary/speciesSummary_averageBMR_overallCOT_different sources.rds")

mod_guig <- lm(log(overallCOT_guig) ~ logMass_kg * source, data = speciesSummary_cot)
summary(mod_guig)
mod_alex <- lm(log(overallCOT_alex) ~ logMass_kg * source, data = speciesSummary_cot)
summary(mod_alex)

coefGuig <- coefficients(mod_guig)
coefAlex <- coefficients(mod_alex)

COT_table <- data.frame(sourceMB=c("Guigueno","Alexander"),
           slope_McKechnie=c(coefGuig["logMass_kg"], coefAlex["logMass_kg"]),
           slope_LondonTropical=c(coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_tropical"],
                                  coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_tropical"]),
           slope_LondonTemperate=c(coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_temperate"],
                                  coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_temperate"]),
           slope_LondonPass=c(coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_pass"],
                                  coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_pass"]))
COT_table

# Export
write.csv(COT_table, "Tables/COT_BMR_MR_effects_suppl.csv", row.names = F)


ggplot(data=speciesSummary_cot, aes(x=logMass_kg, y=log(overallCOT_guig))) +
geom_abline(intercept=coefGuig["(Intercept)"], slope=coefGuig["logMass_kg"], col="black", size=1.2) +
geom_abline(intercept=(coefGuig["(Intercept)"]+coefGuig["sourceLondono_tropical"]), 
            slope=coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_tropical"], col="darkred", linetype="dashed", size=1.2) +
geom_abline(intercept=(coefGuig["(Intercept)"]+coefGuig["sourceLondono_temperate"]), 
            slope=coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_temperate"], col="forestgreen", linetype="dashed", size=1.2) +
geom_abline(intercept=(coefGuig["(Intercept)"]+coefGuig["sourceLondono_pass"]), 
            slope=coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_pass"], col="darkblue", linetype="dashed", size=1.2) +
xlab("log of Body mass (kg)") + ylab("Log of overall COT (J/kg-1 m-1)") +
theme_classic()







