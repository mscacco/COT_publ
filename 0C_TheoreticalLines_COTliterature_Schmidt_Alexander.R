
#_________________________
# THEORETICAL COT LINES
#_________________________
## Digitise theoretical COT lines based on literature

#_______________________________
# Based on SCHMIDT NIELSEN 1972

theoreticalCOTschmidt <- data.frame(species=c("hummingbird","budgerigar","gull",
                                              "Phyllostomus bat","pigeon"), 
                                    bodyMass_g=c(3,35,300,90,384),
                                    MetRateFlight_cal_gh=c(204,105,54,94,58),
                                    speed_kmh=c(50,35,36,17,57.9))
theoreticalCOTschmidt <- theoreticalCOTschmidt[order(theoreticalCOTschmidt$bodyMass_g),]
theoreticalCOTschmidt$bodyMass_Kg <- theoreticalCOTschmidt$bodyMass_g/1000

# Convert in Watt/Kg (J/s*Kg) to compare it to the empirical data from Kyle
# 1 calories / (g * hour) = 1.16222222 W / Kg
theoreticalCOTschmidt$MetRateFlight_W_kg <- 
  theoreticalCOTschmidt$MetRateFlight_cal_gh * 1.16222222
# Calculate metabolic cost for that species with that weight (only in W)
theoreticalCOTschmidt$MetRateFlight_W <- 
  theoreticalCOTschmidt$MetRateFlight_W_kg * theoreticalCOTschmidt$bodyMass_Kg
# Use speed to convert costs to Joule/(Kg*metre) for comparison with Alexander 2003.
theoreticalCOTschmidt$speed_ms <- theoreticalCOTschmidt$speed_kmh / 3.6
theoreticalCOTschmidt$COTflight_J_KgM <- 
  theoreticalCOTschmidt$MetRateFlight_W_kg/theoreticalCOTschmidt$speed_ms


#_________________________
# Based on ALEXANDER 2003

# According to this study, flight cost is calculated as follows and expressed in J kg-1 m-1.
# Plot the costs in Joule/(Kg*m) against the body mass in Kg
theoreticalCOTalex <- data.frame(bodyMass_Kg=seq(0.001,30, by=0.05))
theoreticalCOTalex$COTflight_J_KgM <- 3.6*(theoreticalCOTalex$bodyMass_Kg)^(-0.31)
theoreticalCOTalex$COTswim_J_KgM <- 1.1*(theoreticalCOTalex$bodyMass_Kg)^(-0.38)
theoreticalCOTalex$COTrun_J_KgM <- 10.7*(theoreticalCOTalex$bodyMass_Kg)^(-0.32)

#_____________________________________
## Save and plot theoretical df lines

setwd("...")
save(theoreticalCOTschmidt, theoreticalCOTalex, 
     file="DataFinalSummary/thereticalCOTs_Nielsen-Alexander.RData")

pdf("Plots/finalPlots/theoreticalLines.pdf")
plot(log(COTflight_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTschmidt,
     xlim=c(-6,2.3), ylim=c(0,5),
     xlab="Body mass (Kg)", ylab="Flight cost (J Kg-1 m-1)", 
     main="Cost of transport during flight", type="n")
abline(lm(log(COTflight_J_KgM)~log(bodyMass_Kg), theoreticalCOTschmidt), 
       lty=1, lwd=1.5, col="magenta3")
lines(log(COTflight_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex, 
      lty=2, lwd=1.5, col="forestgreen")
lines(log(COTswim_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex, 
      lty=2, lwd=1.5, col="dodgerblue3")
lines(log(COTrun_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex, 
      lty=2, lwd=1.5, col="gold2")
legend("topright", c("Schmidt-Nielsen 1972, flying","Alexander 2003, flying","Alexander 2003, swimming","Alexander 2003, running"),  
       col=c("magenta3","forestgreen","dodgerblue3","gold2"), bty="n", lty=c(1,2,2,2), lwd=1.5)
dev.off()

#____________________________
# EMPIRICAL FLAPPING METABOLIC RATE (from GUIGUENO et al. 2019)
#____________________________

powerDf_kyle <- read.csv("DataAvailable/EnergyCostEstimation/FlappingPowerValues_Kyle.csv", as.is=T)
powerDf_kyle$Body_mass_kg <- powerDf_kyle$Body_mass/1000
powerDf_kyle$log_bodyMass_Kg <- log(powerDf_kyle$Body_mass_kg)

# cost of gliding relative to BMR for seabirds
gliders <- powerDf_kyle[powerDf_kyle$Flight_mode == "Gliding",]
cbind(gliders$Species, gliders$MetaRate / gliders$BMR)
summary(gliders$MetaRate / gliders$BMR)

# species with sustained flapping
flappers <- powerDf_kyle[powerDf_kyle$Flight_mode == "Flapping",]
summary(flappers$Body_mass[flappers$Taxon == "Bat"])
summary(flappers$Body_mass[flappers$Taxon != "Bat"])
table(flappers$Taxon)
summary(flappers$MetaRate / flappers$BMR)

sub4Emily <- flappers[flappers$Taxon != "Bat",c("Species","Taxon","Flight_mode","Body_mass_kg","MetaRate")]
write.csv(sub4Emily, "Plots/finalPlots/newModelPlot_newRunSwim_phyloModels/toKamiEmily/KylesData_sub4model.csv", row.names = F)

# Model "flapping MR ~ body mass", separately for birds and bats
modBirds <- lm(log(MetaRate)~log_bodyMass_Kg, data=flappers[flappers$Taxon != "Bat",])
summary(modBirds)
modBats <- lm(log(MetaRate)~log_bodyMass_Kg, data=flappers[flappers$Taxon == "Bat",])
summary(modBats)

(betaMassBi <- as.numeric(coefficients(modBirds)["log_bodyMass_Kg"]))
(betaMassBa <- as.numeric(coefficients(modBats)["log_bodyMass_Kg"]))
# a 1% increase in body mass causes about 0.78% increase in metabolic rate in birds and a 0.79% in bats

save(modBirds,modBats,betaMassBi,betaMassBa, file = "DataFinalSummary/flappingModel_KylesData.RData")

