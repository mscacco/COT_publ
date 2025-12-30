
#______________________________________________________________________________________________________
## Segmentation of vedba values in Gs, using MULTIMODE, all tags together
#______________________________________________________________________________________________________

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
#setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")

load("DataFinalSummary/allStudies_allTags_allFlightSegments_binded_birdsBats_thresholdClass_transfGs_March2024_noDupl.RData") #object allStudies_noDups

#_________________
## Some filtering

# Check that there is no missing values in Vedba (none)
table(complete.cases(allStudies_noDups$meanVedba_Gs))

# Check values of Vedba
summary(allStudies_noDups$meanVedba_Gs)
toCheck <- unique(allStudies_noDups$study.name[allStudies_noDups$meanVedba_Gs>6])
summary(allStudies_noDups$meanVedba_Gs[allStudies_noDups$study.name%in%toCheck])
# quantile to remove extreme values, there is a small peak at 2.5 probably due to the max ceiling of the sensor in some tags
hist(allStudies_noDups$meanVedba_Gs, breaks="FD", xlab="Mean Vedba (Gs)")
q <- quantile(allStudies_noDups$meanVedba_Gs, probs=seq(0.9,1,0.0001))
plot(q)
tail(q) #cutting point at around 1.9
q["99.99%"]

# Remove tail, VeDBA > 2 g
allStudies_noDups <- allStudies_noDups[which(allStudies_noDups$meanVedba_Gs <= 2),]
hist(allStudies_noDups$meanVedba_Gs, breaks="FD", xlab="Mean Vedba (Gs)")

# Model Vedba as a function of sampling frequency etc
# vedbaModLog <- lm(meanVedba_Gs ~ acc_sampl_freq_per_axis + acc_burst_duration_s + deviceType, 
#  data=allStudies_noDups)
# summary(vedbaModLog)
# hist(residuals(vedbaModLog))

#____________________
## Filtering species

library(dplyr)

speciesSummary <- allStudies_noDups %>%
  group_by(individual.taxon.canonical.name) %>%
  summarise(meanVedba=mean(meanVedba_Gs),
            Nlocs=n())

speciesSummary %>% arrange(Nlocs) %>% print(n=Inf)

# Explore species with few observations and vedba values of tucans/parrots/hornbills
speciesSummary[speciesSummary$Nlocs < 20,]
speciesSummary[grep(c("Aceros nipalensis|Anthracoceros|Cacatua|Ara|leari"), speciesSummary$individual.taxon.canonical.name),]
table(grep(c("Aceros nipalensis|Anthracoceros|Cacatua|Ara|leari"), allStudies_noDups$individual.taxon.canonical.name, value=T))

# Remove the species identified above from the main dataset
speciesFewObs <- speciesSummary$individual.taxon.canonical.name[speciesSummary$Nlocs < 20]
Tucans <- grep("Aceros nipalensis|Anthracoceros|Cacatua|Ara|leari", speciesSummary$individual.taxon.canonical.name, value=T)
MurreCorm <- grep("Leucocarbo|Uria", speciesSummary$individual.taxon.canonical.name, value=T)
OilBarn <- grep("Steatornis|Tyto", speciesSummary$individual.taxon.canonical.name, value=T)
allStudies_noDups_sub <- allStudies_noDups[!allStudies_noDups$individual.taxon.canonical.name %in% c(speciesFewObs,Tucans,MurreCorm,OilBarn),]

length(unique(allStudies_noDups_sub$individual.taxon.canonical.name))
length(unique(allStudies_noDups$individual.taxon.canonical.name))

# Save the dataset to use in the following steps (after excluding the above species and the VeDBA tail)
saveRDS(allStudies_noDups_sub, file="DataFinalSummary/allStudies_allTags_VedbaGs_March2024_filteredSpecies.rds")

#_____________________
## Vedba segmentation

library(mixR)
library(ggplot2)

finalDf <- readRDS("DataFinalSummary/allStudies_allTags_VedbaGs_March2024_filteredSpecies.rds")

# Find initial values of means from histogram
hist(finalDf$meanVedba_Gs, breaks="FD")
(h1 <- hist(finalDf$meanVedba_Gs[finalDf$meanVedba_Gs < 0.5], breaks=100))
(h2 <- hist(finalDf$meanVedba_Gs[finalDf$meanVedba_Gs > 0.5], breaks=100))

proxyMu1 <- h1$mids[which.max(h1$counts)]
proxyMu2 <- h2$mids[which.max(h2$counts)]

plot(h1); abline(v=proxyMu1, lwd=3)
plot(h2); abline(v=proxyMu2, lwd=3)

# Fit two lognorm binomial distributions, using the points above as first mean approximations
mod <- mixfit(finalDf$meanVedba_Gs, family="lnorm", 
               ncomp=2, mu=c(proxyMu1, proxyMu2))
str(mod)
saveRDS(mod, "DataFinalSummary/mixR_model_VeDBALognormBimodal_final.rds")

plot(mod)
ggsave("Plots/ACCsegmentation/VeDBAsegmentation_allTags_mixR_rawOutput_final.pdf")

# The gamma model fitted better the peak in the distribution of passive flight on the left but worse the distribution of active flight
# The lognormal seems the best fit
# mod2 <- mixfit(finalDf$meanVedba_Gs, family="gamma", 
#               ncomp=2, mu=c(proxyMu1, proxyMu2))
# mod3 <- mixfit(finalDf$meanVedba_Gs, family="weibull", 
#                ncomp=2, mu=c(proxyMu1, proxyMu2))
# mod4 <- mixfit(finalDf$meanVedba_Gs, family="norm", 
#               ncomp=2, mu=c(proxyMu1, proxyMu2), ev=F)
#plot(mod2); plot(mod3); plot(mod4)

mod$mu
mod$sd
str(mod)
# the slot data contains the values of x
table(mod$data == finalDf$meanVedba_Gs)
# the component probability slot contains the probability of each value in x to belong to either distribution
summary(mod$comp.prob[,1])
summary(mod$comp.prob[,2])
length(mod$data) == nrow(mod$comp.prob)

# Extract means from model
mu1 <- mod$mu[1]
mu2 <- mod$mu[2]
# Extract antimode from model
equalProb <- which.min(abs(mod$comp.prob[,1] - mod$comp.prob[,2])) #point of almost equal probability of the two distributions
antimode <- mod$data[equalProb]
c(mu1, antimode, mu2)
# Extract densities
d <- density(mod, at=sort(mod$data)) #sort to allow plotting as line later

#_______________________________
# Add probabilities to the dataset, based on the model

identical(mod$comp.prob[,1], 1 - mod$comp.prob[,2])
summary(round(mod$comp.prob[,1] - (1 - mod$comp.prob[,2]), 3)) # the two components contain 1 - the same information, so we only include one

finalDf$flapping_prob <- mod$comp.prob[,2]

#_______________________________
# Add binary classification to the dataset, using the antimode as threshold

finalDf$flapping_bin <- 0
finalDf$flapping_bin[finalDf$meanVedba_Gs > antimode] <- 1


#_____________
## Some checks per species

library(dplyr)

speciesSummary <- finalDf %>%
  group_by(individual.taxon.canonical.name) %>%
  summarise(meanFlapp_prob=mean(flapping_prob),
            meanFlapp_bin=mean(flapping_bin),
            meanVedba=mean(meanVedba_Gs),
            Nlocs=n())
# Good, flapping prob and bin contain almost the same information
speciesSummary %>% print(n=Inf)

# Save summary per species
write.csv(as.data.frame(speciesSummary), file = "DataFinalSummary/probFlapping_perSpecies_mixR_March2024.csv", row.names = F)

#_______________________________
## Tidy up the final dataframe

colsToKeep <- c("timestamp","location.lat","location.long","height.above.msl","height_gener","height.above.ellipsoid",
                "event.id", "study.id","study.name","individual.id","individual.local.identifier","tag.local.identifier","individual.study.id",
                "individual.tag.id","individual.taxon.canonical.name","BodyMass_value",
                "deviceType","altitudeDiff","timeLag_min","vertSpeed_ms","stepLength_m","groundSpeed_ms","segmentDir","turnAngle",
                "diff_acc_time_s","acc_closest_timestamp","acc_burst_duration_s","acc_sampl_freq_per_axis","n_samples_per_axis",
                "acc_axes","acc_Naxes","flightClass","track_flight_id","meanVedba_Gs","cumVedba_Gs",
                "flapping_bin","flapping_prob","dataSource")
table(colsToKeep %in% names(finalDf))
# check height column:
# height_gener is a generic height column that includes both height.above.msl and height.above.ellipsoid
# missing values in this column correspond in values that are missing in both the other columns
# to reduce data loss we will use this generic column for movebank annotation
# height above ellipsoid and above sea level are not incredibly different so let's pull them together just fore the annotation
table(finalDf$height.above.msl == finalDf$height_gener)
table(finalDf$height.above.ellipsoid == finalDf$height_gener)
table(is.na(finalDf$height_gener))
table(is.na(finalDf$height.above.msl) & is.na(finalDf$height.above.ellipsoid))

# Add empty column if the column name is missing for binding all datasets from all studies
finalDf <- finalDf[,colsToKeep]

# Body mass in kg
names(finalDf)[grep("individual.taxon.canonical.name",names(finalDf))] <- "species"
anyNA(finalDf$BodyMass_value)
finalDf$Body_mass_kg <- finalDf[,"BodyMass_value"]/1000

# Create a column for later merging ENV annotated data with original dataset
finalDf$mergingCol <- 1:nrow(finalDf)

#_______________________________
# Save the final dataset with the new flapping/passive classification

saveRDS(finalDf, file="DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs.rds")


#____________________________________
## Fancy plotting of VeDBA classification

#pdf(file="Plots/ACCsegmentation/VeDBAsegmentation_allTags_mixR_greyHist.pdf")
tiff(file="Plots/ACCsegmentation/VeDBAsegmentation_allTags_mixR_greyHist_final.tiff", res=300, width = 7, height = 7, units = "in")

par(mar = c(5, 4, 4, 4) + 0.3) # Additional space for second y-axis
hist(mod$data, breaks=100, prob=T, xlab="Mean VeDBA (g)", # actual data
     main="Classification of mean VeDBA for all tags (g)", col="grey")
# lines(d$x, d$comp[,1], col="blue", lwd=2) #add densities from model
# lines(d$x, d$comp[,2], col="red", lwd=2)
abline(v=mu1, lty=1, lwd=2, col="blue") #add means of the distributions calculated by the model
abline(v=mu2, lty=1, lwd=2, col="red")
abline(v=antimode, lty=2, lwd=2, col="grey")
par(new = TRUE)
plot(sort(mod$data), mod$comp.prob[order(mod$data),1], col="blue", lwd=2, axes=F, xlab="",ylab="", type="l") # add probabilities calculated by the model
lines(sort(mod$data), mod$comp.prob[order(mod$data),2], col="red", lwd=2)
axis(side = 4, at = pretty(c(0,1)))      # Add second axis
mtext("Probability", side = 4, line = 3)   # Add second axis label

dev.off()

# Same but with colors in the histogram

tiff(file="Plots/ACCsegmentation/VeDBAsegmentation_allTags_mixR_colorHist_final.tiff", res=300, width = 7, height = 7, units = "in")
#pdf(file="Plots/finalPlots/VeDBAsegmentation_allTags_mixR_colorHist_final.pdf", width = 5, height = 5)

par(mar = c(5, 4, 4, 4) + 0.3) # Additional space for second y-axis
h <- hist(finalDf$meanVedba_Gs, breaks="FD", prob=F, main="VeDBA classification for all devices",
          col="white", border="white", xlab="Mean Vedba (g)")
hist(finalDf$meanVedba_Gs[finalDf$flapping_bin==0], breaks=h$breaks, prob=F, add=T, 
     col=alpha("blue", 0.1), border=alpha("blue", 0.15))
hist(finalDf$meanVedba_Gs[finalDf$flapping_bin==1], breaks=h$breaks, prob=F, add=T, 
     col=alpha("red", 0.1), border=alpha("red", 0.15))
legend(x=1.1, y=22000, c("Mean of Active Distribution", "Mean of Passive Distribution","Crossover point"), lty=1, lwd=3, cex=0.8, seg.len = 1,
       col=c("red","blue","darkgrey"), bty="n")

par(new = TRUE) # Add new plot
plot(sort(finalDf$meanVedba_Gs), (1 - sort(finalDf$flapping_prob)), col="blue", cex=0.3, 
     xlab="", ylab="", ylim=c(0,1), axes=F, type="l", lwd=2)
lines(sort(finalDf$meanVedba_Gs), sort(finalDf$flapping_prob), col="red", cex=0.3, lwd=2)
abline(v=mu1, lty=1, lwd=3, col="white")
abline(v=mu1, lty=2, lwd=3, col="blue")
abline(v=mu2, lty=1, lwd=3, col="white")
abline(v=mu2, lty=2, lwd=3, col="red")
abline(v=antimode, lty=1, lwd=3, col="white")
abline(v=antimode, lty=2, lwd=3, col="darkgrey")
axis(side = 4, at = pretty(c(0,1)))      # Add second axis
mtext("Probability", side = 4, line = 3)   # Add second axis label
# legend(x=1.144, y=0.79, legend=c("Passive flight","Active flight"), 
#        col=c(alpha("blue", 0.6), alpha("red", 0.6)), pch=19, bty="n", cex=0.8)

dev.off()


# Plot one classified histogram per study

finalDf <- readRDS("DataFinalSummary/allStudies_allTags_VedbaGs_March2024_filteredSpecies_flappingProbs.rds")

group <- paste(finalDf$study.name, finalDf$individual.taxon.canonical.name, sep="_")
studiesLS <- split(finalDf, group)

lapply(studiesLS, function(df){
  
  studyId <- ifelse(unique(df$dataSource)=="Movebank", unique(df$study.id), unique(df$dataSource))
  fileName <- paste0(studyId,"_",unique(df$individual.taxon.canonical.name),"_",unique(df$deviceType))
  print(fileName)
  
  tiff(file=paste0("Plots/ACCsegmentation/finalSpecies_VeDBAhistograms/study_",fileName,"_VeDBAsegmentation_MixR.tiff"),
       7.5,6, units="in", res=400)
  par(mar = c(5, 4, 4, 4) + 0.3) # Additional space for second y-axis
  h <- hist(df$meanVedba_Gs, breaks="FD", prob=F, main=fileName,
            col="white", border="white", xlab="mean Vedba (g)")
  hist(df$meanVedba_Gs[df$flapping_bin==1], breaks=h$breaks, prob=F, add=T, 
       col=alpha("firebrick", 0.1), border=alpha("firebrick", 0.75))
  hist(df$meanVedba_Gs[df$flapping_bin==0], breaks=h$breaks, prob=F, add=T, 
       col=alpha("blue", 0.1), border=alpha("blue", 0.75))

  par(new = TRUE) # Add new plot
  plot(sort(finalDf$meanVedba_Gs), (1 - sort(finalDf$flapping_prob)), col="firebrick", cex=0.3, 
       xlab="", ylab="", ylim=c(0,1), axes=F, type="l", lwd=2, xlim=range(df$meanVedba_Gs)) #important to plot it in the range of the species' dataset
  lines(sort(finalDf$meanVedba_Gs), sort(finalDf$flapping_prob), col="blue", cex=0.3, lwd=2)
  abline(v=mu1, lty=1, lwd=3, col="white")
  abline(v=mu1, lty=2, lwd=3, col="blue")
  abline(v=mu2, lty=1, lwd=3, col="white")
  abline(v=mu2, lty=2, lwd=3, col="firebrick")
  abline(v=antimode, lty=1, lwd=3, col="white")
  abline(v=antimode, lty=2, lwd=3, col="darkgrey")
  axis(side = 4, at = pretty(c(0,1)))      # Add second axis
  mtext("Probability", side = 4, line = 3)   # Add second axis label
  legend("topright", legend=c("Passive flight","Active flight"),
         col=c(alpha("blue", 0.6), alpha("firebrick", 0.6)), pch=19, cex=0.8)
  dev.off()
})

#________
## Small check with the golden eagles. 
## The histogram shows a firs peak at VeDBA around 0. But still the speeds corresponding to those segments are higher than the set threshold. So we keep them.
names(studiesLS)
df <- studiesLS[[36]]
unique(df$study.id)

summary(df$meanVedba_Gs)
hist(df$meanVedba_Gs[df$meanVedba_Gs < 0.1])
table(df$meanVedba_Gs < 0.04)
summary(df$groundSpeed_ms[df$meanVedba_Gs < 0.04])
summary(df$groundSpeed_ms[df$meanVedba_Gs < 0.01])
df %>% group_by(meanVedba_Gs < 0.01) %>% summarise(medSpeed=median(groundSpeed_ms))
