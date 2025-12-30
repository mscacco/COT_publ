
#setwd("/Users/jeroen/Archive/Research/Projects/Energetics/Birds_Cost of transport")
# setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
setwd("/home/martina/ownCloud/Martina/ProgettiVari/COT/")

#_____________________________________________________
# CLASSIFY SOARERS/FLAPPERS based on prop of flapping
#_____________________________________________________

library(phylolm)
library(geiger)
library(ape)
library(rr2)
library(matrixStats)
library(devtools)
#install_github("JeroenSmaers/evomap")
library(evomap)
#install_github("khabbazian/l1ou")
library(l1ou)
library(dplyr)
library(scales)

#__________________
# Define functions
pGLS.plotGrade<-function (Yvar, Xvar, data, tree, group, model,...) 
{
  dataTemp <- pruneSample(na.omit(data), tree, group)$data
  treeTemp <- pruneSample(dataTemp, tree, group)$tree
  Y <- dataTemp[, which(colnames(dataTemp) == paste(Yvar))]
  X <- dataTemp[, which(colnames(dataTemp) == paste(Xvar))]
  dataGLS <- as.data.frame(cbind(Y, X))
  rownames(dataGLS) <- rownames(dataTemp)
  pGLSTemp <- phylolm(Y ~ X, data=dataGLS, phy=treeTemp,model = model)
  a <- summary(pGLSTemp)$coefficients[1, 1]
  b <- summary(pGLSTemp)$coefficients[2, 1]
  lines(c(min(X), max(X)), c((a + b * min(X)), (a + b * max(X))), 
        ...)
  points(X, Y, ...)
}

# Function to identify outliers using the IQR method
is_outlier <- function(x, threshold = 1.5) {
  # Calculate Q1, Q3, and IQR
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  # Define lower and upper bounds
  lower_bound <- q1 - threshold * iqr
  upper_bound <- q3 + threshold * iqr
  # Identify outliers
  outlier_logic <- x < lower_bound | x > upper_bound
  return(outlier_logic)
}

#___________________________________________________
# Prune tree to species in the dataset (to do once)

# # Import model dataset with COT calculations
# allSegmDfs <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025.rds")
# 
# # Create list of bird species (exclude bats) to download the corresponding trees
# writeLines(unique(allSegmDfs$species_phy), "Phylogeny/Trees_species_2025/speciesList_forPhyloTree_2025.txt")
# 
# # We use this list to download a phylogeny subset of 1000 trees from birdtree.org
# # We then summarise the 1000 trees in one majority-rule consensus (MRC) tree using the phyton library DendroPy (external phyton code)
# 
# # Import the resulting MRC tree
# phyloSum_Er <- read.nexus("Phylogeny/Trees_species_2025/MRCtree_DendroPy_from1000_Ericson_2025.tre")
# 
# is.ultrametric(phyloSum_Er)  # Returns TRUE if the tree is ultrametric
# 
# data <- allSegmDfs
# data$species_phy <- sub(" ","_",allSegmDfs$species_phy)
# 
# # Prune tree and filter data to match species groups
# setdiff(data$species_phy, phyloSum_Er$tip.label)  # Species in data but not in tree
# setdiff(phyloSum_Er$tip.label, data$species_phy)  # Species in tree but not in data
# 
# data <- data[data$species_phy %in% phyloSum_Er$tip.label,] # filter out bats
# pruned_tree <- keep.tip(phyloSum_Er, unique(data$species_phy)) # prune tree to bird species in dataset
# 
# # Save pruned tree for next step
# write.nexus(pruned_tree, file = "Phylogeny/Trees_species_2025/MRCtree_DendroPy_from1000_Ericson_Feb2025_pruned.tre")


#_________________
# Define outliers

# Import model dataset with COT calculations
allSegmDfs <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025.rds")

# individuals' means
perInd <- allSegmDfs %>% group_by(individualID, species_phy) %>%
  summarise(pFlap_geom = exp(mean(log(avg_probFlap), na.rm=T)), # geometric mean
            pFlap_avg = mean(avg_probFlap, na.rm=T), # arithmetic mean
            X1 = unique(log_bodyMass_Kg), #body mass is already log
            n=sum(totNlocs), # number of observations
            duration=sum(segm_timeDuration_h)) # duration in hours

# check distribution of prop. flapping
hist(perInd$pFlap_avg)
quantile(perInd$pFlap_avg, seq(0,1,0.1))
quantile(perInd$pFlap_avg, seq(0,0.01,0.001))

# check number of observations per individual
hist(perInd$n)
quantile(perInd$n, seq(0,1,0.1))
nrow(perInd[perInd$n < 12,])/nrow(perInd) # 106 individuals (about 12%) have less than 12 observations
length(unique(perInd$species_phy)); length(unique(perInd$species_phy[perInd$n > 12])) # no species would get dropped

# Define outlier individuals based on their IQR of pFlap
perInd_ls <- split(perInd, perInd$species_phy)

outliersID <- unlist(sapply(perInd_ls, function(x){
  # which individuals are considered outliers in both average calculations
  out <- which(is_outlier(x$pFlap_avg) & is_outlier(x$pFlap_geom))
  return(x$individualID[out]) # return the individualID
}), recursive = T)

length(outliersID) # 31 individuals considered as outliers
table(perInd$individualID %in% outliersID)

# explore characteristics of outliers
allSegmDfs_out <- allSegmDfs[allSegmDfs$individualID %in% outliersID,]
table(allSegmDfs_out$individualID)
nrow(allSegmDfs_out) # 3993 observations belong to outlier individuals

# Remove outliers from the original dataset
allSegmDfs_noOutliers <- allSegmDfs[!allSegmDfs$individualID %in% outliersID,]

# save dataset without outliers
saveRDS(allSegmDfs_noOutliers, "DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")


#__________________________________________________
# Model mean flapping per species against body mass
# (After removing outliers, recalculate averages for the model)

# Import pruned tree
tree <- read.nexus("Phylogeny/Trees_species_2025/MRCtree_DendroPy_from1000_Ericson_Feb2025_pruned.tre")

# Import data without outlier individuals
allSegmDfs_noOutliers <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")
  
# Variables of interest
data <- allSegmDfs_noOutliers
Y <- "avg_probFlap"; X1 <- "log_bodyMass_Kg"; X2 <- "wingLoading_kgm2"; X3 <- "wingArea_ellipse_cm2"; 
Z <- "species_phy"; individualID <- "individualID"
data <- data[,c(Y, X1, X2, X3, Z, individualID)]
colnames(data) <- c("Y", "X1", "X2", "X3", "Z", "individualID")

# log transform variables (except body mass) before averaging
data$Y <- log(data$Y)
data$X2 <- log(data$X2)
data$X3 <- log(data$X3)

# individuals' means
data <- group_by(data, individualID) %>%
  summarise(Y = mean(Y, na.rm=T), 
            X1 = unique(X1), 
            X2 = mean(X2, na.rm=T),
            X3 = mean(X3, na.rm=T),
            Z = unique(Z))

# species' means
data <- group_by(data, Z) %>%
  summarise(Y = mean(Y),
            X1 = unique(X1),
            X2 = mean(X2),
            X3 = mean(X3)) %>%
  as.data.frame()

# Match tree to data and data to tree
attr(data,"na.action") <- NULL
rownames(data) <- gsub(" ", "_", data$Z) # replace space with _ for matching
setdiff(rownames(data), tree$tip.label)  # Species in data but not in tree
setdiff(tree$tip.label, rownames(data))  # Species in tree but not in data

data_bats <- data[!rownames(data) %in% tree$tip.label,] # extract bats
data <- data[rownames(data) %in% tree$tip.label,] # filter out bats
tree <- keep.tip(tree, rownames(data))
anyNA(data)

# !! Order species in data by tip labels in tree
data <- data[tree$tip.label, ]
identical(rownames(data), tree$tip.label)

# Check correlation of independent variables
cor(data$X1, data$X2) # body mass and wing loading (wing loading = body mass / wing area)
plot(data$X1, data$X2)
cor(data$X1, data$X3) # body mass and wing area (estimated from wing length, secondary length and wing span, see script 4A)
plot(data$X1, data$X3)

# body mass and wing loading have a correlation coef of 0.73
# body mass and wing area of 0.9
# We compare models with one and the other
pGLS_bm <- phylolm(Y ~ X1, data, phy=tree, model="lambda")
pGLS_wl <- phylolm(Y ~ X2, data, phy=tree, model="lambda")
pGLS_wa <- phylolm(Y ~ X3, data, phy=tree, model="lambda")
summary(pGLS_bm)
summary(pGLS_wl)
summary(pGLS_wa)
pGLS_bm$adj.r.squared; pGLS_wl$adj.r.squared; pGLS_wa$adj.r.squared

# Compare models
AIC(pGLS_bm); AIC(pGLS_wl); AIC(pGLS_wa)
# pGLS_BM has similar but the lowest AIC, so this provides the best fit. 
R2(mod = pGLS_bm, mod.r = pGLS_wa, phy = tree)
R2(mod = pGLS_bm, mod.r = pGLS_wl, phy = tree)
# R2_lik can be interpreted as the effect size difference, see https://ives.labs.wisc.edu/pdf/Ives_2019.pdf
# This effect size difference represents how much more variance in one model explains compared to the other, while accounting for phylogenetic relationships
# pGLS_BM explains 2.4% more variance than pGLS_wa, and 5.3% more variance than pGLS_wl 

########### COMBINED MODEL
## pGLS model combining the non-overlapping info between X1 and X2, which are 77% correlated
# residual_X2 <- lm(X2 ~ X1, data = data)$residuals  # Get unique part of X2
# pGLS <- phylolm(Y ~ X1 + residual_X2, data, phy=tree, model="lambda")
# summary(pGLS)
# res <- pGLS$residuals
# plot(Y ~ X1, data=data)
# abline(pGLS)
# lambda <- pGLS$optpar[1]
# Sigma <- vcv(rescale(tree,"lambda",lambda))    #'Sigma' is now our new vcv(tree)
## Excluded it because it does not perform better than using body mass or wing area alone
###########

# We continue using only body mass, as it explained ~=same variance of pFlap as wing area, and more than wing loading
# In addition, body mass information are much easier to find in the literature compared to wing morphology information
data$Independent <- data$X1 # renaming it makes it easier to change it and compare results
pGLS <- phylolm(Y ~ Independent, data, phy=tree, model="lambda")
summary(pGLS)
R2_lik(pGLS)
res <- pGLS$residuals
lambda <- pGLS$optpar[1]
Sigma <- vcv(phytools::rescale(tree,"lambda",lambda))    #'Sigma' is now our new vcv(tree)

# Simple plot
plot(Y ~ Independent, data, pch=19, xlim=c(-1.4,3), ylim=c(-22,1.5),
     xlab="log of Body Mass (kg)", ylab="Mean of log proportion of flapping")
abline(pGLS)
text(Y~Independent, data[data$Z=="Pandion haliaetus",], labels=data$Z[data$Z=="Pandion haliaetus"], pos=4, cex=0.8, offset=0.5, srt=30, font=3)

# Inference of possible grade shifts based on residuals
Model <- adjust_data(tree, res)
eModel <- estimate_shift_configuration(Model$tree, Model$Y, criterion="AICc")
plot(eModel, cex=0.7, label.offset=0.02,edge.ann.cex=0.5,edge.width=2)

# Export plots
pdf("Plots/finalPlots/simplePlot_pFlapBMmodel.pdf", width=5, height=4)
plot(Y ~ Independent, data, pch=19, xlim=c(-1.4,3), ylim=c(-22,1.5),
     xlab="log of Body Mass (kg)", ylab="Mean of log proportion of flapping")
abline(pGLS)
dev.off()
pdf("Plots/finalPlots/simplePlot_gradeshift.pdf", width=7, height=5)
plot(eModel, cex=0.7, label.offset=0.02,edge.ann.cex=0.5,edge.width=2)
dev.off()


#_______________________________________________
# Test hypothesis of species grouping based on Bayou Allometry results
All <- 1:length(tree$tip.label)
Osprey <- which(tree$tip.label=="Pandion_haliaetus")
Telluraves <- getTips(tree, findMRCA(tree, c("Vultur_gryphus","Falco_peregrinus")))
Storks <- getTips(tree, findMRCA(tree, c("Ciconia_ciconia","Mycteria_americana")))

data[c(Storks,setdiff(Telluraves,Osprey)),] #soarers
data[setdiff(All,c(setdiff(Telluraves,Osprey),Storks)),] # flappers


# pANCOVA Test of these possible grade shifts
par(mfrow=c(1,2))
#Intercept
grpI<-rep("XXX",length(rownames(data)))
grpI[All]<-"1"
grpI[c(setdiff(Telluraves,Osprey),Storks)]<-"2"
grpI<-as.factor(grpI)
names(grpI)<-rownames(data)
summary(grpI); length(grpI); length(names(summary(grpI)))
tree$tip.label[which(grpI=="XXX")]
tipCol<-rep("black",length(tree$tip.label))
tipCol[which(grpI=="1")]<-"gold2"
tipCol[which(grpI=="2")]<-"dodgerblue3"
plot(tree,tip.color = tipCol)

#Slope
grpS<-rep("XXX",length(rownames(data)))
grpS[All]<-"1"
grpS[c(setdiff(Telluraves,Osprey),Storks)]<-"2"
grpS<-as.factor(grpS); names(grpS)<-rownames(data)
summary(grpS); length(grpS); length(names(summary(grpS)))
tree$tip.label[which(grpS=="XXX")]
tipCol<-rep("black",length(tree$tip.label))
tipCol[which(grpS=="1")]<-"gold2"
tipCol[which(grpS=="2")]<-"dodgerblue3"
plot(tree,tip.color = tipCol)

Model <- model.matrix(as.formula(Y ~ Independent), data)
Model_SI <- model.matrix(as.formula(Y ~ grpI + grpS:Independent),data)
gls.ancova(Y ~ Independent, Sigma, Model, Model_SI)

Model <- phylolm(as.formula(Y ~ Independent), data, tree, model="lambda")
Model_ML <- phylolm(as.formula(Y ~ grpI + grpS:Independent), data, tree, model="lambda")

Model$aic
Model_ML$aic
Model$aic - Model_ML$aic
R2_lik(Model)
R2_lik(Model_ML)
R2(mod = Model_ML, mod.r = Model, phy = tree)
# The model including two different intercepts and slopes for the two groups
# has a much lower (-34.38) AIC value
# and explains 65.57% more variance than the model with only one intercept and slope 
summary(Model_ML) # slope and intercept of the two regression lines

# Save model with one and two fitted regression lines, for comparison
save(Model_ML, Model, data, tree, file="DataFinalSummary/gradeShiftModel_Feb2025.rdata")

#______________________
# PLOT of grade shifts

# load model and tree 
load("DataFinalSummary/gradeShiftModel_Feb2025.rdata") #objects Model_ML, Model, data, tree

edge_col<-rep("gold2",length(tree$edge.length))
edge_col[getEdges(tree,findMRCA(tree,Storks))]<-"dodgerblue3"
edge_col[getEdges(tree,findMRCA(tree,Telluraves))]<-"dodgerblue3"
edge_col[getEdges(tree,Osprey)]<-"gold2"

identical(data$X1, data$Independent)
birds_bats_abbNames <- rbind(data[,c("Y","X1","Z")], data_bats[,c("Y", "X1","Z")])
birds_bats_abbNames$Z <- paste0(substr(birds_bats_abbNames$Z, 1,1),
                                ". ",
                                sapply(strsplit(birds_bats_abbNames$Z, " "),"[",2))

pdf("Plots/finalPlots/soarersFlappers_gradeShifts_Feb2025.pdf", 7,5)
layout(matrix(c(1,2,2), 1, 3, byrow = TRUE))
plot(tree, edge.col=edge_col, label.offset=5, edge.width=3.5,cex=0.9)
plot(Y~Independent, data, type="n", pch=19, xlim=c(-1.4,3), ylim=c(-22,1.5), #larger x-y range to accommodate bats
     ylab="log(avg Flapping Prob.)", xlab="log(Body mass in Kg)", cex.lab=1.3, cex.axis=1.3)
pGLS.plotGrade("Y","Independent", data[,c("Y","Independent")], tree, model="lambda",group= c(setdiff(Telluraves,Osprey),Storks),col=alpha("dodgerblue3",0.6), cex=1.5, pch=19)
pGLS.plotGrade("Y","Independent", data[,c("Y","Independent")], tree, model="lambda",group= setdiff(All,c(setdiff(Telluraves,Osprey),Storks)),col=alpha("gold2",0.6), cex=1.5, pch=19)
points(Y ~ X1, data_bats, pch=19, col=alpha("black",0.6), cex=1.5)
points(Y ~ X1, data_bats, pch=1, col="black", cex=1.5)
text(Y~X1, birds_bats_abbNames, labels=birds_bats_abbNames$Z, pos=4, cex=1, offset=0.5, srt=30, font=3)
dev.off()


# Add the phylogenetically informed species grouping to the dataset
soarers <- rownames(data[c(setdiff(Telluraves,Osprey),Storks),])
flappers <- rownames(data[setdiff(All,c(setdiff(Telluraves,Osprey),Storks)),])
flappers <- c(flappers, "Eidolon_helvum", "Pteropus_lylei")

species_soarFlap <- rbind(data.frame(species_phy=sub("_"," ",soarers), soarFlap_pgls="soar"),
                          data.frame(species_phy=sub("_"," ",flappers), soarFlap_pgls="flap"))
table(species_soarFlap$soarFlap_pgls)

allSegmDfs_soarFlap <- merge(allSegmDfs_noOutliers, species_soarFlap, by="species_phy", all.x=T)
anyNA(allSegmDfs_soarFlap$soarFlap_pgls)
table(summarise(group_by(allSegmDfs_soarFlap, species_phy), unique(soarFlap_pgls))[,2])

# re-save the original dataset with the added soar/flap categories
saveRDS(allSegmDfs_soarFlap, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")
