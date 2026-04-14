
setwd("...")

#________________________________________________
# S8 - EXPLORE ALTERNATIVE REFERENCE COT LINES
#________________________________________________

# Sequence of possible body masses from our dataset
seq_mass_kg <- round(seq(12, 11000, length.out=100))/1000 # 12 g is the minimum mass in our radar dataset; 11 kg the maximum mass in the our dataset


# COT for all mammalian "athletes" specialised in their locomotory type is the same (although distribution of this COT between maintainance and lcomotion cost are different)
COTmammals_williams <- 10.02 * seq_mass_kg^(-0.31) # from Williams 1999 - cost efficient swimming in marine mammals

#____________________________
# Running/walking
COTrun_alex <- 10.7 * seq_mass_kg^(-0.32) # Alexander 2003 from Taylor 1982

alloSpeed_Heglund_ms <- (5.5 * seq_mass_kg^(0.24)) / 3.6 # allometric scaling of "physiologically similar speed" being running speed at trot-galop gait transition, obtained from HEGLUND, N. C, TAYLOR, C. R. & MCMAHON, T. A. (1974)
#maxSpeed_Garland_ms <- (10^(1.47832 + 0.25892 * (log10(seq_mass_kg)) - 0.06237 * (log10(seq_mass_kg))^2)) / 3.6 # Equation on a dataset of 107 mammals proposed by Garland 1983 suggesting a polynomial relationship for which intermediate-sized animals are the fastest (and suggesting a dynamic similarity between them)
maxSpeed_Garland_small_ms <- 23.3 * seq_mass_kg^(0.225) / 3.6 # Linear relationship from the same paper but for mammals < 300 kg (87 species)

MRrun_taylor_Heglund <- 10.7 * seq_mass_kg^(0.684) * alloSpeed_Heglund_ms + 6.03 * seq_mass_kg^(0.697)# from Taylor 1982, found in Alexander 1998. Original equation in Taylor was expressed in Watts/Kg as: 10.7 Mb^(-0.316) * speed + 6.03 Mb ^(-0.303)
MRrun_taylor_Garland <- 10.7 * seq_mass_kg^(0.684) * maxSpeed_Garland_small_ms + 6.03 * seq_mass_kg^(0.697)# from Taylor 1982, found in Alexander 1998. Original equation in Taylor was expressed in Watts/Kg as: 10.7 Mb^(-0.316) * speed + 6.03 Mb ^(-0.303)

#MRrun_softSub <- 1.13 * MRrun_taylor # White and Yousef 1978

COTrun_taylor_Heglund <- MRrun_taylor_Heglund / (seq_mass_kg * alloSpeed_Heglund_ms)
COTrun_taylor_Garland <- MRrun_taylor_Garland / (seq_mass_kg * maxSpeed_Garland_small_ms)

#COTrun_soft_der <- MRrun_softSub / (seq_mass_kg * walkSpeed)
summary(lm(log(COTrun_taylor_Heglund)~log(seq_mass_kg)))
summary(lm(log(COTrun_taylor_Garland)~log(seq_mass_kg)))

# mass and angle dependent multiplier of COT for steep terrains
terrAngle <- 20
K <- 1 + 0.0550 * seq_mass_kg^(0.222) * terrAngle # from Hasley and White 2017
COTrun_taylor_Garland_steep <- COTrun_taylor_Garland * K
COTrun_taylor_Heglund_steep <- COTrun_taylor_Heglund * K
#COTrun_taylor_Heglund_steep_soft <- (MRrun_taylor_Heglund * 1.3 / (seq_mass_kg * alloSpeed_Heglund_ms)) * K

plot_run <- function() {
  plot(log(COTmammals_williams)~log(seq_mass_kg), col="white", type="l", lwd=1, lty=2, 
       ylab="log COT (J km-1 m-1)", xlab="log Body Mass (Kg)", ylim=c(-2,6))
  lines(log(COTrun_alex)~log(seq_mass_kg), col="gold1", lwd=2)
  lines(log(COTrun_taylor_Heglund)~log(seq_mass_kg), col="#E7BD73", lwd=2, lty=2)
  lines(log(COTrun_taylor_Garland)~log(seq_mass_kg), col="#966919", lwd=2, lty=2)
  lines(log(COTrun_taylor_Heglund_steep)~log(seq_mass_kg), col="#E7BD73", lwd=2, lty=3)
  lines(log(COTrun_taylor_Garland_steep)~log(seq_mass_kg), col="#966919", lwd=2, lty=3)
  #lines(log(COTrun_taylor_Heglund_steep_soft)~log(seq_mass_kg), col="gold3", lwd=1, lty=4)
  legend("bottomleft",
         c("COT run - Alexander 2003/Taylor 1982", 
           "COT run - speeds from Heglund 1974",
           "COT run - speeds from Garland 1983", 
           "COT run - speeds from Heglund 1974 - 20° incline",
           "COT run - speeds from Garland 1983 - 20° incline"),
         lty=c(1,2,2,3,3),
         col=c("gold1", "#E7BD73","#966919","#E7BD73","#966919"), lwd=1, bty="n", cex=1)
}

#___________________
# Flying
COTflight_alex <- (3.6 * seq_mass_kg^(-0.31)) # from Alexander 2003
COTflight_tucker <- 5.2 * seq_mass_kg^(-0.23) # from Tucker 1973 - bird metabolism during flight
#COTflight_schmidt <- 5.40 * seq_mass_kg^(-0.2439) # from Schmidt Nielsen 1972, only 5 values, R2 of 0.3
MRrayner <- 57 * seq_mass_kg^(0.83) #found in Alexander 1998, data from Rayner 1995
maxSpeedsAlexander <- 16 * seq_mass_kg^(0.14) #max range speeds estimated in Alexander 1998, based on data from Pennycuick 1997
soarSpeedAlexander <- 8 * seq_mass_kg^(0.14)
MRalexander_der <- COTflight_alex * maxSpeedsAlexander * seq_mass_kg
# identical to: MRalexander_der2 <- 3.6 * maxSpeedsAlexander * seq_mass_kg^(0.69)
COTalexander_der <- MRrayner / (seq_mass_kg * maxSpeedsAlexander)
BMRnonPass <- 3.6 * seq_mass_kg^(0.72) # from Lasiewski and Dawson 1967; similar to McKechnie 2004 equation, being 3.53 * (mass_kg)^0.669
COTsoar_alex <- (3 * BMRnonPass) / (seq_mass_kg * soarSpeedAlexander)

MRguigueno <- 57.5 * seq_mass_kg^(0.7806)
COTguigueno <- MRguigueno / (seq_mass_kg * maxSpeedsAlexander)
# equivalent to: COTguigueno2 <- as.numeric(exp(predict(modBirds, newdata = data.frame(log_bodyMass_Kg=log(seq_mass_kg))))/(seq_mass_kg * maxSpeedsAlexander))
# COTalexander_der2 <- MRrayner/(seq_mass_kg * 8)
# COTguigueno2 <- MRguigueno/(seq_mass_kg * 8)

eCOT_flap <- 4.89 * seq_mass_kg^(-0.328)
eCOT_soar <- 2.30 * seq_mass_kg^(-0.705)

plot_fly <- function() {
  plot(log(COTmammals_williams)~log(seq_mass_kg), col="white", type="l", lwd=1, lty=2, 
       ylab="log of COT (in J km-1 m-1)", xlab="log of Body Mass (in Kg)", ylim=c(-2,6))
  lines(log(COTflight_alex)~log(seq_mass_kg), col="magenta", lwd=2, lty=1)
  lines(log(COTsoar_alex)~log(seq_mass_kg), col="magenta3", lwd=1, lty=2)
  lines(log(COTflight_tucker)~log(seq_mass_kg), col="darkmagenta", lwd=1, lty=2)
  lines(log(COTguigueno)~log(seq_mass_kg), col="darkorchid", lwd=1, lty=2)
  legend("bottomleft",
         c("COT fly - Alexander 2003/Videler 1993",
           "COT soaring - Alexander 1998",
           "COT fly - Tucker 1973",
           "COT fly - Guigueno 2019"),
         lty=c(1,2,2,2),
         col=c("magenta", "magenta3","darkmagenta","darkorchid"), lwd=1, bty="n", cex=0.6)
}

plot_fly2 <- function() {
  # Flying reference lines
  plot(log(COTmammals_williams)~log(seq_mass_kg), col="white", type="l", lwd=1, lty=2, 
       ylab="log COT (J km-1 m-1)", xlab="log Body Mass (Kg)", ylim=c(-2,6))
  lines(log(COTflight_alex)~log(seq_mass_kg), col="magenta", lwd=2, lty=1)
  lines(log(COTsoar_alex)~log(seq_mass_kg), col="darkorchid", lwd=2, lty=2)
  lines(log(COTflight_tucker)~log(seq_mass_kg), col="#FD3F92", lwd=2, lty=2) 
  lines(log(COTguigueno)~log(seq_mass_kg), col="#9D0759", lwd=2, lty=2)
  # Our eCOT
  lines(log(eCOT_flap)~log(seq_mass_kg), col="black", lwd=2, lty=1)
  lines(log(eCOT_soar)~log(seq_mass_kg), col="grey30", lwd=2, lty=1)
  legend("bottomleft",
         c("COT fly - Alexander 2003/Videler 1993",
           "COT soaring - Alexander 1998",
           "COT fly - Tucker 1973",
           "COT fly - Guigueno 2019",
           "eCOT for flappers - this study", "eCOT for soarers - this study"),
         lty=c(1,2,2,2,1,1),
         col=c("magenta", "darkorchid","#FD3F92","#9D0759","black","grey30"), lwd=1, bty="n", cex=1)
}

#__________________
# Swimming
COTswim_alex <- 1.1 * seq_mass_kg^(-0.38) # Alexander 2003 from Videler 1993
COTswim_brett <- 2.15 * seq_mass_kg^(-0.25) # from Brett 1964
COTmarineMamm_williams <- 7.79 * seq_mass_kg^(-0.29) # data from several sources, equation by Williams 1999 - cost efficient swimming in marine mammals

plot_swim <- function() {
  plot(log(COTmammals_williams)~log(seq_mass_kg), col="white", type="l", lwd=1, lty=2, 
       ylab="log COT (J km-1 m-1)", xlab="log Body Mass (Kg)", ylim=c(-2,6))
  lines(log(COTswim_alex)~log(seq_mass_kg), col="dodgerblue", lwd=2, lty=1)
  lines(log(COTswim_brett)~log(seq_mass_kg), col="#25CEDA", lwd=2, lty=2) # should be same as above but it's not
  lines(log(COTmarineMamm_williams)~log(seq_mass_kg), col="#2338AF", lwd=2, lty=2) # resulting higher than running
  legend("bottomleft",
         c("COT swim - Alexander 2003/Videler 1993", 
           "COT fish - Brett 1964", 
           "COT marine mammals - Wilson 1999"),
         lty=c(1,2,2),
         col=c("dodgerblue", "#25CEDA","#2338AF"), lwd=1, bty="n", cex=1)
}


#_______________________
# Plot all together
#_______________________


pdf("Revision/NewSupplFigures/alternative_COT_lines.pdf", width = 11,height = 8)

layout(matrix(c(1,2,3,4,4,4), nrow=2, byrow=T), widths = c(1,1,1,1.5))
#layout.show(4)
# single locomotion mode panels
plot_run()
plot_fly()
plot_swim()
# One plot with all reference COT lines as 4th panel
# flight
plot(log(COTmammals_williams)~log(seq_mass_kg), col="darkgrey", type="l", lwd=1, lty=2, 
     ylab="log of COT (in J km-1 m-1)", xlab="log of Body Mass (in Kg)", ylim=c(-2,6))
lines(log(COTflight_alex)~log(seq_mass_kg), col="magenta", lwd=2, lty=1)
lines(log(COTsoar_alex)~log(seq_mass_kg), col="magenta3", lwd=1, lty=2)
lines(log(COTflight_tucker)~log(seq_mass_kg), col="darkmagenta", lwd=1, lty=2)
lines(log(COTguigueno)~log(seq_mass_kg), col="darkorchid", lwd=1, lty=2)
# swim
lines(log(COTswim_alex)~log(seq_mass_kg), col="dodgerblue", lwd=2, lty=1)
lines(log(COTswim_brett)~log(seq_mass_kg), col="dodgerblue3", lwd=1, lty=2) # should be same as above but it's not
lines(log(COTmarineMamm_williams)~log(seq_mass_kg), col="dodgerblue4", lwd=1, lty=2) # resulting higher than running
# run
lines(log(COTrun_alex)~log(seq_mass_kg), col="gold1", lwd=2)
lines(log(COTrun_taylor_Heglund)~log(seq_mass_kg), col="gold3", lwd=1, lty=2)
lines(log(COTrun_taylor_Garland)~log(seq_mass_kg), col="gold4", lwd=1, lty=2)
# eCOT
lines(log(eCOT_flap)~log(seq_mass_kg), col="black", lwd=2, lty=1)
lines(log(eCOT_soar)~log(seq_mass_kg), col="grey30", lwd=2, lty=1)
legend(-4.7,6.2,
       c("COT run - Alexander 2003/Taylor 1982", "COT fly - Alexander 2003", "COT swim - Alexander 2003/Videler 1993",
         "eCOT for flappers - this study", "eCOT for soarers - this study"),
       lty=c(1,1,1,1,1),col=c("gold1","magenta","dodgerblue","black","grey30"), lwd=2, bty="n", cex=0.6)
legend(-2,6.2,
       c("COT run - speeds from Heglund 1974","COT run - speeds from Garland 1983",
         "COT fly - Tucker 1973", "COT fly - Guigueno 2019","COT soaring - Alexander 1998", 
         "COT fish - Brett 1964", "COT marine mammals - Wilson 1999",
         "COT run low speed - Alexander 1998","COT run soft sub - Alexander 1998/White Yousef 1978",
         "COT all mammals - Williams 1999"),
       lty=c(2,2,2,2,2,2,2,2),
       col=c("gold3","gold4",
             "darkmagenta","darkorchid","magenta3",
             "dodgerblue3","dodgerblue4",
             "darkgrey"), lwd=1, bty="n", cex=0.6)
dev.off()


pdf("Revision/NewSupplFigures/alternative_COT_lines_onlyABC.pdf", width = 13,height = 5)

layout(matrix(c(1,2,3), nrow=1, byrow=T), widths = c(1,1,1))

plot_run()
plot_fly2()
plot_swim()

dev.off()


### Sanity checks
# par(mfrow=c(1,2))
# plot(MRguigueno~seq_mass_kg, type = "l", lwd=2, col="blue",
#      xlab="Body mass (kg)", ylab="Flapping MR (W)")
# lines(MRrayner~seq_mass_kg, col = "red", lwd=2)
# lines(MRalexander_der~seq_mass_kg, col = "darkgrey", lwd=2)
# legend("topleft",c("COT Alexander 2003 * maxSpeed * Mass", "MR Guigueno 2019", "MR Alexander 1998"),
#        lty=1,col=c("darkgrey","blue","red"), lwd=2,
#        bty="n", cex=0.7)
# 
# plot(log(COTalexander_der)~log(seq_mass_kg), col="red", type="l", ylab="COT (J km-1 m-1)", lwd=2, ylim=c(0,5))
# lines(log(COTflight_alex)~log(seq_mass_kg), col="darkgrey", lwd=2)
# lines(log(COTguigueno)~log(seq_mass_kg), col="blue", lwd=2)
# lines(log(COTswim_alex)~log(seq_mass_kg), col="cyan", lwd=2)
# lines(log(COTrun_alex)~log(seq_mass_kg), col="magenta", lwd=2)
# legend("topleft",c("COT Alexander 2003","MR Guigueno 2019 / (maxSpeed * Mass)", "MR Alexander 1998 / (maxSpeed * Mass)"),
#        lty=1,col=c("darkgrey", "blue", "red"), lwd=2,
#        bty="n", cex=0.7)
###


#________________________________________________
# S9 - EFFECT of BMR estimation on COT estimation
#________________________________________________

#______________________________________________________________
## 1. Plot different calculations of BMR relative to body mass

# Sequence of possible body masses from our dataset
mass_gr <- round(seq(12, 11000, length.out=100)) # 12 g is the minimum mass in our radar dataset; 11 kg the maximum mass in the our dataset

# Different BMR calculations, all in W (i.e. J/s)
BMR_Londono_tropical <- 0.044 * (mass_gr)^0.589
BMR_Londono_temperate <- 0.023 * (mass_gr)^0.729
BMR_Londono_pass <- 0.045 * (mass_gr)^0.627
BMR_Londono_nonPass <- 0.021 * (mass_gr)^0.724 # using mass in kg it would be 3.11 * (mass_kg)^0.724
BMR_McKechnie_birds <- 0.0346 * (mass_gr)^0.669 # same as: 10^(-1.461 + 0.669 * log10(mass_gr)); using mass in kg it would be 3.53 * (mass_kg)^0.669
BMR_Speakman_bats <- exp(1.0895 + 0.744 * log(mass_gr)) * 0.005583

## Plot the different BMR lines against body mass values (use plot in composite below)
plot_BMR <- function(){
  plot(mass_gr/1000, BMR_McKechnie_birds, type = "l", lwd=2, xlab="Body mass (kg)", ylab="BMR (W)" )
  lines(mass_gr/1000, BMR_Speakman_bats, col = "black", lty = 3, lwd=3)
  lines(mass_gr/1000, BMR_Londono_pass, col = "orange", lty = 2, lwd=2)
  lines(mass_gr/1000, BMR_Londono_nonPass, col = "forestgreen", lty = 2, lwd=2)
  lines(mass_gr/1000, BMR_Londono_tropical, col = "magenta4", lty = 2, lwd=2)
  lines(mass_gr/1000, BMR_Londono_temperate, col = "darkgrey", lty = 2, lwd=2)
  legend("bottomright",c("McKechnie 2004 (birds)","Speakman 2003 (bats)","Londono 2015 (passerines)","Londono 2015 (non-passerines)",
                         "Londono 2015 (tropical)","Londono 2015 (temperate)"),
         lty=c(1,3,2,2,2,2),col=c("black","black","orange","forestgreen","magenta4","darkgrey"), lwd=2,
         bty="n", cex=0.7)
}

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

#______________________________________________________________
## 2. Plot effect of BMR calculation on COT calculation
## for different body masses and proportions of flapping

library(scales)


dummyDf <- as.data.frame(cbind(mass_kg=mass_gr/1000, BMR_Londono_nonPass, BMR_Londono_pass, BMR_Londono_temperate, BMR_Londono_tropical,
                               BMR_McKechnie_birds, BMR_Speakman_bats))

# take intermediate and extreme values of prop flap to represent soaring/flapping strategies
propFlap <- c(0.95, 0.5, 0.05) 

dummyDf_prop <- expand.grid(propFlap=propFlap, mass_kg=mass_gr/1000)
dummyDf_prop <- merge(dummyDf_prop, dummyDf, by="mass_kg")

# for speed we don't take all combinations as certain speed-body mass combinations are unrealistic.
# Instead we use the maximum range speeds proposed by Alexander 1998
dummyDf_prop$speeds <- 16 * dummyDf_prop$mass_kg^(0.14) #max range speeds estimated in Alexander 1998, based on data from Pennycuick 1997


# We consider a range of speeds between 3 and 15 m/s (median of the dataset of this study)
# speeds <- c(seq(3, 15, length.out=5), 8) # add 8 m/s as it is the median speed in both datasets
# dummyDf_prop <- do.call(rbind, lapply(speeds, function(s){
#   dummyDf_prop$speed <- s
#   return(dummyDf_prop)}))

head(dummyDf_prop)
summary(dummyDf_prop)

# calculate MR of flapping based on body mass
load("DataFinalSummary/flappingModel_KylesData.RData") #object modBirds, from script 0C
# Guigueno gives us a MR in Watts (J/s) and it tells us that the power is about 57.5 * Mass^(0.7806)
dummyDf_prop$MR_flap_guigueno2019 <- exp(predict(modBirds, newdata = data.frame(log_bodyMass_Kg=log(dummyDf_prop$mass_kg))))

# Alexander gives us a COT in J/(kg m) of 3.6 * Mass^(-0.31)
# to convert to Watts (J/s) we need (3.6 * Mass^(-0.31)) * Mass * Speed
# which simplified becomes 3.6 * Speed * Mass^(0.69). If we use the median speed in our dataset (8 m/s) we obtain a coefficient of 28.8 * Mass^(0.69)
dummyDf_prop$MR_flap_alexander2003 <- 3.6 * dummyDf_prop$speed * (dummyDf_prop$mass_kg)^(0.69)

# use this plot in composite below
plot_MR <- function(){
  plot(MR_flap_guigueno2019~mass_kg, data=dummyDf_prop, type = "l", lwd=2, 
       xlab="Body mass (kg)", ylab="Flapping MR (W)", ylim=c(1,380))
  lines(MR_flap_alexander2003~mass_kg, data=dummyDf_prop, col = "darkgrey", lty = 2, lwd=2)
  legend("bottomright",c("Guigueno 2019","Alexander 2003"),
         lty=c(1,2),col=c("black","darkgrey"), lwd=2,
         bty="n", cex=0.7)
}


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
dummyAlexander$overall_COT_Speakman_bats_Jkgm <- (cost_flap_alexander + cost_pass_Speakman_bats) / dummyAlexander$mass_kg / dummyAlexander$speed

# Bind the two
dummyDf_prop_final <- rbind(dummyGuigueno, dummyAlexander)

sub <- dummyDf_prop_final


pdf("Revision/NewSupplFigures/EffectOnCOT_differentBMR_differentMR_speedAlexander1998.pdf", 10, 10)

layout(matrix(c(1,2,3,4), nrow=2, byrow=T))

plot_BMR()

plot_MR()

plot(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], type="l", 
     xlab="Log body mass (kg)", ylab="Log overall COT (J kg-1 m-1))", ylim=c(-1.5,3.5), col="black", lwd=3, main = "Flapping MR from Alexander 2003")
lines(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Alexander 2003",], col = alpha("black", 0.8), lty = 3, lwd=3)
lines(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Alexander 2003",], col = alpha("black", 0.8), lty = 2, lwd=3)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("orange", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Alexander 2003",], col = alpha("orange", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Alexander 2003",], col = alpha("orange", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("forestgreen", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Alexander 2003",], col = alpha("forestgreen", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Alexander 2003",], col = alpha("forestgreen", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("magenta4", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Alexander 2003",], col = alpha("magenta4", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Alexander 2003",], col = alpha("magenta4", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Alexander 2003",], col = alpha("darkgrey", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Alexander 2003",], col = alpha("darkgrey", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Alexander 2003",], col = alpha("darkgrey", 0.8), lty = 2, lwd=2)
legend("bottomleft",c("McKechnie 2004 (birds)","Londono 2015 (passerines)","Londono 2015 (non-passerines)",
                       "Londono 2015 (tropical)","Londono 2015 (temperate)",
                       "pFlap = 0.05 ('soarers')","pFlap = 0.50","pFlap = 0.95 ('flappers')"),
       lty=c(1,1,1,1,1,1,3,2),col=c("black","orange","forestgreen","magenta4","darkgrey","black","black","black"), lwd=2,
       bty="n", cex=0.7)

plot(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], type="l", 
     xlab="Log body mass (kg)", ylab="Log overall COT (J kg-1 m-1)", ylim=c(-1.5,3.5), col="black", lwd=3, main="Flapping MR from Guigueno 2019")
lines(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Guigueno 2019",], col = alpha("black", 0.8), lty = 3, lwd=3)
lines(log(overall_COT_McKechnie_birds_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Guigueno 2019",], col = alpha("black", 0.8), lty = 2, lwd=3)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("orange", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Guigueno 2019",], col = alpha("orange", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_pass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Guigueno 2019",], col = alpha("orange", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("forestgreen", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Guigueno 2019",], col = alpha("forestgreen", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_nonPass_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Guigueno 2019",], col = alpha("forestgreen", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("magenta4", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Guigueno 2019",], col = alpha("magenta4", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_tropical_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Guigueno 2019",], col = alpha("magenta4", 0.8), lty = 2, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.05 & sub$MRsource == "Guigueno 2019",], col = alpha("darkgrey", 0.8), lty = 1, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.5 & sub$MRsource == "Guigueno 2019",], col = alpha("darkgrey", 0.8), lty = 3, lwd=2)
lines(log(overall_COT_Londono_temperate_Jkgm)~log(mass_kg), data=sub[sub$propFlap == 0.95 & sub$MRsource == "Guigueno 2019",], col = alpha("darkgrey", 0.8), lty = 2, lwd=2)

dev.off()


# isolate pure flapping case to check the difference in slope associated to using these two methods
sub_test <- dummyDf_prop_final[dummyDf_prop_final$propFlap == 0.95,]

lm(log(overall_COT_McKechnie_birds_Jkgm) ~ log(mass_kg),
   data=sub_test[sub_test$MRsource=="Guigueno 2019",])

lm(log(overall_COT_McKechnie_birds_Jkgm) ~ log(mass_kg),
   data=sub_test[sub_test$MRsource=="Alexander 2003",])


# # Reshape df to fit ggplot (different COT estimations in rows instead of columns)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# df_long <- dummyDf_prop_final %>%
#   pivot_longer(cols = starts_with("overall_COT"), 
#                names_to = "BMR_source", 
#                values_to = "overall_COT_Jkgm")
# 
# # summarise COT values for each group
# df_summarized <- df_long %>%
#   group_by(BMR_source, MRsource, log(mass_kg)) %>%
#   summarize(
#     mean_COT = mean(log(overall_COT_Jkgm)),   # Average COT per bodymass per group
#     sd_COT = sd(log(overall_COT_Jkgm)),       # Standard deviation
#     se_COT = sd_COT / sqrt(n()),              # Standard error
#     ci_lower = mean_COT - 1.96 * se_COT,      # Lower 95% CI
#     ci_upper = mean_COT + 1.96 * se_COT       # Upper 95% CI
#   )
# table(df_summarized$BMR_source)
# table(df_summarized$MRsource)
# 
# df_summarized$BMR_source <- gsub("overall_COT_|_Jkgm","",df_summarized$BMR_source)
# 
# # Plot
# ggplot(df_summarized, aes(x = `log(mass_kg)`, y = mean_COT, color = BMR_source, group = BMR_source)) +
#   geom_line(size = 2) +                         # Line plot for each COT group
#   geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = BMR_source), 
#               alpha = 0.05, linetype = "blank") + # Confidence intervals
#   labs(x = "log of Body Mass", y = "log of overall COT (J/(kg * m))", 
#        #title = "COT vs Body Mass with Confidence Intervals",
#        color = "COT Group", fill = "COT Group") +
#   theme_minimal() +
#   facet_wrap(~MRsource) +           # Create facets for each metabolic rate source
#   theme(legend.position = "top")    # Optional: position legend at the top
# ggsave("Plots/finalPlots/EffectOnCOT_differentBMR_differentMR_allPropsAllSpeeds_gg.pdf", width =13,height=7) 
# 

#_______________________________________________
# EXPLORE EFFECT OF CALCULATION IN OUR DATASET!
#_______________________________________________

#__________________________________
## Summary of our species, both gps and radar
## in terms of body mass, proportion of flapping and speed

# Import our data, both gps and radar
gps <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")
radar <- readRDS("DataFinalSummary/RADARdata_finalSummaryDataset_perEcho_COTvariables_WFF-month_echoDurFilter.rds")

# Import kile's MR model
load("DataFinalSummary/flappingModel_KylesData.RData") #object modBirds, from script 0C

# summarise speed prop of flapping and body mass
speciesSummary <- rbind(group_by(gps, species=species) %>% summarise(mass_kg=unique(Body_mass_kg), 
                                                                     avgSpeed=mean(avg_grSpeed_ms),
                                                                     avgPflap=mean(avg_probFlap),
                                                                     soarFlap=unique(soarFlap_pgls),
                                                                     group="Non-passerines"),
                        group_by(radar, species=WFFmonth_species) %>% summarise(mass_kg=unique(meanMass_kg), 
                                                                                avgSpeed=mean(Speed.ms),
                                                                                avgPflap=mean(propPulse),
                                                                                soarFlap="flap",
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
    speciesSummary$BMR[speciesSummary$group=="Non-passerines"] <- 0.021 * (speciesSummary$mass_kg*1000)[speciesSummary$group=="Non-passerines"]^0.724
  }
  speciesSummary$source <- source
  
  costPass <- 2 * speciesSummary$BMR * (1 - speciesSummary$avgPflap)
  costFlap_guig <- (exp(predict(modBirds, newdata = data.frame(log_bodyMass_Kg=speciesSummary$logMass_kg)))) * speciesSummary$avgPflap
  #costFlap_alex <- (3.6*(speciesSummary$mass_kg)^(-0.31)) * speciesSummary$avgPflap
  costFlap_alex <- 3.6 * speciesSummary$avgSpeed * (speciesSummary$mass_kg)^(0.69) * speciesSummary$avgPflap
  
  speciesSummary$overallCOT_guig <- (costPass + costFlap_guig) / speciesSummary$mass_kg / speciesSummary$avgSpeed
  speciesSummary$overallCOT_alex <- (costPass + costFlap_alex) / speciesSummary$mass_kg / speciesSummary$avgSpeed
  return(speciesSummary)
}))
head(speciesSummary_cot)
speciesSummary_cot$source <- factor(speciesSummary_cot$source, levels = c("McKechnie","Londono_tropical","Londono_temperate","Londono_pass"))
saveRDS(speciesSummary_cot, file="DataFinalSummary/speciesSummary_averageBMR_overallCOT_different sources.rds")


# Models
mod_guig_soar <- lm(log(overallCOT_guig) ~ logMass_kg * source, data = speciesSummary_cot[speciesSummary_cot$soarFlap=="soar",])
summary(mod_guig_soar)
mod_alex_soar <- lm(log(overallCOT_alex) ~ logMass_kg * source, data = speciesSummary_cot[speciesSummary_cot$soarFlap=="soar",])
summary(mod_alex_soar)
mod_guig_flap <- lm(log(overallCOT_guig) ~ logMass_kg * source, data = speciesSummary_cot[speciesSummary_cot$soarFlap=="flap",])
summary(mod_guig_flap)
mod_alex_flap <- lm(log(overallCOT_alex) ~ logMass_kg * source, data = speciesSummary_cot[speciesSummary_cot$soarFlap=="flap",])
summary(mod_alex_flap)

# Create a function that extracts the coefficients and creates a table with the slopes for each model
make_COT_table <- function(mod_guig, mod_alex){
  
  coefGuig <- coefficients(mod_guig)
  coefAlex <- coefficients(mod_alex)
  
  data.frame(sourceMB = c("Guigueno","Alexander"),
    slope_McKechnie = c(
      coefGuig["logMass_kg"],
      coefAlex["logMass_kg"]),
    slope_LondonTropical = c(
      coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_tropical"],
      coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_tropical"]),
    slope_LondonTemperate = c(
      coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_temperate"],
      coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_temperate"]),
    slope_LondonPass = c(
      coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_pass"],
      coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_pass"]))
}

# Create two separate tables:
# Soarers
COT_table_soar <- make_COT_table(mod_guig_soar, mod_alex_soar)
# Flappers
COT_table_flap <- make_COT_table(mod_guig_flap, mod_alex_flap)

COT_table_soar
COT_table_flap

# coefGuig <- coefficients(mod_guig)
# coefAlex <- coefficients(mod_alex)
# 
# COT_table <- data.frame(sourceMB=c("Guigueno","Alexander"),
#            slope_McKechnie=c(coefGuig["logMass_kg"], coefAlex["logMass_kg"]),
#            slope_LondonTropical=c(coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_tropical"],
#                                   coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_tropical"]),
#            slope_LondonTemperate=c(coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_temperate"],
#                                   coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_temperate"]),
#            slope_LondonPass=c(coefGuig["logMass_kg"] + coefGuig["logMass_kg:sourceLondono_pass"],
#                                   coefAlex["logMass_kg"] + coefAlex["logMass_kg:sourceLondono_pass"]))
# COT_table


# Export
write.csv(COT_table_soar, "Tables/COT_BMR_effects_soarers.csv", row.names = FALSE)
write.csv(COT_table_flap, "Tables/COT_BMR_effects_flappers.csv", row.names = FALSE)
