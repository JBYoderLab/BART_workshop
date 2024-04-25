# Species distribution modeling with RI BARTs in R
# Jeremy B. Yoder, 19 Apr 2024

# Clears the environment and load key packages
rm(list=ls())

# This script assumes your working directory is the project directory; you may want to change this
# setwd("~/Documents/Active_projects/BART_workshop")

library("tidyverse") 
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

library("embarcadero")
# devtools::install_github("cjcarlson/embarcadero")

#-----------------------------------------------------------
# Load and prepare data

# Joshua tree occurrence records, and pseudoabsences, with environmental data
# (This was generated in the code for topic 01, `BART_SDM.R`)
PA <- read.csv("data/JT_presence-pseudoabsence-envs.txt") %>% filter(!is.na(MAT))
glimpse(PA)

# Visualize the points on a map of the US Southwest
# As you can see, there's quite a bit of heterogeneity in sampling density!
{png("topics/02_RI_models/JT_presence_pseudoabsence_map.png", height=750, width=750)
ggplot() + 
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill="antiquewhite1") + 
  geom_point(data=PA, aes(x=lon, y=lat, color=factor(JT)), size=1) + 
  coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) + 
  theme_bw(base_size=18) + theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), legend.position="none")
}
dev.off()

# Environmental layers we'll use for modeling
envs <- brick("data/Mojave_BioClim-elev-states.grd")

envs
names(envs)

# These are the Bioclim layers at highest resolution (30-second)
plot(envs[["MAT"]]) # Mean annual temperature
plot(envs[["AP"]]) # Annual precipitation

# As a reminder, the full set of bioclim variables are

# BIO1 = Annual Mean Temperature (MAT)
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) (MDR)
# BIO3 = Isothermality (BIO2/BIO7) (×100) (ITH)
# BIO4 = Temperature Seasonality (standard deviation ×100) (TS)
# BIO5 = Max Temperature of Warmest Month (MWaMT)
# BIO6 = Min Temperature of Coldest Month (MCMT)
# BIO7 = Temperature Annual Range (BIO5-BIO6) (TAR)
# BIO8 = Mean Temperature of Wettest Quarter (MWeQT)
# BIO9 = Mean Temperature of Driest Quarter (MDQT)
# BIO10 = Mean Temperature of Warmest Quarter (MWaQT)
# BIO11 = Mean Temperature of Coldest Quarter (MCQT)
# BIO12 = Annual Precipitation (AP)
# BIO13 = Precipitation of Wettest Month (PWeM)
# BIO14 = Precipitation of Driest Month (PDM)
# BIO15 = Precipitation Seasonality (Coefficient of Variation) (PS)
# BIO16 = Precipitation of Wettest Quarter (PWeQ)
# BIO17 = Precipitation of Driest Quarter (PDQ)
# BIO18 = Precipitation of Warmest Quarter (PWaQ)
# BIO19 = Precipitation of Coldest Quarter (PCQ)

# We also have elevation in meters, at the same resolution
plot(envs[["ELEV"]])

# And state boundaries, rendered as numeric values because that's what we can use for modeling
plot(envs[["STATE"]])


#-------------------------------------------------------------------------
# Fit a species distribution model with BARTs

# vector of the Bioclim variables
xnames <- c("MAT", "MDR", "ITH", "TS", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT", "MCQT", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")

# Fit a stepwise BART model. This will be slow!
jtBART.step <- bart.step(y.data=as.numeric(PA[,"JT"]), x.data=PA[,xnames], full=FALSE, quiet=TRUE)

# It may be useful to save the model object for later work. 
# First, we need to "touch" the model state to make sure the trees are saved with the rest of the model.
# This is a `dbarts` feature that is supposed to save storage.
invisible(jtBART.step$fit$state)
# Then write out the model to an Rdata file
write_rds(jtBART.step, file="output/models/jtBART.step.rds") 
jtBART.step <- read_rds(file="output/models/jtBART.step.rds")

summary(jtBART.step)

stepX <- attr(jtBART.step$fit$data@x, "term.labels")
stepX # these are the predictors kept in the stepwise model


#-------------------------------------------------------------------------
# Fit a species distribution model with random-intercept BARTs

# Adding a random-intercept effect to a BART model lets us account for 
# heterogeneity among data sources. Recall that our occurrence records
# are much denser in Nevada --- we can treat state as an RI predictor to
# account for this.

# We can simply add the RI predictor into a model with predictors found by
# stepwise fitting ...
stepX # these are the predictors kept in the stepwise model, loaded above

jtRIBART <- rbart_vi(as.formula(paste(paste('JT', paste(stepX, collapse=' + '), sep = ' ~ '), 'STATE', sep=' - ')), data = PA, group.by = PA[,'STATE'], n.chains = 1, k = 2, power = 2, base = 0.95, keepTrees = TRUE)

invisible(jtRIBART$fit[[1]]$state)
write_rds(jtRIBART, file="output/models/jtRIBART.rds") 
jtRIBART <- read_rds(file="output/models/jtRIBART.rds")

summary(jtRIBART)

# And then we can inspect the estimated RI effects to see whether they 
# make sense to include:
{png("topics/02_RI_models/BART_RI_estimates.png", width=500, height=350)
plot.ri(jtRIBART, temporal=FALSE) + 
scale_x_discrete(labels=c("AZ", "CA", "NV", "UT")) + # relabel numeric coding
labs(title="Random intercept effects of observation year", x="RI term") + theme_bw(base_size=18)
}
dev.off()


# It may be instead that including the RI term changes the effects of other
# enough to change the outcome of stepwise fitting. We can re-do stepwise fitting
# with the RI term included to account for this:

jtRIBART.step <- bart.step(y.data=as.numeric(PA[,"JT"]), x.data=PA[,xnames], ri.data=PA[,"STATE"], full=FALSE, quiet=TRUE) # (This will be slow, again.)

invisible(jtRIBART.step$fit[[1]]$state)
write_rds(jtRIBART.step, file="output/models/jtRIBART.step.rds") 
jtRIBART.step <- read_rds(file="output/models/jtRIBART.step.rds")

# Does this end up with the same predictors as the non-RI stepwise model?
RIstepX <- attr(jtRIBART.step$fit[[1]]$data@x, "term.labels")
RIstepX

# What does the new estimated RI effect look like?
{png("topics/02_RI_models/BART_RI.step_estimates.png", width=500, height=350)
plot.ri(jtRIBART.step, temporal=FALSE) + 
scale_x_discrete(labels=c("AZ", "CA", "NV", "UT")) + # relabel numeric coding
labs(title="Random intercept effects of observation year", x="RI term") + theme_bw(base_size=18)
}
dev.off()

#-------------------------------------------------------------------------
# Make predictions while controlling for the RI term

# We can make predictions based on the marginal effects of the non-RI predictors
pred.ri0 <- predict(jtRIBART.step, envs, ri.data=envs[["STATE"]], ri.name="STATE", splitby=20, ri.pred=FALSE)
# (You have to specify values for the RI predictor, but then we cancel it out with ri.pred=FALSE)


# Create a map of predicted presence
pred.ri0.df <- cbind(coordinates(pred.ri0), as.data.frame(pred.ri0)) %>% rename(prJT = layer, lon=x, lat=y)
glimpse(pred.ri0.df)

{png("topics/02_RI_models/jtRIBART.step_predicted.png", width=1000, height=1000)

ggplot() + 
  geom_tile(data=pred.ri0.df, aes(x=lon, y=lat, fill=prJT)) +
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="white") + 
  scale_fill_distiller(type="seq", palette="YlGn", direction=1, name="Probability of occurrence") + labs(x="Longitude", y="Latitude") + 
  coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
  theme_bw(base_size=24) + theme(legend.position="bottom", legend.key.width=unit(0.05, "npc"))

}
dev.off()

#-------------------------------------------------------------------------
# Compare the RI BART and BART results

# Make predictions based on the original stepwise-fitted model for comparison
pred.step <- predict(jtBART.step, envs[[stepX]], splitby=20)

RIvBART <- pred.ri0 - pred.step

RIvBART.df <- cbind(coordinates(RIvBART), as.data.frame(RIvBART)) %>% rename(prDiff = layer, lon = x, lat = y)
glimpse(RIvBART.df)


{png("topics/02_RI_models/RIBARTvBART_predictions.png", width=1000, height=1000)

ggplot() + 
  geom_tile(data=RIvBART.df, aes(x=lon, y=lat, fill=prDiff)) +
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="black") + 
  scale_fill_distiller(type="seq", palette="PuOr", direction=1, name="RI BART - BART\nprobability of occurrence") + labs(x="Longitude", y="Latitude") + 
  coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
  theme_bw(base_size=24) + theme(legend.position="bottom", legend.key.width=unit(0.05, "npc"))

}
dev.off()




