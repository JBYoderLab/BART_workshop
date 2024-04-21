# Species distribution modeling with BARTs in R
# Jeremy B. Yoder, 19 Apr 2024

# Clears the environment and load key packages
rm(list=ls())

# This script assumes your working directory is the project directory; you may want to change this
# setwd("~/Documents/Active_projects/BART_workshop")

library("tidyverse") 
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

library("dismo")

library("gbm3")
# devtools::install_github("gbm-developers/gbm3")

library("embarcadero")
# devtools::install_github("cjcarlson/embarcadero")

#-----------------------------------------------------------
# Load data

# Joshua tree occurrence records, from published data https://doi.org/10.5061/dryad.6s67t
# ... supplemented with validated iNaturalist records and field notes, and thinned a bit.

jtOcc <- read.csv("data/JT_obs.txt", sep="\t")

glimpse(jtOcc)

# Visualize the points on a map of the US Southwest
# As you can see, there's quite a bit of heterogeneity in sampling density.
{png("topics/01_binomial_SDM/JT_occurrences_map.png", height=750, width=750)
ggplot() + 
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill="antiquewhite1") + 
  geom_point(data=jtOcc, aes(x=lon, y=lat), color="forestgreen", size=1) + 
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


#-----------------------------------------------------------
# Create pseudoabsences to match the occurrences

occ_sf <- st_as_sf(jtOcc, coords=c("lon", "lat"), crs=4326) # coordinates are in degrees
occ_sf <- st_transform(occ_sf, crs=3857) # units in meters

x.5k <- st_union(st_buffer(occ_sf, 5000)) # polygons based on 5km radii around occ_thin
x.40k <- st_union(st_buffer(occ_sf, 40000)) # polygons based on 40km radii
x.donuts <- st_difference(x.40k, x.5k) # difference between above

ggplot() + geom_sf(data=x.donuts) 
# A conservative region for pseudoabsences: Not within 5km of a presence point,
# but not more than 50km from a presence point.

# Now, we draw random points within the "donut" regions --- the same number as 
# we have presence points, for balance
pseudabs <- st_sample(x.donuts, nrow(jtOcc))

# Visualize on our map
{png("topics/01_binomial_SDM/JT_presence_pseudoabsence_map.png", height=750, width=750)
ggplot() + 
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill="antiquewhite1") + 
  geom_sf(data=occ_sf, color="blue", size=1) + 
  geom_sf(data=pseudabs, color="red", size=1) + 
  coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) + 
  theme_bw(base_size=18) + theme(axis.title=element_blank(), panel.background=element_rect(fill="slategray3"), legend.position="none")
}
dev.off()

#-------------------------------------------------------------------------
# Assemble environmental values for our presences and pseudoabsences 

pres_envs <- data.frame(raster::extract(envs, jtOcc[,c("lon","lat")]))
glimpse(pres_envs) # BioClim values for the presence locations

pa_sp <- as_Spatial(st_transform(pseudabs, crs=4326))
pa_envs <- data.frame(raster::extract(envs, pa_sp))
glimpse(pa_envs) # And ditto for the pseudoabsences

# Now, assemble a single data frame with presence/absence coordinates and the
# associated BioClim values:

PA <- rbind(data.frame(lon=jtOcc$lon, lat=jtOcc$lat, JT=1, pres_envs), 
          data.frame(lon=coordinates(pa_sp)[,1], lat=coordinates(pa_sp)[,2], JT=0, pa_envs)) %>% 
          mutate(id = row_number()) %>% filter(!is.na(MAT))

glimpse(PA) # should be ~twice the size of the presence data set

# You may want to write out this data frame for later /read back in 
write.table(PA, "data/JT_presence-pseudoabsence-envs.txt", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
PA <- read.csv("data/JT_presence-pseudoabsence-envs.txt")
glimpse(PA)

#-------------------------------------------------------------------------
# Fit a species distribution model with Boosted Regression Trees

# vector of the Bioclim variables
xnames <- c("MAT", "MDR", "ITH", "TS", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT", "MCQT", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")

# Fit a BRT model with defaults provided by `dismo`
jtBRM <- gbm.step(data=PA, gbm.x=xnames, gbm.y="JT", family = "bernoulli", tree.complexity = 5, learning.rate = 0.05, bag.fraction = 0.5)

jtBRM

# Examine predictor contributions and write out a figure
{png("topics/01_binomial_SDM/jtBRM_predictor_contributions.png", width=500, height=750)
summary(jtBRM)
}
dev.off()

# Try to simplify the model from the full set of 19 bioclim variables plus state ID
jtBRM.simp <- gbm.simplify(jtBRM)

{png("topics/01_binomial_SDM/jtBRM.simp_predictor_contributions.png", width=500, height=750)
summary(jtBRM.simp)
}
dev.off()

# Fit a new stepwise model using the optimal predictors identified (in my run, this is 9, your mileage may vary)
jtBRM2 <- gbm.step(data=PA, gbm.x=jtBRM.simp$pred.list[[9]], gbm.y="JT", family = "bernoulli", tree.complexity = 5, learning.rate = 0.05, bag.fraction = 0.5)

# OR simply jump ahead by using the predictor subset I found in my own analysis:
jtBRM2 <- gbm.step(data=PA, gbm.x=c("MDR", "TS", "MCMT", "TAR", "MWeQT", "MCQT", "AP", "PS", "PDQ", "PWaQ"), gbm.y="JT", family = "bernoulli", tree.complexity = 5, learning.rate = 0.05, bag.fraction = 0.5)

write_rds(jtBRM2, file="output/models/jtBRM.rds") # save the results
jtBRM2 <- read_rds("output/models/jtBRM.rds") # reload later


{png("topics/01_binomial_SDM/jtBRM2_predictor_contributions.png", width=500, height=750)
summary(jtBRM2)
}
dev.off()


# Examine predictor partial effects
{png("topics/01_binomial_SDM/jtBRM_partials.png", width=1000, height=750)
gbm.plot(jtBRM2, n.plots=10, plot.layout=c(3,4), write.title=FALSE)
}
dev.off()


# And, finally, predict to the full range of the Mojave:
jtBRM.pred <- predict(envs, jtBRM2, n.trees=jtBRM2$gbm.call$best.trees, type="response")

jtBRM.pred.df <- cbind(coordinates(jtBRM.pred), as.data.frame(jtBRMpred)) %>% rename(prJT = layer, lon=x, lat=y)
glimpse(jtBRMpred.df)

{png("topics/01_binomial_SDM/jtBRM_predicted.png", width=1000, height=1000)

ggplot() + 
  geom_tile(data=jtBRM.pred.df, aes(x=lon, y=lat, fill=prJT)) +
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="white") + 
  scale_fill_distiller(type="seq", palette="YlGn", direction=1, name="Probability of occurrence") + labs(x="Longitude", y="Latitude") + 
  coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
  theme_bw(base_size=24) + theme(legend.position="bottom", legend.key.width=unit(0.05, "npc"))

}
dev.off()


#-------------------------------------------------------------------------
# Fit a species distribution model with BARTs

# vector of the Bioclim variables
xnames <- c("MAT", "MDR", "ITH", "TS", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT", "MCQT", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")

# Fit a simple BART model with all the predictors
jtBART <- bart(y.train=PA[,"JT"], x.train=PA[,xnames], keeptrees=TRUE)

{png("topics/01_binomial_SDM/jtBART_summary.png", width=750, height=750)

summary(jtBART)

}
dev.off()

# Again, we want to simplify the model from the full set of 19 bioclim variables
# We can do this in stepwise BART model training. This will be slow!
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

# Estimate and visualize predictor partial effects; this will be slow!
jtBART.partials <- partial(jtBART.step, stepX, trace=FALSE, smooth=5) 

write_rds(jtBART.partials, file="output/models/jtBART.step.partials.rds") # save the results
jtBART.partials <- read_rds("output/models/jtBART.step.partials.rds") # reload later

jtBART.partials

# Build a multi-panel partials figure
partvals <- data.frame(predictor=rep(stepX, each=nrow(jtBART.partials[[1]]$data)), do.call("rbind", lapply(jtBART.partials, function(x) x$data))) %>% mutate(predictor=factor(predictor, stepX))


{png("topics/01_binomial_SDM/jtBART.step_partials.png", width=1000, height=1000)

ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, scale="free") + labs(y="Marginal Pr(occ)") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

}
dev.off()


# relative importance of predictors in the model:
varimp(jtBART.step)

# varimp() also generates a plot, if you want:
{png("topics/01_binomial_SDM/jtBART.step_varimp.png", width=500, height=500)
varimp(jtBART.step, plot=TRUE)
}
dev.off()


# And, finally, predict to the full range of the Mojave:

jtBART.step.pred <- predict(jtBART.step, envs[[stepX]], splitby=20)

jtBART.step.pred

jtBART.step.pred.df <- cbind(coordinates(jtBART.step.pred), as.data.frame(jtBART.step.pred)) %>% rename(prJT = layer, lon=x, lat=y)
glimpse(jtBART.step.pred.df)

{png("topics/01_binomial_SDM/jtBART.step_predicted.png", width=1000, height=1000)

ggplot() + 
  geom_tile(data=jtBART.step.pred.df, aes(x=lon, y=lat, fill=prJT)) +
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="white") + 
  scale_fill_distiller(type="seq", palette="YlGn", direction=1, name="Probability of occurrence") + labs(x="Longitude", y="Latitude") + 
  coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
  theme_bw(base_size=24) + theme(legend.position="bottom", legend.key.width=unit(0.05, "npc"))

}
dev.off()

#-------------------------------------------------------------------------
# Compare the BRM and BART results

BRMvBART <- jtBRM.pred - jtBART.step.pred

BRMvBART.df <- cbind(coordinates(BRMvBART), as.data.frame(BRMvBART)) %>% rename(prDiff = layer, lon = x, lat = y)
glimpse(BRMvBART.df)


{png("topics/01_binomial_SDM/BRMvBART_predictions.png", width=1000, height=1000)

ggplot() + 
  geom_tile(data=BRMvBART.df, aes(x=lon, y=lat, fill=prDiff)) +
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="black") + 
  scale_fill_distiller(type="seq", palette="PuOr", direction=1, name="BRM - BART\nprobability of occurrence") + labs(x="Longitude", y="Latitude") + 
  coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
  theme_bw(base_size=24) + theme(legend.position="bottom", legend.key.width=unit(0.05, "npc"))

}
dev.off()




