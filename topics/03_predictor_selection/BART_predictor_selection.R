# Predictor selection and RI models with BARTs
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

library("embarcadero")
# devtools::install_github("cjcarlson/embarcadero")

library("SoftBart")
# devtools::install_github("theodds/SoftBART")

#-----------------------------------------------------------
# Load and prepare data

# Joshua tree occurrence records, and pseudoabsences, with environmental data
# (This was generated in the code for topic 01, SDMs in BART)
PA <- read.csv("data/JT_presence-pseudoabsence-envs.txt") %>% filter(!is.na(MAT))
glimpse(PA)

# Visualize the points on a map of the US Southwest
# As you can see, there's quite a bit of heterogeneity in sampling density.
{png("topics/02_predictor_selection/JT_presence_pseudoabsence_map.png", height=750, width=750)
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

# Fit a simple BART model with all the predictors
jtBART <- bart(y.train=PA[,"JT"], x.train=PA[,xnames], keeptrees=TRUE)

{png("topics/01_binomial_SDM/jtBART_summary.png", width=750, height=750)
summary(jtBART)
}
dev.off()

# Previously, we simplified the model from the full set of 19 bioclim variables using stepwise training.
# We can re-run that, but it'll be slow; if it's already saved, just reload below.
jtBART.step <- bart.step(y.data=as.numeric(PA[,"pres"]), x.data=PA[,xnames], full=FALSE, quiet=TRUE)

# It may be useful to save the model object for later work. 
# First, we need to "touch" the model state to make sure the trees are saved with the rest of the model.
# This is a `dbarts` feature that is supposed to save storage.
invisible(jtBART.step$fit$state)
# Then write out the model to an Rdata file
write_rds(jtBART.step, file="output/models/jtBART.step.rds") 
jtBART.step <- read_rds(file="output/models/jtBART.step.rds")

{png(file="topics/02_predictor_selection/jBART.step_summary.png", width=750, height=750)
summary(jtBART.step)
}
dev.off()

stepX <- attr(jtBART.step$fit$data@x, "term.labels")
stepX # these are the predictors kept in the stepwise model


#-------------------------------------------------------------------------
# Formal predictor selection with BARTs

# Alternatively, use varimp.diag() for predictor selection
# This may be slow!
jtVarimp <- varimp.diag(y.data=as.numeric(PA[,"pres"]), x.data=PA[,xnames])

# a useful fix for a minor issue
jtVarimp$data <- jtVarimp$data %>% mutate(trees = factor(trees, c(10,20,50,100,150,200)))

levels(jtVarimp$data$names) <- xnames

jtVarimp$labels$group <- "Trees"
jtVarimp$labels$colour <- "Trees"

write_rds(jtVarimp, file="output/models/jt_varimp.rds") # save the work
jtVarimp <- read_rds(file="output/models/jt_varimp.rds") # read it back in

# Write out the varimp diagnostic plot
{png(file="topics/02_predictor_selection/jt_varimp_plot.png", width=750, height=500)

jtVarimp + 
theme_bw(base_size=18) +
theme(legend.position="inside", legend.position.inside=c(0.8, 0.7), axis.text.x=element_text(angle=90))  # okay nice

}
dev.off()

# This lets us identify a subset of predictors that are informative in simpler models
varimpX <- c("MAT", "MDR", "ITH", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT")


# Train a new model with these predictors
jtBART.var <- bart(y.train=PA[,"JT"], x.train=PA[,varimpX], keeptrees=TRUE)

{png("topics/02_predictor_selection/jtBART.var_summary.png", width=750, height=750)
summary(jtBART.var)
}
dev.off()

# save the model
invisible(jtBART.var$fit$state)
write_rds(jtBART.var, "output/models/jtBART.var.rds")
jtBART.var <- read_rds("output/models/jtBART.var.rds")


# OR --- which predictors does varimp.diag() suggest that aren't in the stepwise model?
setdiff(varimpX, stepX) 
# And vice versa
setdiff(stepX, varimpX) 
# Which predictors are shared?
intersect(varimpX, stepX) 

# Let's try fitting a model with the predictors selected by both methods
intX <- intersect(varimpX, stepX) 

jtBART.int <- bart(y.train=PA[,"JT"], x.train=PA[,intX], keeptrees=TRUE)

{png("topics/02_predictor_selection/jtBART.int_summary.png", width=750, height=750)

summary(jtBART.int)

}
dev.off()

invisible(jtBART.int$fit$state)
write_rds(jtBART.int, "output/models/jtBART.int.rds")
jtBART.int <- read_rds("output/models/jtBART.int.rds")


#-------------------------------------------------------------------------
# Compare these models

# re-load the partials we estimated for the stepwise model
jtBART.partials <- read_rds("output/models/jtBART.step.partials.rds") # reload later

jtBART.partials

# Build a multi-panel partials figure
step.partvals <- data.frame(predictor=rep(stepX, each=nrow(jtBART.partials[[1]]$data)), do.call("rbind", lapply(jtBART.partials, function(x) x$data))) %>% mutate(predictor=factor(predictor, stepX))

{png("topics/02_predictor_selection/jtBART.step_partials.png", width=1000, height=1000)

ggplot(step.partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, scale="free") + labs(y="Marginal Pr(occurs)") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.01,"npc"))

}
dev.off()

# partials for the varimp.diag predictor set
jtBART.var.partials <- partial(jtBART.var, varimpX, trace=FALSE, smooth=5) 
write_rds(jtBART.var.partials, file="output/models/jtBART.var.partials.rds") # save the results
jtBART.var.partials <- read_rds("output/models/jtBART.var.partials.rds") # reload later

jtBART.var.partials

# Build a multi-panel partials figure
var.partvals <- data.frame(predictor=rep(varimpX, each=nrow(jtBART.var.partials[[1]]$data)), do.call("rbind", lapply(jtBART.var.partials, function(x) x$data))) %>% mutate(predictor=factor(predictor, varimpX))


{png("topics/02_predictor_selection/jtBART.var_partials.png", width=1000, height=1000)

ggplot(var.partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, scale="free") + labs(y="Marginal Pr(occurs)") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.01,"npc"))

}
dev.off()


# partials for the intersection predictor set
jtBART.int.partials <- partial(jtBART.int, intX, trace=FALSE, smooth=5) 
write_rds(jtBART.int.partials, file="output/models/jtBART.int.partials.rds") # save the results
jtBART.int.partials <- read_rds("output/models/jtBART.int.partials.rds") # reload later

jtBART.int.partials

# Build a multi-panel partials figure
int.partvals <- data.frame(predictor=rep(intX, each=nrow(jtBART.int.partials[[1]]$data)), do.call("rbind", lapply(jtBART.int.partials, function(x) x$data))) %>% mutate(predictor=factor(predictor, intX))


{png("topics/02_predictor_selection/jtBART.int_partials.png", width=1000, height=1000)

ggplot(int.partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, scale="free") + labs(y="Marginal Pr(occurs)") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.01,"npc"))

}
dev.off()


# OR alternatively, use the Dirichlet inclusion prior implemented in `SoftBart`
# First, re-fit the full model with softbart:
train <- PA %>% slice_sample(prop=0.8) # use 80% of data for training
test <- PA %>% filter(!id%in%train$id) # use the other 20% for testing

# fit the model
jtDART <- softbart_probit(paste("factor(JT) ~", paste(xnames, collapse="+")), data=train, test_data=test, k=2, opts=Opts(num_burn=5000, num_save=5000))

write_rds(jtDART, "output/models/jtDART.rds")
jtDART <- read_rds("output/models/jtDART.rds")

plot(jtDART$sigma_mu)

variable_selection <- data.frame(varimp=posterior_probs(jtDART)$varimp, post_prob=posterior_probs(jtDART)$post_probs, median_probability_model=posterior_probs(jtDART)$median_probability_model,predictor=xnames)

plot(posterior_probs(jtDART)$post_probs)
print(posterior_probs(jtDART)$median_probability_model)

plot(variable_selection$post_probs, col = ifelse(1:250 < 6, "#386CB0", "#7FC97F"), pch = 20,
     xlab = "Predictor", ylab = "PIP", main = "Variable Selection")

rmse(jtDART$y_hat_test_mean, PA[,xnames])


























# Estimate and visualize predictor partial effects; this will be slow!
jtBART.partials <- partial(jtBART.step, stepX, trace=FALSE, smooth=5) 

write_rds(jtBART.partials, file="topics/01_binomial_SDM/jtBART.step.partials.rds") # save the results
jtBART.partials <- read_rds("topics/01_binomial_SDM/jtBART.step.partials.rds") # reload later

jtBART.partials

# Build a multi-panel partials figure
partvals <- data.frame(predictor=rep(stepX, each=nrow(jtBART.partials[[1]]$data)), do.call("rbind", lapply(jtBART.partials, function(x) x$data))) %>% mutate(predictor=factor(predictor, stepX))


{png("topics/01_binomial_SDM/jtBART.step_partials.png", width=1000, height=1000)

ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, scale="free") + labs(y="Marginal Pr(Flowers)") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

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



#-------------------------------------------------------------------------
# Fit a species distribution model with RI-BARTs

jtBART.RI <- rbart_vi(as.formula(paste(paste('pres', paste(preds, collapse=' + '), sep = ' ~ '), 'STATE', sep=' - ')),
	data = flow,
	group.by = flow[,'STATE'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(jtBART.RI)

# Save again
invisible(jtBART.RI$fit$state)
write_rds(jtBART.RI, file="topics/01_binomial_SDM/jtBART.RI.Rdata") 


# visualize RI estimates --- is there a temporal trend? If so, keep RI, otherwise stick with non-RI model

{cairo_pdf("01_binomial_SDM/BART_RI_estimates.png", width=1000, height=750)

plot.ri(jtBART.RI, temporal=FALSE) + labs(title="Random intercept effects of state", x="State") + theme_bw(base_size=24) + theme(axis.text.x=element_text(angle=0, size=12))

}
dev.off()




# Now we're ready to make a prediction. The `predict()` function can take your
# saved BART sdm and the cropped RasterStack of environmental values, then 
# generate the predicted probability of your species' occurrence at each cell 
# in the RasterStack given the BioClim values assigned to that cell.
pred <- predict(sdm, BClimCrop) 

# The output is another raster layer, with probability of your species' 
# occurrence in each cell:
plot(pred)

# How does that look?

# [WAIT HERE to check in before moving on]

# One possible issue with the SMD you've just made is that it uses all your 
# predictors without evaluating how well they actually predict your species' 
# presence or absence. We can also use utilities in `embarcadero` to perform
# predictor selection, via a sort of model comparison for BARTs.

# The `varimp.diag()` function fits a BUNCH of BART models using different 
# levels of complexity (more or fewer regression trees), with or without each
# predictor variable. Predictors that are more often important in simpler models
# can be understood to have more predictive value.

predSel <- varimp.diag(y.data=PA[,"pres"], x.data=PA[,xnames]) 
# [This may take a while. WAIT HERE to check in before moving on]

# The output of `varimp.diag()` should let you zero in on a subset of predictors
# to use in a new SDM. Let's create a new vector with just those predictors
# named, so we can fit a new BART model with just those predictors.
xnames2 <- c("PDM", "TS", "PWaQ", "PS", "TAR", "MCMT", "ITH", "MDR")

# Now, fit our new BART model
sdm2 <- bart(y.train=PA[,"pres"], x.train=PA[,xnames2], keeptrees=TRUE)

summary(sdm2) # How does this compare to the original model?

# [WAIT HERE to check in before moving on]

# We can use the `varimp()` function to illustrate how each predictor 
# contributes to the model's predictions. 
varimp(sdm2, plots=TRUE)

# And we can use the `partial()` function to summarize the effect of each
# predictor on the probability that the species is present. This may be quite
# slow again, because it's summarizing across all the trees in the BART model.
p <- partial(sdm2, xnames2, trace=FALSE, smooth=5) 

# Calling the object where you've stored the output from `partial()` returns
# a series of plots, one for each predictor. Page through these and you can see
# how each predictor contributes to the model's predicted probability of the 
# species' presence. Do these make biological sense, given the predictors?
p

# In `embarcadero` there's also an option for visualizing "spartials", the 
# spatial partial effects of each predictor. This is also quite slow, but can be
# illuminating once it's finished:

sdmSpart <- spartial(sdm2, BClimCrop, x.vars=xnames2)

# Plotting these spartials shows how each predictor contributes to the model's
# predicted probability of species presence IN SPACE
plot(sdmSpart)

# And then, finally, we can actually make a prediction of species' presence
# from the fitted model, as we did with the first model:
pred2 <- predict(sdm2, BClimCrop) 

# And plot it
plot(pred2)

# How does that compare to the earlier model with all 19 BioClim predictors?
plot(pred2 - pred) 
# This shows places where the new model predicts higher or lower probability of 
# presence than the original 19-predictor model.

