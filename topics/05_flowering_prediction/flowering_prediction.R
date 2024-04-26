# Using BARTs to model Joshua tree mast-flowering
# Jeremy B. Yoder, 20 Apr 2024

rm(list=ls())  # Clears the environment and load key packages

# This script assumes your working directory is the project directory; you may want to change this
# setwd("~/Documents/Active_projects/BART_workshop")

library("tidyverse") 
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

library("embarcadero")
# devtools::install_github("cjcarlson/embarcadero")

#-----------------------------------------------------------
# initial file loading

# flowering/not flowering records for a ~4km raster grid
# biologically-informed candidate predictors, subspecies id'd
flow <- read.csv("data/JT_flowering_obs_climate.csv") %>% filter(!(year==2018.5 & flr==TRUE), year>=2008, year<2023) %>% mutate(year=floor(year)) # drop the late-flowering anomaly

glimpse(flow)

flrsum <- data.frame(table(flow$year, flow$flr)) %>% pivot_wider(names_from=Var2, values_from=Freq) %>% mutate(inc = `TRUE`/(`FALSE`+`TRUE`), Var1=as.numeric(as.character(Var1))) # many more observations in recent years
cor.test(~inc+Var1, data=flrsum, method="sp") # but incidence doesn't have a trend

#-------------------------------------------------------------------------
# Predictor selection

# predictors
xnames <- c("pptY0", "pptY1", "pptY2", "pptY0Y1", "pptY1Y2", "tmaxY0", "tmaxY0Y1", "tminY0", "tminY0Y1", "vpdmaxY0", "vpdmaxY0Y1", "vpdminY0", "vpdminY0Y1") # weather measures for up to two years prior to observation

# variable importance
jotr.varimp <- varimp.diag(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames])
# favors pptY1Y2, pptY0Y1, vpdmaxY0, vpdminY0Y1, tminY0, tmaxY0Y1

write_rds(jotr.varimp, file="output/models/bart.varimp.Jotr.rds")
jotr.varimp <- read_rds("output/models/bart.varimp.Jotr.rds")

# publication-ready figure assembly ...
glimpse(jotr.varimp$data)

jotr.varimp$data <- jotr.varimp$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

levels(jotr.varimp$data$names) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

jotr.varimp$labels$group <- "Trees"
jotr.varimp$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

# the varim.diag() plot
{png(file="topics/05_flowering_prediction/varimp_Jotr_flowering.png", width=750, height=500)

jotr.varimp + scale_x_discrete(label=label_parse) + theme_bw(base_size=18) + theme(legend.position="inside", legend.position.inside=c(0.8, 0.7), axis.text=element_text(size=13), axis.text.x=element_text(size=13, angle=75, hjust=1), legend.text=element_text(size=12), legend.title=element_text(size=13))  # okay nice

}
dev.off()

#-------------------------------------------------------------------------
# Model fitting

# Fitting a model with vars indicated by varimp:
jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

jotr.mod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,jotr.preds], keeptrees=TRUE)

summary(jotr.mod)

invisible(jotr.mod$fit$state)
write_rds(jotr.mod, file="output/models/bart.model.Jotr.rds")
jotr.mod <- read_rds("output/models/bart.model.Jotr.rds")


# Fit a RI-BART model to examine the effect of changing sampling desity over years
jotr.RImod <- rbart_vi(
	as.formula(paste(paste('flr', paste(jotr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = flow,
	group.by = flow[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(jotr.RImod)

{png("topics/05_flowering_prediction/Jotr_flower_RImod_RI_est.png", width=750, height=500)
plot.ri(jotr.RImod, temporal=TRUE) + labs(title="Estimated random intercept effects of observation year", x="Observation year") + theme_bw(base_size=18) + theme(axis.text.x=element_text(angle=0, size=16))

}
dev.off()


#-------------------------------------------------------------------------
# Partial effects and spatial partial effects

p <- partial(jotr.mod, jotr.preds, trace=FALSE, smooth=5) # visualize partials
varimp(jotr.mod)

write_rds(p, file="output/models/bart.model.Jotr.partials.rds")
p <- read_rds("output/models/bart.model.Jotr.partials.rds")

# generate a multi-panel figure of partial effects
partvals <- data.frame(predictor=rep(jotr.preds, each=nrow(p[[1]]$data)), do.call("rbind", lapply(p, function(x) x$data))) %>% mutate(predictor=factor(predictor, jotr.preds))


{png("topics/05_flowering_prediction/Jotr_flowering_partials.png", width=750, height=1000)

ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, scale="free") + labs(y="Marginal Pr(Flowers)") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

}
dev.off()


# Spatial partial effects
# because the predictors vary by year, we need to do this for each year of observation
yr <- 2020

spYr <- spartial(jotr.mod, brick(paste("data/JT_PRISM_predictors/PRISM_derived_predictors_", yr, ".grd", sep="")), x.vars=jotr.preds)

spart.df <- cbind(coordinates(spYr), as.data.frame(spYr)) %>% filter(!is.na(pptY1Y2)) |> rename(lon=x, lat=y) |> pivot_longer(all_of(jotr.preds), names_to="predictor", values_to="prFL") |> mutate(predictor=factor(predictor,jotr.preds))

levels(spart.df$predictor) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

{png(paste("topics/05_flowering_prediction/Jotr_flower_spartials_", yr, ".png", sep=""), width=750, height=1000)

ggplot() + 
	geom_tile(data=spart.df, aes(x=lon, y=lat, fill=prFL)) + 
	geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="white") + 
	facet_wrap("predictor", nrow=3, labeller="label_parsed") + 
	scale_fill_gradient(low="#ffffcc", high="#253494", name="Marginal Pr(Flowers)", breaks=c(0.5,0.55,0.6,0.65), limits=c(0.5,0.675), labels=c("", 0.55, 0.6, 0.65)) + 
	labs(title=paste("Spatial partial effects for", yr)) + 
	coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
	theme_bw(base_size=24) + theme(axis.title=element_blank(), axis.text=element_blank(), legend.position="bottom", panel.background=element_rect(fill="white"), panel.grid=element_blank(), legend.key.width=unit(0.05, "npc"))

}
dev.off()


#-------------------------------------------------------------------------
# Predict flowering in years not observed

yr <- 1995 # try 1994 and 1995

preds <- brick(paste("data/JT_PRISM_predictors/PRISM_derived_predictors_", yr, ".grd", sep=""))

# prediction with the RI predictor (year) removed
pred.yr <- predict(jotr.mod, preds[[jotr.preds]], splitby=20)

writeRaster(pred.yr, paste("topics/05_flowering_prediction/BART_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

pred.df <- cbind(coordinates(pred.yr), as.data.frame(pred.yr)) %>% rename(prFL = layer, lon = x, lat = y)
glimpse(pred.df)

{png(paste("topics/05_flowering_prediction/Jotr_predicted_flowering_", yr, ".png", sep=""), width=1000, height=1000)

ggplot() + 
  geom_tile(data=pred.df, aes(x=lon, y=lat, fill=prFL)) +
  geom_sf(data=ne_states(country = "United States of America", returnclass = "sf"), fill=NA, color="black") + 
  scale_fill_distiller(type="seq", palette="YlGn", direction=1, name="Predicted Pr(Flowers)") + 
  labs(x="Longitude", y="Latitude") + 
  coord_sf(xlim = c(-119, -112.5), ylim = c(33.5, 38), expand = TRUE) +
  theme_bw(base_size=24) + theme(legend.position="bottom", legend.key.width=unit(0.05, "npc"))

}
dev.off()

plot(pred.yr>0.29)








