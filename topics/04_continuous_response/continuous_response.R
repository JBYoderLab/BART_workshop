# Modeling a continuous response with BARTs
# Jeremy B. Yoder, 21 Apr 2024

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

library("SoftBart")
# devtools::install_github("theodds/SoftBART")

#-----------------------------------------------------------
# Load and prepare data

# Burned area of forest fires and weather conditions within 30 min of ignition
#
#  P. Cortez and A. Morais. A Data Mining Approach to Predict Forest Fires using Meteorological Data. 
#  In J. Neves, M. F. Santos and J. Machado Eds., New Trends in Artificial Intelligence, 
#  Proceedings of the 13th EPIA 2007 - Portuguese Conference on Artificial Intelligence, December, 
#  Guimaraes, Portugal, pp. 512-523, 2007. APPIA, ISBN-13 978-989-95618-0-9. 
#  Available at: http://www.dsi.uminho.pt/~pcortez/fires.pdf
# 
#   1. X - x-axis spatial coordinate within the Montesinho park map: 1 to 9
#   2. Y - y-axis spatial coordinate within the Montesinho park map: 2 to 9
#   3. month - month of the year: "jan" to "dec" 
#   4. day - day of the week: "mon" to "sun"
#   5. FFMC - FFMC index from the FWI system: 18.7 to 96.20
#   6. DMC - DMC index from the FWI system: 1.1 to 291.3 
#   7. DC - DC index from the FWI system: 7.9 to 860.6 
#   8. ISI - ISI index from the FWI system: 0.0 to 56.10
#   9. temp - temperature in Celsius degrees: 2.2 to 33.30
#   10. RH - relative humidity in %: 15.0 to 100
#   11. wind - wind speed in km/h: 0.40 to 9.40 
#   12. rain - outside rain in mm/m2 : 0.0 to 6.4 
#   13. area - the burned area of the forest (in ha): 0.00 to 1090.84 

fires <- read.csv("data/forest_fires/forestfires.csv") %>% mutate(logArea=log10(area+1), id=row_number())

glimpse(fires)

#-------------------------------------------------------------------------
# BARTs via embarcadero

# train a dBARTs regression with all the possible predictors
xvars <- c("FFMC", "DMC", "DC", "ISI", "temp", "RH", "wind", "rain")

fireBART <- bart(y.train=fires[,"logArea"], x.train=fires[,xvars], keeptrees=TRUE)

# the summary function, unfortunately, assumes a binary classification model
summary(fireBART)

fireVarimp <- varimp.diag(y.data=fires[,"logArea"], x.data=fires[,xvars])

fireVarimp$data <- fireVarimp$data %>% mutate(trees = factor(trees, c(10,20,50,100,150,200)))

fireVarimp$labels$group <- "Trees"
fireVarimp$labels$colour <- "Trees"

write_rds(fireVarimp, file="output/models/fire_varimp.rds") # save the work
fireVarimp <- read_rds(file="output/models/fire_varimp.rds") # read it back in

# Write out the varimp diagnostic plot
{png(file="topics/04_continuous_response/fire_varimp_plot.png", width=750, height=500)

fireVarimp + 
theme_bw(base_size=18) +
theme(legend.position="inside", legend.position.inside=c(0.8, 0.7))  # okay nice

}
dev.off()

# fit a new BART model with the predictors suggested by varimp.diag()
varX <- c("temp", "rain", "wind")
fireBART2 <- bart(y.train=fires[,"logArea"], x.train=fires[,varX], keeptrees=TRUE)

# examine predictor partial effects (may be slow!)
fireBART2.p <- partial(fireBART2, varX, trace=FALSE, smooth=5)
fireBART2.p

partvals <- data.frame(predictor=rep(varX, each=nrow(fireBART2.p[[1]]$data)), do.call("rbind", lapply(fireBART2.p, function(x) x$data))) %>% mutate(predictor=factor(predictor, varX))


{png("topics/04_continuous_response/fireBART2_partials.png", width=1000, height=500)

ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="slategray2") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=1, scale="free") + labs(y="Marginal area") + theme_bw(base_size=24) + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

}
dev.off()


#-------------------------------------------------------------------------
# SoftBarts regressions

# train a SoftBarts regression with all the possible predictors
xvars <- c("FFMC", "DMC", "DC", "ISI", "temp", "RH", "wind", "rain")

train <- fires %>% slice_sample(prop=0.8) # use 80% of data for training
test <- fires %>% filter(!id%in%train$id) # use the other 20% for testing

# fit the model
fireDART <- softbart_regression(paste("logArea ~", paste(xvars, collapse="+")), data=train, test_data=test, k=2, opts=Opts(num_burn=2000, num_save=1000, num_thin=100)) # this may be slow

write_rds(fireDART, "output/models/fireDART.rds")
fireDART <- read_rds("output/models/fireDART.rds")

plot(fireDART$sigma_mu) # examine MCMC sampling

variable_selection <- data.frame(varimp=posterior_probs(fireDART)$varimp, post_prob=posterior_probs(fireDART)$post_probs,predictor=xvars)

# which predictors have posterior inclusion probability > 0.55? (Arbitrary threshold)
variable_selection
variable_selection %>% filter(post_prob > 0.55)

{png("topics/04_continuous_response/fireDART_variable_selection.png", width=1000, height=500)

ggplot() +
	geom_point(data=variable_selection, aes(x=predictor, y=post_prob, color=post_prob>0.55), size=4) +
	geom_hline(yintercept=0.55, linetype=2) + 
	labs(x = "Predictor", y = "Posterior inclusion prob") +
	theme_bw(base_size=24) + theme(axis.text.x=element_text(angle=75, hjust=1), legend.position="none")

}
dev.off()

# fit a new model with the top predictors
xtop <- variable_selection %>% filter(post_prob > 0.55) %>% .$predictor

fireDART2 <- softbart_regression(paste("logArea ~", paste(xtop, collapse="+")), data=train, test_data=test, k=2, opts=Opts(num_burn=2000, num_save=1000, num_thin=100))

plot(fireDART2$sigma_mu) # examine MCMC sampling


# visualize partial effects of individual predictors
grid_temp <- seq(from = min(fires$temp), to = max(fires$temp), length = 20)
pdf_temp <- partial_dependence_regression(fireDART2, train, "temp", grid_temp)

{png("topics/04_continuous_response/fireDART2_partial_temp.png", width=500, height=500)
ggplot(pdf_temp$pred_df, aes(x = temp, y = mu)) +
	geom_line(stat = "summary", fun = mean) +
	geom_ribbon(stat = "summary", alpha = 0.3, fun.min = function(x) quantile(x, 0.025), fun.max = function(x) quantile(x, 0.975)) + xlab("temp") + ylab(expression(log[10](area+1))) +
	theme_bw(base_size=18)
}
dev.off()


grid_RH <- seq(from = min(fires$RH), to = max(fires$RH), length = 20)
pdf_RH <- partial_dependence_regression(fireDART2, train, "RH", grid_RH)

{png("topics/04_continuous_response/fireDART2_partial_RH.png", width=500, height=500)
ggplot(pdf_RH$pred_df, aes(x = RH, y = mu)) +
	geom_line(stat = "summary", fun = mean) +
	geom_ribbon(stat = "summary", alpha = 0.3, fun.min = function(x) quantile(x, 0.025), fun.max = function(x) quantile(x, 0.975)) + xlab("RH") + ylab(expression(log[10](area+1))) +
	theme_bw(base_size=18)
}
dev.off()

