Working with BARTs in R
=======================


README updated 24 April 2024.


Project description
-------------------

Materials for a workshop demonstrating the uses of Bayesian additive regression tree methods implemented in [`dbarts`](https://cran.r-project.org/web/packages/dbarts) with utilities from [`embarcadero`](https://www.github.com/cjcarlson/embarcadero), as well as [SoftBart](https://www.github.com/theodds/SoftBART). Activities will focus on applications in ecology and biodiversity research.


Contents
--------

- `topics` --- code and slides for worked examples, organized by topic (see below)
- `output` --- materials output by scripts in the worked examples. Saved models amount to more than 1Gb, and are not in the version-controlled repo or posted for download.
- `data` --- data sets needed for the worked examples are more than 100Mb in total, and not posted here, but you can get them from [this Google Drive folder](https://drive.google.com/drive/folders/1JW5iGgF8S2VoYgzjJ5MMXlSP7kWrnto5?usp=share_link).



Module contents
---------------

The code provided should let you dip into any example and take what you want from it. However, this material was developed for a multi-day workshop, and that experience is best reproduced if you work through it in the same order as presented. In some cases, results of computationally intensive (slow) analyses are saved in an earlier module for use in a later module.

Here are the modules and brief descriptions of their contents.

1. Species distribution modeling with BARTs (binary response)
	- Setting up date for an SDM of Joshua trees with presence and pseudo-absence records
	- Comparison of boosted regression tree modeling and BARTs

2. Random intercept BART models	--- Example of the Joshua tree SDM with a random intercept effect to account for uneven sampling effort

3. Predictor selection
	- The `varimp.diag()` predictor selection approach implemented in `embarcadero`
	- Dirichlet predictor selection and posterior inclusion probability assessment from `SoftBart`

4. Modeling a continuous response --- Example of modeling a continuous response (forest fire area) 
	- Continuous-response BART modeling in `dBarts`/`embarcadero`
	- Continuous-response BART regression in `SoftBart`

5. Predicting annual flowering by Joshua trees
	- Modeling Joshua tree mast-flowering as a binary response to annual weather variation
	- Predictor selection
	- Checking for confounding from temporal heterogeneity in sampling with RI-BART
	- Illustrating predictor partial effects and spatial partial effects (spartials)
	- Prediction of flowering in years without direct observation
	
	
