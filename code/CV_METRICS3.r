######################################################
##### Calculate CV Metrics
#######################################################

rm(list=ls())
library(SpatioTemporal)

if(Sys.info()['sysname'] == "Windows"){
	root  <- "P:/"
	group <- "FIXED"
	ii <- 1
} else {
	root <- "~/"
	args <- commandArgs()
	group <- args[3]
	ii <- as.numeric(args[4])
}

######################################################################
#### Load Modified Functions
source(paste(root,"Research/Low_Rank_ST/Code/modified_ST2_functions.r", sep=""))

#### Load Additional Functions
source(paste(root,"Research/Low_Rank_ST/Code/additional_low_rank_functions.r", sep=""))

#### Load Model Objects
load(file=paste(root, "Research/Low_Rank_ST/Output/MESA_Models_", ii, ".RData", sep=""))

#### Load in CV Predictions
load(file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CVPREDICTED_LR_MODELS_", ii,".RData", sep=""))

#### Load in prediction variances
load(file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CVPREDICT_VARIANCE_LR_MODELS_", ii, ".RData", sep=""))

#####################################################################
#### Calc CV Metrics

try(

	cv.metrics <- calc.CV3(object=cv.pred.est.mesa.models, lta=lta, option=tolower(group))

	)

#### Save Results
save(cv.metrics, file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CV_METRICS_LR_MODELS_", ii, ".RData", sep=""))