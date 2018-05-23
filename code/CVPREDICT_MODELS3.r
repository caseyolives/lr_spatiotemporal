######################################################
##### Do AQS/FIXED CV Predictions
#######################################################

rm(list=ls())
library(SpatioTemporal)

if(Sys.info()['sysname'] == "Windows"){
	root  <- "P:/"
	group <- "FIXED"
	ii = 1
	cv=TRUE
} else if (Sys.info()['sysname'] == "Darwin") {
	root <- "/Volumes/hercules.biostat.washington.edu/"
	group="FIXED"
	ii=17
	cv=TRUE
} else {
	root <- "~/"
	args <- commandArgs()
	group <- args[3]
	ii <- as.numeric(args[4])
	cv <- args[5]
}

######################################################################
#### Load Modified Functions
source(paste(root,"Research/Low_Rank_ST/Code/modified_ST2_functions.r", sep=""))

#### Load additional Functions
source(paste(root,"Research/Low_Rank_ST/Code/additional_low_rank_functions.r", sep=""))

#### Load Model Objects
load(file=paste(root, "Research/Low_Rank_ST/Output/MESA_Models_", ii, ".RData", sep=""))

#### Load Estimated Model Objects
if(cv){ # do full CV or base on model fits?
	load(file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CVESTIMATED_LR_MODELS_", ii, ".RData", sep=""))
	x.pred <- cv.est.mesa.models$par.cov
} else {
	load(file=paste(root, "Research/Low_Rank_ST/Output/ESTIMATED_LR_MODELS_", ii, ".RData", sep=""))
	x.pred <- est.mesa.models$res.best$par.cov$par	
}

#### Load CV groups
load(file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CV_GROUPS.RData", sep=""))

# if(ii %in% 32:61){
cv.pred.est.mesa.models <- predictCV.STmodel_opt(object=models, 
													x=x.pred, 
													pred.var=FALSE,
													Ind.cv=Ind.cv, 
													silent=FALSE, 
													type="p")

#### Save Results
save(cv.pred.est.mesa.models, Ind.cv, file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CVPREDICTED_LR_MODELS_", ii, ".RData", sep=""))
