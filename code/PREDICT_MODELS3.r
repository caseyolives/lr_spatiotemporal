######################################################
##### Predicts using various Model objects using new functions
#######################################################

rm(list=ls())
library(SpatioTemporal)

if(Sys.info()['sysname'] == "Windows"){
	root <- "P:/"
	ii=3
} else if (Sys.info()['sysname'] == "Darwin") {
	root <- "/Volumes/hercules.biostat.washington.edu/"
	ii=17
} else {
	root <- "~/"
	args <- commandArgs()
	ii <- as.numeric(args[3]) # which model
	print(ii)
}

######################################################################
#### Load Modified Functions
source(paste(root,"Research/Low_Rank_ST/Code/modified_ST2_functions.r", sep=""))

#### Load additional Functions
source(paste(root,"Research/Low_Rank_ST/Code/additional_low_rank_functions.r", sep=""))

#### Load Model Objects
load(file=paste(root, "Research/Low_Rank_ST/Output/MESA_Models_", ii, ".RData", sep=""))

#### Load Estimated Model Object
load(file=paste(root, "Research/Low_Rank_ST/Output/ESTIMATED_LR_MODELS_", ii, ".RData", sep=""))
x.pred <- est.mesa.models$res.best$par.cov$par

pred.mesa.models <- predict.STmodel_opt(	object=models, 
													STdata=models,
													x=x.pred, 
													pred.var=FALSE, 
													type = "p")

#### Output Results
save(pred.mesa.models, 
		file=paste(root, "Research/Low_Rank_ST/Output/PREDICTED_LR_MODELS_", ii, ".RData", sep=""))