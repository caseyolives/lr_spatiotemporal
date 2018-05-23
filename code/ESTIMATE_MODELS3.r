######################################################
##### Fit various Model objects using new functions
#######################################################

rm(list=ls())
library(SpatioTemporal)

if(Sys.info()['sysname'] == "Windows"){
	root <- "P:/"
	ii=3
} else if (Sys.info()['sysname'] == "Darwin") {
	root <- "/Volumes/hercules.biostat.washington.edu/"
	ii=3
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
load(file=paste(root, "Research/Low_Rank_ST/Output/MESA_Models_", ii,".RData", sep=""))

use.other.mod.init <- TRUE
if(use.other.mod.init){
	mod <- 29
	load(file=paste(root, "Research/Low_Rank_ST/Output/ESTIMATED_LR_MODELS_", mod, ".RData", sep=""))	
	if(nrow(x.init) == nrow(est.mesa.models$res.best$par.cov)){
		x.init <- est.mesa.models$res.best$par.cov$par
	}
}

#### Fit the model
est.sys.time <- system.time(est.mesa.models <- estimate.STmodel_opt(object=models, 
													x=x.init, 
													x.fixed=x.fixed, 
													type = "p"))

#### Output Results
save(est.sys.time, est.mesa.models, 
		file=paste(root, "Research/Low_Rank_ST/Output/ESTIMATED_LR_MODELS_", ii, ".RData", sep=""))