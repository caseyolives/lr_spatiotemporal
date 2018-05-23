######################################################
##### ESTIMATE CV models
#######################################################

rm(list=ls())
library(SpatioTemporal)

if(Sys.info()['sysname'] == "Windows"){
	root <- "P:/"
	group = "HOME"
	ii = 17
} else {
	root <- "~/"
	args <- commandArgs()
	group <- args[3]
	ii <- as.numeric(args[4]) # which model
	print(ii)
}

######################################################################
#### Load Modified Functions
source(paste(root,"Research/Low_Rank_ST/Code/modified_ST2_functions.r", sep=""))

#### Load additional Functions
source(paste(root,"Research/Low_Rank_ST/Code/additional_low_rank_functions.r", sep=""))

#### Load Model Objects
load(file=paste(root, "Research/Low_Rank_ST/Output/MESA_Models_", ii, ".RData", sep=""))

#### Load CV groups
load(file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CV_GROUPS.RData", sep=""))

#### Load Estimated Model Objects to use as starting value
if(ii %in% c(29, 30)){
	x.init.save <- x.init
	load(file=paste(root, "Research/Low_Rank_ST/Output/ESTIMATED_LR_MODELS_", ii, ".RData", sep=""))
	x.init <- cbind(est.mesa.models$res.best$par.cov$par, x.init)
} else if(ii %in% c(32)){
	x.init.save <- x.init
	load(file=paste(root, "Research/Low_Rank_ST/Output/ESTIMATED_LR_MODELS_", ii, ".RData", sep=""))
	x.init <- cbind(est.mesa.models$res.best$par.cov$par, x.init[, -3])
}


#### Fit the model
cv.est.mesa.models <- estimateCV.STmodel_opt( object=models, 
													x=x.init,
													Ind.cv=Ind.cv, 
													silent=FALSE, 
													x.fixed=x.fixed, 
													type = "p")

#### Output Results
save(cv.est.mesa.models, 
		file=paste(root, "Research/Low_Rank_ST/Output/", group, "_CVESTIMATED_LR_MODELS_", ii, ".RData", sep=""))