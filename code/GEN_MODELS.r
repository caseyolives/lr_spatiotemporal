######################################################
##### GENERATE Low-Rank and Full-Rank Model Objects
##### To be evaluated by LL functions and eventually fit
#######################################################

rm(list=ls())
library(SpatioTemporal)

if(Sys.info()['sysname'] == "Windows"){
	root <- "P:/"
} else if (Sys.info()['sysname'] == "Darwin") {
	root <- "/Volumes/hercules.biostat.washington.edu/"
} else {
	root <- "~/"
}

######################################################################
#### Load Modified Functions
source(paste(root,"Research/Low_Rank_ST/Code/modified_ST2_functions.r", sep=""))

#### Load Additional Functions
source(paste(root,"Research/Low_Rank_ST/Code/additional_low_rank_functions.r", sep=""))

#### Load MESA Test Data
# data(mesa.data)

### Read in CA Objects from Assaf - THIS IS THE REV4 OBJECT
load(paste(root, "MESA Data/REV4_NOX_MODELS/Nox4_finalmodel_LA.RData", sep=""))
nox.data <- laxST
nox.model <- laxSTmod
#load(paste(root, "MESA Data/Assaf_Data_Files/noCApredmodels.RData", sep=""))
#nox.model <- nox1CAmodel1

################################################################
################################################################
### Prepare LA ONLY NOx Data for LR Modeling

### Grab Covars
covars <- as.data.frame(nox.data$covars)

### Identify sites in LA
AQS_LA_ID <- covars$ID[(substr(covars$ID, 3, 5) %in% c("037", "059"))]
COMCO_LA_ID <- unique(nox.data$obs$ID[nox.data$obs$date %in% as.Date(c("2006-07-05", "2006-10-25", "2007-01-31"))])
COMCO_LA_ID <- covars$ID[covars$ID %in% COMCO_LA_ID & covars$type=="COMCO"]
HOME_LA_ID <- covars$ID[substr(covars$ID, 2, 2) != "R" & covars$type=="HOME"]
FIXED_LA_ID <- covars$ID[substr(covars$ID, 2, 2) != "R" & covars$type=="FIXED"]
LA_ID <- c(AQS_LA_ID, COMCO_LA_ID, HOME_LA_ID, FIXED_LA_ID)

### Now Prepare LA DATA ONLY ###
nox.data.LA <- createSTdata(obs = nox.data$obs[nox.data$obs$ID %in% LA_ID, c("obs", "date", "ID")], covars=covars[covars$ID %in% LA_ID,])

### Estimate SEOFs using FIXED and AQS data only
sub <- nox.data.LA$covars$ID[nox.data.LA$covars$type%in%c("FIXED", "AQS")]

# Re-estimate trend for LA	
F_LA <- calcSmoothTrends(nox.data.LA, n.basis=2, subset=sub, extra.dates=nox.data.LA$obs$date)	
nox.data.LA$trend <- F_LA$trend
rownames(nox.data.LA$trend) <- nox.data.LA$trend$date

mesa.data <- nox.data.LA

#### Generate Model Objects with nuggets in beta fields and nu field
#### Set RE = FALSE in nu field

#### Optimize over range 
#### M1 = beta=FR
#### M2 = beta=LR with K=100
#### M3 = beta=LR with K=50
#### M4 = beta=LR with K=25
#### M5 = beta=LR with K=10

#### Fixed Range at max dist
#### M6 = beta=FR
#### M7 = beta=LR with K=100
#### M8 = beta=LR with K=50
#### M9 = beta=LR with K=25
#### M10 = beta=LR with K=10

#### Fixed Range at 1/2 dist
#### M11 = beta=FR
#### M12 = beta=LR with K=100
#### M13 = beta=LR with K=50
#### M14 = beta=LR with K=25
#### M15 = beta=LR with K=10

#### Fixed Range at 1/4 dist
#### M16 = beta=FR
#### M17 = beta=LR with K=100
#### M18 = beta=LR with K=50
#### M19 = beta=LR with K=25
#### M20 = beta=LR with K=10

#### Fixed Range at 1/8 dist
#### M21 = beta=FR
#### M22 = beta=LR with K=100
#### M23 = beta=LR with K=50
#### M24 = beta=LR with K=25
#### M25 = beta=LR with K=10

#### TPRS 
#### M26 = beta=FR
#### M27 = beta=LR with K=100
#### M28 = beta=LR with K=50
#### M29 = beta=LR with K=25
#### M30 = beta=LR with K=10

#### M31 = beta = LR with K=0 (ie iid beta fields)

#### Additional Models as of 12/19/2012
#### Fit all of these models without the nugget (except Model 31)
### M32 - M61 are the same as M1 - M30 without nugget in beta field and fit using the optimized LL

### Models without long range covariates
### That is, remove city.hall, log10.caline and m.to.coast from LUR
### Fully optimized range with nugget in beta field
#### M62 = beta=FR
#### M63 = beta=LR with K=100
#### M64 = beta=LR with K=50
#### M65 = beta=LR with K=25
#### M66 = beta=LR with K=10

### Models with knots on a grid
### Fully optimized range with nugget in beta field
#### M67 = beta=FR
#### M68 = beta=LR with K=100
#### M69 = beta=LR with K=50
#### M70 = beta=LR with K=25
#### M71 = beta=LR with K=10

#### Additional LRK models that include lambert_x and lambert_y as covariates
#### Include Nugget in these models
#### M72 = beta=FR
#### M73 = beta=LR with K=100
#### M74 = beta=LR with K=50
#### M75 = beta=LR with K=25



##################################
#### Create Input Parameters #####
##################################

## 12/19/2012
library(sp)
library(grDevices)
hull <- mesa.data$covars[chull(mesa.data$covars[, c("x", "y")]), c("x", "y")]
grid.coords <- expand.grid(seq(min(mesa.data$covars$x), max(mesa.data$covars$x), length=25),
					seq(min(mesa.data$covars$y), max(mesa.data$covars$y), length=25))
ind.inside <- point.in.polygon(grid.coords[, 1], grid.coords[, 2], hull[, 1], hull[, 2])
grid.coords <- grid.coords[ind.inside== 1, ]					
names(grid.coords) <- c("x", "y")
hull <- NULL

STdata=mesa.data
LUR=lapply(nox.model$LUR, function(x) colnames(x)[-1])
LUR[[2]] <- LUR[[3]] <- LUR[[1]]
LUR <- lapply(1:71, function(x) {
						if(x %in% 62:66){
							for(ii in 1:length(LUR)){
							 LUR[[ii]] <- LUR[[ii]][(LUR[[ii]] %in% 
							 	c("log10.m.to.road"))]
							}
							 LUR
						 } else {LUR}
						})
### add in x and y for 72:76
for(ii in 72:75){
	LUR[[ii]] <- lapply(LUR[[1]], function(x) c(x, "x", "y"))
}

ST=NULL
cov.nu=list(covf="exp", nugget=TRUE, random.effect=FALSE)
knot.coords <- vector("list", 75)
for(ii in 67:71){
	knot.coords[[ii]] <- grid.coords
}
strip=FALSE
scale=FALSE 
scale.covars=NULL
covf.beta <- c(rep("exp", 25), rep("tprs", 5), rep("iid", 1), rep("exp", 25), 
				rep("tprs", 5), rep("exp", 10), rep("exp", 4))
nugget.beta <- c(rep(TRUE, 31), rep(FALSE, 30), rep(TRUE, 10), rep(TRUE, 4))
K.beta <- c(rep(c(NA, 100, 50, 25, 10), 6), NA, rep(c(287, 100, 50, 25, 10), 8), 
			c(287, 100, 50, 25))

rho <- max(crossDist(mesa.data$covars[, c("x", "y")]))
rho[2] <- rho[1]/2
rho[3] <- rho[1]/4
rho[4] <- rho[1]/8
rho <- log(rho)

rho <- c(rep(NA, 5), rep(rho, each=5), rep(NA, 6), rep(NA, 5), rep(rho, each=5), 
			rep(NA, 15), rep(NA, 4))

#########################
##### Gen Models ####
#########################

#for(iii in 1:71){
for(iii in c(1, 72:75)){
	### Set BETA Model
	cov.beta = list(covf=covf.beta[iii], nugget=nugget.beta[iii], K=K.beta[iii])

	### Set locations
	locations=list(coords=c("x","y"), 
				long.lat=c("longitude", "latitude"), 
				coords.beta=NULL, 
				coords.nu=NULL, 
				others="type", 
				knot.coords = knot.coords[[iii]])	
	
	### LUR
	LUR.list <- LUR[[iii]]
		
	### Create model object
	models <- createSTmodel3(STdata=mesa.data, 
								LUR=LUR.list, 
								ST=NULL,
								cov.beta=cov.beta,
								cov.nu=cov.nu,
								locations=locations,
								strip=FALSE,
								scale=FALSE, 
								scale.covars=NULL)
	
	### Generate Initial Values - ensure the same initial values for all models
	if(iii == 1){ # full model with exponential everywhere
		set.seed(1983)
		n.par <- loglikeSTdim3(models)$nparam.cov
		x.init <- data.frame(	x1=rep(2, n.par), 
								x2=runif(n.par, 1, 3),
								x3=rep(log(15), n.par)	)
		x.init.save <- x.init 
	} else if (iii %in% 2:25){ # models with exp everywhere ad nuggets in betas
		x.init <- x.init.save
	} else if (iii %in% 26:30) { # tprs with exp nu and nuggets in betas
		x.init <- x.init.save[-seq(1, n.par-3, by=3), ]
	} else if (iii ==31) { # iid beta with exp nu
		x.init <- x.init.save[c(seq(3, n.par-3, by=3), (n.par-2):n.par), ]
	} else if (iii %in% 32:56){ # exp everywhere with no nugget in beta fields
		x.init <- x.init.save[-seq(3, n.par-3, by=3), ]
	} else if (iii %in% 57:61) { # tprs in betas with no nugget with exp nu
		x.init <- x.init.save[c(seq(2, n.par-3, by=3), (n.par-2):n.par), ]
	} else if (iii %in% 62:75){ # models with exp everywhere and nuggets in betas
		x.init <- x.init.save
	}
	
	### Get Fixed Values
	if(is.na(rho[iii])){
		x.fixed=NULL
	} else {
		x.fixed <- 	rep(NA, loglikeSTdim3(models)$nparam.cov)		
		ind.range <- grep(loglikeSTnames3(models, all=FALSE), pattern = "log.range")
		x.fixed[ind.range[-length(ind.range)]] <- rho[iii]
	}
	
	### Save Model output
	save(models, x.init, x.fixed, file=paste(root, "Research/Low_Rank_ST/Output/MESA_Models_", iii, ".RData", sep=""))
}
