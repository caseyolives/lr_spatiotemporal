##########################
### Additional Need Functions

### Casey Olives


## calc.CV is a function which calculates one of 3 types of cross-validated R^2 and RMSE
# The First type is for Cross Validation of AQS and MESA Fixed Sites (option="fixed")
# The second type is for COMCO Snapshot sites (option="comco")
# The third type is for HOME sites - In this option both raw and detrended varieties are provided (option="home")
# home site observations are averages per site (or ID)
# The last is for CV of all observations (option="all", the default)
calc.CV <- function(object, id.type=NULL, option="all", transform=function(x) { return(x) })
{
	# object is predCV model object
	# id.type = model$locations[, c("ID", "type")]
	# option is CV subset option
	# transform is a function used to determine the scale for RMSE and R2
	
    stCheckClass(object, "predCVSTmodel", name = "object")

    if (!is.function(transform)) {
        stop("'transform' should be a function")
    }
	if(is.null(id.type))
	{
		sub <- object$Ind.cv > 0	
	} else {
		object$pred.obs <- merge(object$pred.obs, id.type, by="ID")
		if(option == "all") {
			sub <- (1:nrow(object$pred.obs))
		} else if (option == "fixed") {
			sub <- object$pred.obs$type %in% c("AQS", "FIXED")
		} else if (option == "comco") {
			sub <- object$pred.obs$type %in% c("COMCO")		
		} else {
			sub <- object$pred.obs$type %in% c("HOME")		
		}
	}
	
	if(option == "all" | option == "fixed"){
		
		RMSE <- R2 <- matrix(NA, nr=2, nc=3)		
		rownames(RMSE) <- row.names(R2) <- c("2wk", "lta")
		colnames(RMSE) <- colnames(R2) <- c("EX.mu", "EX.mu.beta", "EX")	
		
		lta <- aggregate(object$pred.obs[sub, c("obs", "EX.mu", "EX.mu.beta", "EX")], 
						list(ID=object$pred.obs$ID[sub]), function(x) mean(transform(x)) )

		RMSE[1, ] <- apply(object$pred.obs[sub, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								sqrt( mean( (transform(x) - transform(object$pred.obs$obs[sub]))^2 )  )
							}
							)
		R2[1, ] <- apply(object$pred.obs[sub, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (transform(x) - transform(object$pred.obs$obs[sub]))^2 )/var(transform(object$pred.obs$obs[sub]))
									)
								)
							}
							)

		RMSE[2, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								sqrt( mean( (x - lta$obs)^2 )  )
							}
							)
		R2[2, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (x - lta$obs)^2 )/var(lta$obs)
										)
									)
							}
							)

	} else if (option == "comco") {
		
		dates <- unique(object$pred.obs$date[sub])		
		RMSE <- t( sapply(dates, 
					function(x) { 
						EX.mu <- sqrt( mean( (	transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) - 
										transform(object$pred.obs$EX.mu[sub][object$pred.obs[sub, "date"] == x]) )^2 )  )
						EX.mu.beta <- sqrt( mean( (	transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) - 
										transform(object$pred.obs$EX.mu.beta[sub][object$pred.obs[sub, "date"] == x]))^2 )  )
						EX <- sqrt( mean( (	transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) - 
										transform(object$pred.obs$EX[sub][object$pred.obs[sub, "date"] == x]) )^2 )  )
						c(EX.mu, EX.mu.beta, EX)
					}
				)
			)
		
		R2 <- t( sapply(dates, 
					function(x) { 
						EX.mu <- mean( 
									(transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) - 
										transform(object$pred.obs$EX.mu[sub][object$pred.obs[sub, "date"] == x]))^2 
									) / 
									var( transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) )
						EX.mu <- max(0, 1-EX.mu)
						EX.mu.beta <- max( c( 0, 1- mean( (transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) - transform(object$pred.obs$EX.mu.beta[sub][object$pred.obs[sub, "date"] == x]))^2 ) / var( transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) )))
						EX <- max( c( 0, 1- mean( (transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) - transform(object$pred.obs$EX[sub][object$pred.obs[sub, "date"] == x]))^2 ) / var( transform(object$pred.obs$obs[sub][object$pred.obs[sub, "date"] == x]) )))
						c(EX.mu, EX.mu.beta, EX)
					}
				)
			)

		rownames(RMSE) <- rownames(R2) <- as.character(dates)
		colnames(RMSE) <- colnames(R2) <- c("EX.mu", "EX.mu.beta", "EX")	
		
	} else if (option == "home") {
		RMSE <- R2 <- matrix(NA, nr=2, nc=3)		
		rownames(RMSE) <- row.names(R2) <- c("raw", "detrended")
		colnames(RMSE) <- colnames(R2) <- c("EX.mu", "EX.mu.beta", "EX")	

		sub <- sub[order(object$pred.obs$date)]
		object$pred.obs <- object$pred.obs[order(object$pred.obs$date), ]
		## mean of AQS and Fixed sites on that date
		object$pred.obs$average <- 	unlist(
										tapply(
											(1:nrow(object$pred.obs)), 
												object$pred.obs$date, 
											function(x) {
												rep(mean( object$pred.obs$obs[x][object$pred.obs$type[x] %in% c("AQS", "FIXED")]), length(x))
													}
											)
										)

		object$pred.obs$average[!sub] <- NA
	
		home.ave <- aggregate(object$pred.obs[sub, c("obs", "EX.mu", "EX.mu.beta", "EX", "average")], 
								list(ID=object$pred.obs$ID[sub]), 
								function(x) {
									mean(transform(x))
								}
							)	
		home.ave$obs.detrend <- home.ave$obs - home.ave$average
		
		print(summary(home.ave))
		
		RMSE[1, ] <- apply(home.ave[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								sqrt( mean( ( x - home.ave$obs)^2 )  )
							}
							)
		R2[1, ] <- apply(home.ave[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (x - home.ave$obs)^2 )/var(home.ave$obs)
										)
									)
							}
							)
		R2[2, ] <- apply(home.ave[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (x - home.ave$obs)^2 )/var(home.ave$obs.detrend)
										)
									)
							}
							)		
	}
	
	out <- list(RMSE=RMSE, R2=R2)
	return(out)
}

##########################
### Additional Need Functions

### Casey Olives
## calc.CV2 slight modification to summary.predCVSTmodel
calc.CV2 <- function(object, pred.naive = NULL, 
			by.date = FALSE, p = 0.95,
         	transform = function(x) { return(x) }, LTA = FALSE, 
         	 option="all", detrend=FALSE, id.type=NULL, ...){
	out <- SpatioTemporal:::summary.predCVSTmodel(object=object, pred.naive = pred.naive, 
			by.date = by.date, p = p, transform = transform, LTA = LTA, option=option, ...)
	if(detrend){
		if(is.null(id.type))
		{
			print("Require Site Type (arg id.type) to Detrend Home Obs")
		} else {
			object$pred.obs <- merge(object$pred.obs, id.type, by="ID")	
			object$pred.obs <- object$pred.obs[order(object$pred.obs$date), ]
			## mean of AQS and Fixed sites on that date
			object$pred.obs$average <- 	unlist(
											tapply(
												(1:nrow(object$pred.obs)), 
													object$pred.obs$date, 
												function(x) {
													rep(mean( object$pred.obs$obs[x][object$pred.obs$type[x] %in% c("AQS", "FIXED")]), length(x))
														}
												)
											)
			home.ave <- aggregate(object$pred.obs[, c("obs", "EX.mu", "EX.mu.beta", "EX", "average")], 
									list(ID=object$pred.obs$ID), 
									function(x) {
										mean(transform(x), na.rm=TRUE)
									}
								)	
			home.ave$obs.detrend <- home.ave$obs - home.ave$average
			tmp <-  1-out$RMSE["obs",]^2/var(home.ave$obs.detrend) 
			out$R2 =rbind(out$R2, ifelse(tmp <= 0, 0, tmp))
			
			rownames(out$R2)[nrow(out$R2)] <- "detrend"
			
		}		
	}
	return(out)
}

### block.sqrt
# This takes square root of matrix blocks in square diagonal matrix
# by computing svd of each block (UDt(V)) and then taking sqrt of diagonal elements of D
block.invsqrt <- function(X, n.blocks=nrow(X), block.sizes=rep(1, nrow(X)), type=rep("exp", n.blocks), inv=TRUE){
	ind <- 1
	for(ll in 1:n.blocks){
		temp.ind <- ind:(ind+block.sizes[ll]-1)
		if(type[ll] != "iid"){
			if(block.sizes[ll] == 1){
				if(inv){
					X[temp.ind,temp.ind] <- 1/sqrt(X[temp.ind,temp.ind])
				} else {
					X[temp.ind,temp.ind] <- sqrt(X[temp.ind,temp.ind])				
				}
			} else {
				temp <- svd(X[temp.ind,temp.ind])
				temp$d <- sqrt(temp$d)
				if(inv){
					X[temp.ind,temp.ind] <- temp$u%*%diag(1/temp$d)%*%t(temp$v)
				} else {
					X[temp.ind,temp.ind] <- temp$u%*%diag(temp$d)%*%t(temp$v)
				}
			}
		}
		ind <- ind + block.sizes[ll]
	}
	return(X)
}

### Function for generating TPRS basis functions
eta.func <- function(m, d, r){
	if(d%%2 == 0){
		out <- (-1)^(m+1+d/2)/(2^(2*m-1)*pi^(d/2)*factorial(m-1)*factorial(m-d/2))*r^(2*m-d)*log(r)
	} else {
		out <- gamma(d/2-m)/(2^(2*m)*pi^(d/2)*factorial(m-1))*r^(2*m-d)
	}
	
	out[is.na(out)] <- 0
	out
}

blockMult2<- function(X, Y, n.blocks=1, block.sizes1=nrow(X), block.sizes2=ncol(X), transpose=FALSE){
	cnt1 <- cnt2 <- 1	
	if(!transpose){
		out <- matrix(0, nr=sum(block.sizes1), nc=ncol(Y))
	} else {
		out <- matrix(0, nr=sum(block.sizes2), nc=ncol(Y))
	}
	for(iii in 1:n.blocks){
		inds1 <- cnt1:(cnt1+block.sizes1[iii]-1)
		inds2 <- cnt2:(cnt2+block.sizes2[iii]-1)		
		if(!transpose){
			if(!is.list(X)){
				out[inds1, ] <- X[inds1, inds2]%*%Y[inds2, , drop=FALSE]
			} else {
				out[inds1, ] <- X[[iii]]%*%Y[inds2, , drop=FALSE]		
			}
		} else {
			if(!is.list(X)){
				out[inds2, ] <- t(X[inds1, inds2])%*%Y[inds1, , drop=FALSE]
			} else {
				out[inds2, ] <- t(X[[iii]])%*%Y[inds1, , drop=FALSE]		
			}		
		}
		cnt1 <- cnt1 + block.sizes1[iii]
		cnt2 <- cnt2 + block.sizes2[iii]		
	}
	return(out)	
}

tZXZ <- function(X, Y, n.blocks=1, block.sizes1=nrow(X), block.sizes2=ncol(X), transpose=FALSE){
	if(!transpose){
		out <- matrix(0, nr=sum(block.sizes1), nc=sum(block.sizes1))
	} else {
		out <- matrix(0, nr=sum(block.sizes2), nc=sum(block.sizes1))
	}
	cnt1 <-  1	
	for(iii in 1:n.blocks){
		inds1 <- cnt1:(cnt1+block.sizes1[iii]-1)
		jjj <- 1
		cnt2 <- 1
		while(jjj <= iii){
			inds2 <- cnt2:(cnt2+block.sizes2[jjj]-1)		
			if(!transpose){
				if(!is.list(X)){
					out[inds1, inds1] <- X[inds1, inds2, drop=FALSE]%*%Y[inds2, inds2, drop=FALSE]%*%t(X[inds1, inds2, drop=FALSE])
				} else {
					out[inds1, inds1] <- X[[iii]]%*%Y[inds2, inds2, drop=FALSE]%*%t(X[[iii]])		
				}
			} else {
				if(!is.list(X)){
					out[inds2, inds2] <- t(X[inds1, inds2, drop=FALSE])%*%Y[inds1, inds1, drop=FALSE]%*%X[inds1, inds2, drop=FALSE]
				} else {
					out[inds2, inds2] <- t(X[[iii]])%*%Y[inds1, inds1, drop=FALSE]%*%X[[iii]]		
				}		
			}
			cnt2 <- cnt2 + block.sizes2[jjj]		
			jjj <- jjj+1
		}
		cnt1 <- cnt1 + block.sizes1[iii]

	}
	return(out)	
}
#########################
#########################

makeSigmaBTPRS <- function (pars, coords1, coords2=coords1, basis.coords=coords1,
							m=2, fit.tprs=TRUE, out.ZB=FALSE, K=nrow(basis.coords), pen.form=FALSE) {
    ## Just constructs one block at a time
    ## pars is just a scalar, the sill
    ## m = parameter in TPRS Smooth (m=2 is for cubic splines)
    ## fit.tprs = TRUE indicates that the output is for ESTIMATION TPRS PARAMETERS, 
    ##				i.e. SigmaB = U_k D_k W_k (W_k' D_k W_k)^{-1} W_k' D_k' U_k'
    ##			= FALSE indicates that the output is for PREDICTION USING TPRS
    ##				i.e. SigmaB = E U_k W_k (W_k' D_k W_k)^{-1} W_k' U_k' E'
    ## where E is dim nrow(coords2) by nrow(coords1)
    ## U_k is dim nrow(coords1) by nrow(coordsK)
    ## W_k is dim nrow(coordsK) by nrow(coordsK)-M
    ## gen.new means we generate for a new set of observations.
    ## i.e. when fit.tprs=FALSE and gen.new=TRUE then output
    ##	E2 <- eta.func(r=crossDist(coords2, coords1), d=2, m=m)
	##	SigmaTPRS = E2%*%Uk%*%Zk%*%O%*%t(O)%*%t(Zk)%*%t(Uk)%*%t(E2)*pars
	##	Otherwise, output
	##	SigmaTPRS = E%*%Uk%*%Zk%*%O%*%t(O)%*%t(Zk)%*%t(Uk)%*%t(E2)*pars
    ## Note that the eigen decompostion is ALWAYS based on coords1
	## NOTE: Assumes that T is in LUR already, based on createSTmodel3
	
	if(fit.tprs & nrow(coords1) != nrow(coords2)){
		print("When fitting TPRS, require coords1 == coords2")
		stop
	}
	
	T <- cbind(1, basis.coords)
	colnames(T) <- c("Int", "x.beta", "y.beta")
	if(m > 2){
		for(jj in 2:(m-1)){
			temp.T <- coords1^jj
			colnames(temp.T) <- c(paste("x.beta", jj, sep=""), paste("y.beta", jj, sep=""))
			T <- cbind(T, temp.T)
			rm(temp.T)
		}
	}
	T <- as.matrix(T)
	M <- choose(m+2-1, 2) #d=2
	E <- eta.func(r=crossDist(basis.coords), d=2, m=m)
	Ek <- eigen(E)
	Uk <- Ek$vectors[, 1:K]
	Dk <- diag(Ek$values[1:K])
	tUkT <- t(Uk)%*%T
	QR <- qr(tUkT)
	Zk <- qr.Q(QR, complete=TRUE)[, (M+1):K]
	O <- t(Zk)%*%Dk%*%Zk
	O <- svd(O)
	O$d <- sqrt(O$d)
	O <- solve(O$u%*%diag(O$d)%*%t(O$v))
	if(!pen.form){			
		if(fit.tprs){
			SigmaTPRS = Uk%*%Dk%*%Zk%*%O
			if(!out.ZB){
				SigmaTPRS = SigmaTPRS%*%t(SigmaTPRS)*pars
			}
		} else {
			E1 <- eta.func(r=crossDist(coords1, basis.coords), d=2, m=m)
			E2 <- eta.func(r=crossDist(coords2, basis.coords), d=2, m=m)
			SigmaTPRS = E1%*%Uk%*%Zk%*%O		
			if(!out.ZB){
				SigmaTPRS = SigmaTPRS%*%O%*%t(Zk)%*%t(Uk)%*%t(E2)*pars
			}
		}
		return(SigmaTPRS)    
	} else {
		zb = Uk%*%Dk%*%Zk
		pen = t(Zk)%*%Dk%*%Zk
		return(list(zb=zb, pen=pen))
	}			
}

### makeSigmaB3 
# Adapt makeSigmaB to allow for low rank kriging model
makeSigmaB3 <- function (pars, coords1, coords2=coords1, coordsK=NULL, basis.coords=coords1, 
				type = "exp", nugget = 0, symmetry = dim(coords1)[1] ==  dim(coords2)[1], 
				ind2.to.1 = 1:dim(coords2)[1], m=2, fit.tprs=TRUE) {
    ## m = parameter in TPRS Smooth (m=2 is for cubic splines)
    ## fit.tprs = TRUE indicates that the output is for ESTIMATION TPRS PARAMETERS, 
    ##				i.e. SigmaB = U_k D_k W_k (W_k' D_k W_k)^{-1} W_k' D_k' U_k'
    ##			= FALSE indicates that the output is for PREDICTION USING TPRS
    ##				i.e. SigmaB = E U_k W_k (W_k' D_k W_k)^{-1} W_k' U_k' E'
    ## where E is dim nrow(coords2) by nrow(coords1)
    ## U_k is dim nrow(coords1) by nrow(coordsK)
    ## W_k is dim nrow(coordsK) by nrow(coordsK)-M
    ## Note that the eigen decompostion is based on basis.coords, typically coords1
    
	if(is.null(coordsK))
	{
		if(all(type != "tprs"))
		{
			SigmaB=makeSigmaB(pars, dist=crossDist(coords1, coords2), 
					type=type, nugget=0, symmetry=symmetry, ind2.to.1=ind2.to.1)
		} else {
			SigmaB=matrix(0, nr=nrow(coords1)*length(type), nc=nrow(coords2)*length(type))
			ind1=ind2=1
			for(ii in 1:length(type))
			{
				inds1 = ind1:(ind1+nrow(coords1)-1)
				inds2 = ind2:(ind2+nrow(coords2)-1)
				if(type[ii] != "tprs"){
					SigmaB[inds1, inds2] = makeSigmaB(pars[[ii]], 
											dist=crossDist(coords1, coords2), 
											type=type[ii], nugget=0, 
											symmetry=symmetry, ind2.to.1=ind2.to.1)
				} else {
					SigmaB[inds1, inds2] = makeSigmaBTPRS(pars[[ii]][parsCovFuns3(type=type[ii]) == 'sill'], 
						coords1=coords1, coords2=coords2, m=m, fit.tprs=fit.tprs, basis.coords=basis.coords)
				}
				ind1=ind1+nrow(coords1)
				ind2=ind2+nrow(coords2)						
			} # end for loop
		} # end if(all(type != "tprs"))
	} else { # if(!is.null(coordsK)
		temp.pars1 <- lapply(1:length(pars), function(x) { if(type[x] != "iid"){ out = pars[[x]]; out[parsCovFuns3(type[x]) == 'sill'] = sqrt(out[parsCovFuns3(type[x]) == 'sill']) } else { out= numeric(0)} ; return(out) })
		temp.pars2 <- lapply(1:length(pars), function(x) { if(type[x] != "iid"){ out = pars[[x]]; out[parsCovFuns3(type[x]) == 'sill'] = 1 } else { out= numeric(0)} ; return(out) })
		if(all(type != "tprs")){
			ZB1 <- makeSigmaB(temp.pars1, 
					dist = crossDist(coords1, coordsK), 
					type = type, nugget = 0) 
			ZB2 <- makeSigmaB(temp.pars1, 
					dist = crossDist(coords2, coordsK), 
					type = type, nugget = 0) 
			O <- makeSigmaB(temp.pars2, 
					dist = crossDist(coordsK, coordsK), 
					type = type, nugget = 0)
			O <- block.invsqrt(O, n.blocks=length(type), 
					block.sizes=rep(nrow(coordsK), length(type)), type=type)
			
			SigmaB1 <- t( blockMult( t(O), t(ZB1), n.blocks=length(type), 
					block.sizes=rep(nrow(coordsK), length(type)) ) )	
			SigmaB2 <- blockMult( t(O), t(ZB2), n.blocks=length(type), 
					block.sizes=rep(nrow(coordsK), length(type)) ) 

			SigmaB <- SigmaB1 %*% SigmaB2
			rm(temp.pars1, temp.pars2, SigmaB1, SigmaB2)
		} else { # if (!all(type != "tprs"))
			SigmaB=matrix(0, nr=nrow(coords1)*length(type), nc=nrow(coords2)*length(type))
			ind1=ind2=1
			for(ii in 1:length(type)){
				inds1 = ind1:(ind1+nrow(coords1)-1)
				inds2 = ind2:(ind2+nrow(coords2)-1)
				if(type[ii] != "tprs"){
					ZB1 <- makeSigmaB(temp.pars1[[ii]], 
							dist = crossDist(coords1, coordsK), 
							type = type[ii], nugget = 0) 
					ZB2 <- makeSigmaB(temp.pars1[[ii]], 
							dist = crossDist(coords2, coordsK), 
							type = type[ii], nugget = 0) 
					O <- makeSigmaB(temp.pars2[[ii]], 
							dist = crossDist(coordsK, coordsK), 
							type = type[ii], nugget = 0)
 					O <- block.invsqrt(O, n.blocks=length(type[ii]), 
 							block.sizes=rep(nrow(coordsK), length(type[ii])), type=type[ii])
					SigmaB[inds1, inds2] = ZB1%*%O%*%t(O)%*%t(ZB2)
				} else {
					SigmaB[inds1, inds2] = makeSigmaBTPRS(pars[[ii]][parsCovFuns3(type=type[ii]) == 'sill'], 
						coords1=coords1, coords2=coords2, m=m, 
						fit.tprs=fit.tprs, basis.coords=basis.coords, K=nrow(coordsK))
				}
				ind1=ind1+nrow(coords1)
				ind2=ind2+nrow(coords2)						
			} # end for loop				
		} # end all(type != "tprs")
		
		### Add Nuggets
		if(any(nugget != 0)){
			temp.pars <- lapply(1:length(type), function(x) numeric(0))
			SigmaB = SigmaB + makeSigmaB(temp.pars, dist=crossDist(coords1, coords2), 
							type=rep("iid", length(type)),
							nugget=nugget, symmetry=symmetry, ind2.to.1=ind2.to.1) 
			rm(temp.pars)
		} # end if(any(nugget != 0)) 
	} # end if (is.null(coordsK)
	return(SigmaB)
}

### new function to make ZB
makeZB <- function (pars, coords1, coords2=coords1, coordsK=NULL, basis.coords=coords1, 
				type = "exp", symmetry = dim(coords1)[1] ==  dim(coords2)[1], 
				ind2.to.1 = 1:dim(coords2)[1], m=2, fit.tprs=TRUE) {
    ## m = parameter in TPRS Smooth (m=2 is for cubic splines)
    ## fit.tprs = TRUE indicates that the output is for ESTIMATION TPRS PARAMETERS, 
    ##				i.e. SigmaB = U_k D_k W_k (W_k' D_k W_k)^{-1} W_k' D_k' U_k'
    ##			= FALSE indicates that the output is for PREDICTION USING TPRS
    ##				i.e. SigmaB = E U_k W_k (W_k' D_k W_k)^{-1} W_k' U_k' E'
    ## where E is dim nrow(coords2) by nrow(coords1)
    ## U_k is dim nrow(coords1) by nrow(coordsK)
    ## W_k is dim nrow(coordsK) by nrow(coordsK)-M
    ## Note that the eigen decompostion is based on basis.coords, typically coords1
    if(any(type == "iid")){
    	stop("iid beta fields do not use ZB for construction")
    }
    
    # Set sills to one to get 
    pars <- lapply(pars, function(z) {z['sill'] = 1; return(z)})
	if(is.null(coordsK))
	{
		if(all(type != "tprs"))
		{
			SigmaB=makeSigmaB(pars, dist=crossDist(coords1, coords2), 
					type=type, nugget=0, symmetry=symmetry, ind2.to.1=ind2.to.1)
		} else {
			SigmaB=matrix(0, nr=nrow(coords1)*length(type), nc=nrow(coords2)*length(type))
			ind1=ind2=1
			for(ii in 1:length(type))
			{
				inds1 = ind1:(ind1+nrow(coords1)-1)
				inds2 = ind2:(ind2+nrow(coords2)-1)
				if(type[ii] != "tprs"){
					SigmaB[inds1, inds2] = makeSigmaB(pars[[ii]], 
											dist=crossDist(coords1, coords2), 
											type=type[ii], nugget=0, 
											symmetry=symmetry, ind2.to.1=ind2.to.1)
				} else {
					SigmaB[inds1, inds2] = makeSigmaBTPRS(pars[[ii]][parsCovFuns3(type=type[ii]) == 'sill'], 
						coords1=coords1, coords2=coords2, m=m, fit.tprs=fit.tprs, basis.coords=basis.coords, 
						out.ZB=TRUE)
				}
				ind1=ind1+nrow(coords1)
				ind2=ind2+nrow(coords2)						
			} # end for loop
		} # end if(all(type != "tprs"))
	} else { # if(!is.null(coordsK)
		if(all(type != "tprs")){
	
			ZB1 <- makeSigmaB(pars, 
					dist = crossDist(coords1, coordsK), 
					type = type, nugget = 0) 
			O <- makeSigmaB(pars, 
					dist = crossDist(coordsK, coordsK), 
					type = type, nugget = 0)
			O <- block.invsqrt(O, n.blocks=length(type), 
					block.sizes=rep(nrow(coordsK), length(type)), type=type)
			
			SigmaB <- t( blockMult( t(O), t(ZB1), n.blocks=length(type), 
					block.sizes=rep(nrow(coordsK), length(type)) ) )	

		} else { # if (!all(type != "tprs"))
			num.tprs <- sum(type == "tprs")
			M <- choose(m+2-1, 2)
			SigmaB=matrix(0, nr=nrow(coords1)*length(type), nc=nrow(coordsK)*(length(type) - num.tprs) + (nrow(coordsK) - M)*num.tprs)
			ind1=ind2=1
			for(ii in 1:length(type)){
				inds1 = ind1:(ind1+nrow(coords1)-1)
				if(type[ii] != "tprs"){
					inds2 = ind2:(ind2+nrow(coordsK)-1)
				} else {
					inds2 = ind2:(ind2+nrow(coordsK)-M-1)
				}
				if(type[ii] != "tprs"){
					ZB1 <- makeSigmaB(pars[[ii]], 
							dist = crossDist(coords1, coordsK), 
							type = type[ii], nugget = 0) 
					O <- makeSigmaB(pars[[ii]], 
							dist = crossDist(coordsK, coordsK), 
							type = type[ii], nugget = 0)
					O <- block.invsqrt(O, n.blocks=length(type[ii]), 
							block.sizes=rep(nrow(coordsK), length(type[ii])), type=type[ii])
					
					SigmaB[inds1, inds2] = ZB1%*%O
				} else {
					SigmaB[inds1, inds2] = makeSigmaBTPRS(pars[[ii]][parsCovFuns3(type=type[ii]) == 'sill'], 
						coords1=coords1, coords2=coords2, 
						m=m, fit.tprs=fit.tprs, basis.coords=basis.coords, out.ZB=TRUE, K=nrow(coordsK))
				}
				ind1=ind1+nrow(coords1)
				if(type[ii] != "tprs"){
					ind2=ind2+nrow(coordsK)						
				} else {
					ind2=ind2+nrow(coordsK)-M						
				}
			} # end for loop				
		} # end all(type != "tprs")
	} # end if (is.null(coordsK)
	return(SigmaB)
}

################
### New loglike to take advantage of computational form without nuggets in the beta field

loglikeST_opt <- function (x = NULL, STmodel, type = "p", x.fixed = NULL) 
{
    stCheckClass(STmodel, "STmodel", name = "STmodel")
    type <- tolower(type)
    x <- SpatioTemporal:::stCheckLoglikeIn(x, x.fixed, type)
    if (is.null(x) || any(is.na(x))) {
        return(loglikeSTnames3(STmodel, all = (type == "f")))
    }
    dimensions <- loglikeSTdim3(STmodel)
    if ((type == "f" && length(x) != dimensions$nparam) || (type != 
        "f" && length(x) != dimensions$nparam.cov)) {
        stop("Number of parameters, length(x), incorrect.")
    }
    tmp <- loglikeSTgetPars3(x, STmodel)
    if (type == "f") {
        gamma <- tmp$gamma
        alpha <- tmp$alpha
    }
    cov.pars.beta <- tmp$cov.beta
    cov.pars.nu <- tmp$cov.nu
    Y <- STmodel$obs$obs
    if(any(STmodel$cov.beta$nugget)){
    	if(!all(STmodel$cov.beta$nugget)){
	    	stop("loglikeST_opt only works when either all or no beta fields have nuggets present")
    	}
    }
    
    if(!all(STmodel$cov.beta != "iid") & !all(STmodel$cov.beta$nugget)){
    	stop("Require either spatial smoothing of all beta fields OR nuggets in all beta fields (or both). Cannot specify no spatial smoothing without a nugget")
    }
    
    if(any(STmodel$cov.beta$covf == "tprs")){
    	m=2
    	STmodel$cov.beta$K <- STmodel$cov.beta$K - choose(m+2-1, 2)*(STmodel$cov.beta$covf=="tprs")
    }
    if (type == "f") {
		stop("FULL Likelihood (type = 'f') not implemented for low rank versions")
    } else if (type == "r") {
		stop("REML Likelihood (type = 'r') not implemented for low rank versions")
	}
    
    sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = STmodel$D.nu, 
        type = STmodel$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = STmodel$nt, 
        ind1 = STmodel$obs$idx)
    sigma.nu <- try(makeCholBlock(sigma.nu, n.blocks = dimensions$T, 
        block.sizes = STmodel$nt), silent = TRUE)
    if (class(sigma.nu) == "try-error") {
        return(-.Machine$double.xmax)
    }
    l <- - sumLogDiag(sigma.nu)
    sigma.nu <- invCholBlock(sigma.nu, block.sizes = STmodel$nt)

	i.sR.Y <- blockMult(sigma.nu, Y, block.sizes = STmodel$nt)
	F.i.sR.F <- calc.tFXF(STmodel$F, sigma.nu, STmodel$obs$idx, 
							block.sizes = STmodel$nt, n.loc = dimensions$n.obs)
	
	if(all(STmodel$cov.beta$covf != "iid")){
		sigma.B <- rep(sqrt(sapply(cov.pars.beta$pars, function(y) return(y['sill']))), 
						each=STmodel$cov.beta$K[1])
		l <- l-sum(log(sigma.B))
		sigma.B <- diag(1/sigma.B^2)
	
		ZB <- makeZB(cov.pars.beta$pars, coords1=STmodel$locations[, c("x.beta", "y.beta")],
		coordsK=STmodel$locations.list$knot.coords[STmodel$cov.beta$knots[[1]],], 
		type = STmodel$cov.beta$covf)

		if(all(STmodel$cov.beta$nugget)){
			sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)
			l <- l - sum(log(sigma.P))
			sigma.P.Y <- F.i.sR.F + diag(1/sigma.P^2)
			
			sigma.P.Y <- try(makeCholBlock(sigma.P.Y, n.blocks=1),silent = TRUE)
			if (class(sigma.P.Y) == "try-error") {
				return(-.Machine$double.xmax)
			}
			l <- l - sumLogDiag(sigma.P.Y)
	
			tF.i.sR <- calc.tFX(STmodel$F, sigma.nu, STmodel$obs$idx, n.loc=dimensions$n.obs)
			tF.i.sR.Y <- tF.i.sR %*% Y
			i.sP.Y.tF.i.sR.Y <- solveTriBlock(sigma.P.Y, tF.i.sR.Y, transpose = TRUE)
			i.sP.Y.tF.i.sR.Y <- solveTriBlock(sigma.P.Y, i.sP.Y.tF.i.sR.Y, transpose = FALSE)
	
			i.sP.Y.tF.i.sR.F <- solveTriBlock(sigma.P.Y, F.i.sR.F, transpose = TRUE)
	
			i.sR.Y <- -t(tF.i.sR) %*% i.sP.Y.tF.i.sR.Y + i.sR.Y		
			F.i.sR.F <- -t(i.sP.Y.tF.i.sR.F) %*% i.sP.Y.tF.i.sR.F + F.i.sR.F			
		} 

		F.i.sR.Y <- calc.tFX(STmodel$F, i.sR.Y, STmodel$obs$idx, 
			n.loc = dimensions$n.obs)
		F.i.sR.F.X <- t(blockMult2(STmodel$LUR, F.i.sR.F, n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=dimensions$p, transpose=TRUE))
	
		sigma.B.Y <- blockMult2(ZB, F.i.sR.F, n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=STmodel$cov.beta$K, transpose=TRUE)
		sigma.B.Y <- blockMult2(ZB, t(sigma.B.Y), n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=STmodel$cov.beta$K, transpose=TRUE)
		sigma.B.Y <- sigma.B.Y + sigma.B
		sigma.B.Y <- try(makeCholBlock(sigma.B.Y, n.blocks = 1), 
			silent = TRUE)
		if (class(sigma.B.Y) == "try-error") {
			return(-.Machine$double.xmax)
		}
		l <- l - sumLogDiag(sigma.B.Y)
		
		ZB.F.i.sR.F.X <- blockMult2(ZB, F.i.sR.F.X, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=STmodel$cov.beta$K, transpose=TRUE)
		sigma.B.Y.ZB.F.i.sR.F.X <- solveTriBlock(sigma.B.Y, ZB.F.i.sR.F.X, transpose = TRUE)
		i.sigma.alpha.Y <- -t(sigma.B.Y.ZB.F.i.sR.F.X)%*%sigma.B.Y.ZB.F.i.sR.F.X +
							blockMult2(STmodel$LUR, F.i.sR.F.X, n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=dimensions$p, transpose=TRUE)
		i.sigma.alpha.Y <- try(makeCholBlock(i.sigma.alpha.Y), 
			silent = TRUE)
		if (class(i.sigma.alpha.Y) == "try-error") {
			return(-.Machine$double.xmax)
		}
		Z.F.i.sR.Y <- blockMult2(ZB, F.i.sR.Y, n.blocks=dimensions$m, 
								block.sizes1=rep(dimensions$n, dimensions$m), 
								block.sizes2=STmodel$cov.beta$K, transpose=TRUE)
		Y.hat <- solveTriBlock(sigma.B.Y, Z.F.i.sR.Y, transpose = TRUE)
		Y.hat.2 <- t(sigma.B.Y.ZB.F.i.sR.F.X)%*%Y.hat
		Y.hat.2 <- -Y.hat.2 + blockMult2(STmodel$LUR, F.i.sR.Y, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=dimensions$p, transpose=TRUE)
		Y.hat.2 <- solveTriBlock(i.sigma.alpha.Y, Y.hat.2, transpose = TRUE)		
	
	} else { # all fields are iid here (i.e. with nugget
		sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)
		l <- l - sum(log(sigma.P))
		sigma.P.Y <- F.i.sR.F + diag(1/sigma.P^2)
		sigma.P.Y <- try(makeCholBlock(sigma.P.Y, n.blocks=1),silent = TRUE)
		if (class(sigma.P.Y) == "try-error") {
			return(-.Machine$double.xmax)
		}
		l <- l - sumLogDiag(sigma.P.Y)

		F.i.sR.Y <- calc.tFX(STmodel$F, i.sR.Y, STmodel$obs$idx, 
			n.loc = dimensions$n.obs)
		F.i.sR.F.X <- t(blockMult2(STmodel$LUR, F.i.sR.F, n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=dimensions$p, transpose=TRUE))
		sigma.P.Y.F.i.sR.F.X <- solveTriBlock(sigma.P.Y, F.i.sR.F.X, transpose = TRUE)
		i.sigma.alpha.Y <- -t(sigma.P.Y.F.i.sR.F.X)%*%sigma.P.Y.F.i.sR.F.X +
							blockMult2(STmodel$LUR, F.i.sR.F.X, n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=dimensions$p, transpose=TRUE)
		i.sigma.alpha.Y <- try(makeCholBlock(i.sigma.alpha.Y), 
			silent = TRUE)
		if (class(i.sigma.alpha.Y) == "try-error") {
			return(-.Machine$double.xmax)
		}
		Y.hat <- solveTriBlock(sigma.P.Y, F.i.sR.Y, transpose = TRUE)
		Y.hat.2 <- t(sigma.P.Y.F.i.sR.F.X)%*%Y.hat
		Y.hat.2 <- -Y.hat.2 + blockMult2(STmodel$LUR, F.i.sR.Y, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=dimensions$p, transpose=TRUE)
		Y.hat.2 <- solveTriBlock(i.sigma.alpha.Y, Y.hat.2, transpose = TRUE)		
	}
	Y.sigma.hat.Y <- dotProd(Y, i.sR.Y) - norm2(Y.hat) - norm2(Y.hat.2)
	l <- l - Y.sigma.hat.Y/2
# 	if (dimensions$L != 0) {
# 		M.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.M
# 		M.hat.2 <- solveTriBlock(i.sigma.alpha.Y, M.hat.2, 
# 			transpose = TRUE)
# 		Y.sigma.hat.M <- (t(Y) %*% i.sR.M - t(Y.hat) %*% 
# 			M.hat - t(Y.hat.2) %*% M.hat.2)
# 		Y.sigma.hat.M <- t(Y.sigma.hat.M)
# 		M.sigma.hat.M <- (t(STmodel$ST) %*% i.sR.M - t(M.hat) %*% 
# 			M.hat - t(M.hat.2) %*% M.hat.2)
# 		M.sigma.hat.M <- try(makeCholBlock(M.sigma.hat.M), 
# 			silent = TRUE)
# 		if (class(M.sigma.hat.M) == "try-error") {
# 			return(-.Machine$double.xmax)
# 		}
# 		if (type == "r") {
# 			l <- l - sumLogDiag(M.sigma.hat.M)
# 		}
# 		Y.sigma.hat.M <- solveTriBlock(M.sigma.hat.M, Y.sigma.hat.M, 
# 			transpose = TRUE)
# 		l <- l + norm2(Y.sigma.hat.M)/2
# 	}

    l <- as.double(l)
    if (!is.finite(l)) 
        l <- -.Machine$double.xmax
    return(l)
}

loglikeSTGrad_opt <- function (x, STmodel, type = "p", x.fixed = NULL, h = 0.001, 
    diff.type = 0) {
    func <- function(x0) {
        loglikeST_opt(x0, STmodel, type, x.fixed)
    }
    df <- genGradient(x, func, h = h, diff.type = diff.type)
    return(df)
}

estimate.STmodel_opt <- function (object, x, x.fixed = NULL, type = "p", h = 0.001, diff.type = 1, 
    hessian.all = FALSE, lower = -15, upper = 15, method = "L-BFGS-B", 
    control = list(trace = 3, maxit = 1000), ...) {
    
    stCheckClass(object, "STmodel", name = "object")
    dimensions <- loglikeSTdim3(object)
    type <- tolower(type)
    SpatioTemporal:::stCheckType(type)
    x <- as.matrix(x)
    if (length(x) == 1) {
        x <- matrix(seq(-5, 5, length.out = x), dimensions$nparam.cov, 
            x, byrow = TRUE)
    }
    tmp <- stCheckX3(x, x.fixed, dimensions, type, object)
    x.all <- tmp$x.all
    x <- tmp$x
    x.fixed <- tmp$x.fixed
    control <- defaultList(control, eval(formals(estimate.STmodel3)$control))
    loglikeSTGrad.loc <- function(x, STmodel, type, x.fixed) {
        	loglikeSTGrad_opt(x, STmodel, type, x.fixed, h = h, diff.type = diff.type)
    	}
    res <- as.list(rep(NA, dim(x)[2]))
    conv <- rep(FALSE, dim(x)[2])
    value <- rep(NA, dim(x)[2])
    err <- tryCatch(loglikeST_opt(x[, 1], STmodel = object, type = type, 
        x.fixed = x.fixed), silent = TRUE)
    if (inherits(err, "try-error")) {
        stop(paste("log-likelihood fails at first starting point with:\n", 
            err[[1]]))
    }
    control$fnscale <- -1
    for (i in 1:dim(x)[2]) {
        if (control$trace != 0) {
            message(paste("Optimisation using starting value ", 
                i, "/", dim(x)[2], sep = ""))
        }
        try(res[[i]] <- optim(x[, i], loglikeST_opt, gr = loglikeSTGrad.loc, 
            STmodel = object, type = type, x.fixed = x.fixed, 
            method = method, control = control, hessian = TRUE, 
            lower = lower, upper = upper), silent = TRUE)
        if (all(!is.na(res[[i]]))) {
            conv[i] <- (res[[i]]$convergence == 0 && all(eigen(res[[i]]$hessian)$value < 
                -1e-10))
            value[i] <- res[[i]]$value
            res[[i]]$conv <- conv[i]
            res[[i]]$par.cov <- data.frame(par = double(dimensions$nparam.cov), 
                sd = NA, fixed = double(dimensions$nparam.cov), 
                init = double(dimensions$nparam.cov), tstat = double(dimensions$nparam.cov))
            res[[i]]$par.all <- data.frame(par = double(dimensions$nparam), 
                sd = NA, fixed = double(dimensions$nparam), init = double(dimensions$nparam), 
                tstat = double(dimensions$nparam))
#             suppressWarnings(par.sd <- sqrt(-diag(solve(res[[i]]$hessian))))
			if(res[[i]]$conv){
	            par.sd <- try(sqrt(-diag(solve(res[[i]]$hessian))), silent=TRUE)
	        } else {
	        	par.sd <- rep(NA, nrow(res[[i]]$hessian))
	        }
            print("bug1")
            print(par.sd)
            if (type != "f") {
                par.type <- "par.cov"
            } else {
                par.type <- "par.all"
            }
            res[[i]][[par.type]]$init <- x.all[, i]
            res[[i]][[par.type]]$par <- res[[i]][[par.type]]$fixed <- x.fixed
            res[[i]][[par.type]]$par[is.na(x.fixed)] <- res[[i]]$par
            res[[i]][[par.type]]$sd[is.na(x.fixed)] <- par.sd
            if (type != "f") {
                 tmp <- predict.STmodel_opt(object, res[[i]]$par.cov$par, 
                   only.pars = TRUE, pred.var = FALSE, type = type)$pars
                res[[i]]$par.all$par <- c(tmp$gamma.E, tmp$alpha.E, 
                  res[[i]]$par.cov$par)
                N.reg <- length(tmp$gamma.E) + length(tmp$alpha.E)
                res[[i]]$par.all$sd <- c(rep(NA, N.reg), res[[i]]$par.cov$sd)
                res[[i]]$par.all$init <- c(rep(NA, N.reg), x.all[, 
                  i])
                res[[i]]$par.all$fixed <- c(rep(NA, N.reg), x.fixed)
            }
            else {
                I <- (dimensions$nparam - dimensions$nparam.cov + 
                  1):dimensions$nparam
                res[[i]]$par.cov <- res[[i]]$par.all[I, , drop = FALSE]
            }
            res[[i]]$par.cov$tstat <- res[[i]]$par.cov$par/res[[i]]$par.cov$sd
            res[[i]]$par.all$tstat <- res[[i]]$par.all$par/res[[i]]$par.all$sd
            rownames(res[[i]]$par.all) <- loglikeSTnames3(object, 
                all = TRUE)
            rownames(res[[i]]$par.cov) <- loglikeSTnames3(object, 
                all = FALSE)
            if (type != "f") {
                names(res[[i]]$par) <- loglikeSTnames3(object, 
                  all = FALSE)[is.na(x.fixed)]
            }
            else {
                names(res[[i]]$par) <- loglikeSTnames3(object, 
                  all = TRUE)[is.na(x.fixed)]
            }
        }
    }
    if (all(is.na(res))) {
        stop("All optimisations failed, consider trying different starting values.")
    }
    status <- data.frame(value = value, convergence = logical(length(res)), 
        conv = (conv == 1))
    par.cov <- matrix(NA, dimensions$nparam.cov, length(res))
    par.all <- matrix(NA, dimensions$nparam, length(res))
    for (i in 1:length(res)) {
        if (all(!is.na(res[[i]]))) {
            status$convergence[i] <- res[[i]]$convergence == 
                0
            par.cov[, i] <- res[[i]]$par.cov$par
            par.all[, i] <- res[[i]]$par.all$par
        }
    }
    rownames(par.all) <- loglikeSTnames3(object, all = TRUE)
    rownames(par.cov) <- loglikeSTnames3(object, all = FALSE)
    Ind.overall <- which.max(value)
    if (any(conv == TRUE)) {
        value[!conv] <- NA
    }
    Ind <- which.max(value)
    res.best <- res[[Ind]]
    summary <- list(status = status, par.all = par.all, par.cov = par.cov, 
        x.fixed = x.fixed)
    if (hessian.all == TRUE) {
        if (type != "f") {
            x.fixed <- res.best$par.all$fixed
            x <- res.best$par.all$par[is.na(x.fixed)]
            res.best$hessian.all <- loglikeSTHessian(x, object, 
                type = "f", x.fixed = x.fixed, h = h)
            suppressWarnings(par.sd <- sqrt(-diag(solve(res.best$hessian.all))))
            res.best$par.all$sd <- NA
            res.best$par.all$sd[is.na(x.fixed)] <- par.sd
            res.best$par.all$tstat <- res.best$par.all$par/res.best$par.all$sd
        }
        else {
            res.best$hessian.all <- res.best$hessian
        }
    }
    out <- list(res.best = res.best, res.all = res, summary = summary)
    class(out) <- "estimateSTmodel"
    return(out)
}

estimateCV.STmodel_opt <- function (object, x, Ind.cv, control = list(trace = 3), ...) {

    stCheckClass(object, "STmodel", name = "object")
    Ind.cv <- SpatioTemporal:::stCheckInternalCV(Ind.cv)
    control <- defaultList(control, eval(formals(estimate.STmodel_opt)$control))
    Ind.cv <- as.matrix(Ind.cv)
    if (dim(Ind.cv)[2] == 1) {
        N.CV.sets <- max(Ind.cv, na.rm = TRUE)
    }else {
        N.CV.sets <- dim(Ind.cv)[2]
    }
    dimensions <- loglikeSTdim3(object)
    res <- list()
    for (i in 1:N.CV.sets) {
        if (dim(Ind.cv)[2] == 1) {
            Ind.current <- Ind.cv == i
        }
        else {
            Ind.current <- as.logical(Ind.cv[, i])
        }
        object.aux <- dropObservations3(object, Ind.current)
        if (control$trace != 0) {
            message(paste("Estimation of CV-set ", i, "/", N.CV.sets, 
                sep = ""))
        }
        res[[i]] <- estimate.STmodel_opt(object.aux, x, ...)
    }
    status <- data.frame(value = sapply(res, function(x) {
        x$res.best$value
    }), convergence = sapply(res, function(x) {
        x$res.best$convergence == 0
    }), conv = sapply(res, function(x) {
        x$res.best$conv
    }))
    tmp <- sapply(res, function(x) {
        range(-eigen(x$res.best$hessian)$value)
    })
    status$eigen.min <- tmp[1, ]
    if (is.null(res[[1]]$res.best$hessian.all)) {
        status$eigen.all.min <- NA
    }
    else {
        tmp <- sapply(res, function(x) {
            range(-eigen(x$res.best$hessian.all)$value)
        })
        status$eigen.all.min <- tmp[1, ]
    }
    par.cov <- matrix(NA, dim(res[[1]]$res.best$par.cov)[1], 
        N.CV.sets)
    par.cov.sd <- matrix(NA, dim(res[[1]]$res.best$par.cov)[1], 
        N.CV.sets)
    par.all <- matrix(NA, dim(res[[1]]$res.best$par.all)[1], 
        N.CV.sets)
    par.all.sd <- matrix(NA, dim(res[[1]]$res.best$par.all)[1], 
        N.CV.sets)
    res.all <- list()
    for (i in 1:N.CV.sets) {
        par.cov[, i] <- res[[i]]$res.best$par.cov$par
        par.cov.sd[, i] <- res[[i]]$res.best$par.cov$sd
        par.all[, i] <- res[[i]]$res.best$par.all$par
        par.all.sd[, i] <- res[[i]]$res.best$par.all$sd
        res.all[[i]] <- list(res.best = res[[i]]$res.best, summary = res[[i]]$summary)
    }
    rownames(par.cov) <- rownames(res[[1]]$res.best$par.cov)
    rownames(par.cov.sd) <- rownames(res[[1]]$res.best$par.cov)
    rownames(par.all) <- rownames(res[[1]]$res.best$par.all)
    rownames(par.all.sd) <- rownames(res[[1]]$res.best$par.all)
    out <- list(par.cov = par.cov, par.cov.sd = par.cov.sd, par.all = par.all, 
        par.all.sd = par.all.sd, res.all = res.all, status = status, 
        Ind.cv = Ind.cv, x.fixed = res.all[[1]]$summary$x.fixed)
    class(out) <- "estCVSTmodel"
    return(out)
}

predict.STmodel_opt <- function (object, x, STdata = NULL, Nmax = 1000, only.pars = FALSE, 
    nugget.unobs = 0, only.obs = FALSE, pred.var = TRUE, pred.covar = FALSE, 
    beta.covar = FALSE, combine.data = FALSE, type = "p", ...) {
 
    stCheckClass(object, "STmodel", name = "object")
    if (!is.null(STdata)) {
        stCheckClass(STdata, c("STdata", "STmodel"), name = "STdata")
    }
    type <- tolower(type)
    SpatioTemporal:::stCheckType(type)
    if (inherits(x, "estimateSTmodel")) {
        x <- coef(x, "all")$par
    }
    dimensions <- loglikeSTdim3(object)
    if (type == "f" && length(x) != dimensions$nparam) {
        stop(paste("type=f, requires", dimensions$nparam, "parameters but length(x) =", 
            length(x)))
    }
    if (type != "f" && length(x) == dimensions$nparam) {
        x <- x[(dimensions$nparam - dimensions$nparam.cov + 1):dimensions$nparam]
    }
    if (type != "f" && length(x) != dimensions$nparam.cov) {
        stop(paste("type!=f, requires", dimensions$nparam.cov, 
            "parameters but length(x) =", length(x)))
    }
    if (only.pars && type == "f") {
        warning("only.pars=TRUE and type=(f)ull only returns KNOWN parameters.", 
            immediate. = TRUE)
    }
    if (only.pars) {
        only.obs <- FALSE
        combine.data <- FALSE
        pred.covar <- FALSE
    }
    if (only.obs && is.null(STdata)) {
        stop("only.obs=TRUE requires STdata.")
    }
    if (combine.data && is.null(STdata)) {
        warning("No data to combine with; predicting for 'object'", 
            immediate. = TRUE)
        combine.data <- FALSE
    }
    if (only.obs && combine.data) {
        warning("only.obs=TRUE implies combine.data=FALSE.", 
            immediate. = TRUE)
        combine.data <- FALSE
    }
    if (pred.covar && !pred.var) {
        warning("pred.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
        pred.var <- TRUE
    }
    if (pred.covar && only.obs) {
        warning("only.obs=TRUE implies pred.covar=FALSE.", immediate. = TRUE)
        pred.covar <- FALSE
    }
    if (beta.covar && !pred.var) {
        warning("beta.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
        pred.var <- TRUE
    }
    if (is.null(STdata)) {
        STdata <- object
    } else {
        if (combine.data) {
            STdata <- c(object, STdata)
        } else if (!inherits(STdata, "STmodel")) {
            if (is.null(STdata$trend)) {
                warning("STdata lacking a trend element, using trend from object.", 
                  immediate. = TRUE)
                STdata$trend <- object$trend
            }
            STdata <- createSTmodel3(STdata, LUR = object$LUR.list, 
                ST = object$ST.list, cov.beta = object$cov.beta, 
                cov.nu = object$cov.nu, locations = object$locations.list, 
                scale = !is.null(object$scale.covars), scale.covars = object$scale.covars)
        } else {
            SpatioTemporal:::areSTmodelsConsistent(object, STdata, "STdata")
        }
    }
    
    if(any(object$cov.beta$covf == "tprs")){
    	m=2
    	object$cov.beta$K <- object$cov.beta$K - choose(m+2-1, 2)*(object$cov.beta$covf=="tprs")
    	STdata$cov.beta$K <- STdata$cov.beta$K - choose(m+2-1, 2)*(STdata$cov.beta$covf=="tprs")
    }

    tmp <- loglikeSTgetPars3(x, object)
    if (type == "f") {
        gamma <- tmp$gamma
        alpha <- tmp$alpha
    }
    cov.pars.beta <- tmp$cov.beta
    cov.pars.nu <- tmp$cov.nu
    nugget.unobs <- SpatioTemporal:::internalFixNuggetUnobs(nugget.unobs, STdata, 
        cov.pars.nu$nugget)
    Xtilde <- calc.FX(object$F, object$LUR, object$obs$idx)
    if (dimensions$L != 0) {
        Xtilde <- cbind(object$ST, Xtilde)
    }
    i.sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = object$D.nu, 
        type = object$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = object$nt, 
        ind1 = object$obs$idx)
    i.sigma.nu <- makeCholBlock(i.sigma.nu, block.sizes = object$nt)
    i.sigma.nu <- invCholBlock(i.sigma.nu, block.sizes = object$nt)

    tF.iS.F <- calc.tFXF(object$F, i.sigma.nu, object$obs$idx, 
        block.sizes = object$nt, n.loc = dimensions$n.obs)
        
	if(all(object$cov.beta$covf != "iid")){ #spatial smoothing
		i.sigma.B <- kronecker(diag(unlist(lapply(cov.pars.beta$pars, function(y) return(y['sill'])))),
								diag(object$cov.beta$K[1]))
		i.sigma.B <- try(makeCholBlock(i.sigma.B, n.blocks = dimensions$m), 
			silent = TRUE)
		if (class(i.sigma.B) == "try-error") {
			return(-.Machine$double.xmax)
		}
		i.sigma.B <- invCholBlock(i.sigma.B, n.blocks = dimensions$m)
		
		ZB <- makeZB(cov.pars.beta$pars, coords1=object$locations[, c("x.beta", "y.beta")],
		coordsK=object$locations.list$knot.coords[object$cov.beta$knots[[1]],], 
		type = object$cov.beta$covf)
		
		if(all(object$cov.beta$nugget)){ # with nugget
			i.sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)	
			i.sigma.P.Y <- tF.iS.F + diag(1/i.sigma.P^2)
			i.sigma.P.Y <- makeCholBlock(i.sigma.P.Y, n.blocks=1)
		
			tF.i.sR <- calc.tFX(object$F, i.sigma.nu, object$obs$idx, n.loc=dimensions$n.obs)	
			i.sP.Y.tF.i.sR <- solveTriBlock(i.sigma.P.Y, tF.i.sR, transpose = TRUE)
			i.sP.Y.tF.iS.F <- solveTriBlock(i.sigma.P.Y, tF.iS.F, transpose = TRUE)
		
			i.sigma.nu <- - t(i.sP.Y.tF.i.sR) %*% i.sP.Y.tF.i.sR + i.sigma.nu
			tF.iS.F <- - t(i.sP.Y.tF.iS.F)%*%i.sP.Y.tF.iS.F + tF.iS.F			
		}
		
		R.i.sigma.B.Y <- blockMult2(ZB, tF.iS.F, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
		R.i.sigma.B.Y <- blockMult2(ZB, t(R.i.sigma.B.Y), n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
		
		R.i.sigma.B.Y <- R.i.sigma.B.Y + i.sigma.B
		R.i.sigma.B.Y <- makeCholBlock(R.i.sigma.B.Y)
		if (type == "f") {
			gamma.E <- c(tmp$gamma)
			alpha.E <- unlist(tmp$alpha)
			gamma.V <- matrix(0, length(gamma.E), length(gamma.E))
			alpha.V <- matrix(0, length(alpha.E), length(alpha.E))
			gamma.alpha.C <- matrix(0, length(gamma.E), length(alpha.E))
			gamma.alpha <- as.matrix(c(gamma.E, alpha.E))
		}   else {
#			iS.nu.X <- blockMult(i.sigma.nu, Xtilde, block.sizes = object$nt)
			iS.nu.X <- i.sigma.nu %*% Xtilde
			tF.iS.nu.X <- calc.tFX(object$F, iS.nu.X, object$obs$idx, 
				n.loc = dimensions$n.obs)
			tZB.tF.iS.nu.X <- blockMult2(ZB, tF.iS.nu.X, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
			iSBY.tZB.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, tZB.tF.iS.nu.X, 
				transpose = TRUE)
			iSBY.tZB.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, iSBY.tZB.tF.iS.X, 
				transpose = FALSE)
			i.XSX <- t(Xtilde) %*% iS.nu.X
			i.XSX <- i.XSX - t(tZB.tF.iS.nu.X) %*% iSBY.tZB.tF.iS.X
#			iS.nu.Y <- blockMult(i.sigma.nu, object$obs$obs, block.sizes = object$nt)
			iS.nu.Y <- i.sigma.nu %*% object$obs$obs
			tF.iS.nu.Y <- calc.tFX(object$F, iS.nu.Y, object$obs$idx, 
				n.loc = dimensions$n.obs)
			tZB.tF.iS.nu.Y <- blockMult2(ZB, tF.iS.nu.Y, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
			tmp2 <- object$obs$obs %*% iS.nu.X
			tmp2 <- tmp2 - t(tZB.tF.iS.nu.Y) %*% iSBY.tZB.tF.iS.X
		}
    } else { # No Spatial Smoothing, only nuggets
        iS.nu.X <- blockMult(i.sigma.nu, Xtilde, block.sizes = object$nt)
        tF.iS.nu.X <- calc.tFX(object$F, iS.nu.X, object$obs$idx, 
            n.loc = dimensions$n.obs)
		iS.nu.Y <- blockMult(i.sigma.nu, object$obs$obs, block.sizes = object$nt)
		tF.iS.nu.Y <- calc.tFX(object$F, iS.nu.Y, object$obs$idx, 
			n.loc = dimensions$n.obs)

		i.sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)	
		i.sigma.P.Y <- tF.iS.F + diag(1/i.sigma.P^2)
		i.sigma.P.Y <- makeCholBlock(i.sigma.P.Y, n.blocks=1)

		i.sP.Y.tF.iS.nu.X <- solveTriBlock(i.sigma.P.Y, tF.iS.nu.X, transpose = TRUE)
		i.sP.Y.tF.iS.nu.X <- solveTriBlock(i.sigma.P.Y, i.sP.Y.tF.iS.nu.X, transpose = FALSE)
		
		i.XSX <- t(Xtilde) %*% iS.nu.X
		i.XSX <- i.XSX - t(tF.iS.nu.X) %*% i.sP.Y.tF.iS.nu.X

        tmp2 <- object$obs$obs %*% iS.nu.X
        tmp2 <- tmp2 - t(tF.iS.nu.Y) %*% i.sP.Y.tF.iS.nu.X
    }
	gamma.alpha <- solve(i.XSX, t(tmp2))
	if (dimensions$L != 0) {
		gamma.E <- c(gamma.alpha[1:dimensions$L])
		alpha.E <- c(gamma.alpha[-(1:dimensions$L)])
	}    else {
		gamma.E <- double(0)
		alpha.E <- c(gamma.alpha)
	}
	i.XSX <- solve(i.XSX)
	if (dimensions$L != 0) {
		gamma.V <- i.XSX[1:dimensions$L, 1:dimensions$L, 
			drop = FALSE]
		alpha.V <- i.XSX[-c(1:dimensions$L), -c(1:dimensions$L), 
			drop = FALSE]
		gamma.alpha.C <- i.XSX[1:dimensions$L, -c(1:dimensions$L), 
			drop = FALSE]
	}      else {
		gamma.V <- matrix(0, 0, 0)
		alpha.V <- i.XSX
		gamma.alpha.C <- matrix(0, length(gamma.E), length(alpha.E))
	}

    gamma.E <- as.matrix(gamma.E)
    alpha.E <- as.matrix(alpha.E)
    names.tmp <- loglikeSTnames3(object, TRUE)
    names.tmp <- names.tmp[1:(dimensions$nparam - dimensions$nparam.cov)]
    if (dimensions$L != 0) {
        rownames(gamma.E) <- names.tmp[1:dimensions$L]
        colnames(gamma.V) <- rownames(gamma.V) <- rownames(gamma.E)
        rownames(gamma.alpha.C) <- rownames(gamma.E)
        rownames(alpha.E) <- names.tmp[-(1:dimensions$L)]
    }   else {
        rownames(alpha.E) <- names.tmp
    }
    colnames(alpha.V) <- rownames(alpha.V) <- rownames(alpha.E)
    colnames(gamma.alpha.C) <- rownames(alpha.E)
    out <- list()
    class(out) <- "predictobject"
    out$opts <- list(only.pars = only.pars, nugget.unobs = nugget.unobs, 
        only.obs = only.obs, pred.var = pred.var, pred.covar = pred.covar, 
        beta.covar = beta.covar, combine.data = combine.data, 
        type = type)
    out$pars <- list(gamma.E = gamma.E, alpha.E = alpha.E, gamma.V = gamma.V, 
        alpha.V = alpha.V, gamma.alpha.C = gamma.alpha.C)
    if (out$opts$only.pars) {
        return(out)
    }
    rm(gamma.E, alpha.E, gamma.V, alpha.V, gamma.alpha.C, names.tmp)
    rm(only.pars, nugget.unobs, only.obs, pred.var, pred.covar, 
        beta.covar, combine.data, type)
    C.minus.mu <- object$obs$obs - (Xtilde %*% gamma.alpha)
# 	iS.nu.C <- blockMult(i.sigma.nu, C.minus.mu, block.sizes = object$nt)
	iS.nu.C <- i.sigma.nu %*% C.minus.mu
	tF.iS.nu.C <- calc.tFX(object$F, iS.nu.C, object$obs$idx, 
		n.loc = dimensions$n.obs)
	tF.iS <- calc.tFX(object$F, i.sigma.nu, object$obs$idx, n.loc = dimensions$n.obs)
	if(all(object$cov.beta$covf != "iid")){ #spatial smoothing
		tZB.tF.iS.nu.C <- blockMult2(ZB, tF.iS.nu.C, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
		iSBY.tZB.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, tZB.tF.iS.nu.C, 
			transpose = TRUE)
		iSBY.tZB.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, iSBY.tZB.tF.iS.C, 
			transpose = FALSE)
		tZB.tF.iS <- blockMult2(ZB, tF.iS, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
		obs.mu <- iS.nu.C - t(tZB.tF.iS) %*% iSBY.tZB.tF.iS.C
		rm(C.minus.mu, iS.nu.C, tZB.tF.iS.nu.C, iSBY.tZB.tF.iS.C)
	} else {
		i.sP.Y.tF.iS.C <- solveTriBlock(i.sigma.P.Y, tF.iS.nu.C, 
			transpose = TRUE)
		i.sP.Y.tF.iS.C <- solveTriBlock(i.sigma.P.Y, i.sP.Y.tF.iS.C, 
			transpose = FALSE)
		obs.mu <- iS.nu.C - t(tF.iS) %*% i.sP.Y.tF.iS.C
	}
    if (out$opts$only.obs) {
        idx.unobs <- STdata$obs$idx
        T1 <- STdata$obs$date
        F <- STdata$F
        if (dimensions$L != 0) {
            ST.unobs <- STdata$ST
        }
        N.unobs <- length(unique(idx.unobs))
        T.unobs <- length(unique(T1))
    }   else {
        N.unobs <- dim(STdata$locations)[1]
        T.unobs <- dim(STdata$trend)[1]
        idx.unobs <- rep(1:N.unobs, each = T.unobs)
        T1 <- rep(STdata$trend$date, N.unobs)
        F <- matrix(1, (T.unobs * N.unobs), dimensions$m)
        for (i in (1:dimensions$m)) {
            if (colnames(STdata$F)[i] != "const") {
                F[, i] <- rep(STdata$trend[, colnames(STdata$F)[i]], 
                  N.unobs)
            }
        }
        if (dimensions$L != 0) {
            ST.unobs <- matrix(STdata$ST.all, (T.unobs * N.unobs), 
                dimensions$L)
        }
        if (out$opts$pred.covar) {
            Nmax <- T.unobs
        }
    }
    date.all <- sort(unique(c(object$trend$date, STdata$trend$date)))
    nt.unobs <- nt.obs <- double(length(date.all))
    for (i in c(1:length(date.all))) {
        nt.obs[i] <- sum(object$obs$date == date.all[i])
    }
    Xtilde.unobs <- calc.FX(F, STdata$LUR.all, idx.unobs)
    if (dimensions$L != 0) {
        Xtilde.unobs <- cbind(ST.unobs, Xtilde.unobs)
    }
    out$EX <- as.matrix(Xtilde.unobs %*% gamma.alpha)
    if (!out$opts$only.obs) {
        dim(out$EX) <- c(T.unobs, N.unobs)
    }
    out$EX.mu <- out$EX
    if (out$opts$pred.var) {
        out$VX <- matrix(NA, dim(out$EX)[1], dim(out$EX)[2])
        out$VX.pred <- matrix(NA, dim(out$EX)[1], dim(out$EX)[2])
    }
    if (out$opts$pred.covar) {
        out$VX.full <- list()
    }
    alpha <- vector("list", dimensions$m)
    offset <- 1
    for (i in 1:length(alpha)) {
        alpha[[i]] <- out$pars$alpha.E[offset:sum(dimensions$p[1:i])]
        offset <- sum(dimensions$p[1:i]) + 1
    }
    out$beta <- list()
    out$beta$mu <- calc.mu.B(STdata$LUR.all, alpha)
    rm(alpha, gamma.alpha)
    loc.unobs.nu <- STdata$locations[, c("x.nu", "y.nu"), drop = FALSE]
    loc.unobs.beta <- STdata$locations[, c("x.beta", "y.beta"), 
        drop = FALSE]
    I.obs <- match(colnames(object$D.nu), object$locations$ID)
    loc.obs.nu <- object$locations[I.obs, c("x.nu", "y.nu"), 
        drop = FALSE]
    loc.obs.beta <- object$locations[I.obs, c("x.beta", "y.beta"), 
        drop = FALSE]
    Ind.2.1 <- match(object$locations$ID, STdata$locations$ID, 
        nomatch = 0)
    loc.knots.beta <- object$locations.list$knot.coords[object$cov.beta$knots[[1]], ]
# 	sigma.B.C <- makeSigmaB(cov.pars.beta$pars, dist = crossDist(loc.unobs.beta, 
#         loc.obs.beta), type = object$cov.beta$covf, nugget = cov.pars.beta$nugget, 
#         ind2.to.1 = Ind.2.1)
 	sigma.B.C <- makeSigmaB3(cov.pars.beta$pars, 
	coords1=loc.unobs.beta, coordsK=loc.knots.beta, coords2=loc.obs.beta, basis.coords=loc.obs.beta,
	type = object$cov.beta$covf, nugget = cov.pars.beta$nugget, ind2.to.1 = Ind.2.1, fit.tprs=FALSE)
 
    out$beta$EX <- (c(out$beta$mu) + (sigma.B.C %*% calc.tFX(object$F, 
        obs.mu, object$obs$idx, n.loc = dimensions$n.obs)))
    dim(out$beta$EX) <- dim(out$beta$mu)
    dimnames(out$beta$EX) <- dimnames(out$beta$mu)
    if (out$opts$pred.var) {
    print("Sorry. pred.var not implemented for low rank models without nugget")
#         Sby.iSb.Sou <- i.sigma.B %*% t(sigma.B.C)
#         Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, 
#             transpose = TRUE)
#         Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, 
#             transpose = FALSE)
# 		if(any(object$cov.beta$covf == "tprs")){
# 			print("sorry, beta.covar not implemented with tprs smooth")
# 			out$beta$VX <- 0*out$beta$EX
# 		} else {
#         	if (out$opts$beta.covar) {
# # 			sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = crossDist(loc.unobs.beta), 
# # 					type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
# 				sigma.B.uu <- makeSigmaB3(cov.pars.beta$pars, 
# 					coords1=loc.unobs.beta, coordsK=loc.knots.beta, basis.coords=loc.obs.beta,
# 					type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
# 				tmp <- sigma.B.uu - sigma.B.C %*% tF.iS.F %*% Sby.iSb.Sou
# 				out$beta$VX.full <- list()
# 				for (i in 1:dim(out$beta$EX)[2]) {
# 					Ind <- (1:dim(out$beta$EX)[1]) + (i - 1) * dim(out$beta$EX)[1]
# 					out$beta$VX.full[[i]] <- tmp[Ind, Ind, drop = FALSE]
# 					rownames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
# 					colnames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
# 				}
# 				names(out$beta$VX.full) <- colnames(out$beta$EX)
# 				tmp <- diag(tmp)
#         	} else {
# 			   sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = matrix(0, 
# 					 1, 1), type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
# 				sigma.B.uu <- matrix(diag(sigma.B.uu), ncol = dim(sigma.B.uu)[1], 
# 					nrow = dim(loc.unobs.beta)[1], byrow = TRUE)
# 				tmp <- c(sigma.B.uu) - rowSums(sigma.B.C * t(tF.iS.F %*% 
# 					Sby.iSb.Sou))
# 			}
# 			out$beta$VX <- matrix(tmp, ncol = dim(out$beta$EX)[2])
# 			dimnames(out$beta$VX) <- dimnames(out$beta$EX)
# 			rm(tmp, Sby.iSb.Sou, tF.iS.F, sigma.B.uu)
#         }
    }
    cross.D.nu <- crossDist(loc.unobs.nu, loc.obs.nu)
    for (i in 1:ceiling(length(out$EX)/Nmax)) {
        Ind <- c((1 + (i - 1) * Nmax):min(i * Nmax, length(out$EX)))
        T1.Ind <- T1[Ind]
        for (j in c(1:length(date.all))) {
            nt.unobs[j] <- sum(T1.Ind == date.all[j])
        }
        T1.order <- order(T1.Ind)
        sigma.B.full.C <- calc.FXtF2(F[Ind, , drop = FALSE], 
            sigma.B.C, loc.ind = idx.unobs[Ind], F2 = object$F, 
            loc.ind2 = object$obs$idx)
        sigma.nu.C <- makeSigmaNu(cov.pars.nu$pars, dist = cross.D.nu, 
            type = object$cov.nu$covf, nugget = 0, random.effect = cov.pars.nu$random.effect, 
            ind1 = (idx.unobs[Ind])[T1.order], ind2 = object$obs$idx, 
            blocks1 = nt.unobs, blocks2 = nt.obs, ind2.to.1 = Ind.2.1)
        sigma.nu.C[T1.order, ] <- sigma.nu.C
        sigma.nu.C <- sigma.nu.C + sigma.B.full.C
        out$EX[Ind] <- out$EX.mu[Ind] + sigma.nu.C %*% obs.mu
        if (out$opts$pred.var) {
	    	print("Sorry, pred.var not implemented for low rank models with no nugget")

  #           I.loc <- sort(unique(idx.unobs[Ind]))
#             unobs.D.beta <- crossDist(loc.unobs.beta[I.loc, , 
#                 drop = FALSE])
# #             sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = unobs.D.beta, 
# #                 type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
#             sigma.B.uu <- makeSigmaB3(cov.pars.beta$pars, 
#             				coords1=loc.unobs.beta[I.loc, , drop=FALSE], 
#             				coordsK=loc.knots.beta,
#             				basis.coords=loc.obs.beta,
#                 			type = object$cov.beta$covf, 
#                 			nugget = cov.pars.beta$nugget, fit.tprs=FALSE)
#             V.uu <- calc.FXtF2(F[Ind, , drop = FALSE], sigma.B.uu, 
#                 loc.ind = idx.unobs[Ind] - min(I.loc) + 1)
#             unobs.D.nu <- crossDist(loc.unobs.nu[I.loc, , drop = FALSE])
#             V.uu <- V.uu + makeSigmaNu(cov.pars.nu$pars, dist = unobs.D.nu, 
#                 type = object$cov.nu$covf, nugget = 0, random.effect = cov.pars.nu$random.effect, 
#                 ind1 = idx.unobs[Ind] - min(I.loc) + 1, blocks1 = nt.unobs)
#             t.sigma.nu.C <- t(sigma.nu.C)
#             iS.Sou <- blockMult(i.sigma.nu, t.sigma.nu.C, block.sizes = object$nt)
#             tF.iS.Sou <- tF.iS %*% t.sigma.nu.C
#             Sby.tF.iS.Sou <- solveTriBlock(R.i.sigma.B.Y, tF.iS.Sou, 
#                 transpose = TRUE)
#             if (out$opts$pred.covar) {
#                 tmp <- -sigma.nu.C %*% iS.Sou + t(Sby.tF.iS.Sou) %*% 
#                   Sby.tF.iS.Sou
#                 V.cond <- V.cond.0 <- V.uu + tmp
#                 diag(V.cond) <- diag(V.cond) + out$opts$nugget.unobs[idx.unobs[Ind]]
#             }       else {
#                 tmp <- (-colSums(t(sigma.nu.C) * iS.Sou) + colSums(Sby.tF.iS.Sou * 
#                   Sby.tF.iS.Sou))
#                 V.cond <- V.cond.0 <- diag(V.uu) + tmp
#                 V.cond <- V.cond + out$opts$nugget.unobs[idx.unobs[Ind]]
#             }
#             if (out$opts$type == "r") {
#                 stop("NOT YET IMPLEMENTED. Use type=\"p\" instead.")
#             }
#             if (out$opts$pred.covar) {
#                 out$VX[Ind] <- diag(V.cond.0)
#                 out$VX.pred[Ind] <- diag(V.cond)
#                 out$VX.full[[i]] <- V.cond.0
#             }      else {
#                 out$VX[Ind] <- V.cond.0
#                 out$VX.pred[Ind] <- V.cond
#             }
         }
    }
    if (out$opts$pred.var) {
    	print("Sorry, pred.var not implemented for low rank models with no nugget")
#         out$VX <- pmax(out$VX, 0)
#         out$VX.pred <- pmax(out$VX.pred, 0)
    }
    EX.beta.list <- list()
    for (i in 1:dim(out$beta$EX)[2]) {
        EX.beta.list[[i]] <- out$beta$EX[, i, drop = FALSE]
    }
    out$EX.mu.beta <- calc.FX(F, EX.beta.list, idx.unobs)
    out$EX.mu.beta <- matrix(rowSums(out$EX.mu.beta))
    if (dimensions$L != 0) {
        out$EX.mu.beta <- out$EX.mu.beta + (ST.unobs %*% out$pars$gamma.E)
    }
    if (!out$opts$only.obs) {
        dim(out$EX.mu.beta) <- c(T.unobs, N.unobs)
    }
    if (!out$opts$only.obs) {
        colnames(out$EX) <- STdata$locations$ID
        rownames(out$EX) <- as.character(STdata$trend$date)
        dimnames(out$EX.mu.beta) <- dimnames(out$EX.mu) <- dimnames(out$EX)
        if (out$opts$pred.var) {
#             dimnames(out$VX) <- dimnames(out$VX.pred) <- dimnames(out$EX)
        }
        if (out$opts$pred.covar) {
#             names(out$VX.full) <- colnames(out$EX)
#             for (i in 1:length(out$VX.full)) colnames(out$VX.full[[i]]) <- rownames(out$VX.full[[i]]) <- rownames(out$EX)
        }
    }
    if (length(STdata$obs$obs) != 0) {
        if (out$opts$only.obs) {
            I <- 1:length(out$EX)
        }      else {
            I <- (match(STdata$obs$ID, colnames(out$EX)) - 1) * 
                dim(out$EX)[1] + match(STdata$obs$date, STdata$trend$date)
        }
        out$I <- data.frame(I = I, date = STdata$obs$date, ID = STdata$obs$ID, 
            stringsAsFactors = FALSE)
    }
    return(out)
}

predictCV.STmodel_opt <- function (object, x, Ind.cv, ..., silent = TRUE) 
{
    stCheckClass(object, "STmodel", name = "object")
    Ind.cv <- SpatioTemporal:::stCheckInternalCV(Ind.cv)
    Ind.cv <- as.matrix(Ind.cv)
    if (dim(Ind.cv)[2] == 1) {
        N.CV.sets <- max(Ind.cv, na.rm = TRUE)
    }else {
        stop("Some observation(s) are left out in several Cv-groups.")
    }
    if (inherits(x, "estCVSTmodel")) {
        x <- coef(x, "all")
    } else if (inherits(x, "estimateSTmodel")) {
        x <- coef(x, "all")$par
    }
    x <- as.matrix(x)
    if (dim(x)[2] == 1) {
        x <- matrix(x, length(x), N.CV.sets)
    } else if (dim(x)[2] != N.CV.sets) {
        stop("Number of parameters does not match the number of cv-sets.")
    }
    pred <- pred_vars <- list()
    for (i in 1:N.CV.sets) {
        if (!silent) 
            message(sprintf("Predicting cv-set %d/%d", i, N.CV.sets))
        if (dim(Ind.cv)[2] == 1) {
            Ind.current <- Ind.cv == i
        }  else {
            Ind.current <- as.logical(Ind.cv[, i])
        }
        object.obs <- dropObservations3(object, Ind.current)
        object.pred <- dropObservations3(object, !Ind.current)
        nugget.unobs <- loglikeSTgetPars3(x[, i], object)$cov.nu$nugget
        nugget.unobs <- nugget.unobs[object.pred$locations$ID, 
            , drop = FALSE]
        pred[[i]] <- predict.STmodel_opt(object.obs, x[, i], STdata = object.pred, 
            nugget.unobs = nugget.unobs, only.pars = FALSE, combine.data = FALSE, 
            ...)
        ### Add in pred_var
        if(pred[[i]]$opts$pred.var){
        	pred_vars[[i]] <- predict.var(object.obs, x[, i], STdata = object.pred, 
            nugget.unobs = nugget.unobs, only.pars = FALSE, combine.data = FALSE, 
            ...)
        }
    }
    out <- list()
    class(out) <- "predCVSTmodel"
    out$opts <- pred[[1]]$opts
    if (!is.null(out$opts$nugget.unobs)) {
        out$opts$nugget.unobs <- unlist(sapply(pred, function(x) {
            x$opts$nugget.unobs
        }))
        names(out$opts$nugget.unobs) <- unlist(sapply(pred, function(x) {
            rownames(x$opts$nugget.unobs)
        }))
    }
    out$Ind.cv <- Ind.cv
    out$pred.obs <- object$obs[, c("obs", "date", "ID"), drop = FALSE]
    out$pred.obs$EX.mu <- NA
    out$pred.obs$EX.mu.beta <- NA
    out$pred.obs$EX <- NA
    if (out$opts$pred.var) {
         out$pred.obs$VX <- NA
         out$pred.obs$VX.pred <- NA
    }
    if(out$opts$pred.covar) {
    	out$VX.full <- list()
    }
    for (i in 1:N.CV.sets) {
        Ind.current <- Ind.cv == i
        I <- pred[[i]]$I$I
        out$pred.obs[Ind.current, c("EX.mu", "EX.mu.beta", "EX")] <- cbind(pred[[i]]$EX.mu[I], 
            pred[[i]]$EX.mu.beta[I], pred[[i]]$EX[I])
        if (out$opts$pred.var) {
             out$pred.obs[Ind.current, c("VX", "VX.pred")] <- cbind(pred_vars[[i]]$VX[I], 
                 pred_vars[[i]]$VX.pred[I])
        }
        if(out$opts$pred.covar) {
        	out$VX.full <- pred_vars[[i]]$VX.full
        }
    }
    out$pred.obs$res <- out$pred.obs$obs - out$pred.obs$EX
    if (out$opts$pred.var) {
         out$pred.obs$res.norm <- out$pred.obs$res/sqrt(out$pred.obs$VX.pred)
    }
    if (out$opts$only.obs) {
        any.duplicates <- any(duplicated(unlist(sapply(pred, 
            function(x) {
                unique(x$I$ID)
            }))))
    }
    else {
        any.duplicates <- any(duplicated(unlist(sapply(pred, 
            function(x) {
                colnames(x$EX)
            }))))
    }
    if (!any.duplicates) {
        out$pred.all <- list()
        out$pred.all$EX.mu <- matrix(NA, dim(object$trend)[1], 
            dim(object$locations)[1])
        colnames(out$pred.all$EX.mu) <- object$locations$ID
        rownames(out$pred.all$EX.mu) <- as.character(object$trend$date)
        out$pred.all$EX <- out$pred.all$EX.mu.beta <- out$pred.all$EX.mu
        if (out$opts$pred.var) {
             out$pred.all$VX <- out$pred.all$VX.pred <- out$pred.all$EX
        }
        out$pred.all$beta <- list()
        out$pred.all$beta$mu <- matrix(NA, dim(object$locations)[1], 
            length(object$LUR))
        colnames(out$pred.all$beta$mu) <- names(object$LUR)
        rownames(out$pred.all$beta$mu) <- object$locations$ID
        out$pred.all$beta$EX <- out$pred.all$beta$mu
        if (out$opts$pred.var) {
             out$pred.all$beta$VX <- out$pred.all$beta$EX
        }
        for (i in 1:N.CV.sets) {
            if (out$opts$only.obs) {
                ID.names <- rownames(pred[[i]]$beta$EX)
                I <- ((match(pred[[i]]$I$ID, object$locations$ID) - 
                  1) * dim(out$pred.all$EX)[1] + match(as.character(pred[[i]]$I$date), 
                  rownames(out$pred.all$EX)))
            }
            else {
                ID.names <- colnames(pred[[i]]$EX)
                I <- (rep(match(ID.names, object$locations$ID) - 
                  1, each = dim(out$pred.all$EX)[1]) * dim(out$pred.all$EX)[1] + 
                  rep(match(rownames(pred[[i]]$EX), rownames(out$pred.all$EX)), 
                    length(ID.names)))
            }
            out$pred.all$EX.mu[I] <- pred[[i]]$EX.mu
            out$pred.all$EX.mu.beta[I] <- pred[[i]]$EX.mu.beta
            out$pred.all$EX[I] <- pred[[i]]$EX
            if (out$opts$pred.var) {
                 out$pred.all$VX[I] <- pred_vars[[i]]$VX
                 out$pred.all$VX.pred[I] <- pred_vars[[i]]$VX.pred
            }
            out$pred.all$beta$mu[ID.names, ] <- pred[[i]]$beta$mu
            out$pred.all$beta$EX[ID.names, ] <- pred[[i]]$beta$EX
            if (out$opts$pred.var) {
#                 out$pred.all$beta$VX[ID.names, ] <- pred[[i]]$beta$VX
            }
        }
    }
    else {
        out$pred.all.by.cv <- vector("list", N.CV.sets)
        for (i in 1:N.CV.sets) {
            if (out$opts$only.obs) {
                out$pred.all.by.cv[[i]] <- list()
                EX.mu <- createDataMatrix(obs = pred[[i]]$EX.mu, 
                  date = pred[[i]]$I$date, ID = pred[[i]]$I$ID)
                EX.mu.beta <- createDataMatrix(obs = pred[[i]]$EX.mu.beta, 
                  date = pred[[i]]$I$date, ID = pred[[i]]$I$ID)
                EX <- createDataMatrix(obs = pred[[i]]$EX, date = pred[[i]]$I$date, 
                  ID = pred[[i]]$I$ID)
                if (out$opts$pred.var) {
                   VX <- createDataMatrix(obs = pred_vars[[i]]$VX, 
                     date = pred[[i]]$I$date, ID = pred[[i]]$I$ID)
                   VX.pred <- createDataMatrix(obs = pred_vars[[i]]$VX.pred, 
                     date = pred[[i]]$I$date, ID = pred[[i]]$I$ID)
                }
                out$pred.all.by.cv[[i]]$EX.mu <- matrix(NA, dim(object$trend)[1], 
                  dim(EX.mu)[2])
                colnames(out$pred.all.by.cv[[i]]$EX.mu) <- colnames(EX.mu)
                rownames(out$pred.all.by.cv[[i]]$EX.mu) <- as.character(object$trend$date)
                out$pred.all.by.cv[[i]]$EX.mu.beta <- out$pred.all.by.cv[[i]]$EX.mu
                out$pred.all.by.cv[[i]]$EX <- out$pred.all.by.cv[[i]]$EX.mu
                if (out$opts$pred.var) {
                   out$pred.all.by.cv[[i]]$VX <- out$pred.all.by.cv[[i]]$EX
                   out$pred.all.by.cv[[i]]$VX.pred <- out$pred.all.by.cv[[i]]$EX
                }
                out$pred.all.by.cv[[i]]$EX.mu[rownames(EX.mu), 
                  ] <- EX.mu
                out$pred.all.by.cv[[i]]$EX.mu.beta[rownames(EX.mu.beta), 
                  ] <- EX.mu.beta
                out$pred.all.by.cv[[i]]$EX[rownames(EX), ] <- EX
                if (out$opts$pred.var) {
                   out$pred.all.by.cv[[i]]$VX[rownames(VX), ] <- VX
                   out$pred.all.by.cv[[i]]$VX.pred[rownames(VX.pred), 
                     ] <- VX.pred
                }
                out$pred.all.by.cv[[i]]$beta <- pred[[i]]$beta
            }
            else {
                if (out$opts$pred.var) {
                   pick.names <- c("EX.mu", "EX.mu.beta", "EX", 
                     "VX", "VX.pred", "beta")
                }
                else {
                  pick.names <- c("EX.mu", "EX.mu.beta", "EX", 
                    "beta")
                }
                out$pred.all.by.cv[[i]] <- pred[[i]][pick.names]
            }
        }
    }
    return(out)
}

calc.STdeterminant <- function (x = NULL, STmodel, x.fixed=NULL, type="p") 
{
    stCheckClass(STmodel, "STmodel", name = "STmodel")
    type <- tolower(type)
    x <- SpatioTemporal:::stCheckLoglikeIn(x, x.fixed, type)
    if (is.null(x) || any(is.na(x))) {
        return(loglikeSTnames3(STmodel, all = (type == "f")))
    }
    dimensions <- loglikeSTdim3(STmodel)
    if ((type == "f" && length(x) != dimensions$nparam) || (type != 
        "f" && length(x) != dimensions$nparam.cov)) {
        stop("Number of parameters, length(x), incorrect.")
    }
    tmp <- loglikeSTgetPars3(x, STmodel)

    cov.pars.beta <- tmp$cov.beta
    cov.pars.nu <- tmp$cov.nu

    if(any(STmodel$cov.beta$nugget)){
    	if(!all(STmodel$cov.beta$nugget)){
	    	stop("Require either all or no beta fields have nuggets present")
    	}
    }
    
    if(!all(STmodel$cov.beta != "iid") & !all(STmodel$cov.beta$nugget)){
    	stop("Require either spatial smoothing of all beta fields OR nuggets in all beta fields (or both). Cannot specify no spatial smoothing without a nugget")
    }
    
    if(any(STmodel$cov.beta$covf == "tprs")){
    	m=2
    	STmodel$cov.beta$K <- STmodel$cov.beta$K - choose(m+2-1, 2)*(STmodel$cov.beta$covf=="tprs")
    }
    if (type == "f") {
		stop("FULL Likelihood (type = 'f') not implemented for low rank versions")
    } else if (type == "r") {
		stop("REML Likelihood (type = 'r') not implemented for low rank versions")
	}
    
    sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = STmodel$D.nu, 
        type = STmodel$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = STmodel$nt, 
        ind1 = STmodel$obs$idx)
    sigma.nu <- try(makeCholBlock(sigma.nu, n.blocks = dimensions$T, 
        block.sizes = STmodel$nt), silent = TRUE)
    if (class(sigma.nu) == "try-error") {
        return(-.Machine$double.xmax)
    }
    l <- - sumLogDiag(sigma.nu)
    sigma.nu <- invCholBlock(sigma.nu, block.sizes = STmodel$nt)

	F.i.sR.F <- calc.tFXF(STmodel$F, sigma.nu, STmodel$obs$idx, 
							block.sizes = STmodel$nt, n.loc = dimensions$n.obs)
    
	if(all(STmodel$cov.beta$covf != "iid")){
		sigma.B <- rep(sqrt(sapply(cov.pars.beta$pars, function(y) return(y['sill']))), 
						each=STmodel$cov.beta$K[1])
		l <- l-sum(log(sigma.B))
		sigma.B <- diag(1/sigma.B^2)
	
		ZB <- makeZB(cov.pars.beta$pars, coords1=STmodel$locations[, c("x.beta", "y.beta")],
		coordsK=STmodel$locations.list$knot.coords[STmodel$cov.beta$knots[[1]],], 
		type = STmodel$cov.beta$covf)

		if(all(STmodel$cov.beta$nugget)){
			sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)
			l <- l - sum(log(sigma.P))
	
			sigma.P.Y <- F.i.sR.F + diag(1/sigma.P^2)
			sigma.P.Y <- try(makeCholBlock(sigma.P.Y, n.blocks=1),silent = TRUE)
			if (class(sigma.P.Y) == "try-error") {
				return(-.Machine$double.xmax)
			}
			l <- l - sumLogDiag(sigma.P.Y)
	
			i.sP.Y.tF.i.sR.F <- solveTriBlock(sigma.P.Y, F.i.sR.F, transpose = TRUE)	
			F.i.sR.F <- -t(i.sP.Y.tF.i.sR.F)%*%i.sP.Y.tF.i.sR.F + F.i.sR.F			
		} 
		sigma.B.Y <- blockMult2(ZB, F.i.sR.F, n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=STmodel$cov.beta$K, transpose=TRUE)
		sigma.B.Y <- blockMult2(ZB, t(sigma.B.Y), n.blocks=dimensions$m,
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=STmodel$cov.beta$K, transpose=TRUE)
		sigma.B.Y <- sigma.B.Y + sigma.B
		sigma.B.Y <- try(makeCholBlock(sigma.B.Y, n.blocks = 1), 
			silent = TRUE)
		if (class(sigma.B.Y) == "try-error") {
			return(-.Machine$double.xmax)
		}
		l <- l - sumLogDiag(sigma.B.Y)
	} else {
		sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)
		l <- l - sum(log(sigma.P))
		sigma.P.Y <- F.i.sR.F + diag(1/sigma.P^2)
		sigma.P.Y <- try(makeCholBlock(sigma.P.Y, n.blocks=1),silent = TRUE)
		if (class(sigma.P.Y) == "try-error") {
			return(-.Machine$double.xmax)
		}
		l <- l - sumLogDiag(sigma.P.Y)	
	}
    l <- as.double(l)
    if (!is.finite(l)) 
        l <- -.Machine$double.xmax
    return(l)
}

calc.STdeterminantNaive <- function (x = NULL, STmodel, x.fixed=NULL, type="p") 
{
    stCheckClass(STmodel, "STmodel", name = "STmodel")
    type <- tolower(type)
    x <- SpatioTemporal:::stCheckLoglikeIn(x, x.fixed, type)
    if (is.null(x) || any(is.na(x))) {
        return(loglikeSTnames3(STmodel, all = (type == "f")))
    }
    dimensions <- loglikeSTdim3(STmodel)
    if ((type == "f" && length(x) != dimensions$nparam) || (type != 
        "f" && length(x) != dimensions$nparam.cov)) {
        stop("Number of parameters, length(x), incorrect.")
    }
    tmp <- loglikeSTgetPars3(x, STmodel)

    cov.pars.beta <- tmp$cov.beta
    cov.pars.nu <- tmp$cov.nu

    sigma.B.full <- makeSigmaB3(cov.pars.beta$pars, coords1=STmodel$locations[, c("x.beta", "y.beta")],
	coordsK=STmodel$locations.list$knot.coords[STmodel$cov.beta$knots[[1]],], 
	type = STmodel$cov.beta$covf, nugget = cov.pars.beta$nugget)

    sigma.B.full <- calc.FXtF2(STmodel$F, sigma.B.full, loc.ind = STmodel$obs$idx)
    sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = STmodel$D.nu, 
        type = STmodel$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = STmodel$nt, 
        ind1 = STmodel$obs$idx)
    sigma.nu <- sigma.nu + sigma.B.full
    sigma.nu <- try(makeCholBlock(sigma.nu, n.blocks = 1), silent = TRUE)
    if (class(sigma.nu) == "try-error") {
        return(-.Machine$double.xmax)
    }
    l <- -sumLogDiag(sigma.nu)
    
    l <- as.double(l)
    if (!is.finite(l)) 
        l <- -.Machine$double.xmax
    return(l)

}

##################
##################
#### New functions for AOAS Revision
#### Add a function to calculate prediction variances

predict.var <- function (object, x, STdata = NULL, Nmax = 1000, only.pars = FALSE, 
    nugget.unobs = 0, only.obs = FALSE, pred.var = TRUE, pred.covar = FALSE, 
    beta.covar = FALSE, combine.data = FALSE, type = "p", ...) {
 
    stCheckClass(object, "STmodel", name = "object")
    if (!is.null(STdata)) {
        stCheckClass(STdata, c("STdata", "STmodel"), name = "STdata")
    }
    type <- tolower(type)
    SpatioTemporal:::stCheckType(type)
    if (inherits(x, "estimateSTmodel")) {
        x <- coef(x, "all")$par
    }
    dimensions <- loglikeSTdim3(object)
    if (type == "f" && length(x) != dimensions$nparam) {
        stop(paste("type=f, requires", dimensions$nparam, "parameters but length(x) =", 
            length(x)))
    }
    if (type != "f" && length(x) == dimensions$nparam) {
        x <- x[(dimensions$nparam - dimensions$nparam.cov + 1):dimensions$nparam]
    }
    if (type != "f" && length(x) != dimensions$nparam.cov) {
        stop(paste("type!=f, requires", dimensions$nparam.cov, 
            "parameters but length(x) =", length(x)))
    }
    if (only.pars && type == "f") {
        warning("only.pars=TRUE and type=(f)ull only returns KNOWN parameters.", 
            immediate. = TRUE)
    }
    if (only.pars) {
        only.obs <- FALSE
        combine.data <- FALSE
        pred.covar <- FALSE
    }
    if (only.obs && is.null(STdata)) {
        stop("only.obs=TRUE requires STdata.")
    }
    if (combine.data && is.null(STdata)) {
        warning("No data to combine with; predicting for 'object'", 
            immediate. = TRUE)
        combine.data <- FALSE
    }
    if (only.obs && combine.data) {
        warning("only.obs=TRUE implies combine.data=FALSE.", 
            immediate. = TRUE)
        combine.data <- FALSE
    }
    if (pred.covar && !pred.var) {
        warning("pred.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
        pred.var <- TRUE
    }
    if (pred.covar && only.obs) {
        warning("only.obs=TRUE implies pred.covar=FALSE.", immediate. = TRUE)
        pred.covar <- FALSE
    }
    if (beta.covar && !pred.var) {
        warning("beta.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
        pred.var <- TRUE
    }
    if (is.null(STdata)) {
        STdata <- object
    } else {
        if (combine.data) {
            STdata <- c(object, STdata)
        } else if (!inherits(STdata, "STmodel")) {
            if (is.null(STdata$trend)) {
                warning("STdata lacking a trend element, using trend from object.", 
                  immediate. = TRUE)
                STdata$trend <- object$trend
            }
            STdata <- createSTmodel3(STdata, LUR = object$LUR.list, 
                ST = object$ST.list, cov.beta = object$cov.beta, 
                cov.nu = object$cov.nu, locations = object$locations.list, 
                scale = !is.null(object$scale.covars), scale.covars = object$scale.covars)
        } else {
            SpatioTemporal:::areSTmodelsConsistent(object, STdata, "STdata")
        }
    }
    
    if(any(object$cov.beta$covf == "tprs")){
    	m=2
    	object$cov.beta$K <- object$cov.beta$K - choose(m+2-1, 2)*(object$cov.beta$covf=="tprs")
    	STdata$cov.beta$K <- STdata$cov.beta$K - choose(m+2-1, 2)*(STdata$cov.beta$covf=="tprs")
    }

    tmp <- loglikeSTgetPars3(x, object)
    if (type == "f") {
        gamma <- tmp$gamma
        alpha <- tmp$alpha
    }

    Xtilde <- calc.FX(object$F, object$LUR, object$obs$idx)
	
    cov.pars.beta <- tmp$cov.beta
    cov.pars.nu <- tmp$cov.nu
    nugget.unobs <- SpatioTemporal:::internalFixNuggetUnobs(nugget.unobs, STdata, 
        cov.pars.nu$nugget)
    i.sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = object$D.nu, 
        type = object$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = object$nt, 
        ind1 = object$obs$idx)
    i.sigma.nu <- makeCholBlock(i.sigma.nu, block.sizes = object$nt)
    i.sigma.nu <- invCholBlock(i.sigma.nu, block.sizes = object$nt)

    tF.iS.F <- calc.tFXF(object$F, i.sigma.nu, object$obs$idx, 
        block.sizes = object$nt, n.loc = dimensions$n.obs)
 
	if(all(object$cov.beta$covf != "iid")){ #spatial smoothing
		i.sigma.B <- kronecker(diag(unlist(lapply(cov.pars.beta$pars, function(y) return(y['sill'])))),
								diag(object$cov.beta$K[1]))
		i.sigma.B <- try(makeCholBlock(i.sigma.B, n.blocks = dimensions$m), 
			silent = TRUE)
		if (class(i.sigma.B) == "try-error") {
			return(-.Machine$double.xmax)
		}
		i.sigma.B <- invCholBlock(i.sigma.B, n.blocks = dimensions$m)
		
		ZB <- makeZB(cov.pars.beta$pars, 
		coords1=object$locations[, c("x.beta", "y.beta")],
		coordsK=object$locations.list$knot.coords[object$cov.beta$knots[[1]],], 
		type = object$cov.beta$covf)
		
		if(all(object$cov.beta$nugget)){ # with nugget
			i.sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)	
			i.sigma.P.Y <- tF.iS.F + diag(1/i.sigma.P^2)
			i.sigma.P.Y <- makeCholBlock(i.sigma.P.Y, n.blocks=1)
		
			tF.i.sR <- calc.tFX(object$F, i.sigma.nu, object$obs$idx, n.loc=dimensions$n.obs)	
			i.sP.Y.tF.i.sR <- solveTriBlock(i.sigma.P.Y, tF.i.sR, transpose = TRUE)
			i.sP.Y.tF.iS.F <- solveTriBlock(i.sigma.P.Y, tF.iS.F, transpose = TRUE)
		
			i.sigma.nu <- - t(i.sP.Y.tF.i.sR) %*% i.sP.Y.tF.i.sR + i.sigma.nu
			tF.iS.F <- - t(i.sP.Y.tF.iS.F)%*%i.sP.Y.tF.iS.F + tF.iS.F			
		}
		
		R.i.sigma.B.Y <- blockMult2(ZB, tF.iS.F, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
		R.i.sigma.B.Y <- blockMult2(ZB, t(R.i.sigma.B.Y), n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
		
		R.i.sigma.B.Y <- R.i.sigma.B.Y + i.sigma.B
		R.i.sigma.B.Y <- makeCholBlock(R.i.sigma.B.Y)

		# This is i.sigma.b + tZb.tF.i.sigma.P.Y.F.Zb (cholesky decomp)
		if (type == "f") {
			print("type = 'f' not implemented")
			stop
		}   else {
			
			t.F.iS.nu <- calc.tFX(object$F, i.sigma.nu, object$obs$idx, 
				n.loc = dimensions$n.obs)
			tZB.tF.iS.nu <- blockMult2(ZB, t.F.iS.nu, n.blocks=dimensions$m, 
									block.sizes1=rep(dimensions$n, dimensions$m), 
									block.sizes2=object$cov.beta$K, transpose=TRUE)
			iSBY.tZB.tF.iS <- solveTriBlock(R.i.sigma.B.Y, tZB.tF.iS.nu, 
				transpose = TRUE)
			iSBY.tZB.tF.iS <- t(iSBY.tZB.tF.iS)%*%iSBY.tZB.tF.iS				
		}
		
	    iS.X <- i.sigma.nu%*%Xtilde
	    i.XSX <- as.matrix(t(Xtilde) %*% iS.X)
		i.XSX <- i.XSX - t(Xtilde) %*% iSBY.tZB.tF.iS %*% Xtilde 
		i.XSX <- solve(i.XSX)

		iSoo.Xtilde <- as.matrix(iS.X - iSBY.tZB.tF.iS %*% Xtilde )
		
    } else { # No Spatial Smoothing, only nuggets

		i.sigma.P <- rep(sqrt(cov.pars.beta$nugget), each=dimensions$n)	
		i.sigma.P.Y <- tF.iS.F + diag(1/i.sigma.P^2)
		i.sigma.P.Y <- makeCholBlock(i.sigma.P.Y, n.blocks=1)
		t.F.iS.nu <- calc.tFX(object$F, i.sigma.nu, object$obs$idx, 
				n.loc = dimensions$n.obs)
		iSBY.tF.iS <- solveTriBlock(i.sigma.P.Y, tF.iS.nu, 
				transpose = TRUE)
		iSBY.tF.iS <- t(iSBY.tF.iS)%*%iSBY.tF.iS				

	    iS.X <- i.sigma.nu%*%Xtilde
	    i.XSX <- as.matrix(t(Xtilde) %*% iS.X)
		i.XSX <- i.XSX - t(Xtilde) %*% iSBY.tF.iS %*% Xtilde 
		i.XSX <- solve(i.XSX)
		
		iSoo.Xtilde <- as.matrix(iS.X - iSBY.tF.iS %*% Xtilde )
    } 
        
    out <- list()
    class(out) <- "predictobject"
    out$opts <- list(only.pars = only.pars, nugget.unobs = nugget.unobs, 
        only.obs = only.obs, pred.var = pred.var, pred.covar = pred.covar, 
        beta.covar = beta.covar, combine.data = combine.data, 
        type = type)
    if (out$opts$only.pars) {
        return(out)
    }
    rm(only.pars, nugget.unobs, only.obs, pred.var, pred.covar, 
        beta.covar, combine.data, type)
        
    if (out$opts$only.obs) {
        idx.unobs <- STdata$obs$idx
        T1 <- STdata$obs$date
        F <- STdata$F
        if (dimensions$L != 0) {
            ST.unobs <- STdata$ST
        }
        N.unobs <- length(unique(idx.unobs))
        T.unobs <- length(unique(T1))
    }   else {
        N.unobs <- dim(STdata$locations)[1]
        T.unobs <- dim(STdata$trend)[1]
        idx.unobs <- rep(1:N.unobs, each = T.unobs)
        T1 <- rep(STdata$trend$date, N.unobs)
        F <- matrix(1, (T.unobs * N.unobs), dimensions$m)
        for (i in (1:dimensions$m)) {
            if (colnames(STdata$F)[i] != "const") {
                F[, i] <- rep(STdata$trend[, colnames(STdata$F)[i]], 
                  N.unobs)
            }
        }
        if (dimensions$L != 0) {
            ST.unobs <- matrix(STdata$ST.all, (T.unobs * N.unobs), 
                dimensions$L)
        }
        if (out$opts$pred.covar) {
            Nmax <- T.unobs
        }
    }
    date.all <- sort(unique(c(object$trend$date, STdata$trend$date)))
    nt.unobs <- nt.obs <- double(length(date.all))
    for (i in c(1:length(date.all))) {
        nt.obs[i] <- sum(object$obs$date == date.all[i])
    }
    if (out$opts$pred.var) {
        if (!out$opts$only.obs) {
	        out$VX <- out$VX.pred <- out$VX.pred.reml <- out$VX.reml <- matrix(NA, T.unobs, N.unobs)
        } else {
	        out$VX <- out$VX.pred <- out$VX.pred.reml <- out$VX.reml <- matrix(NA, ncol(F), 1)
        }
    }
    if (out$opts$pred.covar) {
        out$VX.full <- list()
        out$VX.pred.full <- list()        
        out$VX.pred.reml.full <- list()        
        out$VX.reml.full <- list()        
    }
	
	Xtilde.unobs <- calc.FX(F, STdata$LUR.all, idx.unobs)

    loc.unobs.nu <- STdata$locations[, c("x.nu", "y.nu"), drop = FALSE]
    loc.unobs.beta <- STdata$locations[, c("x.beta", "y.beta"), 
        drop = FALSE]
    I.obs <- match(colnames(object$D.nu), object$locations$ID)
    loc.obs.nu <- object$locations[I.obs, c("x.nu", "y.nu"), 
        drop = FALSE]
    loc.obs.beta <- object$locations[I.obs, c("x.beta", "y.beta"), 
        drop = FALSE]
    Ind.2.1 <- match(object$locations$ID, STdata$locations$ID, 
        nomatch = 0)
    loc.knots.beta <- object$locations.list$knot.coords[object$cov.beta$knots[[1]], ]
    
    ### Calculate Cross covariance in beta/nugget fields
 	sigma.B.C <- makeSigmaB3(cov.pars.beta$pars, 
		coords1=loc.unobs.beta, 
		coordsK=loc.knots.beta, 
		coords2=loc.obs.beta, 
		basis.coords=loc.obs.beta, 
		type = object$cov.beta$covf, 
		nugget = cov.pars.beta$nugget, 
		ind2.to.1 = Ind.2.1, fit.tprs=FALSE)

    cross.D.nu <- crossDist(loc.unobs.nu, loc.obs.nu)
    for (i in 1:ceiling(nrow(F)/Nmax)) {
        Ind <- c((1 + (i - 1) * Nmax):min(i * Nmax, nrow(F)))
        T1.Ind <- T1[Ind]
        for (j in c(1:length(date.all))) {
            nt.unobs[j] <- sum(T1.Ind == date.all[j])
        }
        T1.order <- order(T1.Ind)
        sigma.B.full.C <- calc.FXtF2(F[Ind, , drop = FALSE], 
            sigma.B.C, loc.ind = idx.unobs[Ind], F2 = object$F, 
            loc.ind2 = object$obs$idx)
        sigma.nu.C <- makeSigmaNu(cov.pars.nu$pars, dist = cross.D.nu, 
            type = object$cov.nu$covf, nugget = 0, 
            random.effect = cov.pars.nu$random.effect, 
            ind1 = (idx.unobs[Ind])[T1.order], ind2 = object$obs$idx, 
            blocks1 = nt.unobs, blocks2 = nt.obs, ind2.to.1 = Ind.2.1)
        sigma.nu.C[T1.order, ] <- sigma.nu.C
        sigma.nu.C <- sigma.nu.C + sigma.B.full.C
        if (out$opts$pred.var) {
			I.loc <- sort(unique(idx.unobs[Ind]))
			unobs.D.beta <- crossDist(loc.unobs.beta[I.loc, , drop = FALSE])
			sigma.B.uu <- makeSigmaB3(cov.pars.beta$pars, 
					coords1=loc.unobs.beta[I.loc, , drop=FALSE], 
					coordsK=loc.knots.beta,
					basis.coords=loc.obs.beta,
					type = object$cov.beta$covf, 
					nugget = cov.pars.beta$nugget, fit.tprs=FALSE)
#  			sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = crossDist(loc.unobs.beta), 
#  					type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)

			V.uu <- calc.FXtF2(F[Ind, , drop = FALSE], sigma.B.uu, 
					loc.ind = idx.unobs[Ind] - min(I.loc) + 1)
			unobs.D.nu <- crossDist(loc.unobs.nu[I.loc, , drop = FALSE])
			V.uu <- V.uu + makeSigmaNu(cov.pars.nu$pars, dist = unobs.D.nu, 
					type = object$cov.nu$covf, nugget = 0, 
					random.effect = cov.pars.nu$random.effect, 
					ind1 = idx.unobs[Ind] - min(I.loc) + 1, blocks1 = nt.unobs)
			t.sigma.nu.C <- t(sigma.nu.C)
			iS.Sou <- i.sigma.nu %*% t.sigma.nu.C
			tF.iS.Sou <- t.F.iS.nu %*% t.sigma.nu.C
			if(all(object$cov.beta$covf != "iid")){ #spatial smoothing
				tZB.tF.iS.Sou <- blockMult2(ZB, tF.iS.Sou, n.blocks=dimensions$m, 
										block.sizes1=rep(dimensions$n, dimensions$m), 
										block.sizes2=object$cov.beta$K, transpose=TRUE)
				Sby.tZB.tF.iS.Sou <- solveTriBlock(R.i.sigma.B.Y, tZB.tF.iS.Sou, 
									transpose = TRUE)
			} else { # only nuggets
				Sby.tZB.tF.iS.Sou <- solveTriBlock(i.sigma.P.Y, tF.iS.Sou, 
									transpose = TRUE)			
			}
			
			tmp <- Xtilde.unobs[Ind, , drop = FALSE] - sigma.nu.C %*% iSoo.Xtilde
			
			if (out$opts$pred.covar) {
				V.REML <- (tmp %*% i.XSX) %*% t(tmp)
				tmp <- -sigma.nu.C %*% iS.Sou + t(Sby.tZB.tF.iS.Sou) %*% 
							Sby.tZB.tF.iS.Sou
				V.cond <- V.cond.0 <- V.uu + tmp
				diag(V.cond) <- diag(V.cond) + out$opts$nugget.unobs[idx.unobs[Ind]]
				V.REML <- V.cond + V.REML
				V.REML.0 <- V.cond.0 + V.REML
			}  else {
				V.REML <- rowSums( (tmp %*% i.XSX) %*% t(tmp) )
				tmp <- (-colSums(t(sigma.nu.C) * iS.Sou) + colSums(Sby.tZB.tF.iS.Sou * 
				Sby.tZB.tF.iS.Sou))
				V.cond <- V.cond.0 <- diag(V.uu) + tmp
				V.cond <- V.cond + out$opts$nugget.unobs[idx.unobs[Ind]]
				V.REML <- V.cond + V.REML
				V.REML.0 <- V.cond.0 + V.REML
			}
			if (out$opts$type == "r") {
				stop("NOT YET IMPLEMENTED. Use type=\"p\" instead.")
			}
			if (out$opts$pred.covar) {
				out$VX[Ind] <- diag(V.cond.0)
				out$VX.reml[Ind] <- diag(V.REML.0)
				out$VX.pred[Ind] <- diag(V.cond)
				out$VX.pred.reml[Ind] <- diag(V.REML)
				
				out$VX.full[[i]] <- V.cond.0
				out$VX.reml.full[[i]] <- V.REML.0
				out$VX.pred.full[[i]] <- V.cond
				out$VX.pred.reml.full[[i]] <- V.REML
			} else {
				out$VX[Ind] <- V.cond.0
				out$VX.reml[ind] <- V.REML.0
				out$VX.pred[Ind] <- V.cond
				out$VX.pred.reml[ind] <- V.REML
			}
		}
    }
    
    return(out)
}

calc.CV3 <- function(object, lta=NULL, option="all")
{
	# object is a data.frame two.wk observations and CV predictions/variances
	# lta.object is a data.frame with long-term averages and predictions/variances
	# option is CV subset option
	
    if (!is.function(exp)) {
        stop("'exp' should be a function")
    }

	if(is.null(lta) & option != "comco"){
		print("LTA missing. Cannot calculate LTA Metrics")
	}
	
	if(option == "fixed2")
		option = "fixed"
	
	if(option == "all" | option == "fixed"){
				
		RMSE <- R2 <- COVERAGE <- matrix(NA, nr=2, nc=3)		
		rownames(RMSE) <- row.names(R2) <- row.names(COVERAGE) <- c("2wk", "lta")
		colnames(RMSE) <- colnames(R2) <-  colnames(COVERAGE) <- c("EX.mu", "EX.mu.beta", "EX")	
		

		RMSE[1, ] <- apply(object[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								sqrt( mean( (exp(x) - exp(object$obs))^2 )  )
							}
							)
		R2[1, ] <- apply(object[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (exp(x) - exp(object$obs))^2 )/
												var(exp(object$obs))
									)
								)
							}
							)
							
		COVERAGE[1,3] <- mean(object$coverage)

		RMSE[2, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								sqrt( mean( (x - lta$obs)^2 )  )
							}
							)
		R2[2, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (x - lta$obs)^2 )/var(lta$obs)
										)
									)
							}
							)
		COVERAGE[2,3] <- mean(lta$coverage)
							
	} else if (option == "comco") {
		
		dates <- unique(object$date)		
		RMSE <- t( sapply(dates, 
					function(x) { 
						EX.mu <- sqrt( mean( (	exp(object$obs[object[, "date"] == x]) - 
										exp(object$EX.mu[object[, "date"] == x]) )^2 )  )
						EX.mu.beta <- sqrt( mean( (	exp(object$obs[object[, "date"] == x]) - 
										exp(object$EX.mu.beta[object[, "date"] == x]))^2 )  )
						EX <- sqrt( mean( (	exp(object$obs[object[, "date"] == x]) - 
										exp(object$EX[object[, "date"] == x]) )^2 )  )
						c(EX.mu, EX.mu.beta, EX)
					}
				)
			)
		
		R2 <- t( sapply(dates, 
					function(x) { 
						EX.mu <- mean( 
									(exp(object$obs[object[, "date"] == x]) - 
										exp(object$EX.mu[object[, "date"] == x]))^2 
									) / 
									var( exp(object$obs[object[, "date"] == x]) )
						EX.mu <- max(0, 1-EX.mu)
						EX.mu.beta <- max( c( 0, 1- mean( (exp(object$obs[object[, "date"] == x]) - exp(object$EX.mu.beta[object[, "date"] == x]))^2 ) / var( exp(object$obs[object[, "date"] == x]) )))
						EX <- max( c( 0, 1- mean( (exp(object$obs[object[, "date"] == x]) - exp(object$EX[object[, "date"] == x]))^2 ) / var( exp(object$obs[object[, "date"] == x]) )))
						c(EX.mu, EX.mu.beta, EX)
					}
				)
			)
		
		tmp.coverage <- tapply(object$coverage, object$date, mean) 
		COVERAGE <- matrix(NA, nr=dim(R2)[1], nc=dim(R2)[2])
		COVERAGE[, 3] <- tmp.coverage[as.character(dates)]
		rownames(RMSE) <- rownames(R2) <- rownames(COVERAGE) <- as.character(dates)
		colnames(RMSE) <- colnames(R2)<- colnames(COVERAGE) <- c("EX.mu", "EX.mu.beta", "EX")	
		
	} else if (option == "home") {
	
		RMSE <- R2 <- COVERAGE <- matrix(NA, nr=2, nc=3)		
		rownames(RMSE) <- row.names(R2) <- row.names(COVERAGE) <- c("raw", "detrended")
		colnames(RMSE) <- colnames(R2) <- colnames(COVERAGE) <- c("EX.mu", "EX.mu.beta", "EX")	

		RMSE[1, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								sqrt( mean( ( x - lta$obs)^2 )  )
							}
							)
		R2[1, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (x - lta$obs)^2 )/var(lta$obs)
										)
									)
							}
							)
		COVERAGE[1,3] <- mean(lta$coverage)
		
		R2[2, ] <- apply(lta[, c("EX.mu", "EX.mu.beta", "EX")], 2, 
							function(x) {
								max(c(	0, 
										1 - mean( (x - lta$obs)^2 )/var(lta$obs.detrend)
										)
									)
							}
							)		
	}
	
	out <- list(RMSE=RMSE, R2=R2, COVERAGE=COVERAGE)
	return(out)
}

#### Function to Calculate the effective DF

calc.pd <- function (object, x, type = "p", mean.center=FALSE, ...) {
 
    stCheckClass(object, "STmodel", name = "object")
    type <- tolower(type)
    SpatioTemporal:::stCheckType(type)
    if (inherits(x, "estimateSTmodel")) {
        x <- coef(x, "all")$par
    }
    dimensions <- loglikeSTdim3(object)
    if (type == "f" && length(x) != dimensions$nparam) {
        stop(paste("type=f, requires", dimensions$nparam, "parameters but length(x) =", 
            length(x)))
    }
    if (type != "f" && length(x) == dimensions$nparam) {
        x <- x[(dimensions$nparam - dimensions$nparam.cov + 1):dimensions$nparam]
    }
    if (type != "f" && length(x) != dimensions$nparam.cov) {
        stop(paste("type!=f, requires", dimensions$nparam.cov, 
            "parameters but length(x) =", length(x)))
    }

#     if(any(object$cov.beta$covf == "tprs")){
#     	m=2
#     	object$cov.beta$K <- object$cov.beta$K - choose(m+2-1, 2)*(object$cov.beta$covf=="tprs")
#     }
# 
    tmp <- loglikeSTgetPars3(x, object)
    if (type == "f") {
        gamma <- tmp$gamma
        alpha <- tmp$alpha
    }
    cov.pars.beta <- tmp$cov.beta
    cov.pars.nu <- tmp$cov.nu

	if(mean.center){
		object$LUR <- lapply(object$LUR, function(y) {
							cbind(1, apply(y[, -1], 2, function(yy) yy - mean(yy)))
							})
		object$locations[, c("x.beta", "y.beta")] <- apply(object$locations[, c("x.beta", "y.beta")], 2, 
					function(yy) { yy - mean(yy)})
		
		object$locations.list$knot.coords <- apply(object$locations.list$knot.coords, 2, 
					function(yy){ yy-mean(yy[object$cov.beta$knots[[1]]])})
	}
	
    Xtilde <- lapply(1:ncol(object$F), 
    			function(y) { 
    				calc.FX(object$F[, y, drop=FALSE], list(object$LUR[[y]]), object$obs$idx)
    			})
	D <- lapply(1:length(dimensions$p), function(y) diag(rep(0, dimensions$p[y])))

    i.sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = object$D.nu, 
        type = object$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = object$nt, 
        ind1 = object$obs$idx)
    
    i.sigma.nu <- makeCholBlock(i.sigma.nu, block.sizes = object$nt)
    i.sigma.nu <- invCholBlock(i.sigma.nu, block.sizes = object$nt)

	if(all(object$cov.beta$covf != "iid")){ #spatial smoothing
								
		#### CODE THAT CALCULATES ZB AND ALSO Z, J for TPRS
		ZB.list <- lapply(1:dimensions$m, function(y) {
					if(object$cov.beta$covf[y] != "tprs"){
						makeZB(list(cov.pars.beta$pars[[y]]), 
						coords1=object$locations[, c("x.beta", "y.beta")],
						coordsK=object$locations.list$knot.coords[object$cov.beta$knots[[1]],], 
						type = object$cov.beta$covf[y])
				 	} else {
						makeSigmaBTPRS(cov.pars.beta$pars[[y]], 
						coords1=object$locations[, c("x.beta", "y.beta")],
						K=object$cov.beta$K[y], 
						pen.form=TRUE)$zb
					}
		})
		
		D.list <- lapply(1:dimensions$m, function(y) {
					if(object$cov.beta$covf[y] != "tprs"){
						diag(rep(1/cov.pars.beta$pars[[y]]['sill'], object$cov.beta$K[y]))
					} else {
						makeSigmaBTPRS(cov.pars.beta$pars[[y]], 
						coords1=object$locations[, c("x.beta", "y.beta")],
						K=object$cov.beta$K[y], 
						pen.form=TRUE)$pen*1/cov.pars.beta$pars[[y]]['sill']
					}
		})
		
		Xtilde <- c(Xtilde, 
					lapply(1:ncol(object$F), 
						function(y) {
							calc.FX(F=object$F[, y, drop=FALSE], LUR=list(ZB.list[[y]]), 
										loc.ind=object$obs$idx)
						}))
		D <- c(D, D.list)
		rm(D.list, ZB.list)
		
		if(all(object$cov.beta$nugget)){ # with nugget
			Xtilde <- c(Xtilde, lapply(1:ncol(object$F), 
						function(y){ 
							as.matrix(expandF(object$F[, y, drop=FALSE], object$obs$idx))
						}))
			D <- c(D, lapply(1:length(cov.pars.beta$nugget), 
					function(y) diag(rep(1/cov.pars.beta$nugget[y], dimensions$n))))
		}
	
    } else { # No Spatial Smoothing, only nuggets
			Xtilde <- c(Xtilde, lapply(1:ncol(object$F), 
						function(y){ 
							as.matrix(expandF(object$F[, y, drop=FALSE], object$obs$idx))
						}))
			D <- c(D, lapply(1:length(cov.pars.beta$nugget), 
					function(y) diag(rep(1/cov.pars.beta$nugget[y], dimensions$n))))
   } 
	
	########
	########
	
	mat <- lapply(Xtilde, function(y) t(y)%*%i.sigma.nu)
 	tX.iS.X <- tX.iS.X.D <- matrix(NA, nrow=sum(sapply(Xtilde, ncol)), ncol=sum(sapply(Xtilde, ncol)))		
 	for(ii in 1:length(mat)){
 		for(jj in 1:length(mat)){
 			tmp <- mat[[ii]]%*%t(mat[[jj]])
 			if(ii == 1){
 				rows <- 1:nrow(mat[[ii]])
 			} else {
 				rows <- (1+sum(sapply(mat, nrow)[1:(ii-1)])):sum(sapply(mat, nrow)[1:(ii)])
 			}
 			if(jj == 1){
 				cols <- 1:nrow(mat[[jj]])
 			} else {
 				cols <- (1+sum(sapply(mat, nrow)[1:(jj-1)])):sum(sapply(mat, nrow)[1:(jj)])
 			}
 			tX.iS.X[rows, cols] <- tX.iS.X.D[rows, cols] <- tmp
 			if(ii == jj)
 				tX.iS.X.D[rows, cols] <- tX.iS.X.D[rows, cols] + D[[ii]]
 			rm(tmp)
 		}		
 	}
 	rm(ii, jj, mat)
 	
 	tX.iS.X <- try(solve(tX.iS.X.D, tX.iS.X))
 	if(class(tX.iS.X) == "try-error"){
 		pD2 <- NA
 	} else {
	 	pD2 <- sum(diag(tX.iS.X))		
	}
	rm(tX.iS.X, tX.iS.X.D)
		
	if(object$cov.nu$covf != "iid"){		
		## Create ZV basis functions    
		ZV <- makeZB(list(cov.pars.nu$pars), 
		coords1=object$locations[, c("x.beta", "y.beta")], 
		type = object$cov.nu$covf)
		ZV <- ZV[object$obs$idx, ]
		
		ZV.list <- lapply(1:length(object$nt), function(y) {
							inds <- object$obs$date != unique(object$obs$date)[y]
							tmp <- ZV
							tmp[inds,] <- 0
							tmp[, object$obs$idx[!inds], drop=FALSE]
		})
		#Xtilde <- c(Xtilde, list(ZV.list))
		Xtilde <- c(Xtilde, ZV.list)
		#D <- c(D, diag(rep(1/cov.pars.nu$pars['sill'], nrow(object$obs))))	
		D <- c(D, lapply(object$nt, function(y) {
					if(y == 1){
						as.matrix(1/cov.pars.nu$pars['sill'], nr=1, nc=1)
					} else {
						diag(rep(1/cov.pars.nu$pars['sill'], y))
					}
		}))
	}

	########
	########
	
	tX.iS.X <- tX.iS.X.D <- matrix(NA, nrow=sum(unlist(sapply(Xtilde, ncol))), ncol=sum(unlist(sapply(Xtilde, ncol))))		
	for(ii in 1:length(Xtilde)){
		for(jj in 1:length(Xtilde)){
 			tmp <- t(Xtilde[[ii]])%*%Xtilde[[jj]]
			if(ii == 1){
				rows <- 1:ncol(Xtilde[[ii]])
			} else {
				rows <- (1+sum(unlist(sapply(Xtilde, ncol))[1:(ii-1)])):sum(unlist(sapply(Xtilde, ncol))[1:(ii)])
			}
			if(jj == 1){
				cols <- 1:ncol(Xtilde[[jj]])
			} else {
				cols <- (1+sum(unlist(sapply(Xtilde, ncol))[1:(jj-1)])):sum(unlist(sapply(Xtilde, ncol))[1:(jj)])
			}
			tX.iS.X[rows, cols] <- tX.iS.X.D[rows, cols] <- tmp
			if(ii == jj)
				tX.iS.X.D[rows, cols] <- tX.iS.X.D[rows, cols] + D[[ii]]*cov.pars.nu$nugget[1]
			rm(tmp)
		}		
	}
	rm(ii, jj)
	
	tX.iS.X <- try(solve(tX.iS.X.D, tX.iS.X))
	if(class(tX.iS.X) == "try-error"){
		pD1 <- NA
	} else {
		pD1 <- sum(diag(tX.iS.X))		
	}
	rm(tX.iS.X, tX.iS.X.D)
	
#	pD1 <- NA
	########
	########
	
	out <- list(pD1=pD1, pD2=pD2, pD3=length(x))
    return(out)
}
