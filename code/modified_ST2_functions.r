#####################################################################
##### Casey Olives
##### 07-05-2012

#####################################################################
## This file contains mofified ST2 functions which incorporate
## low-rank modeling work

#####################################################################
### Modified updateCovf3 function
## Use stCheckFields3 function
## Added lines which replicate the number of knots K in nu if only one number provided
## Added data.frame(s) of knot locations to cov.beta
updateCovf3 <- function (STmodel, cov.beta = STmodel$cov.beta, cov.nu = STmodel$cov.nu, select.knots=TRUE) 
{
	require("fields")
    stCheckClass(STmodel, "STmodel", name = "STmodel")
    stCheckFields3(cov.beta, c("covf", "nugget", "K"), name = "cov.beta")
    stCheckFields3(cov.nu, c("covf", "nugget", "random.effect"), 
        name = "cov.nu")
    nt <- dim(STmodel$trend)[2]
    nT <- length(unique(STmodel$obs$date))
	
    for (i in 1:length(cov.beta)) {
        if (length(cov.beta[[i]]) == 1) {
           	cov.beta[[i]] <- rep(cov.beta[[i]], nt)
        } else if (length(cov.beta[[i]]) != nt) {
            stop(sprintf("Covariance specification, %s, for beta has %d!=%d elements.", 
                names(cov.beta)[i], length(cov.beta[[i]]), nt))
        }
    }
	
	if(select.knots)
	{
		## Add in knot locations 
		## if !is.na(cov.beta$K) then select knot locations, else set knot locations at monitoring locations 
		## if!is.na(cov.nu$K) then select knot locations, else set knot locations at monitoring locations at time t
		## Output is the id for the row in STmodel$locations.list$knot.coords
		cov.beta$knots <- list()
		cov.beta$K <- ifelse(is.na(cov.beta$K), nrow(STmodel$locations), cov.beta$K)
		if(cov.beta$K[1] < nrow(STmodel$locations.list$knot.coords))
		{
			tmp <- sort(cover.design(STmodel$locations.list$knot.coords, cov.beta$K[1])$best.id	)
		} else {
			tmp <- 1:nrow(STmodel$locations.list$knot.coords)
		}
		for(i in 1:length(cov.beta$K)){
			cov.beta$knots[[i]] <- tmp
		}
	} else {
		cov.beta$knots <- STmodel$cov.beta$knots
	}# end select.knots
	
    tmp <- unlist(lapply(parsCovFuns(c(cov.beta$covf, cov.nu$covf)), is.null))
    tmp[names(tmp) == "tprs"] <- FALSE
    if (any(tmp)) {
        stop(paste("Unknown covariance specification(s):", paste(names(tmp)[tmp], 
            collapse = ", ")))
    }
    if (!is.logical(cov.nu$random.effect)) {
        stop("'random.effect' specification must be 'logical'.")
    }
    if (!is.logical(cov.beta$nugget)) {
        stop("'nugget' specification for beta-fields must be 'logical'.")
    }
    if (any(cov.beta$covf == "iid" & !cov.beta$nugget)) {
        stop("beta-fields: Covariance model iid requires a nugget.")
    }
    covars <- STmodel$locations
    Ind <- STmodel$locations$ID %in% unique(STmodel$obs$ID)
    covars <- covars[Ind, , drop = FALSE]
    if (is.logical(cov.nu$nugget)) {
        if (cov.nu$nugget[1] == TRUE) {
            cov.nu$nugget <- as.formula("~1", env = .GlobalEnv)
        }
        else {
            cov.nu$nugget <- FALSE
            cov.nu$nugget.matrix <- matrix(NA, dim(covars)[1], 
                0)
        }
    }
    if (class(cov.nu$nugget) == "formula" || is.character(cov.nu$nugget)) {
        if (is.character(cov.nu$nugget)) {
            cov.nu$nugget <- as.formula(paste("~", paste(cov.nu$nugget, 
                collapse = "+")), env = .GlobalEnv)
        }
        covars.tmp <- model.frame(cov.nu$nugget, covars, drop.unused.levels = TRUE)
        tmp <- tryCatch(model.matrix(cov.nu$nugget, covars.tmp), 
            error = function(e) {
                tmp <- model.matrix(cov.nu$nugget, covars)
                test <- apply(tmp, 2, function(x) {
                  length(unique(x))
                })
                test <- test != 1 | names(test) == "(Intercept)"
                tmp <- tmp[, test, drop = FALSE]
                return(tmp)
            })
        cov.nu$nugget.matrix <- as.matrix(tmp)
    } else if (!is.logical(cov.nu$nugget)) {
        stop(paste("Unknown specification of nugget for nu-field:", 
            cov.nu$nugget))
    }
    if (dim(cov.nu$nugget.matrix)[2] == 0 && cov.nu$covf == "iid") {
        stop("nu-field: Covariance model iid requires a nugget.")
    }
    rownames(cov.nu$nugget.matrix) <- covars$ID
    STmodel$cov.beta <- cov.beta
    STmodel$cov.nu <- cov.nu
    return(STmodel)
}

### Modified stCheckFields3 function
## No Changes to this function
stCheckFields3 <- function (x, what, name = "Object") 
{
    if (!all(what %in% names(x))) {
        stop(paste("Field(s):", paste(what[!(what %in% names(x))], 
            collapse = ", "), "- missing from", name))
    }
    return(invisible())
}

### Modified stCheckX function
## Small change that fixed bug with more than one fixed value
stCheckX3 <- function (x, x.fixed, dimensions, type, object) 
{
    if (dim(x)[1] != dimensions$nparam && dim(x)[1] != dimensions$nparam.cov) {
        stop(paste("dim(x)[1] must be", dimensions$nparam, "or", 
            dimensions$nparam.cov))
    }
    if (!is.null(x.fixed) && (!is.vector(x.fixed) || !is.numeric(x.fixed))) {
        stop("'x.fixed' must be a numeric vector.")
    }
    else if (is.null(x.fixed)) {
        x.fixed <- rep(NA, dim(x)[1])
    }
    else if (length(x.fixed) != dimensions$nparam && length(x.fixed) != 
        dimensions$nparam.cov) {
        stop(paste("length(x.fixed) must be", dimensions$nparam, 
            "or", dimensions$nparam.cov))
    }
    if (dim(x)[1] == dimensions$nparam && type != "f") {
        I <- (dimensions$nparam - dimensions$nparam.cov + 1):dimensions$nparam
        if (length(x.fixed) == dimensions$nparam) {
            x.fixed <- x.fixed[I]
        }
        x <- x[I, , drop = FALSE]
    }
    else if (dim(x)[1] == dimensions$nparam.cov && type == "f") {
        x.old <- x
        x <- matrix(NA, dimensions$nparam, dim(x)[2])
        for (i in 1:dim(x)[2]) {
            tmp <- predict.STmodel(object, x.old[, i], only.pars = TRUE, 
                type = "p")$pars
            x[, i] <- c(tmp$gamma.E, tmp$alpha.E, x.old[, i])
        }
        if (length(x.fixed) == dimensions$nparam.cov) {
            x.fixed <- c(rep(NA, dimensions$nparam - dimensions$nparam.cov), 
                x.fixed)
        }
    }
    if (type != "f") {
        names(x.fixed) <- loglikeSTnames3(object, all = FALSE)
    }
    else {
        names(x.fixed) <- loglikeSTnames3(object, all = TRUE)
    }
    x.all <- x
    x <- x[is.na(x.fixed), , drop = FALSE]
    x.all[!is.na(x.fixed), ] <- matrix(x.fixed[!is.na(x.fixed)], 
        nr = sum(!is.na(x.fixed)), ncol = dim(x.all)[2])
    return(list(x.all = x.all, x = x, x.fixed = x.fixed))
}

### Modified createSTmodelInternalDistance function
## Also calculates distance matrix for low rank representations of beta and nu fields
## to be stored in STmodel$locations.list.
## DZ.knots is n \times number of candidate knot locations
## DO.knots is number of candidate knot locations \times number of candidate knot locs
createSTmodelInternalDistance3 <- function (STmodel) {
    if (dim(STmodel$obs)[1] != 0) {
        I.idx <- unique(STmodel$obs$idx)
        I.idx <- 1:max(I.idx)
        STmodel$D.nu <- crossDist(STmodel$locations[I.idx, c("x.nu", 
            "y.nu"), drop = FALSE])
        colnames(STmodel$D.nu) <- rownames(STmodel$D.nu) <- STmodel$locations$ID[I.idx]
        STmodel$D.beta <- crossDist(STmodel$locations[I.idx, 
            c("x.beta", "y.beta"), drop = FALSE])
        colnames(STmodel$D.beta) <- rownames(STmodel$D.beta) <- STmodel$locations$ID[I.idx]
        ## add distance matrix for knot.coords
        if(!is.null(STmodel$locations.list$knot.coords)){
        	STmodel$locations.list$DZ.knots <- crossDist(STmodel$locations[I.idx, c("x", 
            "y"), drop = FALSE], STmodel$locations.list$knot.coords[, , drop = FALSE])
            rownames(STmodel$locations.list$DZ.knots) <- STmodel$locations$ID[I.idx]
            
            STmodel$locations.list$DO.knots <- crossDist(STmodel$locations.list$knot.coords[, , drop = FALSE])
        }
        dates <- sort(unique(STmodel$obs$date))
        STmodel$nt <- double(length(dates))
        for (i in c(1:length(dates))) {
            STmodel$nt[i] <- sum(STmodel$obs$date == dates[i])
        }
    }
    return(STmodel)
}

#### Print ST model
## Now outputs information on number of knots used in Low-Rank Kriging
print.STmodel3 <- function (x, type = x$locations$type, ...) {
    stCheckClass(x, "STmodel", name = "x")
    x$covars <- x$locations
    SpatioTemporal:::commonPrintST(x, "STmodel", 1)
    cat("Models for the beta-fields are:\n")
    print(x$LUR.list)
    if (!is.null(x$scale.covars)) {
        cat("Covariates have been scaled.\n\n")
    }
    if (length(x$ST.list) == 0) {
        cat("No spatio-temporal covariates.\n")
    }
    else {
        cat(sprintf("%d spatio-temporal covariate(s):\n", length(x$ST.list)))
        print(x$ST.list)
    }
    cat("\n")
    cat("Covariance model for the beta-field(s):\n")
    cat(paste("\tCovariance type(s):", paste(x$cov.beta$covf, 
        collapse = ", "), "\n"))
    cat(paste("\tNugget:", paste(c("No", "Yes")[x$cov.beta$nugget + 
        1], collapse = ", "), "\n"))
    cat(paste("\tNumber of Basis Functions for Low-Rank Model:", paste(x$cov.beta$K, 
        collapse = ", "), "\n"))    
    cat("Covariance model for the nu-field(s):\n")
    cat(paste("\tCovariance type:", x$cov.nu$covf, "\n"))
    if (dim(x$cov.nu$nugget.matrix)[2] == 0) {
        cat("\tNugget: No\n")
    }
    else {
        cat(paste("\tNugget:", paste(x$cov.nu$nugget, collapse = ""), 
            "\n"))
    }
    cat(paste("\tRandom effect:", c("No", "Yes")[x$cov.nu$random.effect + 
        1], "\n"))
    cat(paste("\tNumber of Basis Functions for Low-Rank Model:", paste(x$cov.nu$K, 
        collapse = ", "), "\n"))    
         
    SpatioTemporal:::commonPrintST(x, "STmodel", 2, type)
    return(invisible())
}

#### Modified to Allow blocks that  not equally sized
checkDimInternal3 <- function (X) 
{
    dim.X <- sapply(X, dim)

    dim <- list(n = dim.X[1, ], m = length(X), p = dim.X[2, 
        ])
    return(dim)
}

#### Modified parsCovFuns to allow for tprs
parsCovFuns3 <- function(type = namesCovFuns(), list = FALSE){
	out <- parsCovFuns(type, list)
	for(x in 1:length(type)){
		out[[x]] <- unlist(ifelse(type[x] == "tprs", "sill", out[x]))
	}
	return(out)
}

#### Modified loglikeSTdim to allow for tprs
loglikeSTdim3 <- function(STmodel){
   stCheckClass(STmodel, "STmodel", name = "STmodel")
    out <- list()
    out$T <- length(STmodel$nt)
    out$m <- dim(STmodel$F)[2]
    out$n <- dim(STmodel$LUR.all[[1]])[1]
    out$n.obs <- dim(STmodel$LUR[[1]])[1]
    out$p <- sapply(STmodel$LUR, function(x) {
        dim(x)[2]
    })
    out$L <- length(STmodel$ST.list)
    out$npars.beta.covf <- sapply(parsCovFuns3(STmodel$cov.beta$covf, 
        list = TRUE), length)
    out$npars.beta.tot <- (out$npars.beta.covf + STmodel$cov.beta$nugget)
    out$npars.nu.covf <- length(parsCovFuns3(STmodel$cov.nu$covf, 
        list = FALSE))
    out$npars.nu.tot <- out$npars.nu.covf + STmodel$cov.nu$random.effect
    out$npars.nu.tot <- out$npars.nu.tot + dim(STmodel$cov.nu$nugget.matrix)[2]
    out$nparam.cov <- sum(out$npars.beta.tot) + out$npars.nu.tot
    out$nparam <- out$nparam.cov + sum(out$p) + out$L
    return(out)
}

#### Modified loglikeSTgetPars to allow for tprs
loglikeSTgetPars3 <- function(x, STmodel){

   stCheckClass(STmodel, "STmodel", name = "STmodel")
    dimensions <- loglikeSTdim3(STmodel)
    if (length(x) != dimensions$nparam && length(x) != dimensions$nparam.cov) {
        stop(paste("Parameter vector 'x' is", length(x), "should be", 
            dimensions$nparam.cov, "or", dimensions$nparam))
    }
    Ind <- 0
    if (length(x) == dimensions$nparam) {
        if (dimensions$L != 0) {
            gamma <- matrix(x[Ind + (1:dimensions$L)], ncol = 1)
            Ind <- Ind + dimensions$L
        }
        else {
            gamma <- double(0)
        }
        alpha <- vector("list", dimensions$m)
        for (i in 1:dimensions$m) {
            alpha[[i]] <- matrix(x[Ind + (1:dimensions$p[i])], 
                ncol = 1)
            Ind <- Ind + dimensions$p[i]
        }
    }
    else {
        alpha <- double(0)
        gamma <- double(0)
    }
    cov.beta <- list(pars = vector("list", dimensions$m), nugget = double(dimensions$m))
    for (i in 1:dimensions$m) {
        if (dimensions$npars.beta.covf[i] != 0) {
            cov.beta$pars[[i]] <- exp(x[Ind + (1:dimensions$npars.beta.covf[i])])
            Ind <- Ind + dimensions$npars.beta.covf[i]
            names(cov.beta$pars[[i]]) <- parsCovFuns3(STmodel$cov.beta$covf[i])
        }
        else {
            cov.beta$pars[[i]] <- numeric(0)
        }
        if (STmodel$cov.beta$nugget[i]) {
            cov.beta$nugget[i] <- exp(x[Ind + 1])
            Ind <- Ind + 1
        }
        else {
            cov.beta$nugget[i] <- 0
        }
    }
    cov.nu <- list()
    if (dimensions$npars.nu.covf != 0) {
        cov.nu$pars <- exp(x[Ind + (1:dimensions$npars.nu.covf)])
        Ind <- Ind + dimensions$npars.nu.covf
         names(cov.nu$pars) <- parsCovFuns3(STmodel$cov.nu$covf)
    }
    else {
        cov.nu$pars <- numeric(0)
    }
    if (dim(STmodel$cov.nu$nugget.matrix)[2] == 0) {
        cov.nu$nugget <- matrix(0, dim(STmodel$cov.nu$nugget.matrix)[1], 
            1)
    }
    else {
        tmp <- matrix(x[Ind + (1:dim(STmodel$cov.nu$nugget.matrix)[2])], 
            ncol = 1)
        cov.nu$nugget <- exp(STmodel$cov.nu$nugget.matrix %*% 
            tmp)
        Ind <- Ind + dim(STmodel$cov.nu$nugget.matrix)[2]
    }
    rownames(cov.nu$nugget) <- rownames(STmodel$cov.nu$nugget.matrix)
    if (STmodel$cov.nu$random.effect) {
        cov.nu$random.effect <- exp(x[Ind + 1])
    }
    else {
        cov.nu$random.effect <- 0
    }
    return(list(gamma = gamma, alpha = alpha, cov.beta = cov.beta, 
        cov.nu = cov.nu))

}

#### Modified loglikeSTnames to allow for tprs
loglikeSTnames3 <- function (STmodel, all = TRUE) 
{
    stCheckClass(STmodel, "STmodel", name = "STmodel")
    dimensions <- loglikeSTdim3(STmodel)
    out <- NULL
    if (all) {
        if (dimensions$L != 0) {
            out <- c(out, paste("gamma", colnames(STmodel$ST), 
                sep = "."))
        }
        tmp <- lapply(STmodel$LUR, colnames)
        for (i in 1:length(tmp)) {
            tmp[[i]] <- paste("alpha", names(STmodel$LUR)[i], 
                tmp[[i]], sep = ".")
        }
        out <- c(out, unlist(tmp))
    }
    tmp <- parsCovFuns3(STmodel$cov.beta$covf, list = TRUE)
    for (i in 1:length(tmp)) {
        tmp[[i]] <- c(tmp[[i]], switch(STmodel$cov.beta$nugget[i], 
            "nugget"))
        tmp[[i]] <- paste("log", tmp[[i]], names(STmodel$LUR)[i], 
            names(tmp)[i], sep = ".")
    }
    out <- c(out, unlist(tmp))
    if (dim(STmodel$cov.nu$nugget.matrix)[2] == 0) {
        tmp <- c(parsCovFuns3(STmodel$cov.nu$covf, list = FALSE), 
            switch(STmodel$cov.nu$random.effect, "random.effect"))
    }
    else {
        tmp <- c(parsCovFuns3(STmodel$cov.nu$covf, list = FALSE), 
            paste("nugget", colnames(STmodel$cov.nu$nugget.matrix), 
                sep = "."), switch(STmodel$cov.nu$random.effect, 
                "random.effect"))
    }
    tmp <- paste("nu", "log", tmp, STmodel$cov.nu$covf, sep = ".")
    out <- c(out, tmp)
    names(out) <- NULL
    return(out)
}

### Add T to LUR for TPRS Smooth here, not later
### Also added knot.coords and K in cov.beta
createSTmodel3 <- function(STdata, 
							LUR=NULL, 
							ST=NULL, 
							cov.beta=list(covf="exp", nugget=FALSE, K=NA),
                          	cov.nu=list(covf="exp", nugget=TRUE, random.effect=FALSE),
                          	locations=list(coords=c("x","y"), long.lat=NULL, coords.beta=NULL, coords.nu=NULL, others=NULL, knot.coords=NULL),
                         	strip=FALSE, 
                         	scale=FALSE, 
                         	scale.covars=NULL){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")

  ##Default values for several of the inputs
  cov.beta <- defaultList(cov.beta, eval(formals(createSTmodel3)$cov.beta) )
  cov.nu <- defaultList(cov.nu, eval(formals(createSTmodel3)$cov.nu) )
  locations <- defaultList(locations, eval(formals(createSTmodel3)$locations ))  
  
  ##add any nu-field nugget covariates to the locations
  if( class(cov.nu$nugget)=="formula" ){
    nugget.names <- all.vars(cov.nu$nugget)
  }else if( is.character(cov.nu$nugget) ){
    nugget.names <- cov.nu$nugget
  }else{
    nugget.names <- NULL
  }
  ##add the names needed by nugget to locations$others
  locations$others <- c(locations$others, nugget.names)

  ## Either provide knot coordinates or use monitoring locations
  if(is.null(locations$knot.coords)){
	  locations$knot.coords <- STdata$covars[, locations$coords]
  }
  ##create output object, start with a copy
  STmodel <- STdata
  ##and assign both classes
  class(STmodel) <- c("STmodel","STdata")

  ##first add trend if missing
  if( is.null(STmodel$trend) ){
    warning("No smooth trend in 'STdata', assuming only an intercept.")
    STmodel <- updateSTdataTrend(STmodel, n.basis=0)
  }
  ##sort temporal trend
  STmodel$trend <- STmodel$trend[order(STmodel$trend$date),,drop=FALSE]
  ##sort SpatioTemporal covariate to match time-points
  if( !is.null(STmodel$SpatioTemporal) ){
    STmodel$SpatioTemporal <-
      STmodel$SpatioTemporal[as.character(STmodel$trend$date),,,drop=FALSE]
  }
  
  ##should we drop some covariates?
  ID.unique <- sort(unique(STmodel$obs$ID))
  IND <- STmodel$covars$ID %in% ID.unique
  if( strip ){
    STmodel$covars <- STmodel$covars[IND,,drop=FALSE]
  }else{
    ##if not dropping at least order sites with observed locations first
    STmodel$covars <- rbind(STmodel$covars[IND,,drop=FALSE],
                           STmodel$covars[!IND,,drop=FALSE])
  }
  ##sort SpatioTemporal covariate to match locations
  if( !is.null(STmodel$SpatioTemporal) ){
    STmodel$SpatioTemporal <- STmodel$SpatioTemporal[,STmodel$covars$ID,,drop=FALSE]
  }
  
  ##create a seperate locations data.frame
  STmodel$locations.list <- locations
  STmodel$locations <- processLocation(STmodel, locations)
  
  ##add idx variables for the observations
  STmodel$obs$idx <- match(STmodel$obs$ID, STmodel$locations$ID)
  ##sort the observations
  STmodel$obs <- STmodel$obs[order(STmodel$obs$date, STmodel$obs$idx),,drop=FALSE]

  ##Covariance specification, default values
  ##process covariance specification
  STmodel <- updateCovf3(STmodel, cov.beta, cov.nu)
  
  ##should covariates be scaled
  if( scale ){
    if( is.null(scale.covars) ){
      scale.covars <- unlist(lapply(STmodel$covars,is.factor))
      ##mean and sd
      suppressWarnings( mean.covars <- unlist(lapply(STmodel$covars,mean)) )
      suppressWarnings( sd.covars <- unlist(lapply(STmodel$covars,sd)) )
      ##which fields do not have mean/sd
      scale.covars <- (is.na(mean.covars) | is.na(sd.covars) | scale.covars)
      mean.covars[scale.covars] <- NA
      sd.covars[scale.covars] <- NA
      ##save for future use
      STmodel$scale.covars <- list(mean=mean.covars, sd=sd.covars)
    }else{
      STmodel$scale.covars <- scale.covars
    }
  }
  
	##process LUR
	# Check for TPRS SMooths and add x, y to LUR if not already there
	if(any(cov.beta$covf == "tprs")){
		STmodel$covars$x.beta = STmodel$locations$x.beta
		STmodel$covars$y.beta = STmodel$locations$y.beta  	
	}

	for(iii in 1:length(STmodel$cov.beta$covf))
	{
		if(STmodel$cov.beta$covf[iii]=="tprs"){
			if(class(LUR[[iii]]) == "formula"){
				LUR[[iii]] <- as.formula(paste("~", paste(c(as.character(LUR[[iii]])[2], "x.beta", "y.beta"), collapse="+")))
			}
			if (length(LUR[[iii]]) == 0) {
					LUR[[iii]] <- as.formula("~1+x.beta+y.beta", env = .GlobalEnv)
					next
			}
			if (is.numeric(LUR[[iii]])) {
				if (any(LUR[[iii]] > dim(STdata$covars)[2])) {
					warning(sprintf("Unable to find LUR with indeces larger than %d", 
					  dim(STdata$covars)[2]))
					LUR[[iii]] <- LUR[[iii]][LUR[[iii]] <= dim(STdata$covars)[2]]
					if (length(LUR[[iii]]) == 0) {
					  LUR[[iii]] <- as.formula("~1+x.beta+y.beta", env = .GlobalEnv)
					  next
					}
				}
				LUR[[iii]] <- names(STdata$covars)[LUR[[iii]]]
			}
			if (is.character(LUR[[iii]])) {
				LUR[[iii]] <- as.formula(paste("~", paste(c(unique(LUR[[iii]]), "x.beta", "y.beta"), 
					collapse = "+")), env = .GlobalEnv)
			} 
			
# 			temp.name <- names(LUR)[iii] 
# 			LUR[[iii]] = c(LUR[[iii]], "x.beta", "y.beta")
# 			names(LUR)[[iii]] = temp.name
		}
	}  	
	names(LUR) <- c("const", names(STdata$trend)[-which(names(STdata$trend) == "date")])

  STmodel$LUR.list <- processLUR(STmodel, LUR)
  ##...and ST specification
  STmodel$ST.list <- processST(STmodel, ST)
  ##create covariate matrices
  STmodel <- SpatioTemporal:::createLUR(STmodel, STmodel$LUR.list)
  ##create spatio-temporal covariate(s) for observations, ST=M and ST.all
  STmodel <- SpatioTemporal:::createST(STmodel, STmodel$ST.list)
  
  ##create temporal trends for observations, F
  ##match times with observations
  STmodel$F <- STmodel$trend[ match(STmodel$obs$date, STmodel$trend$date),,drop=FALSE]
  ##add intercept
  STmodel$F <- cbind(rep(1,dim(STmodel$F)[1]), STmodel$F)
  names(STmodel$F)[1] <- "const"
  ##drop date column
  STmodel$F <- STmodel$F[,-which(names(STmodel$F)=="date"),drop=FALSE]
  ##cast to matrix
  STmodel$F <- data.matrix(STmodel$F)
  rownames(STmodel$F) <- as.character(STmodel$obs$date)
  
  ##test for colocated monitoring sites, beta-fields
  I.idx <- unique(STmodel$obs$idx)
  if( anyDuplicated(STmodel$locations[I.idx,c("x.beta","y.beta")]) ){
    ##if so, require nugget
    if( any(!STmodel$cov.beta$nugget) ){
      stop("Colocated monitoring sites, nugget is required for beta-fields.")
    }
  }
  if( anyDuplicated(STmodel$locations[I.idx,c("x.nu","y.nu")]) ){
    ##if so, require nugget
    if( dim(STmodel$cov.nu$nugget.matrix)[2]==0 ){
      stop("Colocated monitoring sites, nugget is required for nu-fields.")
    }
  }
  ##compute distance matrix and nt
  STmodel <-  createSTmodelInternalDistance3(STmodel)
  
  ##drop SpatioTemporal, has been replaced by other things
  STmodel$SpatioTemporal <- NULL
  STmodel$covars <- NULL
  ##and asign one class
  class(STmodel) <- "STmodel"
  ##return
  return( STmodel )
}##function createSTmodel3

#### Optimized Log-Likelihood for low rank BETA only
loglikeST3 <- function (x = NULL, STmodel, type = "p", x.fixed = NULL) 
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
    if (type == "f") {
        mu.B <- matrix(calc.mu.B(STmodel$LUR, alpha), ncol = 1)
        if (dimensions$L != 0) {
            Y <- Y - STmodel$ST %*% gamma
        }
    }
 ##############################
 ##############################
	sigma.B <- makeSigmaB3(cov.pars.beta$pars, coords1=STmodel$locations[, c("x.beta", "y.beta")],
	coordsK=STmodel$locations.list$knot.coords[STmodel$cov.beta$knots[[1]],], 
	type = STmodel$cov.beta$covf, nugget = cov.pars.beta$nugget)
    sigma.B <- try(makeCholBlock(sigma.B, n.blocks = dimensions$m), 
        silent = TRUE)
    if (class(sigma.B) == "try-error") {
        return(-.Machine$double.xmax)
    }
    l <- -sumLogDiag(sigma.B)
    sigma.B <- invCholBlock(sigma.B, n.blocks = dimensions$m)
	
	sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = STmodel$D.nu, 
	type = STmodel$cov.nu$covf, nugget = cov.pars.nu$nugget, 
	random.effect = cov.pars.nu$random.effect, blocks1 = STmodel$nt, 
	ind1 = STmodel$obs$idx)
	sigma.nu <- try(makeCholBlock(sigma.nu, n.blocks = dimensions$T, 
	block.sizes = STmodel$nt), silent = TRUE)
	if (class(sigma.nu) == "try-error") {
	return(-.Machine$double.xmax)
	}
	l <- l - sumLogDiag(sigma.nu)
	sigma.nu <- invCholBlock(sigma.nu, block.sizes = STmodel$nt)
			
 ################################
 ################################
 
    i.sR.Y <- blockMult(sigma.nu, Y, block.sizes = STmodel$nt)
    if (type == "f") {
        i.sB.mu.B <- blockMult(sigma.B, mu.B, n.blocks = dimensions$m)
        l <- l - dotProd(i.sB.mu.B, mu.B)/2
        l <- l - dotProd(i.sR.Y, Y)/2
    } else {
        iS.X <- calc.iS.X(STmodel$LUR, sigma.B)
        if (dimensions$L != 0) {
            i.sR.M <- blockMult(sigma.nu, STmodel$ST, block.sizes = STmodel$nt)
        }
        F.i.sR.Y <- calc.tFX(STmodel$F, i.sR.Y, STmodel$obs$idx, 
            n.loc = dimensions$n.obs)
        if (dimensions$L != 0) {
            F.i.sR.M <- calc.tFX(STmodel$F, i.sR.M, STmodel$obs$idx, 
                n.loc = dimensions$n.obs)
        }
    }
    sigma.B.Y <- calc.tFXF(STmodel$F, sigma.nu, STmodel$obs$idx, 
        block.sizes = STmodel$nt, n.loc = dimensions$n.obs)
    sigma.B.Y <- sigma.B.Y + sigma.B
    sigma.B.Y <- try(makeCholBlock(sigma.B.Y, n.blocks = 1), 
        silent = TRUE)
    if (class(sigma.B.Y) == "try-error") {
        return(-.Machine$double.xmax)
    }
    l <- l - sumLogDiag(sigma.B.Y)
    if (type == "f") {
        mu.B.Y <- i.sB.mu.B + calc.tFX(STmodel$F, i.sR.Y, STmodel$obs$idx, 
            n.loc = dimensions$n.obs)
        mu.B.Y <- solveTriBlock(sigma.B.Y, mu.B.Y, transpose = TRUE)
        l <- l + norm2(mu.B.Y)/2
    } else {
        Y.hat <- solveTriBlock(sigma.B.Y, F.i.sR.Y, transpose = TRUE)
        if (dimensions$L != 0) {
            M.hat <- solveTriBlock(sigma.B.Y, F.i.sR.M, transpose = TRUE)
        }
        sigma.B.Y <- invCholBlock(sigma.B.Y)
        sigma.B.Y.iS.X <- sigma.B.Y %*% iS.X
        i.sigma.alpha.Y <- -t(iS.X) %*% sigma.B.Y.iS.X + calc.X.iS.X(STmodel$LUR, 
            iS.X)
        i.sigma.alpha.Y <- try(makeCholBlock(i.sigma.alpha.Y), 
            silent = TRUE)
        if (class(i.sigma.alpha.Y) == "try-error") {
            return(-.Machine$double.xmax)
        }
        if (type == "r") {
            l <- l - sumLogDiag(i.sigma.alpha.Y)
        }
        Y.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.Y
        Y.hat.2 <- solveTriBlock(i.sigma.alpha.Y, Y.hat.2, transpose = TRUE)
        Y.sigma.hat.Y <- dotProd(Y, i.sR.Y) - norm2(Y.hat) - 
            norm2(Y.hat.2)
        l <- l - Y.sigma.hat.Y/2
        if (dimensions$L != 0) {
            M.hat.2 <- t(sigma.B.Y.iS.X) %*% F.i.sR.M
            M.hat.2 <- solveTriBlock(i.sigma.alpha.Y, M.hat.2, 
                transpose = TRUE)
            Y.sigma.hat.M <- (t(Y) %*% i.sR.M - t(Y.hat) %*% 
                M.hat - t(Y.hat.2) %*% M.hat.2)
            Y.sigma.hat.M <- t(Y.sigma.hat.M)
            M.sigma.hat.M <- (t(STmodel$ST) %*% i.sR.M - t(M.hat) %*% 
                M.hat - t(M.hat.2) %*% M.hat.2)
           M.sigma.hat.M <- try(makeCholBlock(M.sigma.hat.M), 
               silent = TRUE)
            if (class(M.sigma.hat.M) == "try-error") {
    	        return(-.Machine$double.xmax)
            }
            if (type == "r") {
               l <- l - sumLogDiag(M.sigma.hat.M)
            }
           Y.sigma.hat.M <- solveTriBlock(M.sigma.hat.M, Y.sigma.hat.M, 
               transpose = TRUE)
        }
    }
    l <- as.double(l)
    if (!is.finite(l)) 
        l <- -.Machine$double.xmax
    return(l)
}

## Modified to allow for use with new optimized LL function
loglikeSTGrad3 <- function (x, STmodel, type = "p", x.fixed = NULL, h = 0.001, 
    diff.type = 0) {
    func <- function(x0) {
        loglikeST3(x0, STmodel, type, x.fixed)
    }
    df <- genGradient(x, func, h = h, diff.type = diff.type)
    return(df)
}

#### Estimation for low rank BETA only
estimate.STmodel3 <- function (object, x, x.fixed = NULL, type = "p", h = 0.001, diff.type = 1, 
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
        	loglikeSTGrad3(x, STmodel, type, x.fixed, h = h, diff.type = diff.type)
    	}
    res <- as.list(rep(NA, dim(x)[2]))
    conv <- rep(FALSE, dim(x)[2])
    value <- rep(NA, dim(x)[2])
    err <- tryCatch(loglikeST3(x[, 1], STmodel = object, type = type, 
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
        try(res[[i]] <- optim(x[, i], loglikeST3, gr = loglikeSTGrad.loc, 
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
            suppressWarnings(par.sd <- sqrt(-diag(solve(res[[i]]$hessian))))
            if (type != "f") {
                par.type <- "par.cov"
            }
            else {
                par.type <- "par.all"
            }
            res[[i]][[par.type]]$init <- x.all[, i]
            res[[i]][[par.type]]$par <- res[[i]][[par.type]]$fixed <- x.fixed
            res[[i]][[par.type]]$par[is.na(x.fixed)] <- res[[i]]$par
            res[[i]][[par.type]]$sd[is.na(x.fixed)] <- par.sd
            if (type != "f") {
                tmp <- predict.STmodel3(object, res[[i]]$par.cov$par, 
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

## predict.STmodel3 does not try to do anything fancy with nu
predict.STmodel3 <- function (object, x, STdata = NULL, Nmax = 1000, only.pars = FALSE, 
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
# 	i.sigma.B <- makeSigmaB(cov.pars.beta$pars, dist = object$D.beta, 
#         type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
	i.sigma.B <- makeSigmaB3(cov.pars.beta$pars, coords1=object$locations[, c("x.beta", "y.beta")],
			coordsK=object$locations.list$knot.coords[object$cov.beta$knots[[1]],], 
			type = object$cov.beta$covf, nugget = cov.pars.beta$nugget, fit.tprs=FALSE)
    i.sigma.B <- makeCholBlock(i.sigma.B, n.blocks = dimensions$m)
    i.sigma.B <- invCholBlock(i.sigma.B, n.blocks = dimensions$m)

    i.sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = object$D.nu, 
        type = object$cov.nu$covf, nugget = cov.pars.nu$nugget, 
        random.effect = cov.pars.nu$random.effect, blocks1 = object$nt, 
        ind1 = object$obs$idx)
    i.sigma.nu <- makeCholBlock(i.sigma.nu, block.sizes = object$nt)
    i.sigma.nu <- invCholBlock(i.sigma.nu, block.sizes = object$nt)
    
    tF.iS.F <- calc.tFXF(object$F, i.sigma.nu, object$obs$idx, 
        block.sizes = object$nt, n.loc = dimensions$n.obs)
    R.i.sigma.B.Y <- tF.iS.F + i.sigma.B
    R.i.sigma.B.Y <- makeCholBlock(R.i.sigma.B.Y)
    if (type == "f") {
        gamma.E <- c(tmp$gamma)
        alpha.E <- unlist(tmp$alpha)
        gamma.V <- matrix(0, length(gamma.E), length(gamma.E))
        alpha.V <- matrix(0, length(alpha.E), length(alpha.E))
        gamma.alpha.C <- matrix(0, length(gamma.E), length(alpha.E))
        gamma.alpha <- as.matrix(c(gamma.E, alpha.E))
    }   else {
        iS.nu.X <- blockMult(i.sigma.nu, Xtilde, block.sizes = object$nt)
        tF.iS.nu.X <- calc.tFX(object$F, iS.nu.X, object$obs$idx, 
            n.loc = dimensions$n.obs)
        iSBY.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, tF.iS.nu.X, 
            transpose = TRUE)
        iSBY.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, iSBY.tF.iS.X, 
            transpose = FALSE)
        i.XSX <- t(Xtilde) %*% iS.nu.X
        i.XSX <- i.XSX - t(tF.iS.nu.X) %*% iSBY.tF.iS.X
        iS.nu.Y <- blockMult(i.sigma.nu, object$obs$obs, block.sizes = object$nt)
        tF.iS.nu.Y <- calc.tFX(object$F, iS.nu.Y, object$obs$idx, 
            n.loc = dimensions$n.obs)
        tmp2 <- object$obs$obs %*% iS.nu.X
        tmp2 <- tmp2 - t(tF.iS.nu.Y) %*% iSBY.tF.iS.X
        gamma.alpha <- solve(i.XSX, t(tmp2))
        rm(iS.nu.X, tF.iS.nu.X, iSBY.tF.iS.X, iS.nu.Y, tF.iS.nu.Y, 
            tmp2)
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
    class(out) <- "predictSTmodel"
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
    iS.nu.C <- blockMult(i.sigma.nu, C.minus.mu, block.sizes = object$nt)
    tF.iS.nu.C <- calc.tFX(object$F, iS.nu.C, object$obs$idx, 
        n.loc = dimensions$n.obs)
    iSBY.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, tF.iS.nu.C, 
        transpose = TRUE)
    iSBY.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, iSBY.tF.iS.C, 
        transpose = FALSE)
    tF.iS <- calc.tFX(object$F, i.sigma.nu, object$obs$idx, n.loc = dimensions$n.obs)
    obs.mu <- iS.nu.C - t(tF.iS) %*% iSBY.tF.iS.C
    rm(C.minus.mu, iS.nu.C, tF.iS.nu.C, iSBY.tF.iS.C)
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
        Sby.iSb.Sou <- i.sigma.B %*% t(sigma.B.C)
        Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, 
            transpose = TRUE)
        Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, 
            transpose = FALSE)
		if(any(object$cov.beta$covf == "tprs")){
			print("sorry, beta.covar not implemented with tprs smooth")
			out$beta$VX <- 0*out$beta$EX
		} else {
        	if (out$opts$beta.covar) {
# 				sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = crossDist(loc.unobs.beta), 
# 					type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
				sigma.B.uu <- makeSigmaB3(cov.pars.beta$pars, 
					coords1=loc.unobs.beta, coordsK=loc.knots.beta, basis.coords=loc.obs.beta,
					type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
				tmp <- sigma.B.uu - sigma.B.C %*% tF.iS.F %*% Sby.iSb.Sou
				out$beta$VX.full <- list()
				for (i in 1:dim(out$beta$EX)[2]) {
					Ind <- (1:dim(out$beta$EX)[1]) + (i - 1) * dim(out$beta$EX)[1]
					out$beta$VX.full[[i]] <- tmp[Ind, Ind, drop = FALSE]
					rownames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
					colnames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
				}
				names(out$beta$VX.full) <- colnames(out$beta$EX)
				tmp <- diag(tmp)
        	} else {
			   sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = matrix(0, 
					 1, 1), type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
				sigma.B.uu <- matrix(diag(sigma.B.uu), ncol = dim(sigma.B.uu)[1], 
					nrow = dim(loc.unobs.beta)[1], byrow = TRUE)
				tmp <- c(sigma.B.uu) - rowSums(sigma.B.C * t(tF.iS.F %*% 
					Sby.iSb.Sou))
			}
			out$beta$VX <- matrix(tmp, ncol = dim(out$beta$EX)[2])
			dimnames(out$beta$VX) <- dimnames(out$beta$EX)
			rm(tmp, Sby.iSb.Sou, tF.iS.F, sigma.B.uu)
        }
    }
    rm(i.sigma.B)
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
            I.loc <- sort(unique(idx.unobs[Ind]))
            unobs.D.beta <- crossDist(loc.unobs.beta[I.loc, , 
                drop = FALSE])
#             sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = unobs.D.beta, 
#                 type = object$cov.beta$covf, nugget = cov.pars.beta$nugget)
            sigma.B.uu <- makeSigmaB3(cov.pars.beta$pars, 
            				coords1=loc.unobs.beta[I.loc, , drop=FALSE], 
            				coordsK=loc.knots.beta,
            				basis.coords=loc.obs.beta,
                			type = object$cov.beta$covf, 
                			nugget = cov.pars.beta$nugget, fit.tprs=FALSE)
            V.uu <- calc.FXtF2(F[Ind, , drop = FALSE], sigma.B.uu, 
                loc.ind = idx.unobs[Ind] - min(I.loc) + 1)
            unobs.D.nu <- crossDist(loc.unobs.nu[I.loc, , drop = FALSE])
            V.uu <- V.uu + makeSigmaNu(cov.pars.nu$pars, dist = unobs.D.nu, 
                type = object$cov.nu$covf, nugget = 0, random.effect = cov.pars.nu$random.effect, 
                ind1 = idx.unobs[Ind] - min(I.loc) + 1, blocks1 = nt.unobs)
            t.sigma.nu.C <- t(sigma.nu.C)
            iS.Sou <- blockMult(i.sigma.nu, t.sigma.nu.C, block.sizes = object$nt)
            tF.iS.Sou <- tF.iS %*% t.sigma.nu.C
            Sby.tF.iS.Sou <- solveTriBlock(R.i.sigma.B.Y, tF.iS.Sou, 
                transpose = TRUE)
            if (out$opts$pred.covar) {
                tmp <- -sigma.nu.C %*% iS.Sou + t(Sby.tF.iS.Sou) %*% 
                  Sby.tF.iS.Sou
                V.cond <- V.cond.0 <- V.uu + tmp
                diag(V.cond) <- diag(V.cond) + out$opts$nugget.unobs[idx.unobs[Ind]]
            }       else {
                tmp <- (-colSums(t(sigma.nu.C) * iS.Sou) + colSums(Sby.tF.iS.Sou * 
                  Sby.tF.iS.Sou))
                V.cond <- V.cond.0 <- diag(V.uu) + tmp
                V.cond <- V.cond + out$opts$nugget.unobs[idx.unobs[Ind]]
            }
            if (out$opts$type == "r") {
                stop("NOT YET IMPLEMENTED. Use type=\"p\" instead.")
            }
            if (out$opts$pred.covar) {
                out$VX[Ind] <- diag(V.cond.0)
                out$VX.pred[Ind] <- diag(V.cond)
                out$VX.full[[i]] <- V.cond.0
            }      else {
                out$VX[Ind] <- V.cond.0
                out$VX.pred[Ind] <- V.cond
            }
        }
    }
    if (out$opts$pred.var) {
        out$VX <- pmax(out$VX, 0)
        out$VX.pred <- pmax(out$VX.pred, 0)
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
            dimnames(out$VX) <- dimnames(out$VX.pred) <- dimnames(out$EX)
        }
        if (out$opts$pred.covar) {
            names(out$VX.full) <- colnames(out$EX)
            for (i in 1:length(out$VX.full)) colnames(out$VX.full[[i]]) <- rownames(out$VX.full[[i]]) <- rownames(out$EX)
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

## dropObservations3 just does what dropObservations does but calls updateCovf3 and createSTmodelInternalDistance3
dropObservations3 <- function (STmodel, Ind.cv) {
    stCheckClass(STmodel, "STmodel", name = "STmodel")
    N <- dim(STmodel$obs)[1]
    if ((is.vector(Ind.cv) && length(Ind.cv) != N) || (!is.vector(Ind.cv) && 
        any(dim(Ind.cv) != c(N, 1)))) {
        stop("Length of Ind.cv must match dim(STmodel$obs)[1]")
    }
    STmodel$obs <- STmodel$obs[!Ind.cv, , drop = FALSE]
    STmodel$F <- STmodel$F[!Ind.cv, , drop = FALSE]
    if (length(STmodel$ST) != 0) {
        STmodel$ST <- STmodel$ST[!Ind.cv, , drop = FALSE]
    }
    IND <- STmodel$locations$ID %in% unique(STmodel$obs$ID)
    STmodel$locations <- STmodel$locations[IND, , drop = FALSE]
    if (length(STmodel$ST.all) != 0) {
        STmodel$ST.all <- STmodel$ST.all[, IND, , drop = FALSE]
    }
    for (i in 1:length(STmodel$LUR)) {
        IND <- rownames(STmodel$LUR[[i]]) %in% unique(STmodel$obs$ID)
        STmodel$LUR[[i]] <- STmodel$LUR[[i]][IND, , drop = FALSE]
    }
    STmodel$LUR.all <- STmodel$LUR
    STmodel$obs$idx <- match(STmodel$obs$ID, STmodel$locations$ID)
    STmodel <- createSTmodelInternalDistance3(STmodel)
    STmodel <- updateCovf3(STmodel, select.knots=FALSE)
    if(nrow(STmodel$locations) < STmodel$cov.beta$K[1]){
    	print("Fewer locations than knots. Resetting K so that model is full rank (ie n = k).")
    	STmodel$cov.beta$K <- rep(nrow(STmodel$locations), length(STmodel$cov.beta$K))
#	    STmodel <- updateCovf3(STmodel, select.knots=TRUE)    	
		tmp <- as.numeric(rownames(STmodel$locations))	
		for(ii in 1:length(STmodel$cov.beta$K)){
			STmodel$cov.beta$knots[[ii]] <- tmp
		}
    }
    return(STmodel)
}

### estimateCV.STmodel3 for low rank BETA only
estimateCV.STmodel3 <- function (object, x, Ind.cv, control = list(trace = 3), ...) {

    stCheckClass(object, "STmodel", name = "object")
    Ind.cv <- SpatioTemporal:::stCheckInternalCV(Ind.cv)
    control <- defaultList(control, eval(formals(estimate.STmodel3)$control))
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
        res[[i]] <- estimate.STmodel3(object.aux, x, ...)
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

### predictCV.STmodel3 for low rank BETA only
predictCV.STmodel3 <- function (object, x, Ind.cv, ..., silent = TRUE) 
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
    pred <- list()
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
        pred[[i]] <- predict.STmodel3(object.obs, x[, i], STdata = object.pred, 
            nugget.unobs = nugget.unobs, only.pars = FALSE, combine.data = FALSE, 
            ...)
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
    for (i in 1:N.CV.sets) {
        Ind.current <- Ind.cv == i
        I <- pred[[i]]$I$I
        out$pred.obs[Ind.current, c("EX.mu", "EX.mu.beta", "EX")] <- cbind(pred[[i]]$EX.mu[I], 
            pred[[i]]$EX.mu.beta[I], pred[[i]]$EX[I])
        if (out$opts$pred.var) {
            out$pred.obs[Ind.current, c("VX", "VX.pred")] <- cbind(pred[[i]]$VX[I], 
                pred[[i]]$VX.pred[I])
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
                out$pred.all$VX[I] <- pred[[i]]$VX
                out$pred.all$VX.pred[I] <- pred[[i]]$VX.pred
            }
            out$pred.all$beta$mu[ID.names, ] <- pred[[i]]$beta$mu
            out$pred.all$beta$EX[ID.names, ] <- pred[[i]]$beta$EX
            if (out$opts$pred.var) {
                out$pred.all$beta$VX[ID.names, ] <- pred[[i]]$beta$VX
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
                  VX <- createDataMatrix(obs = pred[[i]]$VX, 
                    date = pred[[i]]$I$date, ID = pred[[i]]$I$ID)
                  VX.pred <- createDataMatrix(obs = pred[[i]]$VX.pred, 
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

### Naive log likelihood for models with and without nugget
loglikeSTnaive3 <- function (x = NULL, STmodel, type = "p", x.fixed = NULL) {
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
    if (type == "f") {
        mu.B <- calc.mu.B(STmodel$LUR, alpha)
        mean.val <- 0
        for (i in 1:dimensions$m) {
            mean.val <- mean.val + mu.B[STmodel$obs$idx, i] * 
                STmodel$F[, i]
        }
        Y <- Y - mean.val
        if (dimensions$L != 0) {
            Y <- Y - STmodel$ST %*% gamma
        }
    }
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
    Y <- solveTriBlock(sigma.nu, Y, transpose = TRUE)
    l <- l - norm2(Y)/2
    if (type != "f") {
        Ftmp <- calc.FX(STmodel$F, STmodel$LUR, STmodel$obs$idx)
        if (dimensions$L != 0) {
            Ftmp <- cbind(Ftmp, STmodel$ST)
        }
        Ftmp <- solveTriBlock(sigma.nu, Ftmp, transpose = TRUE)
        FY <- t(Ftmp) %*% Y
        sigma.alt <- t(Ftmp) %*% Ftmp
        sigma.alt <- try(makeCholBlock(sigma.alt, n.blocks = 1), 
            silent = TRUE)
        if (class(sigma.alt) == "try-error") {
            return(-.Machine$double.xmax)
        }
        if (type == "r") {
            l <- l - sumLogDiag(sigma.alt)
        }
        FY <- solveTriBlock(sigma.alt, FY, transpose = TRUE)
        l <- l + norm2(FY)/2
    }
    l <- as.double(l)
    if (!is.finite(l)) 
        l <- -.Machine$double.xmax
    return(l)
}
