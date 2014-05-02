###########################################################################################################
###########################################################################################################
## 
## Name:		rags2ridges
## Authors:		Carel F.W. Peeters & Wessel N. van Wieringen
##			Molecular Biostatistics Unit
##			Dept. of Epidemiology & Biostatistics
##			VU University medical center
##			Amsterdam, the Netherlands
## Email:		cf.peeters@vumc.nl; w.vanwieringen@vumc.nl
## 
## Version:		1.2
## Last Update:	01/05/2014
## Description:	Ridge estimation, with supporting functions, for high-dimensional precision matrices
## Publication:	Wessel N. van Wieringen & Carel F.W. Peeters (2014)
##			"Ridge Estimation of Inverse Covariance Matrices for High-Dimensional Data"
##			arXiv:1403.0904 [stat.ME]. 
##
###########################################################################################################
###########################################################################################################




##---------------------------------------------------------------------------------------------------------
## 
## Hidden support functions
##
##---------------------------------------------------------------------------------------------------------

.trace <- function(M){
	#####################################################################################################
	# - Internal function to compute the trace of a matrix
	# - Faster support function (as opposed to 'matrix.trace') when input M is already forced to 'matrix'
	# - M > matrix input
	#####################################################################################################

	return(sum(diag(M)))
}



.is.int <- function(x, tolerance = .Machine$double.eps){
	#####################################################################################################
	# - Logical function that checks if a number is an integer within machine precision
	# - x         > input number
	# - tolerance > tolerance threshold for determining integer quality
	#####################################################################################################

	abs(x - round(x)) < tolerance
}



.LL <- function(S, P){
	#####################################################################################################
	# - Function that computes the value of the (negative) log-likelihood
	# - S > sample covariance matrix
	# - P > precision matrix (possibly regularized inverse of covariance or correlation matrix)
	#####################################################################################################

    	LL <- -log(det(P)) + .trace(S %*% P) 
    	return(LL)
}



.FrobeniusLoss <- function(O, P){
	#####################################################################################################
	# - Function computing Frobenius loss
	# - O > Estimated (possibly regularized) precision matrix
	# - P > True (population) precision matrix
	#####################################################################################################

	return(sum(abs(O - P)^2))
}



.QuadraticLoss <- function(O, C){
	#####################################################################################################
	# - Function computing Quadratic loss
	# - O > Estimated (possibly regularized) precision matrix
	# - C > True (population) covariance matrix
	#####################################################################################################

	return((sum(abs((O %*% C - diag(ncol(O))))^2)))
}




##---------------------------------------------------------------------------------------------------------
## 
## Support functions
##
##---------------------------------------------------------------------------------------------------------

symm <- function(M){
	#####################################################################################################
	# - Large objects that are symmetric sometimes fail to be recognized as such by R due to
	#   rounding under machine precision. This function symmetrizes for computational purposes
	#   matrices that are symmetric in numeric ideality
	# - M > symmetric (in numeric ideality) square matrix
	#####################################################################################################

	# Dependencies
	# require("base")

	if (!is.matrix(M)){
		stop("M should be a matrix")
	}
	else if (nrow(M) != ncol(M)){
		stop("M should be a square matrix")
	}
	else {
		# Symmetrize
		Msym <- (M + t(M))/2

		# Return
		return(Msym)
	}
}



adjacentMat <- function(M, diag = FALSE){
	#####################################################################################################
	# - Function that transforms a real matrix into an adjacency matrix
	# - Intended use: Turn sparsified precision matrix into an adjacency matrix for undirected graph
	# - M    > (sparsified precision) matrix
	# - diag > logical indicating if the diagonal elements should be retained
	#####################################################################################################

	# Dependencies
	# require("base")

	if (!is.matrix(M)){
		stop("M should be a matrix")
	}
	else if (nrow(M) != ncol(M)){
		stop("M should be square matrix")
	}
	else {
		# Create adjacency matrix
		AM <- M
		AM[AM != 0] <- 1
		diag(AM) <- 0

		if (diag){
			diag(AM) <- 1
		}

		# Return
		return(AM)
	}
}



covML <- function(Y){
	#####################################################################################################
	# - function that gives the maximum likelihood estimate of the covariance matrix
	# - Y   > (raw) data matrix, assumed to have variables in columns
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("stats")

	if (!is.matrix(Y)){
		stop("Input (Y) should be a matrix")
	}
	else {
		Ys  <- scale(Y, center = TRUE, scale = FALSE)
		Sml <- (t(Ys) %*% Ys)/nrow(Ys)
		return(Sml)
	}
}



evaluateS <- function(S, verbose = TRUE){
	#####################################################################################################
	# - Function evualuating various properties of an input matrix
	# - Intended use is to evaluate the various properties of what is assumed to be a covariance matrix
	# - Another use is to evaluate the various properties of a (regularized) precision matrix
	# - S       > sample covariance/correlation matrix or (regularized) precision matrix
	# - verbose > logical indicating if output should be printed on screen
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("stats")
	
	if (!is.matrix(S)){
		stop("S should be a matrix")
	} 
	else if (nrow(S) != ncol(S)){
		stop("S should be a square matrix")
	}
	else if (class(verbose) != "logical"){
		stop("Input (verbose) is of wrong class")
	} 
	else {
		Sproperties <- list()

		# Is S symmetric?
		Sproperties$symm <- isSymmetric(S)

		# Are eigenvalues S real and positive?
		evs                   <- eigen(S)$values 
		Sproperties$realEigen <- all(Im(evs) == 0) 
		Sproperties$posEigen  <- all(evs > 0)

		# Is S diagonally dominant?
		Sproperties$dd <- all(abs(cov2cor(S)[upper.tri(S)]) < 1)

		# Trace and determinant S
		Sproperties$trace <- sum(diag(S))
		Sproperties$det   <- det(S)

		# Spectral condition number S
		Sproperties$condNumber <- max(evs) / min(evs)	

		if (verbose){
			cat("Properties of S:", "\n")
			cat("----------------------------------------", "\n")
			cat(paste("-> symmetric : ", Sproperties$symm, sep=""), "\n")
			cat(paste("-> eigenvalues real : ", Sproperties$realEigen, sep=""), "\n")
			cat(paste("-> eigenvalues > 0 : ", Sproperties$posEigen, sep=""), "\n")
			cat(paste("-> diag. dominant : ", Sproperties$dd, sep=""), "\n")
	 		cat("", "\n")
			cat(paste("-> trace : ", round(Sproperties$trace, 5), sep=""), "\n")
			cat(paste("-> determinant : ", round(Sproperties$det, 5), sep=""), "\n")
			cat(paste("-> cond. number : ", round(Sproperties$condNumber, 5), sep=""), "\n")
			cat("----------------------------------------", "\n")
		}

		# Return
		return(Sproperties)
	}
}



pcor <- function(P, pc = TRUE){
	#####################################################################################################
	# - Function computing partial correlation/standardized precision matrix from a precision matrix
	# - P  > precision matrix (possibly regularized inverse of covariance or correlation matrix)
	# - pc > logical indicating if the partial correlation matrix should be computed
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("stats")

	if (!is.matrix(P)){
		stop("P should be a matrix")
	}
	else if (!isSymmetric(P)){
		stop("P should be a symmetric matrix")
	}
	else {
		# Compute partial correlation matrix
		if (pc){
			P       <- -P
			diag(P) <- -diag(P)
			Pcor    <- cov2cor(P)
			return(Pcor)
		} 

		# Compute standardized precision matrix
		else {
			SP <- cov2cor(P)
			return(SP)
		}
	}
}




##---------------------------------------------------------------------------------------------------------
## 
## Function for Ridge Estimators of the Precision matrix
##
##---------------------------------------------------------------------------------------------------------

ridgeS <- function(S, lambda, type = "Alt", target = diag(1/diag(S))){
	#####################################################################################################
	# - Function that calculates Ridge estimators of a precision matrix
	# - S       > sample covariance matrix
	# - lambda  > penalty parameter (choose in accordance with type of Ridge estimator)
	# - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
	# - Alt     > van Wieringen-Peeters alternative ridge estimator of a precision matrix
	# - ArchI   > Archetypal I ridge estimator of a precision matrix
	# - ArchII  > Archetypal II ridge estimator of a precision matrix
	# - target  > target (precision terms) for Type I estimators, default = diag(1/diag(S))
	#
	# - NOTES:
	# - All types are rotation equivariant
	# - When type = "Alt" and target is p.d., one obtains the van Wieringen-Peeters type I estimator
	# - When type = "Alt" and target is null-matrix, one obtains the van Wieringen-Peeters type II est.
	# - When target is not the null-matrix it is expected to be p.d. for the vWP type I estimator
	# - The target is always expected to be p.d. in case of the archetypal I estimator
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("expm")

	if (!isSymmetric(S)){
		stop("S should be a covariance matrix")
	}
	else if (lambda <= 0){
		stop("lambda should be positive")
	}
	else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
		stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
	}
	else{
		# Calculate Ridge estimator
		# Alternative estimator
		if (type == "Alt"){
			if (!isSymmetric(target)){
				stop("Shrinkage target should be symmetric")
			} else if (dim(target)[1] != dim(S)[1]){
				stop("S and target should be of the same dimension")
			} else if (!identical(target, mat.or.vec(nrow(S),ncol(S))) & any(eigen(target, only.values = T)$values <= 0)){
				stop("When target is not a null-matrix it should be p.d.")
			} else {
				P_Alt <- solve((S - lambda * target)/2 + sqrtm(((S - lambda * target) %*% (S - lambda * target))/4 + lambda * diag(nrow(S))))
				return(P_Alt)
			}
		}

		# Archetypal I
		if (type == "ArchI"){
			if (lambda > 1){
				stop("lambda should be in (0,1] for this type of Ridge estimator")
			} else if (!isSymmetric(target)){
				stop("Shrinkage target should be symmetric")
			} else if (dim(target)[1] != dim(S)[1]){
				stop("S and target should be of the same dimension")
			} else if (any(eigen(target, only.values = T)$values <= 0)){
				stop("Target should be p.d.")
			} else {
				P_ArchI <- solve((1-lambda) * S + lambda * solve(target))
				return(P_ArchI)
			}
		}

		# Archetypal II
		if (type == "ArchII"){
			P_ArchII <- solve(S + lambda * diag(nrow(S)))
			return(P_ArchII)
		}
	}
}




##---------------------------------------------------------------------------------------------------------
## 
## Functions for Penalty Parameter selection
##
##---------------------------------------------------------------------------------------------------------

optPenaltyCV <- function(Y, lambdaMin, lambdaMax, step, type = "Alt", target = diag(1/diag(covML(Y))), 
                         targetScale = TRUE, output = "light", graph = TRUE, verbose = TRUE){ 
	#####################################################################################################
	# - Function that selects the optimal penalty parameter by leave-one-out cross-validation
	# - Y           > (raw) Data matrix, variables in columns
	# - lambdaMin   > minimum value penalty parameter (dependent on 'type')
	# - lambdaMax   > maximum value penalty parameter (dependent on 'type')
	# - step        > determines the coarseness in searching the grid [lambdaMin, lambdaMax]
	# - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
	# - target      > target (precision terms) for Type I estimators, default = diag(1/diag(covML(Y)))
	# - targetScale > logical, if TRUE: default target diag(1/diag(S)) is used, where S is based on 
	#			CV sample; otherwise target is fixed for all CV samples and given by 'target'
	# - output      > must be one of {"all", "light"}, default = "light"
	# - graph       > Optional argument for visualization optimal penalty selection, default = TRUE
	# - verbose     > logical indicating if intermediate output should be printed on screen
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("stats")
	# require("graphics")

	if (class(verbose) != "logical"){
		stop("Input (verbose) is of wrong class")
	} 
	if (verbose){
		cat("Perform input checks...", "\n")
	}
	if (!is.matrix(Y)){
		stop("Y should be a matrix")
	} 
	else if (class(lambdaMin) != "numeric"){
		stop("Input (lambdaMin) is of wrong class")
	} 
	else if (length(lambdaMin) != 1){
		stop("lambdaMin must be a scalar")
	} 
	else if (lambdaMin <= 0){
		stop("lambdaMin must be positive")
	} 
	else if (class(lambdaMax) != "numeric"){
		stop("Input (lambdaMax) is of wrong class")
	} 
	else if (length(lambdaMax) != 1){
		stop("lambdaMax must be a scalar")
	} 
	else if (lambdaMax <= lambdaMin){
		stop("lambdaMax must be larger than lambdaMin")
	}
	else if (class(step) != "numeric"){
		stop("Input (step) is of wrong class")
	}
	else if (!.is.int(step)){
		stop("step should be integer")
	}
	else if (step <= 0){
		stop("step should be a positive integer")
	}
	else if (class(targetScale) != "logical"){
		stop("Input (targetScale) is of wrong class")
	}
	else if (!(output %in% c("all", "light"))){
		stop("output should be one of {'all', 'light'}")
	}
	else if (class(graph) != "logical"){
		stop("Input (graph) is of wrong class")
	}
	else {
		# Set preliminaries
		LLs     <- numeric()
		lambdas <- seq(lambdaMin, lambdaMax, len = step)

		# Calculate CV scores
		if (verbose){cat("Calculating cross-validated negative log-likelihoods...", "\n")}
		for (k in 1:length(lambdas)){
			slh <- numeric()
        		for (i in 1:nrow(Y)){
				S   <- covML(Y[-i,])
				if (targetScale)  {slh <- c(slh, .LL(t(Y[i,,drop = F]) %*% Y[i,,drop = F], ridgeS(S, lambdas[k], type = type)))}
				if (!targetScale) {slh <- c(slh, .LL(t(Y[i,,drop = F]) %*% Y[i,,drop = F], ridgeS(S, lambdas[k], type = type, target = target)))}
			}
			
			LLs <- c(LLs, mean(slh))
			if (verbose){cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")}
		}

		# Visualization
		optLambda <- min(lambdas[which(LLs == min(LLs))])
		if (graph){
			plot(log(lambdas), type = "l", log(LLs), xlab = "log(penalty value)", ylab = "log(cross-validated neg. log-likelihood)", main = type)
			abline(h = log(min(LLs)), v = log(optLambda), col = "red")
			legend("topleft", legend = c(paste("min. CV LL: ", round(min(LLs),3), sep = ""), 
			paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
		}

		# Return
		S <- covML(Y)
		if (output == "all"){
			if (targetScale) {return(list(optLambda = optLambda, optPrec = ridgeS(S, optLambda, type = type), lambdas = lambdas, LLs = LLs))}
			if (!targetScale){return(list(optLambda = optLambda, optPrec = ridgeS(S, optLambda, type = type, target = target), lambdas = lambdas, LLs = LLs))}
		}
		if (output == "light"){
			if (targetScale) {return(list(optLambda = optLambda, optPrec = ridgeS(S, optLambda, type = type)))}
			if (!targetScale){return(list(optLambda = optLambda, optPrec = ridgeS(S, optLambda, type = type, target = target)))}
		}
	}
}   



optPenalty.aLOOCV <- function(Y, lambdaMin, lambdaMax, step, type = "Alt", target = diag(1/diag(covML(Y))), 
                              output = "light", graph = TRUE, verbose = TRUE){
	#####################################################################################################
	# - Function that selects the optimal penalty parameter by approximate leave-one-out cross-validation
	# - Y           > (raw) Data matrix, variables in columns
	# - lambdaMin   > minimum value penalty parameter (dependent on 'type')
	# - lambdaMax   > maximum value penalty parameter (dependent on 'type')
	# - step        > determines the coarseness in searching the grid [lambdaMin, lambdaMax]
	# - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
	# - target      > target (precision terms) for Type I estimators, default = diag(1/diag(covML(Y)))
	# - output      > must be one of {"all", "light"}, default = "light"
	# - graph       > Optional argument for visualization optimal penalty selection, default = TRUE
	# - verbose     > logical indicating if intermediate output should be printed on screen
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("graphics")

	if (class(verbose) != "logical"){
		stop("Input (verbose) is of wrong class")
	} 
	if (verbose){
		cat("Perform input checks...", "\n")
	}
	if (!is.matrix(Y)){
		stop("Y should be a matrix")
	} 
	else if (class(lambdaMin) != "numeric"){
		stop("Input (lambdaMin) is of wrong class")
	} 
	else if (length(lambdaMin) != 1){
		stop("lambdaMin must be a scalar")
	} 
	else if (lambdaMin <= 0){
		stop("lambdaMin must be positive")
	} 
	else if (class(lambdaMax) != "numeric"){
		stop("Input (lambdaMax) is of wrong class")
	} 
	else if (length(lambdaMax) != 1){
		stop("lambdaMax must be a scalar")
	} 
	else if (lambdaMax <= lambdaMin){
		stop("lambdaMax must be larger than lambdaMin")
	}
	else if (class(step) != "numeric"){
		stop("Input (step) is of wrong class")
	}
	else if (!.is.int(step)){
		stop("step should be integer")
	}
	else if (step <= 0){
		stop("step should be a positive integer")
	}
	else if (!(output %in% c("all", "light"))){
		stop("output should be one of {'all', 'light'}")
	}
	else if (class(graph) != "logical"){
		stop("Input (graph) is of wrong class")
	}
	else {
		# Set preliminaries
		S       <- covML(Y)
		n       <- nrow(Y)
		lambdas <- seq(lambdaMin, lambdaMax, len = step)
		aLOOCVs <- numeric()

		# Calculate approximate LOOCV scores
		if (verbose){cat("Calculating approximate LOOCV scores...", "\n")}
		for (k in 1:length(lambdas)){
			P    <- ridgeS(S, lambdas[k], type = type, target = target)
			nLL  <- .5 * .LL(S, P)
			isum <- numeric()

        		for (i in 1:n){
				S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
				isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
			}

			aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
			if (verbose){cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")}
		}

		# Visualization
		optLambda <- min(lambdas[which(aLOOCVs == min(aLOOCVs))])
		if (graph){
			plot(log(lambdas), type = "l", log(aLOOCVs), xlab = "log(penalty value)", ylab = "log(Approximate LOOCV neg. log-likelihood)", main = type)
			abline(h = log(min(aLOOCVs)), v = log(optLambda), col = "red")
			legend("topleft", legend = c(paste("min. approx. LOOCV LL: ", round(min(aLOOCVs),3), sep = ""), 
			paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
		}

		# Return
		if (output == "all"){
			return(list(optLambda = optLambda, optPrec = ridgeS(S, optLambda, type = type, target = target), lambdas = lambdas, aLOOCVs = aLOOCVs))
		}
		if (output == "light"){
			return(list(optLambda = optLambda, optPrec = ridgeS(S, optLambda, type = type, target = target)))
		}
	}
}



conditionNumber <- function(S, lambdaMin, lambdaMax, step, type = "Alt", target = diag(1/diag(S)),
				    verticle = FALSE, value = 1, verbose = TRUE){
	#####################################################################################################
	# - Function that visualizes the spectral condition number against the regularization parameter
	# - Can be used to heuristically determine the (minimal) value of the penalty parameter
	# - The ridges are rotation equivariant, meaning they work by shrinking the eigenvalues
	# - Maximum shrinkage implies that all eigenvalues will be equal
	# - Ratio of maximum and minimum eigenvalue of P can then function as a heuristic 
	# - It's point of stabilization can give an acceptable value for the penalty
	# - The ratio boils down to the (spectral) condition number of a matrix
	# - S         > sample covariance/correlation matrix
	# - lambdaMin > minimum value penalty parameter (dependent on 'type')
	# - lambdaMax > maximum value penalty parameter (dependent on 'type')
	# - step      > determines the coarseness in searching the grid [lambdaMin, lambdaMax]
	# - type      > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
	# - target    > target (precision terms) for Type I estimators, default = diag(1/diag(S))
	# - verticle  > optional argument for visualization verticle line in graph output, default = FALSE
	#               Can be used to indicate the value of, e.g., the optimal penalty as indicated by some
	#               routine. Can be used to assess if this optimal penalty will lead to a 
	#               well-conditioned estimate
	# - value     > indicates constant on which to base verticle line when verticle = TRUE. Default = 1
	# - verbose   > logical indicating if intermediate output should be printed on screen
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("graphics")

	if (class(verbose) != "logical"){
		stop("Input (verbose) is of wrong class")
	} 
	if (verbose){
		cat("Perform input checks...", "\n")
	}
	if (!is.matrix(S)){
		stop("S should be a matrix")
	} 
	else if (class(lambdaMin) != "numeric"){
		stop("Input (lambdaMin) is of wrong class")
	} 
	else if (length(lambdaMin) != 1){
		stop("lambdaMin must be a scalar")
	} 
	else if (lambdaMin <= 0){
		stop("lambdaMin must be positive")
	} 
	else if (class(lambdaMax) != "numeric"){
		stop("Input (lambdaMax) is of wrong class")
	} 
	else if (length(lambdaMax) != 1){
		stop("lambdaMax must be a scalar")
	} 
	else if (lambdaMax <= lambdaMin){
		stop("lambdaMax must be larger than lambdaMin")
	}
	else if (class(step) != "numeric"){
		stop("Input (step) is of wrong class")
	}
	else if (!.is.int(step)){
		stop("step should be integer")
	}
	else if (step <= 0){
		stop("step should be a positive integer")
	}
	else if (class(verticle) != "logical"){
		stop("Input (verticle) is of wrong class")
	}
	else {
		# Set preliminaries
		lambdas <- seq(lambdaMin, lambdaMax, len = step)
		condNR  <- numeric()

		# Calculate spectral condition number ridge estimate on lambda grid
		if (verbose){cat("Calculating spectral condition numbers...", "\n")}
		for (k in 1:length(lambdas)){
			P         <- ridgeS(S, lambdas[k], type = type, target = target)
			condNR[k] <- as.numeric((eigen(P, only.values = T)$values)[1]/(eigen(P, only.values = T)$values)[nrow(P)])
			if (verbose){cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")}
		}

		# Visualization
		plot(log(lambdas), type = "l", condNR, xlab = "log(penalty value)", ylab = "spectral condition number", main = type)
		if (verticle){
			if (class(value) != "numeric"){
				stop("Input (value) is of wrong class")
			} else if (length(value) != 1){
				stop("Input (value) must be a scalar")
			} else if (value <= 0){
				stop("Input (value) must be positive")
			} else {
				abline(v = log(value), col = "red")
			}
		}
	}
}




##---------------------------------------------------------------------------------------------------------
## 
## Test for Vanishing Partial Correlations
##
##---------------------------------------------------------------------------------------------------------

sparsify <- function(P, type = c("threshold", "localFDR"), threshold = .25, FDRcut = .8, verbose = TRUE){
	#####################################################################################################
	# - Function that sparsifies/determines support of a partial correlation matrix
	# - Support can be determined ad-hoc by thresholding or by local FDRs 
	# - Local FDR operates on the nonredundant non-diagonal elements of a partial correlation matrix
	# - Function is to some extent a wrapper around certain 'GeneNet' and 'fdrtool' functions
	# - P         > (possibly shrunken) precision matrix
	# - type      > signifies type of sparsification: either thresholding or local FDR testing
	# - threshold > cut-off for partial correlation elements selection based on thresholding
	#               Only when type = 'threshold'. Default = .25
	# - FDRcut    > cut-off for partial correlation element selection based on local FDR
	#               Only when type = 'localFDR'. Default = .8
	# - verbose   > logical indicating if intermediate output should be printed on screen
	#               Only when type = 'localFDR'. Default = TRUE
	#####################################################################################################

	# Dependencies
	# require("base")
	# require("stats")
	# require("corpcor")
	# require("longitudinal")
	# require("fdrtool")
	# require("GeneNet")

	if (!is.matrix(P)){
		stop("P should be a matrix")
	}
	else if (!isSymmetric(P)){
		stop("P should be a symmetric matrix")
	}
	else if (!evaluateS(P, verbose = FALSE)$posEigen){
		stop("P is expected to be positive definite")
	}
	else if (missing(type)){
		stop("Need to specify type of sparsification ('threshold' or 'localFDR')")
	}
	else if (!(type %in% c("threshold", "localFDR"))){
		stop("type should be one of {'threshold', 'localFDR'}")
	}
	else {
		# Obtain partial correlation matrix
		if (sum(diag(P)) == ncol(P)){
			PC <- P
		} else {
      		PC <- pcor(P)
		}

		# Obtain sparsified matrix
		if (type == "threshold"){
			if (class(threshold) != "numeric"){
				stop("Input (threshold) is of wrong class")
			} else if (length(threshold) != 1){
				stop("threshold must be a scalar")
			} else if (threshold <= 0){
				stop("threshold must be positive")
			} else {
				PC0 <- PC
				PC0[!(abs(PC0) >= threshold)] <- 0
			}
		}

		if (type == "localFDR"){
			if (class(FDRcut) != "numeric"){
				stop("Input (FDRcut) is of wrong class")
			} else if (length(FDRcut) != 1){
				stop("FDRcut must be a scalar")
			} else if (FDRcut < 0 | FDRcut > 1){
				stop("FDRcut must be in the interval [0,1]")
			} else if (class(verbose) != "logical"){
				stop("Input (verbose) is of wrong class")
			} else {
				if (verbose){
					PCtest <- ggm.test.edges(PC, fdr = TRUE, direct = FALSE, plot = TRUE)
					PCnet  <- extract.network(PCtest, cutoff.ggm = FDRcut); print(PCnet)
				} else {
					PCtest <- ggm.test.edges(PC, fdr = TRUE, direct = FALSE, plot = FALSE)
					PCnet  <- extract.network(PCtest, cutoff.ggm = FDRcut)
				}
				EdgeSet <- PCnet[,2:3]
				PC0     <- diag(nrow(PC))
				for(k in 1:dim(EdgeSet)[1]){
					PC0[EdgeSet[k,1],EdgeSet[k,2]] = PC0[EdgeSet[k,2],EdgeSet[k,1]] <- PC[EdgeSet[k,1],EdgeSet[k,2]]
				} 
			}
		}

		# Return
		colnames(PC0) = rownames(PC0) <- colnames(P)
		return(PC0)
	}
}




##---------------------------------------------------------------------------------------------------------
## 
## Functions for Loss/Entropy Evaluation
##
##---------------------------------------------------------------------------------------------------------

loss <- function(E, T, precision = TRUE, type = c("frobenius", "quadratic")){
	#####################################################################################################
	# - Function evualuating various loss functions on the precision
	# - E         > Estimated (possibly regularized) precision matrix
	# - T         > True (population) covariance or precision matrix
	# - precision > Logical indicating if T is a precision matrix (when TRUE)
	# - type      > character indicating which loss function is to be used
	#####################################################################################################

	if (!is.matrix(E)){
		stop("Input (E) is of wrong class")
	} 
	else if (!isSymmetric(E)){
		stop("E should be a symmetric matrix")
	}
	else if (!is.matrix(T)){
		stop("Input (T) is of wrong class")
	} 
	else if (!isSymmetric(T)){
		stop("T should be a symmetric matrix")
	}
	else if (dim(E)[1] != dim(T)[1]){
		stop("E and T should be of the same dimension")
	}
	else if (class(precision) != "logical"){
		stop("Input (precision) is of wrong class")
	}
	else if (missing(type)){
		stop("Need to specify loss type ('frobenius' or 'quadratic')")
	}
	else if (!(type %in% c("frobenius", "quadratic"))){
		stop("type should be one of {'frobenius', 'quadratic'}")
	}
	else {
		# Frobenius loss
		if (type == "frobenius"){
			if (precision)  {loss <- .FrobeniusLoss(E, T)}
			if (!precision) {loss <- .FrobeniusLoss(E, solve(T))}
		}

		# Quadratic loss
		if (type == "quadratic"){
			if (precision)  {loss <- .QuadraticLoss(E, solve(T))}
			if (!precision) {loss <- .QuadraticLoss(E, T)}
		}

		# Return
		return(loss)
	}
}



KLdiv <- function(Mtest, Mref, Stest, Sref, symmetric = FALSE){
	#####################################################################################################
	# - Function that calculates the Kullback-Leibler divergence between two normal distributions
	# - Mtest     > mean vector approximating m.v. normal distribution
	# - Mref      > mean vector 'true'/reference m.v. normal distribution
	# - Stest     > covariance matrix approximating m.v. normal distribution
	# - Sref      > covariance matrix 'true'/reference m.v. normal distribution
	# - symmetric > logical indicating if original symmetric version of KL div. should be calculated
	#####################################################################################################

	# Dependencies
	# require("base")

	if (class(Mtest) != "numeric"){
		stop("Input (Mtest) is of wrong class")
	} 
	else if (class(Mref) != "numeric"){
		stop("Input (Mref) is of wrong class")
	} 
	else if (length(Mtest) != length(Mref)){
		stop("Mtest and Mref should be of same length")
	} 
	else if (!is.matrix(Stest)){
		stop("Input (Stest) is of wrong class")
	} 
	else if (!is.matrix(Sref)){
		stop("Input (Sref) is of wrong class")
	} 
	else if (!isSymmetric(Stest)){
		stop("Stest should be symmetric")
	}
	else if (!isSymmetric(Sref)){
		stop("Sref should be symmetric")
	}
	else if (dim(Stest)[1] != length(Mtest)){
		stop("Column and row dimension of Stest should correspond to length Mtest")
	}
	else if (dim(Sref)[1] != length(Mref)){
		stop("Column and row dimension of Sref should correspond to length Mref")
	}
	else if (class(symmetric) != "logical"){
		stop("Input (symmetric) is of wrong class")
	} 
	else {
		# Evaluate KL divergence
		KLd <- (sum(diag(solve(Stest) %*% Sref)) + t(Mtest - Mref) %*% solve(Stest) %*% (Mtest - Mref) 
			  - nrow(Sref) - log(det(Sref)) + log(det(Stest)))/2

		# Evaluate (original) symmetric version KL divergence
		if (symmetric){
			KLd <- KLd + (sum(diag(solve(Sref) %*% Stest)) + t(Mref - Mtest) %*% solve(Sref) %*% (Mref - Mtest) 
				 	  - nrow(Sref) - log(det(Stest)) + log(det(Sref)))/2
		}

	# Return
	return(as.numeric(KLd))
	}
}




##---------------------------------------------------------------------------------------------------------
## 
## Functions for Visualization
##
##---------------------------------------------------------------------------------------------------------

if(getRversion() >= "2.15.1") utils::globalVariables(c("X1", "X2", "value"))

edgeHeat <- function(M, lowColor = "blue", highColor = "red", textsize = 10, diag = TRUE, legend = TRUE, main = ""){
	#####################################################################################################
	# - function that visualizes precision matrix as a heatmap
	# - can be used to assess (visually) the performance of set of graphical modeling techniques
	# - M         > Precision matrix
	# - lowColor  > determines color scale in the negative range, default = "blue"
	# - highColor > determines color scale in the positive range, default = "red"
	# - textsize  > set textsize row and column labels, default = 10
	# - diag      > logical determining treatment diagonal elements M. If FALSE, then the diagonal
	#               elements are given the midscale color of white; only when M is a square matrix
	# - legend    > optional inclusion of color legend, default = TRUE
	# - main      > character specifying the main title, default = ""
	#####################################################################################################

	# Dependencies
	#require("ggplot2")
	#require("reshape")

	if (!is.matrix(M)){
		stop("Supply 'M' as matrix")
	}
	else if (class(lowColor) != "character"){
		stop("Input (lowColor) is of wrong class")
	}
	else if (length(lowColor) != 1){
		stop("Length lowColor must be one")
	}
	else if (class(highColor) != "character"){
		stop("Input (highColor) is of wrong class")
	}
	else if (length(highColor) != 1){
		stop("Length highColor must be one")
	}
	else if (class(textsize) != "numeric"){
			stop("Input (textsize) is of wrong class")
	}
	else if (length(textsize) != 1){
			stop("Length textsize must be one")
	}
	else if (textsize <= 0){
			stop("textsize must be positive")
	} 
	else if (class(diag) != "logical"){
		stop("Input (diag) is of wrong class")
	}
	else if (class(legend) != "logical"){
		stop("Input (legend) is of wrong class")
	}
	else if (class(main) != "character"){
		stop("Input (main) is of wrong class")
	}
	else {
		# Put matrix in data format
		if (nrow(M) == ncol(M) & !diag) {diag(M) <- 0}
		Mmelt    <- melt(M)
		Mmelt$X1 <- factor(as.character(Mmelt$X1), levels = unique(Mmelt$X1), ordered = TRUE)
		Mmelt$X2 <- factor(as.character(Mmelt$X2), levels = unique(Mmelt$X2), ordered = TRUE)
	
		# Visualize
		if (legend){
			ggplot(Mmelt, aes(X2, X1, fill = value)) + geom_tile() + 
         	 	 	 scale_fill_gradient2("", low = lowColor,  mid = "white", high = highColor, midpoint = 0) +
		 	 	 theme(axis.ticks = element_blank()) +
				 theme(axis.text.y = element_text(size = textsize)) +
				 theme(axis.text.x = element_text(angle = -90, vjust = .5, size = textsize)) +
	   	 	 	 xlab(" ") + ylab(" ") +
         		 	 ylim(rev(levels(Mmelt$X1))) +
				 ggtitle(main)
		} else {
			ggplot(Mmelt, aes(X2, X1, fill = value)) + geom_tile() + 
         	 	 	 scale_fill_gradient2("", low = lowColor,  mid = "white", high = highColor, midpoint = 0) +
		 	 	 theme(axis.ticks = element_blank()) +
				 theme(axis.text.y = element_text(size = textsize)) +
				 theme(axis.text.x = element_text(angle = -90, vjust = .5, size = textsize)) +
	   	 	 	 xlab(" ") + ylab(" ") +
         		 	 ylim(rev(levels(Mmelt$X1))) +
				 ggtitle(main) +
			 	 theme(legend.position = "none")
		}
	}
}



Ugraph <- function(M, type = c("plain", "fancy", "weighted"), lay = layout.circle, Vsize = 15, Vcex = 1, 
			 Vcolor = "orangered", VBcolor = "darkred", VLcolor = "black", prune = FALSE, legend = FALSE, 
			 label = "", Lcex = 1.3, PTcex = 4, cut = .5, scale = 10, pEcolor = "black", nEcolor = "grey", 
			 main = ""){
	#####################################################################################################
	# - Function that visualizes the sparsified precision matrix as an undirected graph
	# - Function is partly a wrapper around certain 'igraph' functions
	# - M       > (Possibly sparsified) precision matrix
	# - type    > graph type: 'plain' gives plain undirected graph. 'fancy' gives undirected graph in
	#		  which dashed lines indicate negative partial correlations while solid lines indicate
	#		  positive partial correlations, and in which black lines indicate strong edges. 'weighted'
	#             gives an undirected graph in which edge thickness indicates the strenght of the partial
	#             correlations. Grey lines then indicate negative partial correlations while black lines
	#             represent positive partial correlations.
	# - lay     > determines layout of the graph. All layouts in 'layout{igraph}' are accepted.
	#		  Default = layout.circle, giving circular layout
	# - Vsize   > gives vertex size, default = 15
	# - Vcex    > gives size vertex labels, default = 1
	# - Vcolor  > gives vertex color, default = "orangered"
	# - VBcolor > gives color of the vertex border, default = "darkred"
	# - VLcolor > gives color of the vertex labels, default = "black"
	# - prune   > logical indicating if vertices of degree 0 should be removed
	# - legend  > optional inclusion of color legend, default = FALSE
	# - label   > character label for the endogenous variables, default = ""; only when legend = TRUE
	# - Lcex    > scaling legend box, default = 1.3; only when legend = TRUE
	# - PTcex   > scaling node in legend box, default = 4; only when legend = TRUE
	# - cut     > cut-off for indication of strong edge, default = .5; only when type = "fancy"
	# - scale   > scale factor for visualizing strenght of edges, default = 10; only when type = "weighted"
	# - pEcolor > gives edge color for edges tied to positive precision elements, default = "black"; only
	#             when type = "weighted"
	# - nEcolor > gives edge color for edges tied to negative precision elements, default = "grey"; only
	#             when type = "weighted"
	# - main    > character specifying heading figure, default = "" 
	#####################################################################################################

	# Dependencies
	# require("igraph")
	# require("reshape")

	if (!is.matrix(M)){
		stop("M should be a matrix")
	}
	else if (nrow(M) != ncol(M)){
		stop("M should be square matrix")
	}
	else if (missing(type)){
		stop("Need to specify graph type ('plain' or 'fancy' or 'weighted')")
	}
	else if (!(type %in% c("plain", "fancy", "weighted"))){
		stop("type should be one of {'plain', 'fancy', 'weighted'}")
	}
	else if (class(Vsize) != "numeric"){
		stop("Input (Vsize) is of wrong class")
	}
	else if (length(Vsize) != 1){
		stop("Length Vsize must be one")
	}
	else if (Vsize <= 0){
		stop("Vsize must be positive")
	}
	else if (class(Vcex) != "numeric"){
		stop("Input (Vcex) is of wrong class")
	}
	else if (length(Vcex) != 1){
		stop("Length Vcex must be one")
	}
	else if (Vcex <= 0){
		stop("Vcex must be positive")
	}
	else if (class(Vcolor) != "character"){
		stop("Input (Vcolor) is of wrong class")
	}
	else if (length(Vcolor) != 1){
		stop("Length Vcolor must be one")
	}
	else if (class(VBcolor) != "character"){
		stop("Input (VBcolor) is of wrong class")
	}
	else if (length(VBcolor) != 1){
		stop("Length VBcolor must be one")
	}
	else if (class(VLcolor) != "character"){
		stop("Input (VLcolor) is of wrong class")
	}
	else if (length(VLcolor) != 1){
		stop("Length VLcolor must be one")
	}
	else if (class(prune) != "logical"){
		stop("Input (prune) is of wrong class")
	}
	else if (class(legend) != "logical"){
		stop("Input (legend) is of wrong class")
	}
	else if (class(main) != "character"){
		stop("Input (main) is of wrong class")
	}
	else {
		# Preliminaries
		AM <- adjacentMat(M)
		GA <- graph.adjacency(AM, mode = "undirected")
		if (prune){GA <- delete.vertices(GA, which(degree(GA) < 1))}

		# Plain graph
		if (type == "plain"){
			plot(GA, layout = lay, vertex.size = Vsize, vertex.label.family = "sans", vertex.label.cex = Vcex, 
			     vertex.color = Vcolor, vertex.frame.color = VBcolor, vertex.label.color = VLcolor, main = main)
		}

		# Fancy graph
		if (type == "fancy"){
			if (class(cut) != "numeric"){
				stop("Input (cut) is of wrong class")
			} else if (length(cut) != 1){
				stop("Length cut must be one")
			} else if (cut <= 0){
				stop("cut must be positive")
			} else {
				Names <- colnames(M)
				colnames(M) = rownames(M) <- seq(1, ncol(M), by = 1)
				Mmelt <- melt(M)
				Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
				Mmelt <- Mmelt[Mmelt$value != 0,]
				E(GA)$weight <- Mmelt$value
				E(GA)$color  <- "grey"
				E(GA)[E(GA)$weight < 0]$style <- "dashed"
				E(GA)[E(GA)$weight > 0]$style <- "solid"		
				E(GA)[abs(E(GA)$weight) > cut]$color <- "black"
				plot(GA, layout = lay, vertex.size = Vsize, vertex.label.family = "sans", vertex.label.cex = Vcex, 
				     vertex.color = Vcolor, vertex.frame.color = VBcolor, vertex.label.color = VLcolor, 
				     edge.color = E(GA)$color, edge.lty = E(GA)$style, main = main)
			}
		}

		# Weighted graph
		if (type == "weighted"){
			if (class(scale) != "numeric"){
				stop("Input (scale) is of wrong class")
			} else if (length(scale) != 1){
				stop("Length scale must be one")
			} else if (scale <= 0){
				stop("scale must be positive")
			} else if (class(pEcolor) != "character"){
				stop("Input (pEcolor) is of wrong class")
			} else if (length(pEcolor) != 1){
				stop("Length pEcolor must be one")
			} else if (class(nEcolor) != "character"){
				stop("Input (nEcolor) is of wrong class")
			} else if (length(nEcolor) != 1){
				stop("Length nEcolor must be one")
			} else {
				Names <- colnames(M)
				colnames(M) = rownames(M) <- seq(1, ncol(M), by = 1)
				Mmelt <- melt(M)
				Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
				Mmelt <- Mmelt[Mmelt$value != 0,]
				E(GA)$weight <- Mmelt$value
				E(GA)[E(GA)$weight < 0]$color <- nEcolor
				E(GA)[E(GA)$weight > 0]$color <- pEcolor
				plot(GA, layout = lay, vertex.size = Vsize, vertex.label.family = "sans", vertex.label.cex = Vcex, 
				     vertex.color = Vcolor, vertex.frame.color = VBcolor, vertex.label.color = VLcolor, 
				     edge.color = E(GA)$color, edge.width = scale*abs(E(GA)$weight), main = main)
			}
		}
				
		# Legend
		if (legend){
			if (class(label) != "character"){
				stop("Input (label) is of wrong class")
			} else if (length(label) != 1){
				stop("Length label must be one")
			} else if (class(Lcex) != "numeric"){
				stop("Input (Lcex) is of wrong class")
			} else if (length(Lcex) != 1){
				stop("Length Lcex must be one")
			} else if (Lcex <= 0){
				stop("Lcex must be positive")
			} else if (class(PTcex) != "numeric"){
				stop("Input (PTcex) is of wrong class")
			} else if (length(PTcex) != 1){
				stop("Length PTcex must be one")
			} else if (PTcex <= 0){
				stop("PTcex must be positive")
			} else{
				legend("bottomright", label, pch = 20, col = Vcolor, cex = Lcex, pt.cex = PTcex)
			}
		}
	}
}




##---------------------------------------------------------------------------------------------------------
## 
## Miscellaneous
##
##---------------------------------------------------------------------------------------------------------

.TwoCents <- function(){
	#####################################################################################################
	# - Unsolicited Advice
	#####################################################################################################

	cat("
            ##########################
            ##########################
            Get ridge or die trying. 
                              - 2Cents 
            ##########################
            ########################## \n")
}




##---------------------------------------------------------------------------------------------------------
## 
## NOTES
##
##---------------------------------------------------------------------------------------------------------

## Updates from version 1.1 to 1.2
#- Inclusion function for ML estimation of the sample covariance matrix: 'covML'
#- Inclusion function for approximate leave-one-out cross-validation: 'optPenalty.aLOOCV'
#- Inclusion function 'conditionNumber' to visualize the spectral condition number over the regularization path
#- Inclusion function 'evaluateS' to evaluate basic properties of a covariance matrix
#- Inclusion function 'KLdiv' that calculates the Kullback-Leibler divergence between two normal distributions
#- Inclusion option to suppress on-screen output in 'sparsify' function
#- Corrected small error in 'optPenaltyCV' function
#- Both 'optPenaltyCV' and 'optPenalty.aLOOCV' now utilize 'covML' instead of 'cov'
#- Default output option in 'optPenaltyCV' (as in 'optPenalty.aLOOCV') is now "light"



