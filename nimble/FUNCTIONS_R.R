#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

##################################################
## A list of R functions for case study I & II  ##
##################################################

##############
## Contents ##
##############
## muSc2muVar
## muVar2alpha
## muSc2alpha
## getM
## TW1
## TW2
## R_dJSMD
## rConditBetaCopula
## dConditBetaCopula
## testCopulaFunctions
## medianMatrix
## SUCH
## DoubleSUCH
## priorUniformU
## sumSqError 
## IniValAlpha 
## subsetOutput
## tailSamples2subsetList
## mProb


#####################
## Transformations ##
#####################
muSc2muVar <- function(muSc) {
    ## Convert mean and variance of beta distribution to canonical parameters
    muSc <- matrix(muSc, ncol=2)
    mu   <- muSc[,1] 
    sc   <- muSc[,2] 
    if (any(sc > 1 | sc <= 0)) {
        stop("Ensure scale is in (0, 1)")
    }
    var <- mu * (1-mu) * sc
    cbind(mu=mu, var=var)
}

muVar2alpha <- function(muVar) {
    ## Convert mean and variance of beta distribution to canonical parameters
    muVar <- matrix(muVar, ncol=2)
    mu    <- muVar[,1]
    var   <- muVar[,2]
    if (any(var > mu*(1-mu)))
        stop("Ensure var <= mu*(1-mu)")
    denom  <- mu * (1 - mu) / var 
    alpha1 <- mu * (denom - 1)
    alpha2 <- denom - 1 - alpha1
    cbind(alpha1=alpha1, alpha2=alpha2)
}

muSc2alpha <- function(muSc) {
    muVar2alpha(muSc2muVar(muSc))
}


###########################################
## Generating and using projectin matrix ##
getM <- function(paras, Dim, family) {
    ## Create a vector to store development level per class
    (res <- Dim-1)
    (dev <- 0:res / res) 
    ## Initialise the matrix
    M  <- diag(0, Dim)  ## creates a matrix of zeros with zeros on the diagonal
    for(i in 1:Dim) M[i,i] <- 1
    ## Parameterise the first column of M (pdfK)
    if (family=="truncNorm") { ## First column of M defined via a truncated (at dev=0) Gaussian
        ## Obtain parameters from input
        muK  = paras[1]
        sigK = paras[2]
        ## Evaluate distribution function
        pdfK <- diff( ptruncnorm(dev, a=0, b=Inf, m=muK, sd=sigK) ) ## The development vector.
    } else if (family=="exp") { ## Define first column of M with exponential CDF
        ## Obtain parameters from input
        muK  = paras[1]
        ## Evaluate distribution function
        pdfK <- diff(pexp(dev, rate=1/muK))
    } else if (family=="beta") { ## Define first column of M with beta distribution
        ## Obtain parameters from input
        mu     <- paras[1]
        sc     <- paras[2]
        Alphas <- muSc2alpha(c(mu, sc))
        ## Evaluate distribution function
        pdfK <- diff(pbeta(dev, Alphas[1], Alphas[2]))
    } 
    ## Parameterise all other Columns 
    ## Generate a lower-triangular matrix with near-identical columns    
    for (j in 1:res) { ## Loop on columns
        M[j:res,j] <- pdfK[1:(Dim-j)]
    }
    ## Columns must sum to 1, so define final row to comply with this constraint
    if (res > 1) {
        for (i in 1:res) {
            if (pdfK[Dim-i] != 0)
                M[Dim, i] <- max(0, 1-sum(M[1:res, i]))
        }            
    } else { ## res==1, Dim==2
        M[2,1] <- max(0, 1-M[1,1])
    }
    ## Assign pdfK to M's attributes
    attributes(M) <- list(dim=attributes(M)$dim, pdfK=as.numeric(pdfK))
    M
}  



TW1 <- function (M, sur=1, Dim, iterMax=365, thresh=1E-6, verbose=FALSE) {## stop condition : threshold
    ## Initialise population vector
    Nini  <- c(1, rep(0, Dim-1))  
    ## Initialise population
    N <- Nini
    sur <- c(rep(sur, Dim-1), 1)
    ## Looping using M  
    End <- FALSE
    tDevCum  <- rep(0, iterMax) ## Cumulative number of individuals completing E stage at each time step
    tDeadCum <- tDevCum         ## Cumulative number of individuals dead before completing E stage at each time step
    iter <- 0
    while (!End) {
        (iter <- iter + 1)
        if (!End) {
            (N <- sur * (M %*% N))
            (tDevCum[iter]  <- N[Dim])
            (tDeadCum[iter] <- 1 - sum(N))
            (Nstill <- sum(N[-Dim]))  ## Everybody who's still developing
            (Condition <- Nstill < thresh) ##  | iter >= iterMax))
            if (Condition) {
                End <- TRUE
                maxIter <- iter
                if (maxIter > iterMax) 
                    warning ("population development too long, check parameters")
                tDevCum  <- tDevCum[1:iter]    ## truncate output vector
                tDeadCum <- tDeadCum[1:iter]   ## truncate output vector                    
            }
        }  
    }      
    ## Probabilities to become fully developed in each time step given the kernel
    (pdfDev  <- diff(c(0, tDevCum)))
    (pdfDead <- diff(c(0, tDeadCum)))
    ## Filter out potential tiny -ves which can occassionaly arise due to rounding error
    pdfDev[pdfDev < 0]   <- 0
    pdfDead[pdfDead < 0] <- 0
    ## Combine vectors and return output
    pdf <- cbind(pdfDev, pdfDead)
    colnames(pdf) <- c("pDev", "pDead")
    return(pdf)
}



    
TW2 <- function (M, sur, Dim, nIter, thresh = 1E-20, verbose=FALSE) {
    ## Initialise population vector
    Nini  <- c(1, rep(0, Dim-1))  ## For egg only traveling wave
    ## Initialise population
    N <- Nini
    sur <- c(rep(sur, Dim-1), 1)
    ## Looping using M  
    End <- FALSE
    tDevCum  <- numeric()  ## Cumulative number of individuals completing  stage at each time step
    tDeadCum <- tDevCum    ## Cumulative number of individuals dead before completing  stage at each time step
    iter <- 0
    while (!End) {
        (iter <- iter + 1)
        if (!End) {
            (N <- sur * (M %*% N))
            (tDevCum[iter]  <- N[Dim])
            (tDeadCum[iter] <- 1 - sum(N))
            (Condition <- iter >= nIter)
            if (Condition) {
                End <- TRUE
            } 
        }     
    }   
    ## Probabilities to become fully developed in each time step given the kernel
    (pdfDev  <- diff(c(0, tDevCum)))
    (pdfDead <- diff(c(0, tDeadCum)))
    ## Filter out potential tiny -ves which can occassionaly arise due to rounding error
    pdfDev[pdfDev < 0]   <- 0
    pdfDead[pdfDead < 0] <- 0
    ## Combine vectors and return output
    pdf <- cbind(pdfDev, pdfDead)
    colnames(pdf) <- c("pDev", "pDead")
    return(pdf)
}


###############################################################################
## R functions for calculating the conditional density of a bivariate copula ##
###############################################################################

## This function takes as input
## q: a vector of qualities drawn from a uniform(0,1), assumed to the be marginal distribution of the first dimension
## rho: the correlation of the standard normals
## betaParams: a vector of two parameters for the marginal beta distribution of the second dimension
##
## A Gaussian copula is used to correlate the values of q and beta2
##
## It returns as output a vector of beta2 values whose marginal (integrating over quality)
##           is the beta distribution with parameters (inputs) betaParams (given as shape1, shape2)
R_rConditBetaCopula <- function(n, q, rho, betaParams) {
    ## Generates beta2 (a developmental increment) conditionally on quality q
    ## x1 and x2 will refer to the standard normal values
    x1 <- qnorm(q) ## set x1 to the qth quantile of a standard normal
    conditionalMean2 <- rho * x1 ## mean of x2 given x1
    conditionalVar2 <- 1-rho*rho ## variance of x2 given x1
    x2 <- rnorm(n, conditionalMean2, sqrt(conditionalVar2))
    u2 <- pnorm(x2) ## marginal cdf of x2
    beta2 <- qbeta(u2, betaParams[1], betaParams[2])
    beta2
}

R_dConditBetaCopula <- function(x, q, rho, betaParams, log = FALSE) {
    ## Returns the conditional density of x (WHO'S MARGINAL IS BETA WITH BETAPARAMS)
    ## given "quality" q (who's marginal is uniform(0,1))
    ## and a Gaussian copula with correlation rho    
    x1 <- qnorm(q)               ## Set x1 to the qth quantile of a standard normal
    conditionalMean2 <- rho * x1 ## Mean of x2 given x1
    conditionalVar2 <- 1-rho*rho ## Variance of x2 given x1
    conditionalSD2 <- sqrt(conditionalVar2)
    x2 <- qnorm(pbeta(x, betaParams[1], betaParams[2]))
    logDensity <- dnorm(x2, conditionalMean2, conditionalSD2, log = TRUE) +
        dbeta(x, betaParams[1], betaParams[2], log = TRUE) - dnorm(x2, log = TRUE)
    if(log) return(logDensity)
    return(exp(logDensity))
}

testCopulaFunctions <- function(q, rho, betaParams = c(3, 2), n = 10000) {
    ## this simulates a sample, makes a histogram, and plots the corresponding density
    ## If the density and the histogram agree, then they seem to be correct
    rbSample <- R_rConditBetaCopula(n, q, rho, betaParams)
    xAxisGrid <- seq(0, 1, by = 0.001)
    conditionalDensity <- R_dConditBetaCopula(xAxisGrid, q, rho, betaParams)
    hist(rbSample, prob = TRUE, xlim=c(0,1))
    points(xAxisGrid, conditionalDensity, type = 'l')
}




## R function to generate random draws from joint sojourn~mortality distribution
R_dJSMD <- function(x, paras,Dim,log) {
    if ( min(paras) <= 0 | max(paras) >= 1) {
        if (log) 
            return(-Inf)
        return (0)
    } else {
        M      <- nf_getM(paras[1:3], Dim=Dim)
        nIter  <- x[1]
        (probTw <- nf_TW2(Dim=Dim, M=M, nIter=nIter))
        (prob   <- probTw[x[1],x[2]])
        if (log) {
            return(log(prob))
        } else {
            return (prob)
        }
    }
}

### ############################################################################
## FUNCTION medianMatrix & SUCH are used to track mcmc outputs in trackData.R ##
##  ############################################################################
## arguments have dim(nLines, 365)

## medianMatrix provides matrix dim(365,3), giving the median of mu, sc, sur (cols) for every day(row)
medianMatrix <- function(mat1, mat2,mat3) {
    medianVec <-  apply ( cbind(mat1,mat2,mat3), 2, function(x) median(x) )
    ## re-arrange
    l <- ncol(mat1)
    ## create matrix of NAs
    medianMat <- matrix(,l,3) ## 3 is the number of arguments, 1 by parameter by stage
    for (i in 1:l)
    medianMat[i, ] <- rbind( medianVec[i], medianVec[i+l],medianVec[i+2*l] )
    return(medianMat)
}


SUCH_ORI <- function(x, y) {
    ## Stochastic Unimodal Cubic Hermite spline 
    require(stats)
    #####################
    ## Test order in x ##
    #####################
    if (any(order(x) != 1:length(x))) {
        stop("x must be orderred in SUCH.")
    }
    ###########################
    ## Test unimodality in y ##
    ###########################
    if (sum((diff(sign(diff(y))))==-2) > 1) {
        stop("Unimodality constraint violated by data input into SUCH.")
    }
    ## 
    (n      <- length(x))          ## Number of data points
    (delta  <- diff(y)/diff(x))    ## Slope of secant lines
    (m      <- vector(length=n))   ## Slope at each data point
    (alpha <- 3*runif(n))          ## Random initial parameter values 
    (beta  <- 3*runif(n))          ## Random initial parameter values 
    (m[1]  <- alpha[1]*delta[1])   ## Slope at first data point
    (peak  <- which(y==max(y)))    ## 
    if (is.element(peak, c(1,n)))  ## If peak is at start or end of vector then assume monotonic trend (no peak) beyond data range
        peak <- -99
    for (int in 2:n) { 
        if (int == peak) {
            ## At mode sample angle uniformly between the two neighbouring secants
            (m[int] <- (tan(runif(1, atan(delta[int]), atan(delta[int-1])))))
        } else {
            if (int == n) {
                ## At final point do as we did at the first point
                (m[int] <- beta[int-1]*delta[int-1])
            } else {
                ## Set slope based on alpha[int]
                (m[int] <- alpha[int]*delta[int])
            }
        }
        ## Derive beta | m
        (beta[int-1] <- m[int] / delta[int-1])
        if (m[int] == 0 & delta[int-1] == 0)
            ## This can occur due ot rounding error on rare occassions
            (beta[int-1] <- 0)
        ## Verify constraint
        if (beta[int-1]^2 > 9 - alpha[int-1]^2 & int != peak) {
            ## If constraint is violated...
            ## ... re-sample beta given alpha 
            (beta[int-1] <- (runif(1, 0, sqrt(9-alpha[int-1]^2))))
            (m[int]      <- beta[int-1]*delta[int-1])
            ## ... and re-set alpha to satisfy constraint 
            (alpha[int]  <- m[int]/delta[int])
            ## Re-verify constraint
            if (beta[int-1]^2 > 9 - alpha[int-1]^2) {
                ## Should not happen
                browser()
            }
        }
    }
    ## Return a function for the fitted spline 
    fittedSplineFunction <- splinefunH(x, y, m=m) 
    ## curve(fittedSplineFunction, add = TRUE, col = rgb(0,0,1,0.1), n = 1001)
    fittedSplineFunction
}


SUCH <- function(x, y) { 
    ## Stochastic Unimodal Cubic Hermite spline 
    require(stats) 
    ## Test order in x 
    if (any(order(x) != 1:length(x))) 
        stop("x must be orderred in SUCH2.") 
    ## Test unimodality in y 
    if (sum((diff(sign(diff(y))))==-2) > 1) 
        stop("Unimodality constraint violated by data input into SUCH2.") 
    ## Set constants
    (n     <- length(x))          ## Number of data points
    (n1    <- n-1)                
    (delta <- diff(y)/diff(x))    ## Slope of secant lines
    ## Identify location of max(y)
    (iMode <- which(y==max(y)))
    if (is.element(iMode, c(1,n))) ## If iMode is at start or end of vector then assume monotonic trend (no iMode) beyond data range
        iMode <- 2 * n * (-1)^(iMode==1)
    ## Set bounds for slice sampling m 
    LBm <- UBm <- 3 * c(delta[1], delta)
    LBm[LBm > 0] <- 0
    UBm[UBm < 0] <- 0
    if (1 < iMode & iMode < n) {
        LBm[iMode] <- delta[iMode]
        UBm[iMode] <- delta[iMode-1]
    } 
    ## Transform to polar angles
    (atanLBm <- atan(LBm))
    (atanUBm <- atan(UBm))
    ## Initial gradients, alphas and betas
    (m     <- tan(runif(n, atanLBm, atanUBm)))
    (alpha <- m[1:n1] / delta)
    (beta  <- m[2:n]  / delta)
    ## Reset gradients via rejection sampling
    (a2b2 <- alpha^2 + beta^2)
    while(any(a2b2 > 9)) {
        ## Identify which interval with overshoot to modify
        (kk <- which(a2b2==max(a2b2))) ## a2b2[kk]
        ## Identify which end of interval kk to modify 
        if (abs(beta[kk]) > abs(alpha[kk]))
            kk <- kk+1 
        ## Stepping in 
        if (kk < iMode) {
            atanUBm[kk] <- atan(m[kk])
        } else if (kk > iMode) {
            atanLBm[kk] <- atan(m[kk])
        } else {       
            if (m[kk] < 0) {
                atanLBm[kk] <- atan(m[kk])
            } else {   
                atanUBm[kk] <- atan(m[kk])
            }          
        }              
        ## Sample      
        (m[kk] <- tan(runif(1, atanLBm[kk], atanUBm[kk])))
        (alpha <- m[1:n1] / delta)
        (beta  <- m[2:n]  / delta)
        ## Reset gradients via rejection sampling
        (a2b2 <- alpha^2 + beta^2) 
    }
    ## Return a function for the fitted spline 
    fittedSplineFunction <- splinefunH(x, y, m=m) 
    ## curve(fittedSplineFunction, add = TRUE, col = rgb(0,0,1,0.1), n = 1001)
    fittedSplineFunction
}


DoubleSUCH <- function(fOri, gOri, TempsOri, RangeTempsNew) {
    ## Fit SUCH to two data sets, reject proposals where curves cross and return the two functions in a list
    ## fOri          - value of q01 at each observed temperature
    ## gOri          - value of q99 at each observed temperature
    ## TempsOri      - each observed temperature
    ## RangeTempsNew - range over which to extrapolate
    if (any(gOri <= fOri ))
        stop("Ensure gOri > fOri")
    ## Simple Rejection Sampler
    ConditionSatisfied <- FALSE
    while (!ConditionSatisfied) {
        ## ##############################
        ## Apply SUCH to fOri and gOri ## 
        ## ##############################
        fSUCH <- SUCH(TempsOri, fOri)
        gSUCH <- SUCH(TempsOri, gOri)
        ## ##########################################################################
        ## Next, interpolate/extrapolate at finer resolution using Kruger's method ##
        ## ##########################################################################
        (Temps4extrapolation <- seq(RangeTempsNew[1], RangeTempsNew[2], by=1))
        (fExtrapolation      <- fSUCH(Temps4extrapolation))
        (gExtrapolation      <- gSUCH(Temps4extrapolation))
        ## plot(Temps4extrapolation, fExtrapolation, typ="l")
        ## lines(Temps4extrapolation, gExtrapolation)
        ## browser()
        ## Test Condition
        if (all(gExtrapolation > fExtrapolation))## if (q99 > q01)
            ConditionSatisfied <- TRUE
    }
    list(fSUCH, gSUCH)
}



priorUniformU <- function(xPS, maxNmax=1, maxNmin=0, maxNeither=1) {  ## nInt, diffTemps
    ## Uniform priors with Uni-Exterma shape constraint
    ## xPS       - Vector of parameters at different temperatures on parameter scale
    ##           - e.g (0, 1) when family == "beta"
    if (all(is.finite(xPS))) {
        (diffXPS <- diff(xPS))
        (dsd     <- diff(sign(diffXPS)))        
        (maxima  <- dsd == -2)
        (minima  <- dsd == 2)
        if ( (sum(maxima) > maxNmax) | (sum(minima) > maxNmin) | (sum(minima+maxima) > maxNeither)) { ## A more flexible constraint, allows 1 minima or maxima
            ## Uni-extrema constraint is violated
            logPrior <- -Inf
        } else {
            ## Relatively uninformative alternative
            logPrior <- 0
        }
    } else {
        logPrior <- -Inf
    }
    ## Return log of prior weight
    logPrior
}
priorUniformU_AllParas <- function(ParasPS, res) { ## nInt, diffTemps
    nParas <- ncol(ParasPS)
    if (res!=2) {
        Evaluate <- 1:nParas                ## Evaluate prior for all parameters
    } else {
        if (family == "beta") {
            Evaluate <- c(1,3)              ## Skip the second parameter
        } else if (family == "truncExp") {
            Evaluate <- 2:nParas            ## Skip the first parameter
        }
    }
    LogPriorWeight <- 0
    for (ii in Evaluate)        
        LogPriorWeight <- LogPriorWeight + priorUniformU(ParasPS[,ii])
    ## Return log of prior weight
    return(as.numeric(LogPriorWeight))
}




## ################################################################
## Functions for Nelder-Mead estimation of alpha given quantiles ##
## ################################################################

sumSqError <- function(alphas, target, p=c(0.01, 0.99)) {
    ## Error function for estimating alphas from extreme percentiles
    if (any(alphas<1))
        return(-Inf)
    hat <- qbeta(p=p, alphas[1], alphas[2])
    as.numeric(t(hat - target) %*% (hat - target))
}

IniValAlpha <- function(p, qn) {
    ## Determines some reasonable initial values for alpha given qn and p
    (muIni <- mean(qn))
    (kappaUB <- min(c(muIni/(1+muIni), (1-muIni)/(2-muIni))))
    (vIni <- 0.00001 * muIni * (1-muIni))
    (aIni <- muVar2alpha(c(muIni, vIni)))
}


## #####################################
## Functions for interpolation script ##
## #####################################
subsetOutput <- function(output, parnames) {
    ## Generates a list for MCMC output 
    ## output   : MCMC output matrix with colNames
    ## parnames : names of parameters to include in the list
    namesCol <- colnames(output)
    l        <- length(parnames)
    subset   <- list()
    indx     <- sapply(parnames, function(x, namesCol) {grep(x, namesCol, ignore.case = TRUE)}, namesCol=namesCol)
    for (i in 1:l)
        subset[[i]] <- output[,indx[[i]]]
    return(subset) 
}

tailSamples2subsetList <- function(samples, nLines, parNames) {
    ## samples a matrix
    ## This function applies uses tail to obtain nLines of output and returns a subset defined by parNames
    if (nrow(samples) < nLines) {
        samples <- tail(samples, nrow(samples) - round(nrow(samples)/2))  ## Ad-hoc burn-in for working with running chains
        warning(paste("Target number of samples not obtained in", mcmcDir))
    } else {
        ## Apply tail to completed chain
        samples <- tail(samples, nLines)
    }
    ## Return subset of data as a list
    subsets <- subsetOutput(output=samples, parnames=parNames) 
}



## ###############################################################
## Functions for  taking model probability based on mcmc outpus ##
## ###############################################################

## Function that takes a list composed of matrices (every matrix corresponds to a IPLM model)
## and returns the model probabilities (used in project_culicoides)
mProb <- function(L) { ## L : a list of matrices
    nLines  <- nrow(L[[1]]) 
    nbRes   <- length(L) ## nb of models/resolutions to consider
    if (nbRes==1)
        stop("L must have more than one element to calculate mProb")
    logliks <- matrix(0,nLines,nbRes) 
    for (j in 1:nbRes) 
        logliks[,j] <- rowSums(L[[j]]) ## sum logliks of samples2 by line
    ## Take model probability
    mProb <- t(apply(logliks, 1, function(x) {
        x <- x - max(x)
        x <- exp(x)
        x <- x / sum(x)
        x
    }))
}
