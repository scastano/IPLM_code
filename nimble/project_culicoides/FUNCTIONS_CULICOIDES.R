#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

##################################################
## A list of Nimble functions for Case study II ##
##################################################
library(nimble)

###############
## Contents  ##
###############
## nf_get_quantiles
## nf_get_quantiles_vec
## diff1
## diff2
## dUnimodal
## rUnimodal
## condMuScSurv
## rMuScSurv
## dLogUnif
## rLogUnif
## dObMuSig_sub
## rObMuSig
## nfIsNAN
## nf_isFinite
## nf_TW1_MuSig
## nf_TW1_MuSigSurv
## travellingWave_paras2pDevDead
## travellingWave_doubleParas2pDevDead
## dMaturation
## rMaturation
## dGonotrophicCycle.single
## rGonotrophicCycle.single
## dGonotrophicCycle
## rGonotrophicCycle
## twoColumnMatrix2Vector
## twoColumnMatrix2ShortVector
## myRound
## dPseudoBeta
## rPseudoBeta
## matrix2modelValues

########################################################################
## Functions to extract percentiles from one or more beta             ## 
## distributions parameterised via expected value (mu) and scale (sc) ##
########################################################################
nf_get_quantiles <- nimbleFunction(
    ## Returns quantile of beta distribution given expected value and scale
    run = function(muSc = double(1), P=double(0)) { ## mu, sc
        alpha <- numeric(2)
        alpha[1:2] <- nf_muVar2alpha(nf_muSc2muVar(muSc[1:2]))
        quantileGivenMuSc <- qbeta(P, alpha[1], alpha[2])
        ## Return quantile
        ## return(alpha)
        return(quantileGivenMuSc)
        returnType(double(0))
    }
)


nf_get_quantiles_vec <- nimbleFunction(
    ## Returns quantile of beta distribution given expected value and scale
    run = function(muSc = double(2), P=double(0)) { ## mu, sc
        nRow  <- dim(muSc)[1]
        quant <- numeric(nRow)
        for (i in 1:nRow)
            quant[i] <- nf_get_quantiles(muSc = muSc[i,1:2], P=P)
        ## Return vector of quantiles
        return(quant)
        ## return(quantileGivenMuSc)
        returnType(double(1))
    }
)


##############################################
## First and second differences of a vector ##
##############################################
diff1 <- nimbleFunction( 
    ## Returns first differences of a numeric vector
    run = function (x = double(1)) { ## Some vector 
        nMin1  <- length(x)-1
        if (nMin1 < 1)
            stop("vector x too short")
        diff1x <- numeric(nMin1)
        for (i in 1:nMin1)
            diff1x[i] <- x[i+1] - x[i]
        returnType(double(1))
        return(diff1x)
    }
)


diff2 <- nimbleFunction( 
    ## Returns second differences of a numeric vector
    run = function (x = double(1)) { ## Some vector 
        nMin1  <- length(x)-1
        if (nMin1 < 2)
            stop("vector x too short")
        diff1x <- diff1(x)
        diff2x <- diff1(diff1x)            
        ## 
        returnType(double(1))
        return(diff2x)
    }
)

#########################################################################
## A "distribution" to check that a vector has no more than one maxima ##
#########################################################################
dUnimodal <- nimbleFunction(
    ## Density for unimodal constraint
    ## Returns 1 (or log(1)) if satisfied
    ## Returns 0 (or -Inf) if not  
    run = function(x   = double(1), 
                   nx  = integer(0, default=10), ## Sets vector length in rUnimodal.
                   log = integer(0, default = 0)) {
        nx     <- length(x) ## Override nx
        nxMin1 <- nx - 1
        nxMin2 <- nx - 2        
        if (nx < 3) { 
            Prior    <- 1
            logPrior <- 0
        } else { 
            maxNmax    <- 1
            maxNmin    <- 0
            maxNeither <- 1
            diff1x     <- numeric(nxMin1)
            SignDf     <- numeric(nxMin1)
            diff2x     <- numeric(nxMin2)
            dsd        <- numeric(nxMin2)
            maxima     <- numeric(nxMin2)
            minima     <- numeric(nxMin2)
            ## diff2x[1:nxMin2] <- diff2(x)
            diff1x[1:nxMin1]  <- diff1(x)
            for (i in 1:nxMin1)
                SignDf[i]  <- 2 * (diff1x[i] >= 0) - 1 
            dsd[1:nxMin2]     <- diff1(SignDf[1:nxMin1])
            for (i in 1:nxMin2) {
                (maxima[i] <- dsd[i] == -2)
                (minima[i] <- dsd[i] == 2)
            }
            ## Constraint: allow 1 minima or 1 maxima only
            if ( (sum(maxima) > maxNmax) | (sum(minima) > maxNmin) |
                 (sum(minima)+sum(maxima) > maxNeither)) { 
                ## Uni-extrema constraint is violated
                Prior    <- 0
                logPrior <- -Inf
            } else {
                ## Relatively uninformative alternative
                Prior    <- 1
                logPrior <- 0
            }            
        }
        returnType(double(0))
        if (log) {
            return(logPrior)
        } else {
            return(Prior)
        }
    }
)


rUnimodal <- nimbleFunction(
    ## Density for unimodal constraint
    ## Returns 1 (or log(1)) if satisfied
    ## Returns 0 (or -Inf) if not
    run = function(n  = integer(0, default = 1), ## number of random vectors
                   nx = integer(0, default = 10) ## length of desired vector
                   ) {
        if (n != 1)
            stop("Function only defined for n=1")
        vec        <- numeric(nx, 1/nx)
        idMax      <- rcat(1, prob=vec)
        vec[idMax] <- runif(1,0,1)
        if (idMax != 1) {
            for (i in (idMax-1):1) {
                vec[i] <- runif(1, 0, vec[i+1])
            }
        }
        if (idMax != nx) {
            for (i in (idMax+1):nx) {
                vec[i] <- runif(1, 0, vec[i-1])
            }
        }
        returnType(double(1))        
        return(vec)
    }
)

registerDistributions(
    list(dUnimodal    = list(
             BUGSdist = "dUnimodal(nx)",
             range    = c(-Inf, Inf),
             pqAvail  = FALSE,
             types    = c("value=double(1)", ## 'value' denotes the random variable itself
                          "nx=integer(0)"))
         )
)


################################################################################################
## Tests that kernel percentiles and survival are unimodal w.r.t. temperature & all alpha > 1 ##
################################################################################################
condMuScSurv <- nimbleFunction(
    ## Returns a logical indicating if mu, sc and surv satisy the conditions of unimodality and alpha > 1 w.r.t. temperature
    ## Returns TRUE / FALSE
    run = function(paras = double(2)) { ## Entire parameter matrix: columns (mu, sc, surv) vs. rows (ranked experimental temperatures)
        ## Initialise
        nT     <- dim(paras)[1]
        P01    <- numeric(nT)
        P99    <- numeric(nT)
        alphas <- nimMatrix(0, nrow=nT, ncol=2)
        ## Loop on temperatures
        for (tt in 1:nT) {
            alphas[tt,1:2] <- nf_muVar2alpha(nf_muSc2muVar(paras[tt,1:2]))
            P01[tt] <- qbeta(0.01, alphas[tt,1], alphas[tt,2])
            P99[tt] <- qbeta(0.99, alphas[tt,1], alphas[tt,2])
        }           
        ## Identify if conditions of the constraint are meet
        cond <- 1 == (dUnimodal(P01) * dUnimodal(P99) * dUnimodal(paras[,3]) * (min(alphas)>1))
        if (is.na(cond))
            cond <- FALSE
        return(cond)
        returnType(logical(0))
    }
)



##################################################################
## Generates a random parameter set that satisfies condMuScSurv ##
##################################################################
rMuScSurv <- nimbleFunction(
    ## Generates mu, sc and surv which satisy unimodality and alpha > 1 constraints w.r.t. temperature
    ## Returns paras
    run = function(n   = integer(0, default =  1),
                   nT  = integer(0, default = 10)) { ## Number of temperatures.
        if (n != 1)
            stop("Function only defined for n=1")
        ## INITIALISE     
        mu     <- numeric(nT)
        sc     <- numeric(nT)
        surv   <- numeric(nT)
        paras  <- nimMatrix(0, nrow=nT, ncol=3)
        ## RANDOM NUMBER GENERATION        
        constraintMet <- 0
        while (constraintMet == 0) {        
            mu[1:nT]      <- rUnimodal(1, nT)
            sc[1:nT]      <- mu[1:nT] / 2
            surv[1:nT]    <- rUnimodal(1, nT)
            paras[1:nT,1] <- mu[1:nT]
            paras[1:nT,2] <- sc[1:nT]
            paras[1:nT,3] <- surv[1:nT]
            ## CHECK      
            constraintMet <- condMuScSurv(paras=paras)
        }
        returnType(double(2))
        return(paras)        
    }
)

#######################################################
## Uniform distribution transformed to the log scale ##
#######################################################
dLogUnif <- nimbleFunction (
    ## Returns density of x = log(y)
    ##              where y ~ Unif(L,U)
    run = function(x   = double(0),
                   L   = double(0, default=0.5),
                   U   = double(0, default=100.5),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        if ( U > L ) {
            y = exp(x)  
            probX = (y / (U-L)) * (y > L) * (y < U) 
        } else 
            probX = 0
        if (log)
            return(log(probX))
        return(probX)
    }
)


rLogUnif <- nimbleFunction (
    ## Generates y ~ Unif(L,U)
    ## Returns   x = log(y)
    run = function(n = integer(0, default=1),
                   L = double(0, default=0.5),
                   U = double(0, default=100.5)) {
        returnType(double(0))
        if(n != 1) 
            nimPrint("Warning: rLogUnif only allows n = 1; Using n = 1.\n")
        y <- runif(n=1,min=L,max=U)
        x <- log(y)
        return(x)
    }
)


registerDistributions(list(dLogUnif = list(
                               BUGSdist = "dLogUnif(L,U)",
                               discrete = FALSE,
                               types    = c("value=double(0)", "L=double(0)", "U=double(0)"),
                               pqAvail  = FALSE)))


############################################################
## Likelihood of observed mu and sigma given mu and sigma ##
############################################################
dObMuSig_sub <- substitute(nimbleFunction ( 
    ## Returns density of x = observed mu and sigma
    ## Given the sample size and some estimates of the true values for mu and sigma
    ##              where y ~ Unif(L,U)
    run = function(x   = double(1), ## Row matrix with mu and sigma
                   mu  = double(0), ## Estimated or proposed mu
                   sig = double(0), ## Estimated or proposed sig
                   N   = double(0), ## Number of observations
                   log = integer(0, default = 0)) {
        if (is.na(mu) | is.nan(mu) | !nf_isFinite(mu) | mu < 0 |
            is.na(sig) | is.nan(sig) | !nf_isFinite(sig) | sig < 0)
            logDensityX = -Inf
        else {    
            muOb  <- x[1]
            sigOb <- x[2]
            sig2  <- sig^2
            logDensityX <- (-N/2)*log(2*PI*sig2) - N*(sigOb^2+(muOb-mu)^2)/(2*sig2)
        }
        if (log)
            return(logDensityX)
        return(exp(logDensityX))
        returnType(double(0))
    }
), list(PI = pi))
dObMuSig <- eval(dObMuSig_sub)


rObMuSig <- nimbleFunction (
    ## Simulates and "observed" sample estimate of of mean and standard deviation given "true" parameters and a sample size
    run = function(n   = integer(0, default=1),
                   mu  = double(0), ## Estimated or proposed mu
                   sig = double(0), ## Estimated or proposed sig
                   N   = double(0)  ## Number of observations
                   ) {
        if (sig < 0) 
            stop("Error: standard deviation must be positive.")
        x <- numeric(N)
        for (i in 1:N)
            x[i] <- rnorm(1, mean=mu, sd=sig)
        simMuSig    <- numeric(2)
        simMuSig[1] <- mean(x[1:N])
        simMuSig[2] <- sd(x[1:N])        
        return(simMuSig)
        returnType(double(1))
    }
)

registerDistributions(list(dObMuSig = list(
                           BUGSdist = "dObMuSig(mu, sig, N)",
                           discrete = FALSE,
                           types    = c("value=double(1)", "mu=double(0)", "sig=double(0)", "N=double(0)"),
                           pqAvail  = FALSE)))


#############################
## Test if a scalar is NAN ##
#############################
nfIsNAN <- nimbleFunction(
    run = function (x = double(0)) {
        returnType(double(0))
        if (is.na(x) | is.nan(x))
            return(TRUE) else return (FALSE)
    }
)

################################
## Test if a scalar is finite ##
################################
nf_isFinite <- nimbleFunction(
    run = function (x = double(0)) {
        returnType(double(0))
        if (is.na(x) | is.nan(x) | x == Inf | x == -Inf)
            return(FALSE) else return (TRUE)
    }
)

#####################################################################################
## Function to run travelling wave and return moments of waiting time distribution ##
#####################################################################################
nf_TW1_MuSig <- nimbleFunction(
    ## Based on nf_TW1, returns mean and standard deviation of maturation times (given survival) and p(survive) 
    ## STOPPING RULE: Nstill < thresh | iter == iterMax 
    run = function (paras  = double(1),
                    res    = integer(0),
                    thresh = double(0, default=1E-11)) {
        ## nimPrint("TW1: paras = ", paras)
        ## Initialise
        MuSig <- numeric(2)
        if (res < 1) {
            MuSig[1] = Inf 
            MuSig[2] = Inf 
        } else {
            Dim            <- res + 1
            M              <- matrix(0, Dim, Dim)
            M[1:Dim,1:Dim] <- nf_getM(paras=paras[1:3], res=res)
            if (M[1,1] == 0) {
                iterMax <- Dim
            } else {
                iterMax <- res * ceiling(log(thresh) / log(M[1,1]))
            }
            if (nf_isFinite(iterMax)) {
                N       <- numeric(Dim, 0)
                N[1]    <- 1 
                End     <- FALSE
                CumDev  <- numeric(iterMax) ## Cumulative nb. completing each time step
                iter    <- 0      
                while (!End) {    
                    iter          <- iter + 1
                    MN            <- M %*% N 
                    N[]           <- MN[,1] 
                    CumDev[iter]  <- N[Dim] 
                    Nstill        <- sum(N[1:res]) ## Density still developing
                    End           <- (iter >= iterMax) | (Nstill < thresh)
                }   
                ## Probabilities to become fully developed
                pdfDev        <- numeric(iter)
                pdfDevNorm    <- numeric(iter)
                tDevCumlower  <- 0
                for(i in 1:iter) {
                    tDevCumUpper <- CumDev[i]
                    pdfDev[i]    <- tDevCumUpper - tDevCumlower
                    tDevCumlower <- tDevCumUpper
                }
                ## Filter out tiny -ves which can occassionaly arise due to rounding error
                for (i in 1:iter) {
                    if(pdfDev[i] < 0) {
                        if (pdfDev[i] < -1E-11) ## Shouldn't happen. If it does it will warn about rounding error from some serious bug(s).
                            nimPrint("Warning: pdfDev[", i, "] = ", pdfDev[i], "\n")
                        pdfDev[i]  <- 0
                    }
                    ## Normalise the developing density
                    pdfDevNorm[i] <- pdfDev[i] / CumDev[iter]
                }   
                ## Obtain mu and sig
                oneToIter <- matrix(0,1,iter)
                for (i in 1:iter)
                    oneToIter[1,i] <- i
                muMatTime <- (oneToIter[1,1:iter] %*% pdfDevNorm[1:iter])
                MuSig[1]  <- muMatTime[1,1]
                ##  
                difSqVec <- matrix(0,1,iter)
                for (i in 1:iter)
                    difSqVec[1,i] <- (i - muMatTime[1,1])^2 * pdfDevNorm[i]
                sig2MatTime <- sum(difSqVec[1,1:iter])
                MuSig[2] <- sqrt(sig2MatTime)
            } else {         
                ## iterMax is infinite, so mu and sigma should be too
                MuSig[1] = Inf
                MuSig[2] = Inf 
            }
        }
        return(MuSig)
        returnType(double(1))
    }
)


nf_TW1_MuSigSurv <- nimbleFunction(
    ## Based on nf_TW1_MuSig, returns mean and standard deviation of maturation times (given survival) and p(survive) 
    ## STOPPING RULE: Nstill < thresh | iter == iterMax 
    run = function (paras  = double(1),
                    res    = integer(0),
                    thresh = double(0, default=1E-11),
                    maxLifeSpanInStage = double(0, default=3650) ## Used to prevent ridiculous proposals making MCMC ridiculously slow. But should be large relative to the prior to prevent this constraint having non-negligible effects on the estimated posterior.
                    ) {
        MuSigSurv <- numeric(3)
        MinParas  <- min(paras[1], min(paras[2], paras[3]))
        MaxParas  <- max(paras[1], max(paras[2], paras[3]))
        if (res < 1 | is.na(MinParas) | is.nan(MinParas) | is.na(MaxParas) | is.nan(MaxParas) | MinParas <= 0 | MaxParas >= 1 ) {
            MuSigSurv[1] <- Inf 
            MuSigSurv[2] <- Inf 
            MuSigSurv[3] <- Inf 
        } else {
            Dim            <- res + 1
            M              <- matrix(0, Dim, Dim)
            M[1:Dim,1:Dim] <- nf_getM(paras=paras[1:3], res=res)
            if (M[1,1] == 0) {
                iterMax <- Dim
            } else {
                iterMax <- res * ceiling(log(thresh) / log(M[1,1]))
            }
            if (nf_isFinite(iterMax) & iterMax < maxLifeSpanInStage) {
                N       <- numeric(Dim, 0)
                N[1]    <- 1 
                End     <- FALSE
                CumDev  <- numeric(iterMax) ## Cumulative nb. completing each time step
                CumDead <- numeric(iterMax) ## Cumulative nb. completing each time step
                iter    <- 0      
                while (!End) {    
                    iter          <- iter + 1
                    MN            <- M %*% N 
                    N[]           <- MN[,1] 
                    CumDev[iter]  <- N[Dim] 
                    Nstill        <- sum(N[1:res])               ## Density still developing
                    CumDead[iter] <- max(0, 1 - Nstill - N[Dim]) ## Density dead (cumulative)
                    End           <- (iter >= iterMax) | (Nstill < thresh)
                }   
                ## Probabilities to become fully developed
                pdfDev        <- numeric(iter)
                pdfDevNorm    <- numeric(iter)
                tDevCumlower  <- 0
                for(i in 1:iter) {
                    ## Developing
                    tDevCumUpper <- CumDev[i]
                    pdfDev[i]    <- tDevCumUpper - tDevCumlower
                    tDevCumlower <- tDevCumUpper
                }
                ## Filter out tiny -ves which can occassionaly arise due to rounding error
                for (i in 1:iter) {
                    if (pdfDev[i] < 0) {
                        if (pdfDev[i] < -1E-11) ## Shouldn't happen. If it does it will warn about rounding error from some serious bug(s).
                            nimPrint("WARNING: pdfDev[", i, "] = ", pdfDev[i], "\n")
                        pdfDev[i]  <- 0
                    }
                    ## Normalise the developing density
                    pdfDevNorm[i] <- pdfDev[i] / CumDev[iter]
                }   
                ## Obtain mu, sig & surv
                oneToIter <- matrix(0,1,iter)
                for (i in 1:iter)
                    oneToIter[1,i] <- i
                muMatTime    <- (oneToIter[1,1:iter] %*% pdfDevNorm[1:iter])
                difSqVec     <- matrix(0,1,iter)
                for (i in 1:iter)
                    difSqVec[1,i] <- (i - muMatTime[1,1])^2 * pdfDevNorm[i]
                sig2MatTime  <- sum(difSqVec[1,1:iter])
                MuSigSurv[1] <- muMatTime[1,1]
                MuSigSurv[2] <- sqrt(abs(sig2MatTime))
                MuSigSurv[3] <- 1-CumDead[iter]
            } else {         
                ## iterMax is infinite, so mu and sigma should be too
                MuSigSurv[1] <- Inf
                MuSigSurv[2] <- Inf 
                MuSigSurv[3] <- NA
            }
        }
        return(MuSigSurv)
        returnType(double(1))
    }
)  


###############################
## Travelling Wave Functions ##
###############################
travellingWave_paras2pDevDead <- nimbleFunction( ## stop condition : nIter
    run = function (paras = double(1),
                    res   = integer(0),
                    nIter = integer(0)) {
        pdfDevDead <- matrix(nrow = nIter, ncol = 2)
        MinParas   <- min(paras[1], min(paras[2], paras[3]))
        MaxParas12 <- max(paras[1], paras[2])
        MaxParas   <- max(MaxParas12, paras[3])
        if (res < 1 | is.na(MinParas) | is.nan(MinParas) | is.na(MaxParas) | is.nan(MaxParas) | MinParas <= 0 | MaxParas12 >= 1 | paras[3] > 1) {
            ## Perhaps do nothing, just return pdfDevDead as zero matrix
        } else {
            pdfDev   <- numeric(nIter)
            pdfDead  <- numeric(nIter)
            tDevCum  <- numeric(nIter) ## Cumulative number of individuals completing  stage at each time step
            tDeadCum <- tDevCum        ## Cumulative number of individuals dead before completing  stage at each time step
            Dim      <- res + 1
            N        <- numeric(Dim)
            N[1]     <- 1 
            M        <- matrix(0, Dim, Dim)
            M[,]     <- nf_getM(paras=paras[1:3], res=res)
            ##       
            iter <- 0
            End  <- FALSE
            while (!End) {
                iter <- iter + 1
                if (!End) {
                    MN             <- M %*% N          
                    N[]            <- MN[,1]
                    tDevCum[iter]  <- N[Dim]
                    tDeadCum[iter] <- 1 - sum(N)
                    if (iter == nIter)
                        End <- TRUE
                }      
            }          
            ## Probabilities to become fully developed
            tDevCumlower <- 0
            for(i in 1:nIter) { 
                tDevCumUpper <- tDevCum[i]
                pdfDev[i]    <- max(0, tDevCumUpper - tDevCumlower)
                tDevCumlower <- tDevCumUpper
            }
            ## Probabilities of dying
            tDeadCumlower <- 0
            for(i in 1:nIter) {
                tDeadCumUpper <- tDeadCum[i]
                pdfDead[i]    <- max(0, tDeadCumUpper - tDeadCumlower)
                tDeadCumlower <- tDeadCumUpper
            }
            ## Combine vectors and return output
            pdfDevDead[1:nIter,1] <- pdfDev
            pdfDevDead[1:nIter,2] <- pdfDead
        }
        return(pdfDevDead)
        returnType(double(2))
    }
)




travellingWave_doubleParas2pDevDead <- nimbleFunction( ## stop condition : nIter
    ## Adapts travellingWave_paras2pDevDead to apply to two successive stages
    run = function (paras1 = double(1),
                    paras2 = double(1),
                    res1   = integer(0),
                    res2   = integer(0),
                    nIter  = integer(0)) {
        pdfDevDead <- matrix(nrow = nIter, ncol = 2)
        ## 
        MinParas1  <- min(paras1[1], min(paras1[2], paras1[3]))
        MinParas2  <- min(paras2[1], min(paras2[2], paras2[3]))
        MinParas   <- min(MinParas2, MinParas2)        
        MaxParas12 <- max(max(paras1[1], paras1[2]), max(paras2[1], paras2[2]))
        Max3       <- max(paras1[3],paras2[3])
        MaxParas   <- max(MaxParas12, Max3)
        if (min(res1,res2) < 1 | is.na(MinParas) | is.nan(MinParas) | is.na(MaxParas) |
              is.nan(MaxParas) | MinParas <= 0 | MaxParas12 >= 1 | Max3 > 1) {
            ## Do nothing, just return pdfDevDead as zero matrix
        } else {
            pdfDev   <- numeric(nIter)
            pdfDead  <- numeric(nIter)
            tDevCum  <- numeric(nIter) ## Cumulative nb completing stage / time step
            tDeadCum <- tDevCum        ## Cumulative nb dead before completing stage / step
            ## Matrices 
            Dim1     <- res1 + 1
            Dim2     <- res2 + 1
            Dim12    <- res1 + res2 + 1
            N        <- numeric(Dim12)
            N[1]     <- 1 
            M1       <- matrix(0, Dim1, Dim1)
            M1[,]    <- nf_getM(paras=paras1[1:3], res=res1)
            M2       <- matrix(0, Dim2, Dim2)
            M2[,]    <- nf_getM(paras=paras2[1:3], res=res2)
            M        <- matrix(0, Dim12, Dim12)
            M[1:Dim1,1:Dim1]         <- M1[1:Dim1,1:Dim1]
            M[Dim1:Dim12,Dim1:Dim12] <- M2[1:Dim2,1:Dim2]
            iter <- 0
            End  <- FALSE
            while (!End) {
                iter <- iter + 1
                if (!End) {
                    MN             <- M %*% N          
                    N[]            <- MN[,1]
                    tDevCum[iter]  <- N[Dim12]
                    tDeadCum[iter] <- 1 - sum(N)
                    if (iter == nIter)
                        End <- TRUE
                }      
            }          
            ## Probabilities to become fully developed
            tDevCumlower <- 0
            for(i in 1:nIter) { 
                tDevCumUpper <- tDevCum[i]
                pdfDev[i]    <- max(0, tDevCumUpper - tDevCumlower)
                tDevCumlower <- tDevCumUpper
            }
            ## Probabilities of dying
            tDeadCumlower <- 0
            for(i in 1:nIter) {
                tDeadCumUpper <- tDeadCum[i]
                pdfDead[i]    <- max(0, tDeadCumUpper - tDeadCumlower)
                tDeadCumlower <- tDeadCumUpper
            }
            ## Combine vectors and return output
            pdfDevDead[1:nIter,1] <- pdfDev
            pdfDevDead[1:nIter,2] <- pdfDead
        }
        return(pdfDevDead)
        returnType(double(2))
    }
)

## #######################################################
## A DISTRIBUTION TAYLORED TO THE ABOVE MATURATION DATA ##
## #######################################################
dMaturation <- nimbleFunction (
    ## Returns density of cumulative number of maturations with censoring
    run = function(x     = double(1), ## Cumulative number of maturations. NAs (censored intervals) must be represented by -ve values.
                   sSize = double(0), ## Sample size
                   pDev  = double(1), ## Probabilities to mature each day
                   log   = integer(0, default = 0)) { 
        nSteps <- length(x)
        if (nSteps != length(pDev))
            stop("x and pDev must have same length in dMaturation")
        ## ################
        logProbX  <- 0 
        pInterval <- 0 
        pUsed     <- 0 ## A running total on used probability
        xLastOb   <- 0 ## Nb. maturations at previous observation (time zero initially)
        for (Step in 1:nSteps) { 
            pInterval <- pInterval + pDev[Step]
            ## if ( !is.na(x[Step]) ) {
            if ( x[Step] >= 0 ) { ## NOTE: -VE VALUES REPRESENT MISSING DATA
                pUsed     <- pUsed + pInterval
                xDiff     <- x[Step] - xLastOb
                logProbX  <- logProbX + log(pInterval^xDiff) ## ASSUMES x IS MONOTONIC (not including missing data flags)
                ## Note: log(pInterval^xDiff) more robust than xDiff*log(pInterval) [can give NaN when both arguments are zero]
                pInterval <- 0
                xLastOb   <- x[Step]
            } 
        } 
        ## Account for the dead and right-censored
        pDeadOrCensored <- max(0, 1 - pUsed)
        logProbX        <- logProbX + (sSize-xLastOb) * log(pDeadOrCensored)
        if (sSize < xLastOb)
            stop("Sample size must be larger than cumulative number of maturations.")
        ##
        if (log)
            return(logProbX)
        else 
            return(exp(logProbX))
        returnType(double(0))
    }
)


rMaturation <- nimbleFunction (
    ## Generates x ~ dMaturation(n, pDev)
    run = function(n     = integer(0, default=1), 
                   sSize = double(0), ## Sample size
                   pDev  = double(1)  ## Probabilities to mature each day
                   ) { 
        if(n != 1) 
            nimPrint("Warning: rMaturation only allows n = 1; Using n = 1.\n") 
        nSteps <- length(pDev)
        if (1 < sum(pDev[1:nSteps]))
            stop("Total of pDev cannot be greater than 1")
        if (min(pDev[1:nSteps]) < 0)
            stop("Elements of pDev cannot be negative")
        ## 
        z               <- numeric(nSteps + 1) ## Final element for pDead + pCens
        pVec            <- numeric(nSteps + 1)
        pDeadOrCensored <- max(0, 1 - sum(pDev[1:nSteps]))
        pVec[1:nSteps]  <- pDev[1:nSteps]
        pVec[nSteps+1]  <- pDeadOrCensored
        for (ii in 1:sSize) {
            stepMat    <- rcat(n=1, prob=pVec)
            z[stepMat] <- z[stepMat] + 1 
        }
        ## Clip off the dead & censored
        x <- numeric(nSteps)
        x[1:nSteps] <- z[1:nSteps]
        ## Generate Cumulative 
        cx <- numeric(nSteps)
        cx[1] <- x[1]
        for (ii in 2:nSteps) {
            cx[ii] <- cx[ii-1] + x[ii]
        }
        return(cx)
        returnType(double(1))
    }
)

registerDistributions(list(dMaturation  = list(
                               BUGSdist = "dMaturation(sSize, pDev)",
                               discrete = TRUE,
                               types    = c("value=double(1)", "sSize=double(0)", "pDev=double(1)"),
                               pqAvail  = FALSE)))


## ##########################################################################
## LIKELIHOOD OF GONOTROPHIC CYCLE DATA FROM THE INSECTARIUM | TEMPERATURE ##
## ##########################################################################
dGonotrophicCycle.single <- nimbleFunction (
    run = function(x          = double(1), ## Two element vector c(dayInCycle, eggs), where eggs==0 implies death. 
                   pLayDie    = double(2), ## Two column matrix cbind(pLay, pDie)
                   EFecundity = double(0), ## Expected fecundity
                   log        = integer(0, default = 0)) {
        if (nimDim(pLayDie)[2] != 2)
            stop("pLayDie must have two columns.")
        if (nimDim(pLayDie)[1] < x[1])
            stop("pLayDie must have at least as many rows as x[1]).")
        LL   <-  0         
        if (x[2] > 0) {
            ## Case 1: Female oviposited x[1] days into a gonotrophic cycle
            LL <- LL + log(pLayDie[x[1], 1])
            ## Account for egg number 
            LL <- LL + dpois(x[2], lambda = EFecundity, log=TRUE)
        } else { 
            ## Case 2: Female died x[1] days into a gonotrophic cycle
            LL <- LL + log(pLayDie[x[1], 2])
        }
        ## Filter out NAs
        if (is.na(LL))
            LL <- -Inf
        ## Return log-likelihood
        if (log)
            return(LL)
        else 
            return(exp(LL))
        returnType(double(0))
    }
)


rGonotrophicCycle.single <- nimbleFunction (
    run = function(n          = integer(0, default = 1), 
                   pLayDie    = double(2), ## Two column matrix cbind(pLay, pDie)
                   EFecundity = double(0)  ## Expected fecundity
                   ) {
        ## browser()
        if (n != 1)
            stop("n in rGonotrophicCycle currently restricted to one individual")
        nSteps  <- nimDim(pLayDie)[1] ## How many steps was travelling wave
        SumP    <- sum(pLayDie[1:nSteps,1:2])
        SumDead <- sum(pLayDie[1:nSteps,2])
        if (SumP != 1)
            stop("Ensure probabilities in pLayDie sum to one. Try running travelling wave for longer.")
        if (SumDead == 0)
            stop("Ensure not all elements of pLayDie[,2] are zero")
        dayAndEggs <- numeric(2)
        alive      <- TRUE
        nRows      <- 1
        pVec       <- nimNumeric(length=2*nSteps)
        pVec[1:nSteps]              <- pLayDie[1:nSteps,1]
        pVec[(1+nSteps):(2*nSteps)] <- pLayDie[1:nSteps,2]
        dayLayOrDie <- rcat(1,pVec[1:(2*nSteps)]) 
        day <- dayLayOrDie
        if (day > nSteps) {
            day   <- day - nSteps
            alive <- FALSE
        } 
        dayAndEggs[1] <- day
        if (alive) 
            ## Case 1: Finished gonotrophic cycle before dying
            dayAndEggs[2] <- rpois(1, EFecundity)        
        ## Return matrix
        return(dayAndEggs)
        returnType(double(1))
    }
)


registerDistributions(list(dGonotrophicCycle.single = list(
                               BUGSdist = "dGonotrophicCycle.single(pLayDie, EFecundity)",
                               ## discrete = FALSE,
                               types    = c("value=double(1)", "pLayDie=double(2)", "EFecundity=double(0)"),
                               pqAvail  = FALSE)))


## ##########################################################################
## LIKELIHOOD OF GONOTROPHIC CYCLE DATA FROM THE INSECTARIUM | TEMPERATURE ##
## ##########################################################################
dGonotrophicCycle <- nimbleFunction (
    run = function(x          = double(2), ## Two column matrix cbind(dayInCycle, eggs), where eggs==0 implies death. 
                   pLayDie    = double(2), ## Two column matrix cbind(pLay, pDie)
                   EFecundity = double(0), ## Expected fecundity
                   log        = integer(0, default = 0)) {
        nObs <- nimDim(x)[1] ## Total number of observations
        if (nimDim(pLayDie)[2] != 2)
            stop("pLayDie must have two columns.")
        if (nimDim(pLayDie)[1] < max(x[1:nObs,1]))
            stop("pLayDie must have at least as many rows as max(x[,1]).")
        LL   <-  0 
        for (ii in 1:nObs) {
            ## browser()
            if (x[ii,2] > 0) {
                ## Case 1: Female oviposited x[ii,1] days into a gonotrophic cycle
                LL <- LL + log(pLayDie[x[ii,1], 1])
                ## Account for egg number                
                LL <- LL + dpois(x[ii,2], lambda = EFecundity, log=TRUE)
            } else { 
                ## Case 2: Female died x[ii,1] days into a gonotrophic cycle
                LL <- LL + log(pLayDie[x[ii,1], 2])
            }
        }
        ## Filter out NAs
        if (is.na(LL))
            LL <- -Inf
        ## Return log-likelihood
        if (log)
            return(LL)
        else 
            return(exp(LL))
        returnType(double(0))
    }
)


rGonotrophicCycle <- nimbleFunction (
    run = function(n          = integer(0, default = 1), 
                   pLayDie    = double(2), ## Two column matrix cbind(pLay, pDie)
                   EFecundity = double(0)  ## Expected fecundity
                   ) {
        ## browser()
 ## if (n != 1)
 ##     stop("n in rGonotrophicCycle currently restricted to one individual")
        nSteps  <- nimDim(pLayDie)[1] ## How many steps was travelling wave
        SumP    <- sum(pLayDie[1:nSteps,1:2])
        SumDead <- sum(pLayDie[1:nSteps,2])
        if (SumP != 1)
            stop("Ensure probabilities in pLayDie sum to one")
        if (SumDead == 0)
            stop("Ensure not all elements of pLayDie[,2] are zero")
        alive <- TRUE
        dayAndEggs <- nimMatrix(0, nrow=1, ncol=2)
        nRows <- 1
        pVec <- nimNumeric(length=2*nSteps)
        pVec[1:nSteps]              <- pLayDie[1:nSteps,1]
        pVec[(1+nSteps):(2*nSteps)] <- pLayDie[1:nSteps,2]
        ## browser()        
        while (alive) {
            dayLayOrDie <- rcat(1,pVec[1:(2*nSteps)]) 
            day <- dayLayOrDie
            if (day > nSteps) {
                day   <- day - nSteps
                alive <- FALSE
            }
            dayAndEggs[nRows,1] <- day
            if (alive) {
                ## Case 1: Finished gonotrophic cycle before dying
                dayAndEggs[nRows,2] <- rpois(1, EFecundity)
                nRows <- nRows + 1
                temp  <- dayAndEggs
                setSize(dayAndEggs, nRows, 2)
                dayAndEggs[1:(nRows-1),1:2] <- temp
            } 
        }
        ## Return matrix
        return(dayAndEggs)
        returnType(double(2))
    }
)

registerDistributions(list(dGonotrophicCycle = list(
                               BUGSdist = "dGonotrophicCycle(pLayDie, EFecundity)",
                               discrete = TRUE,
                               types    = c("value=double(2)", "pLayDie=double(2)", "EFecundity=double(0)"),
                               pqAvail  = FALSE)))


####################################################################################
## To rearrange 2 column travelling wave output into a vector for use with dcat() ##
####################################################################################
twoColumnMatrix2Vector <- nimbleFunction(
    run = function (mat = double(2)) {
        if (nimDim(mat)[2] != 2)
            stop("twoColumnMatrix2Vector defined for two column matrixes only")
        nRow             <- nimDim(mat)[1]
        nRow1            <- nRow + 1
        nRow2            <- nRow * 2
        nRow21           <- nRow2 + 1 
        vec              <- nimNumeric(nRow21)
        vec[1:nRow]      <- mat[1:nRow,1]
        vec[nRow1:nRow2] <- mat[1:nRow,2]
        vec[nRow21]      <- max(0, 1 - sum(vec[1:nRow2])) ## p(Right Censor)
        return(vec)
        returnType(double(1))
    }
)

twoColumnMatrix2ShortVector <- nimbleFunction(
    run = function (mat = double(2)) {
        if (nimDim(mat)[2] != 2)
            stop("twoColumnMatrix2Vector defined for two column matrixes only")
        nRow        <- nimDim(mat)[1]
        nRow1       <- nRow + 1
        vec         <- nimNumeric(nRow1)
        vec[1:nRow] <- mat[1:nRow,1]
        vec[nRow1]  <- max(0, 1 - sum(vec[1:nRow])) ## p(Dead) + p(Right Censor)
        return(vec)
        returnType(double(1))
    }
)


myRound <- nimbleFunction (
    ## Rounds a scalar to a given level of precision
    run = function(x = double(0), prec=double(0)) {
        roundX <- floor(0.5 + x * 10^prec) / 10^prec
        return(roundX)
        returnType(double(0))
    }
)


#########################################################################################################
## Likelihood function for pseudo data : a surrogate for the likelihood of PostSS | PreSS and survival ##
#########################################################################################################
dPseudoBeta <- nimbleFunction (
    run = function (x      = double(0), ## Pseudo data - should be 1
                    SSPre  = double(0), ## "pre-mortality" sample size 
                    SSPost = double(0), ## "post-mortality" sample size
                    pSurv  = double(0), ## probability to survive developmental process
                    log    = integer(0, default=0)
                    ) {
        if (x != 1) 
            stop("pseudo data x should be set to 1")
        if (SSPost > SSPre)
            stop("Ensrure SSPost <= SSPre")
        LL <- dbeta(pSurv * x, SSPost + 1, SSPre - SSPost + 1, log=TRUE) ## canonical parameterization
        if (log)
            return(LL)
        else
            return (exp(LL))
        returnType(double(0))
    }
)
 
rPseudoBeta <- nimbleFunction (
    run = function (n      = integer(0, default=1), ## Pseudo data - should be 1
                    SSPre  = double(0), ## "pre-mortality" sample size 
                    SSPost = double(0), ## "post-mortality" sample size
                    pSurv  = double(0)  ## probability to survive developmental process
                    ) {
        if (n != 1)
            stop("pseudo data x should be set to 1")
        if (SSPost > SSPre)
            stop("Ensrure SSPost <= SSPre")
        x <- 1
        return (x)
        returnType(double(0))
    }
)


registerDistributions(list(dPseudoBeta = list(
                              BUGSdist = "dPseudoBeta(SSPre, SSPost, pSurv)",
                              discrete = TRUE,
                              types    = c("value=double(0)", "SSPre=double(0)", "SSPost=double(0)", "pSurv=double(0)"),
                              pqAvail  = FALSE)))

#################################################################
## Function to copy values in a matrix to a modelValues object ##
#################################################################
matrix2modelValues <- nimbleFunction(
    setup = function (model, ModelValues, targetNode) { } ,
    run   = function (aMatrix = double(2), indexVector = double(1)) { 
        nrMatrix <- dim(aMatrix)[1]     
        len      <- dim(indexVector)[1] 
        if (nrMatrix != getsize(ModelValues))
            resize(ModelValues, nrMatrix)
        for (ii in 1:nrMatrix) {
            for (jj in 1:len) {    
                ModelValues[targetNode,ii][indexVector[jj]] <<- aMatrix[ii,indexVector[jj]] 
            }
        }
    } 
)

