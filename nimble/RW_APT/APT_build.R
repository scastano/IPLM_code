##############################################################################################
## DRJP: an adaptation of buildMCMC that generates an adaptive parallel tempering algorithm ##
##############################################################################################

#' Create an MCMC function, from an MCMCconf object
#'
#' Accepts a single required argument, which may be of class MCMCconf, or inherit from class modelBaseClass (a NIMBLE model object).  Returns an MCMC function; see details section.
#'
#' @param conf An object of class MCMCconf that specifies the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{configureMCMC} for details of creating MCMCconf objects.  Alternatively, \code{MCMCconf} may a NIMBLE model object, in which case an MCMC function corresponding to the default MCMC configuration for this model is returned.
#' @param ... Additional arguments to be passed to \code{configureMCMC} if \code{conf} is a NIMBLE model object
#'
#' @author Daniel Turek
#' @export
#' @details
#' Calling buildMCMC(conf) will produce an uncompiled (R) R mcmc function object, say 'Rmcmc'.
#'
#' The uncompiled MCMC function will have arguments:
#'
#' \code{niter}: The number of iterations to run the MCMC.
#'
#' \code{reset}: Boolean specifying whether to reset the model and stored samples.  This will simulate into any stochastic nodes with value NA, propagate values through any deterministic nodes, and calculate all model probabilities. This will also reset the internal stored MCMC samples. Specifying \code{reset=FALSE} allows the MCMC algorithm to continue running from where it left off. Generally, \code{reset=FALSE} should only be used when the MCMC has already been run (default = TRUE).
#'
#' \code{simulateAll}: Boolean specifying whether to simulate into all stochastic nodes.  This will overwrite the current values in all stochastic nodes (default = FALSE).
#'
#' \code{time}: Boolean specifying whether to record runtimes of the individual internal MCMC samplers.  When \code{time=TRUE}, a vector of runtimes (measured in seconds) can be extracted from the MCMC using the method \code{mcmc$getTimes()} (default = FALSE).
#'
#' \code{progressBar}: Boolean specifying whether to display a progress bar during MCMC execution (default = TRUE).
#'
#' Samples corresponding to the \code{monitors} and \code{monitors2} from the MCMCconf are stored into the interval variables \code{mvSamples} and \code{mvSamples2}, respectively.
#' These may be accessed and converted into R matrix objects via:
#' \code{as.matrix(mcmc$mvSamples)}
#' \code{as.matrix(mcmc$mvSamples2)}
#'
#' The uncompiled (R) MCMC function may be compiled to a compiled MCMC object, taking care to compile in the same project as the R model object, using:
#' \code{Cmcmc <- compileNimble(Rmcmc, project=Rmodel)}
#'
#' The compiled function will function identically to the uncompiled object, except acting on the compiled model object.
#' 
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(code)
#' conf <- configureMCMC(Rmodel)
#' Rmcmc <- buildMCMC(conf)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#' Cmcmc$run(10000)
#' samples <- as.matrix(Cmcmc$mvSamples)
#' head(samples)
#' }

## sampler_APT <- nimbleFunctionVirtual(
##     methods = list(
##         reset = function() { }
##     )
## )


#' @rdname samplers
#' @export
sampler_APT <- nimbleFunctionVirtual(
    ## run = function(temperture=double(0, default=1)) {}, 
    methods = list(
        setTemp           = function(temp = double()) {},
        ## turnOffAdaptation = function() {},
        reset             = function() {}
    )
)


buildAPT <- nimbleFunction(
    setup = function(conf,                ## As for buildMCMC 
                     Temps,               ## Vector of temperatures. Typically, lowest temperature should be 1.
                     monitorTmax=TRUE,    ## Logical, save MCMC output for Tmax too.
                     ULT=1E6,             ## Scalar, Upper Limit on Temperatures
                     thinPrintTemps=1,    ## Nb. iterations between prints of temps when adaptTemps==TRUE
                     ...) {
        if(inherits(conf, 'modelBaseClass'))
            conf <- configureMCMC(conf, ...)
        else if(!inherits(conf, 'MCMCconf'))
            stop('conf must either be a nimbleModel or
                  a MCMCconf object (created by configureMCMC(...) )')
        if (missing(ULT)) {
            ULT <- 1E6
            nimPrint("ULT set to ", ULT)
        }
        model              <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved            <- modelValues(model) ## Used to restore model following rejection in MCMC
        nSamplersPerT      <- length(seq_along(conf$samplerConfs)) ## Nb. samplers / temp'
        samplerFunctions   <- nimbleFunctionList(sampler_APT)
        ## Tempering related objects
        Temps        <- as.numeric(sort(Temps)) ## Ensures compiler knows Temps is not integer
        nTemps       <- length(Temps)           ## Number of temperatures for tempering
        invTemps     <- 1/Temps[1:nTemps]
        TempsCurrent <- Temps[1:nTemps]
        TempsOri     <- Temps[1:nTemps]
        nTemps_1     <- nTemps - 1        
        mvTemps      <- modelValues(model, nTemps) ## Stores state of MCMC at each temperature
        tempTraj     <- nimMatrix(0, nrow=1000, ncol=nTemps) ## Stores trajectories of temperatures
        logProbTemps <- numeric(nTemps)
        pSwapMatrix  <- nimMatrix(0, nTemps, nTemps)
        temporary    <- numeric(nTemps)
        accCountSwap <- nimMatrix(0, nTemps, nTemps)
        nimPrint("Initial temperatures:", Temps,"\n")
        ## 
        for (tt in 1:nTemps) {
            for(ss in 1:nSamplersPerT) { 
                ist <- ss + (tt-1) * nSamplersPerT ## index for Sampler & Temperature
                samplerFunctions[[ist]] <-
                    conf$samplerConfs[[ss]]$buildSampler(model=model, mvSaved=mvSaved)
                samplerFunctions[[ist]]$setTemp(temp=Temps[tt])
            }
        }
        ## nSamplers <- nSamplersPerT * nTemps
        ##
        thin              <- conf$thin
        thin2             <- conf$thin2
        mvSamplesConf     <- conf$getMvSamplesConf(1)
        mvSamples2Conf    <- conf$getMvSamplesConf(2)
        monitors          <- processMonitorNames(model, conf$monitors)
        monitors2         <- processMonitorNames(model, conf$monitors2)
        mvSamples         <- modelValues(mvSamplesConf)   ## For storing MCMC output (T=1)
        mvSamples2        <- modelValues(mvSamples2Conf)  ## For storing MCMC output (T=Tmax)
        if (monitorTmax==TRUE) {
            mvSamplesTmax  <- modelValues(mvSamplesConf)  ## For MCMC output (T=Tmax) 
            mvSamples2Tmax <- modelValues(mvSamples2Conf) ## For MCMC output (T=Tmax) too
        }
        samplerTimes      <- c(0,0) ## Establish as a vector
        progressBarLength <- 52     ## Multiples of 4 only
        resetAnyway       <- TRUE
        totalIters        <- 1
    },
    #################################################
    run = function(niter          = integer(),
                   reset          = logical(default=TRUE),
                   resetTempering = logical(default=FALSE),
                   simulateAll    = logical(default=FALSE),
                   time           = logical(default=FALSE),
    ## adaptiveOff = logical(default=FALSE), # Never worked. Boost iteration counter instead?
                   adaptTemps     = logical(default=TRUE),
                   printTemps     = double(default=FALSE),
                   tuneTemper1    = double(default=10), 
                   tuneTemper2    = double(default=1), 
                   progressBar    = logical(default=TRUE)) {
        if(simulateAll)     simulate(model)    ## Default behavior excludes data nodes
        if(resetAnyway==TRUE) { ## Force reset to be TRUE on first usage to avoid segfault.
            reset          <- resetAnyway
            resetTempering <- resetAnyway
            resetAnyway   <<- !reset
            totalIters    <<- 1
        }
        my_initializeModel$run()
        ## if(adaptiveOff)
        ##     for(i in seq_along(samplerFunctions))
        ##         samplerFunctions[[i]]$turnOffAdaptation()
        if(resetTempering) {
            nimPrint("Resetting adaptation rate for tempering") 
            totalIters <<- 1
        }
        accCountSwap[,] <<- 0 * accCountSwap[,] 
        if(reset) {
            nimPrint("Resetting MCMC samplers and initial values set from mvTemps row 1") 
            for (tt in 1:nTemps)
                nimCopy(from = model, to = mvTemps, row = tt, logProb = TRUE)
            for(i in seq_along(samplerFunctions))
                samplerFunctions[[i]]$reset()
            mvSamples_offset  <- 0
            mvSamples2_offset <- 0
            setSize(tempTraj,  niter/thin, nTemps)
            resize(mvSamples,  niter/thin) 
            resize(mvSamples2, niter/thin2)
            if (monitorTmax==TRUE) {
                resize(mvSamplesTmax,  niter/thin)
                resize(mvSamples2Tmax, niter/thin2)
            }
        } else {
            mvSamples_offset  <- getsize(mvSamples)
            mvSamples2_offset <- getsize(mvSamples2)
            setSize(tempTraj,  niter/thin, nTemps)
            resize(mvSamples,  mvSamples_offset  + niter/thin)
            resize(mvSamples2, mvSamples2_offset + niter/thin2)
            if (monitorTmax==TRUE) {
                resize(mvSamplesTmax,  mvSamples_offset  + niter/thin)
                resize(mvSamples2Tmax, mvSamples2_offset + niter/thin2)
            } 
        }
        ##### Monitors & Progress Bar #####
        if(dim(samplerTimes)[1] != length(samplerFunctions))            
            setSize(samplerTimes, length(samplerFunctions)) ## samplerTimes <<- numeric(length(samplerFunctions)) ## A BUG!!!!
        if(niter < progressBarLength+3)
            progressBar <- progressBar & 0  ## avoids compiler warning
        if(progressBar) {
            for(iPB1 in 1:4) {
                cat('|'); for(iPB2 in 1:(progressBarLength/4)) cat('-')
            }
            print('|'); cat('|')
        }
        progressBarIncrement <- niter/(progressBarLength+3)
        progressBarNext      <- progressBarIncrement
        progressBarNextFloor <- floor(progressBarNext)
        ##########################
        ######## SAMPLING ########
        for(iter in 1:niter) {
            ## nimPrint(iter)
            checkInterrupt()
            ## ###############################################################
            ## Random Walk & Adaptation Phase:                              ##
            ## i.e. MCMC at each temperature, with or without time tracking ##
            if(time) { ## time == TRUE
                for(tt in 1:nTemps) {
                    ## Copy state of MCMC at temperature tt to model & mvSaved and continue
                    nimCopy(from=mvTemps, to=model,   row=tt,  logProb = TRUE) 
                    nimCopy(from=model,   to=mvSaved, rowTo=1, logProb = TRUE) 
                    for(ss in 1:nSamplersPerT) {
                        iST <- ss + (tt-1) * nSamplersPerT
                        samplerFunctions[[iST]]$setTemp(temp=Temps[tt])
                        samplerTimes[iST] <<- samplerTimes[iST] +
                            run.time(samplerFunctions[[iST]]$run())
                    }
                    ## Copy state of MCMC at temperature tt back to mvTemps
                    nimCopy(from=model, to=mvTemps, rowTo=tt, logProb = TRUE)
                    logProbTemps[tt] <<- model$getLogProb()
                } 
            } else {   ## time == FALSE 
                for(tt in 1:nTemps) {
                    ## Copy state of MCMC at temperature tt to model & mvSaved and continue
                    nimCopy(from=mvTemps, to=model,   row=tt,  logProb = TRUE) 
                    nimCopy(from=model,   to=mvSaved, rowTo=1, logProb = TRUE) 
                    for(ss in 1:nSamplersPerT) {
                        iST <- ss + (tt-1) * nSamplersPerT
                        ## if (iter>=7300) nimPrint(iST)
                        samplerFunctions[[iST]]$setTemp(temp=Temps[tt])
                        samplerFunctions[[iST]]$run()
                    }
                    ## Copy state of MCMC at temperature tt back to mvTemps
                    nimCopy(from=model, to=mvTemps, rowTo=tt, logProb = TRUE)
                    logProbTemps[tt] <<- model$getLogProb()
                }   
            }       
            ## ################################################
            ## Random Swap Phase: Jumps between temperatures ##
            pSwapCalc() ## Updates pSwapMatrix            
            for (ii in 1:nTemps) { 
                if (iter %% 2 == 0) 
                    iFrom <- ii
                else
                    iFrom <- nTemps + 1 - ii
                iTo <- rcat(n=1, prob=pSwapMatrix[iFrom, 1:nTemps]) ##
                if (iFrom==iTo) {
                    nimPrint("THIS SHOULD NOT HAPPEN!!!") ## Delete this line once tested
                }
                pProp <- pSwapMatrix[iFrom,iTo]
                pRev  <- pSwapMatrix[iTo,iFrom]
                lMHR  <- (invTemps[iTo]-invTemps[iFrom]) * 
                    (logProbTemps[iFrom]-logProbTemps[iTo]) + log(pRev) - log(pProp)
                lu    <- log(runif(1))
                if (lu < lMHR) { ## ACCEPT
                    ## Swap model values
                    nimCopy(from=mvTemps, to=mvSaved, row=iTo,   rowTo=1,     logProb = TRUE)
                    nimCopy(from=mvTemps, to=mvTemps, row=iFrom, rowTo=iTo,   logProb = TRUE)
                    nimCopy(from=mvSaved, to=mvTemps, row=1,     rowTo=iFrom, logProb = TRUE)
                    ## Swap logProbTemps elements
                    temporary[1]        <<- logProbTemps[iTo]
                    logProbTemps[iTo]   <<- logProbTemps[iFrom]
                    logProbTemps[iFrom] <<- temporary[1]
                    ## Swap rows of pSwapMat
                    ## browser()
                    temporary[1:nTemps]         <<- pSwapMatrix[iTo,1:nTemps]
                    pSwapMatrix[iTo,1:nTemps]   <<- pSwapMatrix[iFrom,1:nTemps]                    
                    pSwapMatrix[iFrom,1:nTemps] <<- temporary[1:nTemps]
                    ## Swap cols of pSwapMat
                    temporary[1:nTemps]         <<- pSwapMatrix[1:nTemps,iTo]
                    pSwapMatrix[1:nTemps,iTo]   <<- pSwapMatrix[1:nTemps,iFrom]
                    pSwapMatrix[1:nTemps,iFrom] <<- temporary[1:nTemps]
                    ## Acceptance counter
                    accCountSwap[iFrom, iTo]    <<- accCountSwap[iFrom, iTo] + 1 
                } 
            }
            ## nimPrint(pSwapMatrix)
            ## nimPrint(Temps, "     ", invTemps)
            ###################################
            ## Temperature adaptation scheme ##
            totalIters <<- totalIters + 1
            if (adaptTemps) {
                gammaTSA    <- 1 / ((totalIters/tuneTemper1 + 3) ^ tuneTemper2)  ## gammaTSA <- 1 / ((totalIters/thin + 3) ^ 0.8) 
                TempsCurrent[1:nTemps] <<- Temps
                for (ii in 1:nTemps_1) {
                    accProb          <- min(1, exp( (invTemps[ii+1]-invTemps[ii]) * (logProbTemps[ii]-logProbTemps[ii+1]) ))
                    Temps[ii+1]     <<- Temps[ii] + (TempsCurrent[ii+1] - TempsCurrent[ii]) * exp(gammaTSA*(accProb-0.234))
                    if (Temps[ii+1] > ULT)
                        Temps[ii+1] <<- ULT + ii ## To prevent problems when temperature ladder explodes                    
                } 
                for (ii in 1:nTemps_1) 
                    invTemps[ii+1]  <<- 1 / Temps[ii+1]
                ## nimPrint("iter: ", iter)
                ## nimPrint("logProbTemps: ", logProbTemps)
                ## nimPrint(accCountSwap)
            }            
            ####################################
            ## Dropping unneeded temperatures ##
            ##
            ## Can be done manually for now...
            ##
            #################
            ## MCMC Output ##
            if(iter %% thin  == 0) {
                nimCopy(from = mvTemps, to = mvSamples,
                        row = 1, rowTo = mvSamples_offset + iter/thin,  nodes = monitors)
                tempTraj[iter/thin, 1:nTemps] <<- Temps[1:nTemps]
                if(printTemps == 1.0) {
                    if(totalIters %% thinPrintTemps == 0) {
                        nimPrint(iter, asRow(Temps))  
                        if(!adaptTemps)
                            printTemps <- 0.0
                    }
                }
            }
            if(iter %% thin2 == 0) 
                nimCopy(from = mvTemps, to = mvSamples2, 
                        row = 1, rowTo = mvSamples2_offset + iter/thin2, nodes = monitors2) 
            if (monitorTmax==TRUE) {
                if(iter %% thin  == 0)
                    nimCopy(from = mvTemps, to = mvSamplesTmax, row = nTemps,
                            rowTo = mvSamples_offset  + iter/thin,  nodes = monitors) 
                if(iter %% thin2 == 0)
                    nimCopy(from = mvTemps, to = mvSamples2Tmax, row = nTemps,
                            rowTo = mvSamples2_offset + iter/thin2, nodes = monitors2)
            } 
            ## Progress Bar
            if(progressBar & (iter == progressBarNextFloor)) {
                cat('-')
                progressBarNext <- progressBarNext + progressBarIncrement
                progressBarNextFloor <- floor(progressBarNext)
            }
        }
        if(progressBar) print('|')
        nimCopy(from = mvTemps, to = model, row = 1, logProb = TRUE) ## Ensures model is parameterised with T=1 samples prior to exiting.
    },
    ## #########################################################
    methods = list(
        getTimes = function() {
            returnType(double(1))
            return(samplerTimes)
        }, 
        pSwapCalc = function() {
            ## Fill matrix
            for (ii in 2:nTemps) {
                for (jj in 1:(ii-1)) {
                    pSwapMatrix[ii,jj] <<- exp(-abs(logProbTemps[ii] - logProbTemps[jj]))
                    pSwapMatrix[jj,ii] <<- pSwapMatrix[ii,jj]
                }
            }
            ## Normalise rows            
            for (ii in 1:nTemps) {
                rowSum <- sum(pSwapMatrix[ii,1:nTemps])
                if(rowSum>0) {
                    for (jj in 1:nTemps)
                        pSwapMatrix[ii,jj] <<- pSwapMatrix[ii,jj] / rowSum
                } else {
                    for (jj in 1:nTemps) {
                        pSwapMatrix[ii,jj] <<- (ii!=jj) / nTemps_1
                    }
                }
            }
        },
        mvTemps2model = function(row = double()) {
            ## Useful in R for exploring node values at each temp
            nimCopy(mvTemps, model, row=row, logProb=TRUE)
        }
    ),
    where = getLoadingNamespace()
)


# This is a function that will weed out missing indices from the monitors
processMonitorNames <- function(model, nodes){
	isLogProbName <- grepl('logProb', nodes)
	expandedNodeNames <- model$expandNodeNames(nodes[!isLogProbName])
	origLogProbNames <- nodes[isLogProbName]
	expandedLogProbNames <- character()
	if(length(origLogProbNames) > 0){
		nodeName_fromLogProbName <- gsub('logProb_', '', origLogProbNames)
		expandedLogProbNames <- model$modelDef$nodeName2LogProbName(nodeName_fromLogProbName)
	}
	return( c(expandedNodeNames, expandedLogProbNames) )
}


sampler_RW_tempered <- nimbleFunction(
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        temperPriors  <- control$temperPriors
        if (is.null(temperPriors))
            stop("control$temperPriors unspecified in sampler_RW_tempered")
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target) 
        targetNode     <- model$expandNodeNames(target)
        dependantNodes <- calcNodes[!is.element(calcNodes, targetNode)]
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        optimalAR     <- 0.44
        gamma1        <- 0
        range         <- getDistribution(model$getNodeDistribution(target))$range
        ## checks
        if(length(targetAsScalar) > 1) stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))   stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)      stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
        ## initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is just a hack. Only 1st element is used.
    },  
    run = function() { 
        ## nimPrint(temperature[1])
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) {
            propLogScale <- rnorm(1, mean = 0, sd = scale)
            propValue    <- currentValue * exp(propLogScale)
        } else        
            propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective)
            while(propValue < range[1] | propValue > range[2]) {
                if(propValue < range[1]) propValue <- 2*range[1] - propValue
                if(propValue > range[2]) propValue <- 2*range[2] - propValue }        
        if (temperPriors) {
            model[[target]] <<- propValue 
            logMHR <- calculateDiff(model, calcNodes) / temperature[1] + propLogScale ## Original. Tempers everything.
        } else {
            logPriorWeightOri   <- getLogProb(model, target)
            model[[target]]    <<- propValue
            diffLogProb         <- calculateDiff(model, calcNodes)
            logPriorWeightProp  <- getLogProb(model, target)
            diffLogPriorWeights <- logPriorWeightProp - logPriorWeightOri            
            logMHR <- (diffLogProb - diffLogPriorWeights) / temperature[1] + diffLogPriorWeights + propLogScale ## POSSIBLY BUGGY ???            
        }
        jump <- decide(logMHR)
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },                
    methods = list(   
        setTemp = function(temp=double()) {
            temperature[1] <<- temp
        },            
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            scale         <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            gamma1        <<- 0
        }
    ), where = getLoadingNamespace()
)



sampler_RW_block_tempered <- nimbleFunction(
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        temperPriors   <- control$temperPriors
        if (is.null(temperPriors))
            stop("control$temperPriors unspecified in sampler_RW_block_tempered")
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        targetNodes    <- model$expandNodeNames(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        d <- length(targetAsScalar)
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        ## nested function and function list definitions
        my_setAndCalculateDiff  <- setAndCalculateDiff(model, target)
        my_decideAndJump        <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
        ## checks
        if(class(propCov)      != 'matrix')  stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric') stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov)   == d))        stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))            stop('propCov matrix must be symmetric')
        ## Initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is a hack. Only temperature[1] is used.
    },
    run = function() { 
        propValueVector <- generateProposalVector()
        if (temperPriors) 
            lpMHR <- my_setAndCalculateDiff$run(propValueVector) / temperature[1]
        else {
            logPriorWeightOri   <- getLogProb(model, targetNodes)
            lpMHR               <- my_setAndCalculateDiff$run(propValueVector)
            diffLogPriorWeights <- getLogProb(model, targetNodes) - logPriorWeightOri
            lpMHR               <- (lpMHR - diffLogPriorWeights) / temperature[1] + diffLogPriorWeights
        }        
        jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
        ## nimPrint(temperature[1], " ", jump)
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        setTemp = function(temp=double()) {
            temperature[1] <<- temp
        },
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$gamma1
                    for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                    empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                }
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            scale              <<- scaleOriginal
            propCov            <<- propCovOriginal
            chol_propCov       <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
            timesRan           <<- 0
            timesAccepted      <<- 0
            timesAdapted       <<- 0
            my_calcAdaptationFactor$reset()
        }
    ), where = getLoadingNamespace()
)


sampler_slice_tempered <- nimbleFunction(
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptInterval  <- control$adaptInterval
        width          <- control$sliceWidth
        maxSteps       <- control$sliceMaxSteps
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        ## numeric value generation
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        discrete      <- model$isDiscrete(target)
        ## checks     
        if(length(targetAsScalar) > 1)     stop('cannot use slice sampler on more than one target node')
        ## initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is a hack. Only first element is used.
    },
    run = function() { 
        ## nimPrint(temperature[1])
        u  <- getLogProb(model, calcNodes) / temperature[1] - rexp(1, 1)    # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
        x0 <- model[[target]]    # create random interval (L,R), of width 'width', around current value of target
        L  <- x0 - runif(1, 0, 1) * width
        R  <- L + width
        maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
        maxStepsR <- maxSteps - 1 - maxStepsL
        lp <- setAndCalculateTarget(L) / temperature[1]
        while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
            L <- L - width
            lp <- setAndCalculateTarget(L) / temperature[1]
            maxStepsL <- maxStepsL - 1
        }
        lp <- setAndCalculateTarget(R) / temperature[1]
        while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
            R <- R + width
            lp <- setAndCalculateTarget(R) / temperature[1]
            maxStepsR <- maxStepsR - 1
        }
        x1 <- L + runif(1, 0, 1) * (R - L)
        lp <- setAndCalculateTarget(x1) / temperature[1]
        while(is.nan(lp) | lp < u) {   # must be is.nan()
            if(x1 < x0) {
                L <- x1
            } else {
                R <- x1
            }
            x1 <- L + runif(1, 0, 1) * (R - L)           # sample uniformly from (L,R) until sample is inside of slice (with shrinkage)
            lp <- setAndCalculateTarget(x1) / temperature[1]
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        jumpDist <- abs(x1 - x0)
        if(adaptive)     adaptiveProcedure(jumpDist)
    },
    methods = list(
        setTemp = function(temp = double()) {
            temperature[1] <<- temp
        },
        setAndCalculateTarget = function(value = double()) {
            if(discrete)     value <- floor(value)
            model[[target]] <<- value
            lp <- calculate(model, calcNodes)
            returnType(double())
            return(lp)
        },
        adaptiveProcedure = function(jumpDist = double()) {
            timesRan <<- timesRan + 1
            sumJumps <<- sumJumps + jumpDist   # cumulative (absolute) distance between consecutive values
            if(timesRan %% adaptInterval == 0) {
                adaptFactor <- (3/4) ^ timesAdapted
                meanJump <- sumJumps / timesRan
                width <<- width + (2*meanJump - width) * adaptFactor   # exponentially decaying adaptation of 'width' -> 2 * (avg. jump distance)
                timesAdapted <<- timesAdapted + 1
                timesRan <<- 0
                sumJumps <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            width        <<- widthOriginal
            timesRan     <<- 0
            timesAdapted <<- 0
            sumJumps     <<- 0
        }
    ), where = getLoadingNamespace()
)


sampler_RW_multinomial_tempered <- nimbleFunction( 
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        useTempering  <- control$useTempering
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        targetAllNodes <- unique(model$expandNodeNames(target))
        calcNodes      <- model$getDependencies(target) 
        lTarget        <- length(targetAsScalar)
        Ntotal         <- sum(values(model,target))
        NOverL         <- Ntotal / lTarget
        ## numeric value generation
        Zeros             <- matrix(0, lTarget, lTarget)
        Ones              <- matrix(1, lTarget, lTarget)
        timesRan          <- Zeros
        AcceptRates       <- Zeros
        ScaleShifts       <- Zeros
        totalAdapted      <- Zeros
        timesAccepted     <- Zeros
        ENSwapMatrix      <- Ones
        ENSwapDeltaMatrix <- Ones
        RescaleThreshold  <- 0.2 * Ones
        lpProp  <- 0
        lpRev   <- 0
        Pi      <- pi 
        PiOver2 <- Pi / 2 ## Irrational number prevents recycling becoming degenerate
        u       <- runif(1, 0, Pi)
        ## nested function and function list definitions
        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        my_decideAndJump       <- decideAndJump(model, mvSaved, calcNodes)
        ## checks
        if(model$getNodeDistribution(target) != 'dmulti')   stop('can only use RW_multinomial sampler for multinomial distributions')
        if(length(targetAllNodes) > 1)                      stop('cannot use RW_multinomial sampler on more than one target')
        if(adaptive & adaptInterval < 100)                  stop('adaptInterval < 100 is not recommended for RW_multinomial sampler')
        ## initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is a hack. Only first element is used.
    },
    run = function() {
        for(iFROM in 1:lTarget) {            
            for(iTO in 1:(lTarget-1)) {
                if(u > PiOver2) {                
                    iFrom <- iFROM
                    iTo   <- iTO
                    if (iFrom == iTo)
                        iTo <- lTarget
                    u <<- 2 * (u - PiOver2)   # recycle u
                } else {
                    iFrom <- iTO
                    iTo   <- iFROM
                    if (iFrom == iTo)
                        iFrom <- lTarget
                    u <<- 2 * (PiOver2 - u)   # recycle u
                }
                propValueVector <- generateProposalVector(iFrom, iTo)
                if (useTempering) 
                    lpMHR <- my_setAndCalculateDiff$run(propValueVector) / temperature[1] + lpRev - lpProp
                else
                    lpMHR <- my_setAndCalculateDiff$run(propValueVector) + lpRev - lpProp
                jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## returns lpMHR + 0 - 0 + 0
                if(adaptive)   adaptiveProcedure(jump=jump, iFrom=iFrom, iTo=iTo)
            }
        }
    },
    methods = list(
        setTemp = function(temp=double()) {
            temperature[1] <<- temp
        },
        generateProposalVector = function(iFrom = integer(), iTo = integer()) { 
            propVector <- values(model,target) 
            pSwap      <- min(1, max(1, ENSwapMatrix[iFrom,iTo]) / propVector[iFrom]) 
            nSwap      <- rbinom(n=1,   size=propVector[iFrom], prob=pSwap) 
            lpProp    <<- dbinom(nSwap, size=propVector[iFrom], prob=pSwap, log=TRUE) 
            propVector[iFrom] <- propVector[iFrom] - nSwap 
            propVector[iTo]   <- propVector[iTo]   + nSwap 
            pRevSwap   <- min(1, max(1, ENSwapMatrix[iTo,iFrom]) / (propVector[iTo] + nSwap)) 
            lpRev     <<- dbinom(nSwap, size=propVector[iTo], prob=pRevSwap, log=TRUE) 
            returnType(double(1)) 
            return(propVector) 
        },
        adaptiveProcedure = function(jump=logical(), iFrom=integer(), iTo=integer()) {
            NVector <- values(model,target) 
            timesRan[iFrom, iTo] <<- timesRan[iFrom, iTo] + 1
            if(jump)
                timesAccepted[iFrom, iTo] <<- timesAccepted[iFrom, iTo] + 1
            if (timesRan[iFrom, iTo] %% adaptInterval == 0) {
                totalAdapted[iFrom, iTo] <<- totalAdapted[iFrom, iTo] + 1
                accRate                   <- timesAccepted[iFrom, iTo] / timesRan[iFrom, iTo]
                AcceptRates[iFrom, iTo]  <<- accRate
                if (accRate > 0.5) {
                    ENSwapMatrix[iFrom, iTo] <<-
                        min(Ntotal,
                            ENSwapMatrix[iFrom,iTo] + ENSwapDeltaMatrix[iFrom, iTo] / totalAdapted[iFrom,iTo])
                } else {
                    ENSwapMatrix[iFrom, iTo] <<-
                        max(1,
                            ENSwapMatrix[iFrom,iTo] - ENSwapDeltaMatrix[iFrom,iTo] / totalAdapted[iFrom,iTo])
                } 
                if(accRate<RescaleThreshold[iFrom,iTo] | accRate>(1-RescaleThreshold[iFrom,iTo])) {
                    ## rescale iff ENSwapMatrix[iFrom, iTo] is not set to an upper or lower bound 
                    if (ENSwapMatrix[iFrom, iTo] > 1 & ENSwapMatrix[iFrom, iTo] < Ntotal) {
                        ScaleShifts[iFrom, iTo]       <<- ScaleShifts[iFrom, iTo] + 1 
                        ENSwapDeltaMatrix[iFrom, iTo] <<- min(NOverL, ENSwapDeltaMatrix[iFrom, iTo] * totalAdapted[iFrom,iTo] / 10)
                        ENSwapDeltaMatrix[iTo, iFrom] <<- ENSwapDeltaMatrix[iFrom, iTo] 
                        RescaleThreshold[iFrom,iTo]   <<- 0.2 * 0.95^ScaleShifts[iFrom, iTo]
                    }
                }
                ## lower Bound 
                if(ENSwapMatrix[iFrom, iTo] < 1)
                    ENSwapMatrix[iFrom, iTo] <<- 1                
                ## symmetry in ENSwapMatrix helps maintain good acceptance rates
                ENSwapMatrix[iTo,iFrom]   <<- ENSwapMatrix[iFrom,iTo]
                timesRan[iFrom, iTo]      <<- 0
                timesAccepted[iFrom, iTo] <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            timesRan          <<- Zeros
            AcceptRates       <<- Zeros
            ScaleShifts       <<- Zeros
            totalAdapted      <<- Zeros
            timesAccepted     <<- Zeros
            ENSwapMatrix      <<- Ones
            ENSwapDeltaMatrix <<- Ones
            RescaleThreshold  <<- 0.2 * Ones
        }
    ), where = getLoadingNamespace()
)
