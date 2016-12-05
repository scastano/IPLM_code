#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

####################################################################
## MODEL in NIMBLE. CASE STUDY I (correlated (via copula) stages) ## 
####################################################################
rm(list=ls())

setwd("~/IPLM_paper/nimble")
options(width=200)
library(nimble)
library(latex2exp)
library(coda)
source("FUNCTIONS_R.R")
source("FUNCTIONS_NIMBLE.R")
baseDir <- "~/IPLM_paper/nimble/project_correlatedStages"
setwd(baseDir)

PDF <- FALSE

## MCMC PARAMS 
BLOCK       <- "RhoQ_ResParas"  
SAVE_OUTPUT <- TRUE ## FALSE
NITER       <- 1E3

## ###########################################
## Take arguments from script, if available ##
CA <- commandArgs(TRUE)
if (length(CA)==0) {
    UseScript <- FALSE
} else {
    UseScript <- TRUE
}

if (UseScript) { 
    print(CA)
    qsubID <- as.integer(CA)[2]
    print(qsubID)
} else {
    qsubID <- 123456789
}

## Source to get proposed parameter
source("corr_parameters.R")


####################################################
## Multi-stage model INCLUDING individual quality ##
####################################################
CodeQuality <- nimbleCode ({
    ## Priors for marginal kernel parameters of each stage
    for (j in 1:S) {
        logit(paras[j,1]) ~ dLogitUnif()
        logit(paras[j,2]) ~ dLogitUnif()
        logit(paras[j,3]) ~ dLogitUnif()
        alphas[j,1:2]    <- nf_muVar2alpha( nf_muSc2muVar ( paras[j,1:2] ))
        constraintData[j] ~ dconstraint( min(alphas[j,1:2]) > 1 )
        logit(resCont[j]) ~ dLogitUnif()                 ## Continuous
        res[j]           <- ceiling(ResMax * resCont[j]) ## Discretised
    }
    ## Observed data likelihood
    for (i in 1:N) {
        ## stage 1
        y[i,1:3,1] ~ dJSMD_condit_quality(paras = paras[1,1:3], res = res[1],
                                          CensorTime = CensorTime, 
                                          quality = qual[i], rho = rho)
        for (j in 2:S) {
            ## stage j
            y[i,1:3,j] ~ dJSMD_condit_quality(
                paras = paras[j,1:3], res = res[j], quality = qual[i], rho = rho,
                ## CensoredTime = 1(completed last stage) * (T_censor - T_0)
                CensorTime = (1==y[i,2,j-1])*(CensorTime-sum(y[i,1,1:(j-1)])))            
        }
        ## Prior on individual quality
        logit(qual[i]) ~ dLogitUnif()
    }
    logit(rho) ~ dLogitUnif() 
})

## ASSEMBLE NIMBLE MODEL
## Warning messages can be ignored here
Qmodel <- nimbleModel(CodeQuality, constants=Constants, inits = Inits, data = Data, check= FALSE)

InitsOri <- Inits
rm(Inits, Data) ## To be replaced below following simulation

##################################
## Simulate data from the model ##
system.time(Qmodel$simulate('y', includeData = TRUE))
Qmodel$y
Qmodel$calculate(Qmodel$getDependencies("y")) ## Update the associated log likelihood


#####################################################
## Tabulate outcomes for individuals in each stage ##
## Columns : Developed, Died, interval (1 == observed, 2== dead,  == right censor)
(fy1  <- factor(Qmodel$y[,2,1], levels=c(1,2,9)))
(fy2  <- factor(Qmodel$y[,2,2], levels=c(1,2,9)))
(t1   <- table(fy1)) ## Tabulate Stage 1
(t2   <- table(fy2)) ## Tabulate Stage 2, censored not yet true
t2[3] <- t2[3] - t1[3] - t1[2] ## Remove individuals who did not start stage 2
tab   <- rbind(t1,t2)
colnames(tab) <- c("Dev","Die","Censor") ; rownames(tab) <- c("Stage1", "Stage2")
print(tab)
(fulldev <- tab[2,1])


## ###########################
## REJECT/ACCEPT parameters ##
## "Al" starts being lower than 1
CQuality <- compileNimble(Qmodel)

## Define "quasi-random" lower bound  
(LB <- 35 + (qsubID%%( ((N-35) - 2) )))

nloop        <- 0
resample_alp <- TRUE
resample_sur <- FALSE

while (resample_alp==TRUE | resample_sur==TRUE) {
    nloop <- nloop + 1
    print(paste("while-loop number is", nloop))
    if(resample_alp==TRUE) {
        Nodes <- CQuality$getDependencies(c('logit_paras', 'logit_resCont', 'logit_qual', 'rho', 'y'))
        CQuality$simulate(Nodes, includeData = TRUE, constraintData=c(1,1)) 
        print(CQuality$alphas)
    } else if (resample_sur==TRUE) {
        low    <- which(CQuality$paras[,3]==min(CQuality$paras[,3]))
        lowLim <- min(CQuality$paras[,3])
        if(low==1) {
            CQuality$logit_paras[1,3] <- logit(runif(1,lowLim,1))
            Nodes <- CQuality$getDependencies('logit_paras[1,3]', self=FALSE)
        } else { ## low==2
            CQuality$logit_paras[2,3] <- logit(runif(1,lowLim,1))
            Nodes <- CQuality$getDependencies('logit_paras[2,3]', self=FALSE)
        }
        CQuality$simulate(Nodes, includeData=TRUE)
    }
    ## Optional visual checks
    if (FALSE) { ## TRUE
        CQuality$simulate(Nodes, includeData = TRUE)
        par(mfrow=c(2,1)); Al <- CQuality$alphas
        curve(dbeta(x, Al[1,1], Al[1,2]),0,1,n=1001)
        curve(dbeta(x, Al[2,1], Al[2,2]),0,1,n=1001)
    }
    ## Count dev,die,cens
    (die1 <- sum(CQuality$y[,2,1]==2))
    (cen1 <- sum(CQuality$y[,2,1]==9))
    (die2 <- sum(CQuality$y[,2,2]==2)) 
    ## Remove individuals who did not start stage 2
    (cen2 <- sum(CQuality$y[,2,2]==9) - sum(CQuality$y[,2,1]==9) - sum(CQuality$y[,2,1]==2))
    (fulldev <- sum(CQuality$y[,2,2][CQuality$y[,2,2]==1]))
    print(paste("nb of full developed individuals is ", fulldev))
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    ##
    resample_alp <- ( (fulldev < LB) & ( (cen1+cen2)  > (die1+die2) ) ) | any(CQuality$alphas < 1)
    resample_sur <- ( (fulldev < LB) & ( (cen1+cen2) <= (die1+die2) ) )
}

## Check
if (FALSE) { ## TRUE
    fulldev
    par(mfrow=c(2,1)); Al <- CQuality$alphas
    curve(dbeta(x, Al[1,1], Al[1,2]),0,1,n=1001)
    curve(dbeta(x, Al[2,1], Al[2,2]),0,1,n=1001)
}


###########################
## PREPARE INFO TO TRACK ##

## Update
qual <- CQuality$qual

## To distinguish individual history to use in plots of posteriors of Quality
(fy1  <- factor(CQuality$y[,2,1], levels=c(1,2,9)))
(fy2  <- factor(CQuality$y[,2,2], levels=c(1,2,9)))
## Right Censor indicator
censStg2 <- as.numeric(CQuality$y[,,2][,2]==9 -(CQuality$y[,,1][,2] == 9 | CQuality$y[,,1][,2] == 2))
## Quality-values for individuals of different destinies
(qual_DieStg1  <- qual[which(fy1==2)])      ## Died in stage 1
(qual_CensStg1 <- qual[which(fy1==9)])      ## Censored in stage 1 (rare)
(qual_DieStg2  <- qual[which(fy2==2)])      ## Died in stage 2
(qual_CensStg2 <- qual[which(censStg2==1)]) ## Censored in stage 2
(qual_FullDev  <- qual[which(fy2==1)])      ## Fully developed individuals
## Check total vector length equals N 
(length(qual_DieStg1) + length(qual_CensStg1) + length(qual_CensStg2) + length(qual_DieStg2) + length(qual_FullDev))
## Quality list for individuals of different destinies
qual_history <- list(qual_DieStg1=qual_DieStg1, qual_DieStg2=qual_DieStg2,
                     qual_CensStg1=qual_CensStg1, qual_CensStg2=qual_CensStg2,
                     qual_FullDev=qual_FullDev)         
## Associated indices
(qualOrder_DieStg1  <- which(fy1==2))
(qualOrder_CensStg1 <- which(fy1==9)) 
(qualOrder_DieStg2  <- which(fy2==2))
(qualOrder_CensStg2 <- which(censStg2==1)) 
(qualOrder_FullDev  <- which(fy2==1)) 
## Check total vector length equals N 
(length(qualOrder_DieStg1) + length(qualOrder_CensStg1) + length(qualOrder_DieStg2) + length(qualOrder_CensStg2) +length(qualOrder_FullDev))
## List of indices for individuals of different destinies
qualOrder_history <- list(qualOrder_DieStg1=qualOrder_DieStg1, qualOrder_DieStg2=qualOrder_DieStg2,
                          qualOrder_CensStg1=qualOrder_CensStg1, qualOrder_CensStg2=qualOrder_CensStg2,
                          qualOrder_FullDev=qualOrder_FullDev)         


## Redefine Data and Inits in order to re-build the Nimble model
Data  <- list(y=CQuality$y, constraintData=c(1,1)) 
Inits <- list(paras         = CQuality$paras,
              rho           = CQuality$rho,
              qual          = CQuality$qual,
              logit_resCont = CQuality$logit_resCont,
              resCont       = CQuality$resCont,
              res           = ceiling(ResMax * CQuality$resCont), 
              logit_paras   = CQuality$logit_paras,
              logit_rho     = CQuality$logit_rho,
              logit_qual    = CQuality$logit_qual)


## Rebuild with new data (avoids memory based crash)
Qmodel <- nimbleModel(CodeQuality, constants=Constants, inits = Inits, data = Data, check= FALSE) ## TRUE
CQuality <- compileNimble(Qmodel)


## #########################
## INITIALISE output files #
setwd(baseDir)
source("corr_init_output.R")
nimPrint("Working directory is", getwd())


####################
## CONFIGURE MCMC ##
Mon1 <- c('paras','qual', 'res', 'resCont', 'rho')
Mon2 <- c('logProb_logit_paras', 'logProb_logit_resCont',
          'logProb_logit_qual', 'logProb_logit_rho')

ConfQuality <- configureMCMC(CQuality, monitors=Mon1, monitors2=Mon2)
MCMCQuality <- buildMCMC(ConfQuality)

if (BLOCK=="RhoQ_ResParas") {
    ConfQuality$addSampler(target = c('rho', 'qual'),
                           type = 'RW_block', control = list(adaptInterval = 100))
    ConfQuality$addSampler(target = c('res', 'paras'),
                           type = 'RW_block', control = list(adaptInterval = 100))
}
ConfQuality$printSamplers()


## ###################
## BUILD & COMPILE  ##
MCMCQuality <- buildMCMC(ConfQuality)
ConfQuality$getMonitors() 
cMCMCQuality <- compileNimble(MCMCQuality)


## ##########################
## MCMC - ADAPTIVE BURN-IN ##

print("####################################")
print("starting adaptive burn-in WHILE loop")
print("####################################")
## Redefine thin & thin2 if necessary
THIN               <- 1
cMCMCQuality$thin  <- THIN
cMCMCQuality$thin2 <- THIN 
(meanL             <- CQuality$calculate())
meanL_previous     <- -Inf
ii                 <- 0

while(meanL > meanL_previous + 2) {
    ii <- ii+1
    meanL_previous <- meanL
    print(paste0("iteration nb.", ii, "within while loop. meanL = ", meanL_previous))
    ## ##############
    ## Short run 1 ##
    nIter <- 1E3
    syst1 <- system.time(cMCMCQuality$run(nIter))
    nimPrint("While loop: 1st short run finished. sysT = ", syst1[3])
    samples <- tail(as.matrix(cMCMCQuality$mvSamples), nIter/THIN) 
    dim(samples); head(samples); tail(samples)
    samples2 <- tail(as.matrix(cMCMCQuality$mvSamples2), floor(nIter/THIN)) 
    meanL    <- mean(rowSums(samples2))
    nimPrint("meanL = ", meanL)
    ## ###############
    ## Short run 2 ##
    nIter <- 1E3
    syst2 <- system.time((cMCMCQuality$run(nIter, reset=FALSE))) ## Do not reset the adaptive MCMC
    nimPrint("While loop: 2nd short run finished. sysT = ", syst2[3])
    ## update meanL
    samples2 <- tail(as.matrix(cMCMCQuality$mvSamples2), floor(nIter/THIN)) 
    meanL    <- mean(rowSums(samples2))
    nimPrint("meanL = ", meanL)
    ## Calculate & print ESS 
    samples  <- tail(as.matrix(cMCMCQuality$mvSamples),  floor(nIter/10))
    alphCols <- substring(colnames(samples),1,6)=="alphas"
    sub      <- samples[,!alphCols]
    mc        <- as.mcmc(sub)
    ESS       <- effectiveSize(mc) 
    (ESS      <- ESS[order(ESS)])
    nimPrint(ESS)
}
 

## ######################################################################
## Extract MCMC samples & LogProbs from last round of adaptive burn-in ##
samples     <- tail(as.matrix(cMCMCQuality$mvSamples),     floor(nIter/THIN)) 
samples2    <- tail(as.matrix(cMCMCQuality$mvSamples2),    floor(nIter/THIN))
print("Finished extracting model values samples & samples2")

## Check dimensions match
if( dim(samples)[1] != dim(samples2)[1])
    stop("nrow of samples & samples2 don't match")
print("DIMENSION OF samples & samples2 are OK")


## Optionally, examine trajectories of parameters
if (PDF) {
    mc <- as.mcmc(tail(samples, nIter/THIN))
    pdf(file="Paras-trajectories.pdf")
    plot(mc)
    dev.off()
}

## ###########################################
## ### DIAGNOSTICS - Effective Sample Size ###
print("starting ESS calculation")

## Filter columns
alphCols <- substring(colnames(samples),1,6)=="alphas"
sub      <- samples[,!alphCols]
## Check: these should match
matrix(sub[nrow(sub),1:6], nrow=1)
as.numeric(ilogit(CQuality$logit_paras))

mc         <- as.mcmc(sub)
ESS        <- effectiveSize(mc)
(ESS       <- ESS[order(ESS)])
iCandidate <- 1  
candidate  <- names(ESS[iCandidate])  ## CANDIDATE FOR NEW BLOCK SAMPLERS (to add to the list of samplers above)
(SuggestedThinning <- round(nIter / ESS[iCandidate]))  
nimPrint("SuggestedThinning is ", SuggestedThinning)

## ############################################################################
## If min(ESS) & SuggestedThinning are terrible at this point...             ##
## perhaps return to sampler declarations and add more samplers              ##
## Alternatively pay attention to what Temps and samplesTMax have been doing ##
## ############################################################################

## #################################
## Some checks and optional plots ##
if(PDF) {
    ## These should match
    mv2 <- tail(as.matrix(cMCMCQuality$mvSamples2), nIter/THIN/5)
    logliks <- rowSums(mv2)
    tail(logliks, 1)
    CQuality$calculate()
    pdf(file="LogProb.pdf")
    plot(logliks, typ="l")
    dev.off()
}

##########################
##########################
#### MCMC - INFERENCE ####
##########################
##########################
print("Starting MCMC for inference")
THIN  <- ceiling(SuggestedThinning * 2)  ## Estimated from last short run
nIter <-  NITER * THIN
cMCMCQuality$thin  <- THIN
cMCMCQuality$thin2 <- THIN

sysT <- system.time(cMCMCQuality$run(nIter, reset=FALSE))

nimPrint("MCMC finished. sysT = ", sysT[3])

## Get again samples & some other info to store in output
samples  <- tail(as.matrix(cMCMCQuality$mvSamples),  floor(nIter/THIN))
samples2 <- tail(as.matrix(cMCMCQuality$mvSamples2), floor(nIter/THIN))
alphCols <- substring(colnames(samples),1,6)=="alphas"
sub      <- samples[,!alphCols]
mc       <- as.mcmc(sub, floor(nrow(samples)/2))
ESS      <- effectiveSize(mc) 
nimPrint(ESS)

## GET friendly colnames for samples in csv file
colnames(samples) <- c('mu1', 'sc1', 'mu2', 'sc2', 'sur1', 'sur2',
                       paste0("qual" , 1:N), 'res1', 'res2',
                       'resCont1', 'resCont2', 'rho')

setwd(mcmcDir)
outputfile <- paste0("fullOutput", ".Rdata")
nimPrint("Current working directory is: ", getwd())
nimPrint("MCMC output going to file: ", outputfile)


## ##############
## SAVE OUTPUT ##
if (SAVE_OUTPUT) {
    save(Constants, Inits, Data, THIN, SuggestedThinning,  NITER, fulldev, qualOrder_history,
         qual_history, samples, sysT, BLOCK, InitsOri,
         file=outputfile)
}



############################
## Generate FINISHED file ##
if (TRUE) { ## FALSE   
    finishedFile <- file("FINISHED.txt","w")    
    write("MCMC FINISHED", file=finishedFile, ncolumns=1, sep="\t")
    close(finishedFile)
}

print("%%%%%%% FINISHED %%%%%%%%")

