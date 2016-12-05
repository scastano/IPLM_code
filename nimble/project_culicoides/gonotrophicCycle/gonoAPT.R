#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

#####################################################
## IPLM model and MCMC for stage=Gonotrophic cycle ## 
## Sampler: adaptive parallel tempering algorithm  ##
#####################################################
rm(list=ls())

## Set basic parameters
options(width=500)
## libraries
library(latex2exp)
library(nimble)
library(coda)
PDF <-  FALSE 

baseDir <- "~/IPLM_paper/nimble/project_culicoides/gonotrophicCycle"

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
    iRES   <- as.integer(CA)[1]
    qsubID <- as.integer(CA)[2]
    print(iRES)
    print(qsubID)
} else {
    qsubID <- 123456789
}

## ###########################################
## Source scripts to initialise the process ##
setwd(baseDir)
source("../../FUNCTIONS_R.R")
source("../../FUNCTIONS_NIMBLE.R")
source("../FUNCTIONS_CULICOIDES.R")
source("../../RW_APT/APT_build.R")
source("gonoDataConstantsInitial.R") 

#################################################################################################
## BUGS code for integral projection Lefkovitch matrix (IPLM) model for gonotrophic cycle data ##
CodeGono <- nimbleCode ({
    ## KERNEL PARAMETERS AND TRAVELLING WAVE
    for (tt in 1:nTempsGC) {
        logit(paras[tt,1])  ~ dLogitUnif()  ## Expected value
        logit(paras[tt,2])  ~ dLogitUnif()  ## Scale
        logit(paras[tt,3])  ~ dLogitUnif()  ## Survival
        EFecundity[tt]      ~ dgamma(shape = 0.5 + totalEggs[tt], rate = totalOvipos[tt]) ## Poisson likelihood * Jeffreys prior => Gamma posterior
        ## REPARAMETERISATION 
        alphas[tt,1:2] <- nf_muVar2alpha(nf_muSc2muVar(paras[tt,1:2]))
        P01[tt]        <- qbeta(0.01, alphas[tt,1], alphas[tt,2])
        P99[tt]        <- qbeta(0.99, alphas[tt,1], alphas[tt,2])    
    }
    ## TRAVELLING WAVES | TEMPERATURE 
    pLayDie15[1:nSteps[1], 1:2] <- travellingWave_paras2pDevDead(paras[1, 1:3], res = resGC, nSteps[1]) 
    pLayDie20[1:nSteps[2], 1:2] <- travellingWave_paras2pDevDead(paras[2, 1:3], res = resGC, nSteps[2]) 
    pLayDie25[1:nSteps[3], 1:2] <- travellingWave_paras2pDevDead(paras[3, 1:3], res = resGC, nSteps[3]) 
    ## LIKELIHOOD | TEMPERATURE
    pVec15[1:nSteps21[1]] <- twoColumnMatrix2Vector(pLayDie15[1:nSteps[1], 1:2])
    pVec20[1:nSteps21[2]] <- twoColumnMatrix2Vector(pLayDie20[1:nSteps[2], 1:2])
    pVec25[1:nSteps21[3]] <- twoColumnMatrix2Vector(pLayDie25[1:nSteps[3], 1:2])
    for (ii in 1:nObs[1]) { 
        gonoData15[ii] ~ dcat(pVec15[1:nSteps21[1]]) 
    }
    for (ii in 1:nObs[2]) { 
        gonoData20[ii] ~ dcat(pVec20[1:nSteps21[2]]) 
    }
    for (ii in 1:nObs[3]) { 
        gonoData25[ii] ~ dcat(pVec25[1:nSteps21[3]]) 
    } 
    ## SHAPE CONSTRAINTS
    constraintData ~ dconstraint(1 ==
                                 (1 < min(alphas[1:nTempsGC, 1:2])) *
                                 dUnimodal(EFecundity[1:nTempsGC])  *
                                 dUnimodal(paras[1:nTempsGC, 1])    * 
                                 dUnimodal(paras[1:nTempsGC, 3])    *
                                 dUnimodal(P01[1:nTempsGC])         *
                                 dUnimodal(P99[1:nTempsGC]))
})

## ASSEMBLE NIMBLE MODEL
gonoModel <- nimbleModel(CodeGono, constants=Constants, inits=Inits, data=Data, debug=FALSE, check=FALSE) 

###################
## EXAMINE MODEL ##
gonoModel$getNodeNames()
gonoModel$paras
gonoModel$resGC
gonoModel$logit_paras
gonoModel$EFecundity
gonoModel$getDependencies(c('logit_paras'))

#################################################
## Ensure deterministic nodes are up to date   ##
## Update log-likelihood of all dependent data ##
gonoModel$getLogProb(gonoModel$getDependencies("paras"))   ## Could give NA if some nodes were not initialised
gonoModel$simulate(gonoModel$getDependencies("paras"))     ## Update deterministic nodes
gonoModel$calculate(gonoModel$getDependencies("paras"))    ## Update logProb of all dependent nodes
gonoModel$pLayDie15
gonoModel$pLayDie20
gonoModel$pLayDie25


########################################################
## SET resGC from shell script input here & calculate ## 
if (UseScript)
    gonoModel$resGC <- iRES 

print(gonoModel$calculate())
resGC <- gonoModel$resGC
nimPrint("resGC=", gonoModel$resGC[1])

## Now resGC is fixed at the desired value, initialise the output files
setwd(baseDir)
getwd()
source('initialise_output.R')
nimPrint("Working directory is", getwd())

#################################################
## COMPILE MODEL (required for compiling MCMC) ##
cGonoModel <- compileNimble(gonoModel) 


#######################################################
## Check all logProb monitors are included in output ##
## The following lines should give the same output   ## 
sum(gonoModel$logProb_logit_paras,
    gonoModel$logProb_EFecundity,
    gonoModel$logProb_gonoData15,
    gonoModel$logProb_gonoData20,
    gonoModel$logProb_gonoData25)
gonoModel$calculate() 


#################################################
## CONFIGURE ADAPTIVE PARALLEL TEMPERING (APT) ##
mcmcConfGono <- configureMCMC(gonoModel, nodes=NULL, control=list(temperPriors=FALSE))
## Monitors
mcmcConfGono$getMonitors()
mcmcConfGono$resetMonitors()
mcmcConfGono$addMonitors(c('resGC', 'paras', 'EFecundity', 'P01', 'P99'))
mcmcConfGono$addMonitors2(c('logProb_logit_paras', 'logProb_EFecundity',
                            'logProb_gonoData15', 'logProb_gonoData20', 'logProb_gonoData25')) 

## Add samplers
mcmcConfGono$addSampler(target=c("logit_paras[1:3, 1:3]"), type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("EFecundity[1:3]"),       type='sampler_RW_block_tempered')
mcmcConfGono$addSampler(target=c("logit_paras[1:3, 1:2]"), type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[1:3, 3]"),   type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[1:2, 1:2]"), type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[2:3, 1:2]"), type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("EFecundity[1]"),         type='sampler_RW_tempered') 
mcmcConfGono$addSampler(target=c("EFecundity[2]"),         type='sampler_RW_tempered') 
mcmcConfGono$addSampler(target=c("EFecundity[3]"),         type='sampler_RW_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[1, 3]"),     type='sampler_RW_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[2, 3]"),     type='sampler_RW_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[3, 3]"),     type='sampler_RW_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[1, 1:2]"),   type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[2, 1:2]"),   type='sampler_RW_block_tempered') 
mcmcConfGono$addSampler(target=c("logit_paras[3, 1:2]"),   type='sampler_RW_block_tempered') 

mcmcConfGono$printSamplers()

#########################
## BUILD & COMPILE APT ##
nTemps    <- 11 
mcmcGono  <- buildAPT(mcmcConfGono, Temps=exp(seq(0,log(100),l=nTemps)), monitorTmax=TRUE, ULT=1E6) 
cMcmcGono <- compileNimble(mcmcGono, project=gonoModel) 

#############################
## MCMC - ADAPTIVE BURN-IN ##
print("#########################")
print("Starting adaptive burn-in")
print("#########################")

THIN <- 1
cMcmcGono$thin  <- THIN
cMcmcGono$thin2 <- THIN
TuneTemper      <- c(1, 1)
meanL           <- cGonoModel$calculate() 
meanL_previous  <- -Inf
ii              <- 0

while(meanL > meanL_previous + 2) {
    ii <- ii+1
    meanL_previous <- meanL
    print(paste0("iteration nb.", ii, "within while loop. meanL = ", meanL_previous))
    #################
    ## Short run 1 ##
    nIter <- 1E4; cMcmcGono$thinPrintTemps <- nIter / 10    
    syst1 <- system.time(cMcmcGono$run(nIter,
                                       reset          = TRUE,  ## Resets the adaptive MCMC. Let's proposal distributions change direction if required.
                                       adaptTemps     = FALSE, ## Prevents temperature ladder adaptation (to avoid volatile behaviour when counter is reset)
                                       resetTempering = TRUE,  ## Resets counter used in temperature ladder adaptation
                                       printTemps     = TRUE,  ## Will print once only
                                       tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
    ## Update meanL
    nimPrint("While loop: 1st short run finished. sysT = ", syst1[3])
    samples2       <- tail(as.matrix(cMcmcGono$mvSamples2), floor(nIter/THIN)) ## LogProbs 
    meanL          <- mean(rowSums(samples2))
    nimPrint("meanL = ", meanL)
    #################
    ## Short run 2 ##
    nIter <- 1E4; cMcmcGono$thinPrintTemps <- round(nIter / 10) ## Ensures temps are only printed 10 times
    syst2 <- system.time(cMcmcGono$run(nIter,
                                       reset          = FALSE, ## Do not reset the adaptive MCMC, let adaptation continue as it is
                                       adaptTemps     = TRUE,  ## Allows temperature ladder to adjust 
                                       resetTempering = FALSE, ## Keeps the adjustments modest so avoids volatile behaviour
                                       printTemps     = TRUE,  ## Prevents verbose printing of temperature ladder updates
                                       tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
    nimPrint("While loop: 2nd short run finished. sysT = ", syst2[3])
    ## Update meanL
    samples2 <- tail(as.matrix(cMcmcGono$mvSamples2), floor(nIter/THIN)) ## LogProbs 
    meanL    <- mean(rowSums(samples2))
    ## Calculate & print ESS 
    samples   <- tail(as.matrix(cMcmcGono$mvSamples),  floor(nIter/THIN))
    parasCols <- substring(colnames(samples),1,5)=="paras"
    sub       <- samples[,parasCols]
    mc        <- as.mcmc(sub)
    ESS       <- effectiveSize(mc) 
    (ESS      <- ESS[order(ESS)])
    nimPrint(ESS)     
}
print(paste0("iteration nb.", ii, "within while loop. meanL = ", meanL))

print("accCountSwap")
print(cMcmcGono$accCountSwap)


## ######################################################################
## Extract MCMC samples & LogProbs from last round of adaptive burn-in ##
samples     <- tail(as.matrix(cMcmcGono$mvSamples),     floor(nIter/THIN)) 
samples2    <- tail(as.matrix(cMcmcGono$mvSamples2),    floor(nIter/THIN)) ## LogProbs
print("Finished extracting model values samples & samples2")

## Check dimensions match
if( dim(samples)[1] != dim(samples2)[1])
    stop("nrow of samples & samples2 don't match")
print("DIMENSION OF samples & samples2 are OK")

## Optionally, Plot Temperature Trajectories
if (PDF) {
    print("start plotting first PDF")
    pdf(file=paste0("gono_APT_temperature_trajectory_", nTemps, ".pdf"))
    par(mfrow=n2mfrow(1))
    tt     <- cMcmcGono$tempTraj 
    myCols <- rainbow(nTemps)
    sub    <- 1:nrow(tt)
    plot(sub, tt[sub,1], ylim=range(tail(tt, ceiling(4/5*nrow(tt)))),
         typ="l", col=myCols[1], xlab="nb. iterations", ylab="Temperatures ", main="APT - Temperature trajectories")
    for (i in 2:ncol(tt)) lines(sub, tt[sub,i], col=myCols[i])
    dev.off()
}
print("Finished  plotting first  PDF")

## Optionally Examine trajectories of parameters
if (PDF) {
    mc <- as.mcmc(tail(samples, nIter/THIN))
    pdf(file=paste0("gono_APT_paras-trajectories_nTemps", nTemps, ".pdf"))
    plot(mc)
    dev.off()
}

print("Finished last PDFs")

## ###########################################
## ### DIAGNOSTICS - Effective Sample Size ###
print("starting ESS calculation")

## Filter columns
parasCols <- substring(colnames(samples),1,5)=="paras" ## filter columns
sub       <- samples[,parasCols]
## Check: these should match
matrix(sub[nrow(sub),], nrow=3)
cGonoModel$paras
## 
mc         <- as.mcmc(sub)
ESS        <- effectiveSize(mc)  
(ESS       <- ESS[order(ESS)])
iCandidate <- 1
candidate  <- names(ESS[iCandidate]) ## CANDIDATE FOR NEW BLOCK SAMPLERS (to add to the list of samplers above)
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
    mv2 <- tail(as.matrix(cMcmcGono$mvSamples2), nIter/THIN)
    LP <- rowSums(mv2)
    tail(LP, 1)
    cGonoModel$calculate()
    pdf(file=paste0("gono_APT_LogProb_nTemps", nTemps, ".pdf"))
    plot(LP, typ="l")
    dev.off()
}




##########################
##########################
#### MCMC - INFERENCE ####
##########################
##########################

print("Starting MCMC for inference") 
THIN  <- SuggestedThinning * 3 
nIter <-  1E4 * THIN
cMcmcGono$thin <- cMcmcGono$thin2 <- THIN
sysT <- system.time(cMcmcGono$run(nIter,
                                  reset          = FALSE, ## Keeps proposal distributions as they are (initially)
                                  adaptTemps     = FALSE, ## Fix the temperature ladder
                                  resetTempering = FALSE, ## Do not reset the counter for temperature ladder adaptation
                                  printTemps     = TRUE,  ## Do print the temperature ladder once
                                  tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
nimPrint("MCMC finished. sysT = ", sysT[3])


## Get samples & some info to store in output
samples     <- tail(as.matrix(cMcmcGono$mvSamples),     floor(nIter/THIN)) ## Sampledparameters
samples2    <- tail(as.matrix(cMcmcGono$mvSamples2),    floor(nIter/THIN)) ## LogProbs
parasCols   <- substring(colnames(samples),1,5)=="paras"                   ## filter columns
sub         <- samples[,parasCols]
mc          <- as.mcmc(sub, floor(nrow(samples)/2))
ESS         <- effectiveSize(mc) 
nimPrint(ESS)

## Dimension check
if( dim(samples)[1] != dim(samples2)[1])
    stop("nrow of samples & samples2 don't match")
print("DIMENSION OF samples & samples2 are OK")

## Get friendly colnames for samples in csv file
colnames(samples) <- c(paste0("EFecundity",temps),
                       paste("P01",  temps, sep='_'),
                       paste("P99",  temps, sep='_'),
                       paste0("mu",  temps),
                       paste0("sc",  temps),
                       paste0("surv",temps),
                       "resGC")


## ###########################################
## Set output directory & Save MCMC samples ##
setwd(mcmcDir)
nimPrint("Current working directory is: ", getwd())
nimPrint("MCMC output going to file: "   , mcmcFile)
nimPrint("LogProb output going to file: ", mcmcFile2)

print("START WRITING SAMPLES")
write.table(samples,   file=mcmcFile,   col.names=TRUE, row.names = FALSE, sep=",")
write.table(samples2,  file=mcmcFile2,  col.names=TRUE, row.names = FALSE, sep=",")

print("Finished writing samples to file")


##########################
## Save MCMC parameters ##
write(paste("iterations  =", nIter), file = parametersFile, append = FALSE)
write(paste("thin        =", cMcmcGono$thin), file = parametersFile, append = TRUE)
write(paste("thin2       =", cMcmcGono$thin2), file = parametersFile, append = TRUE)
write(paste("nb. samples =", nrow(samples)),  file = parametersFile, append = TRUE)
write(paste("minimum ESS =", floor(min(ESS))),file = parametersFile, append = TRUE)
##
## Add sampler info
##
sink("mcmc_parameters.txt", append=TRUE)
cat("\n")
STOP <- FALSE
ii   <- 0
while (!STOP) {
    ii  <- ii + 1
    TXT <- try(mcmcConfGono$printSamplers(ii))
    if (class(TXT) == "try-error") 
        STOP <- TRUE
    else 
        cat(TXT) 
}
cat("\n")
ESS
sink()


print("%%%%%%% FINISHED %%%%%%%%")

