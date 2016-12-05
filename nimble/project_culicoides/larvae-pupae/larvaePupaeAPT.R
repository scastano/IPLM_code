#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

####################################################
## IPLM model and MCMC for stage= larvae-pupae    ## 
## Sampler: adaptive parallel tempering algorithm ##
####################################################
rm(list=ls())

## Set basic parameters
options(width=500)
## libraries
library(latex2exp)
library(nimble)
library(coda)
PDF <- FALSE

baseDir <- "~/IPLM_paper/nimble/project_culicoides/larvae-pupae"
setwd(baseDir)
nimPrint("Directory is: ", baseDir)

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
    iResL  <- as.integer(CA)[1]
    iResP  <- as.integer(CA)[2]
    qsubID <- as.integer(CA)[3]
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
source("larvaePupaeDataConstantsInitial.R")

############################################################################################
## BUGS code for integral projection Lefkovitch matrix (IPLM) model for larvae-pupae data ##
LPCode <- nimbleCode ({
    ## KERNEL PARAMETERS AND TRAVELLING WAVE
    for (tt in 1:nTempsL) {
        ## LARVAE
        logit(parasL[tt,1]) ~ dLogitUnif()  ## Expected value
        logit(parasL[tt,2]) ~ dLogitUnif()  ## Scale
        logit(parasL[tt,3]) ~ dLogitUnif()  ## Daily Survival
        alphasL[tt,1:2]    <- nf_muVar2alpha(nf_muSc2muVar(parasL[tt,1:2]))
        P01_L[tt]          <- qbeta(0.01, alphasL[tt,1], alphasL[tt,2])
        P99_L[tt]          <- qbeta(0.99, alphasL[tt,1], alphasL[tt,2])
    }
    for (tt in 1:nTempsP) {
        ## PUPAE
        logit(parasP[tt,1]) ~ dLogitUnif()  ## Expected value
        logit(parasP[tt,2]) ~ dLogitUnif()  ## Scale
        logit(parasP[tt,3]) ~ dLogitUnif()  ## Daily Survival
        muSigSurvP[tt,1:3] <- nf_TW1_MuSigSurv(paras = parasP[tt,1:3], res = resP, thresh = thresh)
        alphasP[tt,1:2]    <- nf_muVar2alpha(nf_muSc2muVar(parasP[tt,1:2]))
        P01_P[tt]          <- qbeta(0.01, alphasP[tt,1], alphasP[tt,2])
        P99_P[tt]          <- qbeta(0.99, alphasP[tt,1], alphasP[tt,2])    
    }
    ## LIKELIHOOD OF MULLENS PUPAE DATA
    for (ii in 1:nTempsMullens) {
        MullensPupaeMuSig[ii,1:2] ~ dObMuSig(mu  = muSigSurvP[MullensRowToTempsIndex[ii],1],
                                             sig = muSigSurvP[MullensRowToTempsIndex[ii],2],
                                             N   = MullensSSPost[ii])
    } 
    ## IMPUTATION AND LIKELIHOOD OF VAUGHAN PUPAE DATA
    VaughanSSPre[1:nTempsVaughan] ~ dmultinom(size = VaughanTSSPre, prob = pVaughanSSPre[1:nTempsVaughan])
    for (ii in 1:nTempsVaughan) { 
        VaughanSSPost[ii]    <- myRound(VaughanSSPre[ii] * VaughanSurvivalObs[ii], 0)
        VaughanPseudoData[ii] ~ dPseudoBeta(VaughanSSPre[ii], VaughanSSPost[ii], muSigSurvP[VaughanRowToTempsIndex[ii],3])
        ## DEVELOPMENT
        VaughanPupaeMuSig[ii,1:2] ~ dObMuSig(mu  = muSigSurvP[VaughanRowToTempsIndex[ii],1],
                                             sig = muSigSurvP[VaughanRowToTempsIndex[ii],2],
                                             N   = VaughanSSPost[ii])
    }
    ## LIKELIHOOD COMBINED LARVAE-PUPAE DATA
    pDevDead1[1:nLP[1],1:2] <- travellingWave_doubleParas2pDevDead(paras1=parasL[LPRowToTempIndex[1],1:3], paras2=parasP[LPRowToTempIndex[1],1:3], res1=resL, res2=resP, nLP[1])
    pDevDead2[1:nLP[2],1:2] <- travellingWave_doubleParas2pDevDead(paras1=parasL[LPRowToTempIndex[2],1:3], paras2=parasP[LPRowToTempIndex[2],1:3], res1=resL, res2=resP, nLP[2])
    pDevDead3[1:nLP[3],1:2] <- travellingWave_doubleParas2pDevDead(paras1=parasL[LPRowToTempIndex[3],1:3], paras2=parasP[LPRowToTempIndex[3],1:3], res1=resL, res2=resP, nLP[3])
    pDevDead4[1:nLP[4],1:2] <- travellingWave_doubleParas2pDevDead(paras1=parasL[LPRowToTempIndex[4],1:3], paras2=parasP[LPRowToTempIndex[4],1:3], res1=resL, res2=resP, nLP[4])
    pDevDead5[1:nLP[5],1:2] <- travellingWave_doubleParas2pDevDead(paras1=parasL[LPRowToTempIndex[5],1:3], paras2=parasP[LPRowToTempIndex[5],1:3], res1=resL, res2=resP, nLP[5])
    pDevDeadOrCensored1[1:nLP1[1]] <- twoColumnMatrix2ShortVector(pDevDead1[1:nLP[1],1:2])
    pDevDeadOrCensored2[1:nLP1[2]] <- twoColumnMatrix2ShortVector(pDevDead2[1:nLP[2],1:2])
    pDevDeadOrCensored3[1:nLP1[3]] <- twoColumnMatrix2ShortVector(pDevDead3[1:nLP[3],1:2])
    pDevDeadOrCensored4[1:nLP1[4]] <- twoColumnMatrix2ShortVector(pDevDead4[1:nLP[4],1:2])
    pDevDeadOrCensored5[1:nLP1[5]] <- twoColumnMatrix2ShortVector(pDevDead5[1:nLP[5],1:2])
    mullensLP17[1:nLP1[1]] ~ dmulti(size=SampSizeLP[1], prob=pDevDeadOrCensored1[1:nLP1[1]])
    mullensLP20[1:nLP1[2]] ~ dmulti(size=SampSizeLP[2], prob=pDevDeadOrCensored2[1:nLP1[2]])
    mullensLP23[1:nLP1[3]] ~ dmulti(size=SampSizeLP[3], prob=pDevDeadOrCensored3[1:nLP1[3]])
    mullensLP27[1:nLP1[4]] ~ dmulti(size=SampSizeLP[4], prob=pDevDeadOrCensored4[1:nLP1[4]])
    mullensLP30[1:nLP1[5]] ~ dmulti(size=SampSizeLP[5], prob=pDevDeadOrCensored5[1:nLP1[5]])
    ## SHAPE CONSTRAINTS
    constraintDataL ~ dconstraint(1 ==
                                  (1 < min(alphasL[1:nTempsL, 1:2])) *
                                  dUnimodal(parasL[1:nTempsL, 1])    *
                                  dUnimodal(parasL[1:nTempsL, 3])    *
                                  dUnimodal(P01_L[1:nTempsL])        *
                                  dUnimodal(P99_L[1:nTempsL]))
    constraintDataP ~ dconstraint(1 ==
                                  (1 < min(alphasP[1:nTempsP, 1:2])) *
                                  dUnimodal(parasP[1:nTempsP, 1])    *
                                  dUnimodal(parasP[1:nTempsP, 3])    *
                                  dUnimodal(P01_P[1:nTempsP])        *
                                  dUnimodal(P99_P[1:nTempsP]))
})

## ASSEMBLE NIMBLE MODEL
## the error message concerning nodes at 17 & 20 is fixed below when sourcing fixed "initialParasGivenRes.R"
LPModel <- nimbleModel(LPCode,
                       constants=Constants,
                       inits=Inits,
                       data=Data, debug=FALSE, check=TRUE)
 

###################
## EXAMINE MODEL ##
LPModel$getNodeNames()
LPModel$resL; LPModel$resP
LPModel$logit_parasL
LPModel$parasL
LPModel$logit_parasP
LPModel$parasP
LPModel$getDependencies(c('resL'))
LPModel$getDependencies(c('resP'))
LPModel$getDependencies(c('logit_parasL'))
LPModel$getDependencies(c('logit_parasP'))
LPModel$calculate()

#################################################
## Ensure deterministic nodes are up to date   ##
## Update log-likelihood of all dependent data ##
LPModel$simulate(LPModel$getDependencies("parasL"))       ## Update deterministic nodes only (logit_parasL not included)
LPModel$simulate(LPModel$getDependencies("parasP"))       ## Update deterministic nodes only (logit_parasP not included)
LPModel$simulate(LPModel$getDependencies("VaughanSSPre")) ## Update deterministic nodes only 

LPModel$calculate(LPModel$getDependencies("parasL"))     ## Update logProb of all dependent nodes
LPModel$calculate(LPModel$getDependencies("resL"))   
LPModel$calculate(LPModel$getDependencies("parasP"))     ## Update logProb of all dependent nodes
LPModel$calculate(LPModel$getDependencies("resP"))   
LPModel$calculate()

## ###########################################################
## SET resL & resP FROM SHELL SCRIPT INPUT HERE & CALCULATE ## 
if (UseScript) {
    LPModel$resL <- iResL
    LPModel$resP <- iResP 
}

resL <- LPModel$resL
resP <- LPModel$resP
nimPrint("resL=", LPModel$resL[1])
nimPrint("resP=", LPModel$resP[1])

###############################################
## Update initial values to match resL and resP
setwd(baseDir)
source("initialParasGivenRes.R") 
(LP1 <- LPModel$calculate())

(LP2 <- sum(LPModel$logProb_logit_parasL) + sum(LPModel$logProb_logit_parasP) + 
     sum(LPModel$logProb_MullensPupaeMuSig) + sum(LPModel$logProb_VaughanPupaeMuSig) + 
     sum(LPModel$logProb_VaughanSSPre) + 
     sum(LPModel$logProb_VaughanPseudoData) + 
     sum(LPModel$logProb_constraintDataL) + sum(LPModel$logProb_constraintDataP) + 
     sum(LPModel$logProb_mullensLP17) + sum(LPModel$logProb_mullensLP20) + sum(LPModel$logProb_mullensLP23) + sum(LPModel$logProb_mullensLP27) + sum(LPModel$logProb_mullensLP30))
nimPrint("LP1 - LP2 = ", LP1-LP2)


## Once resL & resP are defined at the desired value, initialise the output file
setwd(baseDir)
nimPrint("Working directory is", getwd())
source("initialise_output.R")
nimPrint("Working directory is", getwd())


## ##############################################
## COMPILE MODEL (required for compiling MCMC) ##
cLPModel <- compileNimble(LPModel) 

## ##############################################
## CONFIGURE ADAPTIVE PARALLEL TEMPERING (APT) ##
mcmcConfLP <- configureMCMC(LPModel, nodes=NULL, control=list(temperPriors=FALSE))
## Monitors
mcmcConfLP$getMonitors()
mcmcConfLP$resetMonitors()
mcmcConfLP$addMonitors (c('resL', 'resP', 'parasL', 'parasP', 'VaughanSSPre', 'VaughanSSPost', 'P01_L', 'P99_L', 'P01_P', 'P99_P')) 
mcmcConfLP$addMonitors2(c("logProb_logit_parasL", "logProb_logit_parasP",
                          "logProb_MullensPupaeMuSig", "logProb_VaughanPupaeMuSig", "logProb_VaughanSSPre", "logProb_VaughanPseudoData",
                          "logProb_mullensLP17", "logProb_mullensLP20", "logProb_mullensLP23", "logProb_mullensLP27", "logProb_mullensLP30" ))

## Set Samplers
## Round 1: Imputation
mcmcConfLP$addSampler(target="VaughanSSPre", type="sampler_RW_multinomial_tempered", control=list(useTempering=FALSE)) 
## Round 1: Blocks of 14
mcmcConfLP$addSampler(target=c("logit_parasL[1:2,1:3]","logit_parasP[1:2,1:3]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[2:3,1:3]","logit_parasP[2:3,1:3]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[3:4,1:3]","logit_parasP[3:4,1:3]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[4:5,1:3]","logit_parasP[4:5,1:3]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[5,  1:3]","logit_parasP[5:6,1:3]"), type="sampler_RW_block_tempered")
## Round 1: Blocks of 10
mcmcConfLP$addSampler(target=c("logit_parasL[1:2,1:2]","logit_parasP[1:2,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[2:3,1:2]","logit_parasP[2:3,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[3:4,1:2]","logit_parasP[3:4,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[4:5,1:2]","logit_parasP[4:5,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[5,  1:2]","logit_parasP[5:6,1:2]"), type="sampler_RW_block_tempered")
## Round 1: Univariate
mcmcConfLP$addSampler(target="logit_parasL[1,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[2,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[3,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[4,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[5,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[1,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[2,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[3,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[4,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[5,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[1,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[2,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[3,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[4,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasL[5,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[1,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[2,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[3,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[4,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[5,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[6,1]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[1,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[2,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[3,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[4,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[5,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[6,2]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[1,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[2,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[3,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[4,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[5,3]", type="sampler_RW_tempered")
mcmcConfLP$addSampler(target="logit_parasP[6,3]", type="sampler_RW_tempered")

## Round 2
mcmcConfLP$addSampler(target=c("logit_parasP[3:4, 1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[1:2, 1:2]"), type="sampler_RW_block_tempered")

## Round 3
mcmcConfLP$addSampler(target=c("logit_parasL[1:2,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[1:2,1]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[1:2,2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[1,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[2,1:2]"), type="sampler_RW_block_tempered")

## Round 4 
mcmcConfLP$addSampler(target="VaughanSSPre", type="sampler_RW_multinomial_tempered", control=list(useTempering=FALSE)) 
mcmcConfLP$addSampler(target=c("logit_parasP[3:4,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[1:2,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target="VaughanSSPre", type="sampler_RW_multinomial_tempered", control=list(useTempering=FALSE)) 
mcmcConfLP$addSampler(target=c("logit_parasL[5,1:2]"), type="sampler_RW_block_tempered")
mcmcConfLP$addSampler(target=c("logit_parasL[5,1:3]"), type="sampler_RW_block_tempered")  
##
mcmcConfLP$printSamplers()

######################### 
## BUILD & COMPILE APT ##
nTemps  <- 11 ## 
mcmcLP  <- buildAPT(mcmcConfLP, Temps=exp(seq(0,log(100),l=nTemps)), monitorTmax=TRUE, ULT=1E4) 
cMcmcLP <- compileNimble(mcmcLP, project=LPModel) 


#############################
## MCMC - ADAPTIVE BURN-IN ##
print("####################################")
print("starting adaptive burn-in WHILE loop")
print("####################################")

THIN           <- 1
cMcmcLP$thin  <- THIN
cMcmcLP$thin2 <- THIN
TuneTemper     <- c(1, 1)
meanL          <- cLPModel$calculate()
meanL_previous <- -Inf
ii             <- 0

while(meanL > meanL_previous + 2) {
    ii <- ii+1
    meanL_previous <- meanL
    print(paste0("iteration nb.", ii, "within while loop. meanL = ", meanL_previous))
    #################
    ## Short run 1 ##
    nIter <- 1E4; cMcmcLP$thinPrintTemps <- nIter / 10
    syst1 <- system.time(cMcmcLP$run(nIter,
                                      reset          = TRUE,  ## Resets the adaptive MCMC. Let's proposal distributions change direction if required.
                                      adaptTemps     = FALSE, ## Prevents temperature ladder adaptation (to avoid volatile behaviour when counter is reset)
                                      resetTempering = TRUE,  ## Resets counter used in temperature ladder adaptation
                                      printTemps     = TRUE,  ## Will print once only
                                      tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2])) 
    ## Update meanL
    nimPrint("While loop: 1st short run finished. sysT = ", syst1[3])
    samples2 <- tail(as.matrix(cMcmcLP$mvSamples2), floor(nIter/THIN)) ## LogProbs 
    meanL    <- mean(rowSums(samples2))
    nimPrint("meanL = ", meanL)
    #################
    ## Short run 2 ##
    nIter <- 1E4; cMcmcLP$thinPrintTemps <- round(nIter / 10) ## Ensures temps are only printed 10 times
    syst2 <- system.time(cMcmcLP$run(nIter,
                                      reset          = FALSE, ## Do not reset the adaptive MCMC, let adaptation continue as it is
                                      adaptTemps     = TRUE,  ## Allows temperature ladder to adjust 
                                      resetTempering = FALSE, ## Keeps the adjustments modest so avoids volatile behaviour
                                      printTemps     = TRUE,  ## Prevents verbose printing of temperature ladder updates
                                      tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
    nimPrint("While loop: 2nd short run finished. sysT = ", syst2[3])
    ## update meanL
    samples2 <- tail(as.matrix(cMcmcLP$mvSamples2), floor(nIter/THIN)) ## LogProbs 
    meanL    <- mean(rowSums(samples2))
    ## Calculate & print ESS 
    samples   <- tail(as.matrix(cMcmcLP$mvSamples),  floor(nIter/THIN))
    parasCols <- substring(colnames(samples),1,5)=="paras"
    sub       <- samples[,parasCols]
    mc        <- as.mcmc(sub)
    ESS       <- effectiveSize(mc) 
    (ESS      <- ESS[order(ESS)])
    nimPrint(ESS)     
}
print(paste0("iteration nb.", ii, "within while loop. meanL = ", meanL))


## ######################################################################
## Extract MCMC samples & LogProbs from last round of adaptive burn-in ##
samples     <- tail(as.matrix(cMcmcLP$mvSamples),     floor(nIter/THIN)) ## Sampled parameters T=1
samplesTMax <- tail(as.matrix(cMcmcLP$mvSamplesTmax), floor(nIter/THIN)) ## Sampled parameters T=Tmax
samples2    <- tail(as.matrix(cMcmcLP$mvSamples2),    floor(nIter/THIN)) ## LogProbs
print("Finished extracting model values samples & samples2")

## Check dimensions match
if( dim(samples)[1] != dim(samples2)[1])
    stop("nrow of samples & samples2 don't match")
print("DIMENSION OF samples & samples2 are OK")

## Optionally Plot Temperature Trajectories
if (PDF) {
    print("start plotting first PDF")
    pdf(file=paste0("LP_APT_temperature_trajectory_", nTemps, ".pdf"))
    par(mfrow=n2mfrow(1))
    tt     <- cMcmcLP$tempTraj 
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
    pdf(file=paste0("LP_APT_paras-trajectories_nTemps", nTemps, ".pdf"))
    plot(mc)
    dev.off()
}

if(PDF) {
    mcTMax <- as.mcmc(tail(samplesTMax, nIter/THIN))
    pdf(file=paste0("LP_APT_paras-trajectories_TMAX_nTemps", nTemps, ".pdf"))
    plot(mcTMax)
    dev.off()
}
print("Finished last PDFs")


## ###########################################
## ### DIAGNOSTICS - Effective Sample Size ###
print("starting ESS calculation")

## Filter columns
parasCols <- substring(colnames(samples),1,5)=="paras"
sub       <- samples[,parasCols]
## Check: these should match
matrix(sub[nrow(sub),1:15], nrow=5)
cLPModel$parasL
##
matrix(sub[nrow(sub),-(1:15)], nrow=6)
cLPModel$parasP
##
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
    mv2 <- tail(as.matrix(cMcmcLP$mvSamples2), nIter/THIN/5)
    LP <- rowSums(mv2)
    tail(LP, 1)
    cLPModel$calculate()
    pdf(file=paste0("LP_APT_LogProb_nTemps", nTemps, ".pdf"))
    plot(LP, typ="l")
    dev.off()
}


##########################
##########################
#### MCMC - INFERENCE ####
##########################
##########################

print("Starting MCMC for inference")
THIN  <- ceiling(SuggestedThinning * 3)
nIter <-  1E4 * THIN
cMcmcLP$thin <- cMcmcLP$thin2 <- THIN
sysT <- system.time(cMcmcLP$run(nIter,
                                reset          = FALSE, ## Keeps proposal distributions as they are (initially)
                                adaptTemps     = FALSE, ## Fix the temperature ladder
                                resetTempering = FALSE, ## Do not reset the counter for temperature ladder adaptation
                                printTemps     = TRUE,  ## Do print the temperature ladder once
                                tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
nimPrint("MCMC finished. sysT = ", sysT[3])


## Get again samples & some info to store in output
samples     <- tail(as.matrix(cMcmcLP$mvSamples),  floor(nIter/THIN))
samples2    <- tail(as.matrix(cMcmcLP$mvSamples2), floor(nIter/THIN))
samplesTMax <- tail(as.matrix(cMcmcLP$mvSamplesTmax), floor(nIter/THIN))
parasCols   <- substring(colnames(samples),1,5)=="paras" ## filter columns
sub         <- samples[,parasCols]
mc          <- as.mcmc(sub, floor(nrow(samples)/2))
ESS         <- effectiveSize(mc) 
nimPrint(ESS)

## Dimension check
if( dim(samples)[1] != dim(samples2)[1])
    stop("nrow of samples & samples2 don't match")
print("DIMENSION OF samples & samples2 are OK")


if (FALSE) { ## TRUE
    ## CODA
    samples <- as.matrix(cMcmcLP$mvSamples)
    save(samples, file="samples.Rdata")
    ##    
    samples <- tail(samples, nIter / cMcmcLP$thin)
    dim(samples)
    mc      <- as.mcmc(samples)
    summary(mc)
    plot(mc)
    tail(samples)
    table(samples[,"resL"])
    hist(samples[,"resP"])
    crosscorr(mc)
    crosscorr.plot(mc)
    autocorr.plot(mc)
    autocorr(mc, lags=2)
    ##    
    samples2 <- as.matrix(cMcmcLP$mvSamples2)
    dim(samples2)
    samples2 <- tail(samples2, nIter)
    mc2 <- as.mcmc(samples2)
    summary(mc2)
    plot(mc2)
    crosscorr.plot(mc2)
    crosscorr(mc2)
    autocorr.plot(mc2)
}

## Get friendly colnames for samples in csv file
colnames(samples) <- c(paste("P01_L" , temps[MullensRowToTempsIndex], sep='_'),
                       paste("P01_P" , temps, sep='_'),
                       paste("P99_L",  temps[MullensRowToTempsIndex], sep='_'),
                       paste("P99_P",  temps, sep='_'),
                       paste0("VaughanSSPost",temps[VaughanRowToTempsIndex]),
                       paste0("VaughanSSPre" ,temps[VaughanRowToTempsIndex]),
                       paste0("muL"  , temps[MullensRowToTempsIndex]),
                       paste0("scL"  , temps[MullensRowToTempsIndex]),
                       paste0("survL", temps[MullensRowToTempsIndex]),
                       paste0("muP"  , temps),
                       paste0("scP"  , temps),
                       paste0("survP", temps),
                       "resL",
                       "resP")
colnames(samplesTMax) <- colnames(samples)


## ###########################################
## Set output directory & Save MCMC samples ##
setwd(mcmcDir)
nimPrint("Current working directory is: ", getwd())
nimPrint("MCMC output going to file: ", mcmcFile)
nimPrint("LogProb output going to file: ", mcmcFile2)
nimPrint("TMax MCMC output going to file: ", mcmcFile3)

print("START WRITING SAMPLES")
write.table(samples,     file=mcmcFile,   col.names=TRUE, row.names = FALSE, sep=",")
write.table(samples2,    file=mcmcFile2,  col.names=TRUE, row.names = FALSE, sep=",")
write.table(samplesTMax, file=mcmcFile3,  col.names=TRUE, row.names = FALSE, sep=",")

print("Finished writing samples to file")


##########################
## Save MCMC parameters ##
write(paste("iterations  =", nIter),          file = parametersFile, append = FALSE)
write(paste("thin        =", cMcmcLP$thin),   file = parametersFile, append = TRUE)
write(paste("thin2       =", cMcmcLP$thin2),  file = parametersFile, append = TRUE)
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
    TXT <- try(mcmcConfLP$printSamplers(ii))
    if (class(TXT) == "try-error") 
        STOP <- TRUE
    else 
        cat(TXT)    
}
cat("\n")
ESS
sink()


print("%%%%%%% FINISHED %%%%%%%%")

