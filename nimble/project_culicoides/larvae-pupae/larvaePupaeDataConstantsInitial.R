#####################
## MODEL CONSTANTS ##
thresh  <- 1E-6 ## Threshold probability at which right censor is applied. Should be set to a negligable value. 

##########################################
## READ DATA COMBINED PUPAE-LARVAE DATA ##
dataDir <- "~/IPLM_paper/nimble/project_culicoides/larvae-pupae/Data/"
setwd(dataDir)
source("mullens_full_data.R")

mullensLP17 <- rep(0, tail(LP17[,"daysL2A"],1)) ## Vector for storing totalByDay for each successive day 
mullensLP20 <- rep(0, tail(LP20[,"daysL2A"],1))
mullensLP23 <- rep(0, tail(LP23[,"daysL2A"],1))
mullensLP27 <- rep(0, tail(LP27[,"daysL2A"],1))
mullensLP30 <- rep(0, tail(LP30[,"daysL2A"],1))

mullensLP17[LP17[,"daysL2A"]] <- LP17[,"totalByDay"] 
mullensLP20[LP20[,"daysL2A"]] <- LP20[,"totalByDay"]
mullensLP23[LP23[,"daysL2A"]] <- LP23[,"totalByDay"]
mullensLP27[LP27[,"daysL2A"]] <- LP27[,"totalByDay"]
mullensLP30[LP30[,"daysL2A"]] <- LP30[,"totalByDay"]

mullensLP17 <- c(mullensLP17, nDeadOrCensorredLP[1]) 
mullensLP20 <- c(mullensLP20, nDeadOrCensorredLP[2])
mullensLP23 <- c(mullensLP23, nDeadOrCensorredLP[3])
mullensLP27 <- c(mullensLP27, nDeadOrCensorredLP[4])
mullensLP30 <- c(mullensLP30, nDeadOrCensorredLP[5]) 

(nLP  <- as.numeric(rbind(tail(LP17,1), tail(LP20,1), tail(LP23,1), tail(LP27,1), tail(LP30,1))$daysL2A))
(nLP1 <- 1 + as.numeric(rbind(tail(LP17,1), tail(LP20,1), tail(LP23,1), tail(LP27,1), tail(LP30,1))$daysL2A))

length(mullensLP17) == nLP1[1]
length(mullensLP20) == nLP1[2]
length(mullensLP23) == nLP1[3]
length(mullensLP27) == nLP1[4]
length(mullensLP30) == nLP1[5]

(SampSizeLP <- c(sum(mullensLP17), sum(mullensLP20), sum(mullensLP23), sum(mullensLP27), sum(mullensLP30)))

tempsMullensLP <- c(17, 20, 23, 27, 30)

############################
## READ DATA SETS - PUPAE ## 
############################
source("mullens_pupae_data.R")
source("vaughan_pupae_data.R")
(nTempsMullens <- nrow(MullensPupaeData))
(nTempsVaughan <- nrow(VaughanPupaeData))

(tempsVaughanP <- VaughanPupaeData$temp)
(tempsMullensP <- MullensPupaeData$temp)
(tempsLP  <- tempsMullensLP)
(tempsP   <- sort(unique(c(tempsMullensP, tempsVaughanP))))
(nTempsLP <- length(tempsLP))
(nTempsL  <- length(tempsLP)) ## Identical in our data set
(nTempsP  <- length(tempsP))

(temps <- sort(unique(c(MullensPupaeData$temp, VaughanPupaeData$temp, tempsMullensLP))))
nTemps <- length(temps)

VaughanRowToTempsIndex  <- sapply(VaughanPupaeData$temp, function(x, temps) {which(temps == x)}, temps=temps)
VaughanRowToTempsIndex
nTempsVaughan <- nrow(VaughanPupaeData)
##
MullensRowToTempsIndex  <- sapply(MullensPupaeData$temp, function(x, temps) {which(temps == x)}, temps=temps)
MullensRowToTempsIndex 
nTempsMullens <- nrow(MullensPupaeData)
##
LPRowToTempIndex  <- sapply(tempsMullensLP, function(x, temps) {which(temps == x)}, temps=temps)
LPRowToTempIndex 
nTempsMullensLP <- length(LPRowToTempIndex )

################################
## INITIAL PARAMETERS - PUPAE ##
################################

## FROM A PUPAE ONLY INITIAL RUN
(parasP <- matrix(c(0.205395480,0.322141067,0.441154901,0.667875746,0.916065387,0.861617058,
                    0.007594670,0.005569199,0.019513633,0.092439909,0.055390690,0.120462453,
                    0.160106866,0.945772702,0.954589788,0.974832493,0.953589167,0.894751528), 6,3))
resP <- 44

## FROM A PUPAE-FIXED PRE-RUN
parasL   <- matrix(c(0.096192456, 0.102633083, 0.107934355, 0.093000000, 0.092000000, 
                     0.003698174, 0.007284897, 0.011969644, 0.080554204, 0.083454855, 
                     0.999896007, 0.992328182, 0.985853358, 0.984616774, 0.975381396), 5, 3)
resL <- 4


#######################################################################
## Check there are no zeros after day 1 in initial pDevDead matrices ##
#######################################################################
nf_muVar2alpha(nf_muSc2muVar(parasP[1,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasP[2,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasP[3,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasP[4,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasP[5,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasP[6,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasL[1,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasL[2,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasL[3,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasL[4,1:2]))
nf_muVar2alpha(nf_muSc2muVar(parasL[5,1:2]))
travellingWave_doubleParas2pDevDead(parasL[LPRowToTempIndex[1],1:3], parasP[LPRowToTempIndex[1],1:3], resL, resP, nLP[1]); nLP1
travellingWave_doubleParas2pDevDead(parasL[LPRowToTempIndex[2],1:3], parasP[LPRowToTempIndex[2],1:3], resL, resP, nLP[2]); nLP1
travellingWave_doubleParas2pDevDead(parasL[LPRowToTempIndex[3],1:3], parasP[LPRowToTempIndex[3],1:3], resL, resP, nLP[3]); nLP1
travellingWave_doubleParas2pDevDead(parasL[LPRowToTempIndex[4],1:3], parasP[LPRowToTempIndex[4],1:3], resL, resP, nLP[4]); nLP1
travellingWave_doubleParas2pDevDead(parasL[LPRowToTempIndex[5],1:3], parasP[LPRowToTempIndex[5],1:3], resL, resP, nLP[5]); nLP1



###################################
## Initial Value Checks : Alphas ##
###################################
for (ii in 1:nTempsL) {
    print( alphasL <- nf_muVar2alpha( nf_muSc2muVar ( parasL[ii,1:2] ))) ## Initial alphasL
    if (min(alphasL) <= 1)
        stop("Warning: Ensure initial alphasL are > 1")
}
for (ii in 1:nTempsP) {
    print( alphasP <- nf_muVar2alpha( nf_muSc2muVar ( parasP[ii,1:2] ))) ## Initial alphasP
    if (min(alphasP) <= 1)
        stop("Warning: Ensure initial alphasL are > 1")
}

###############################################
## Initial Values for muSigSurv and pDevDead ## 
###############################################
muSigSurvP <- 0 * parasP
colnames(muSigSurvP) <- c("mu","sig","surv")
for (tt in 1:nTempsP) {
    muSigSurvP[tt,1:3] <- nf_TW1_MuSigSurv(paras = parasP[tt,1:3], res = resP, thresh = thresh)
}
muSigSurvP


##########################################
## TEST UNIMODALITY OF KERNEL QUANTILES ##
##########################################
if (FALSE) { ## TRUE
    (P01 <- nf_get_quantiles_vec(parasP[,1:2], P=0.01))
    (P99 <- nf_get_quantiles_vec(parasP[,1:2], P=0.99))
    plot(1:nTemps, P99, ylim=0:1, typ="l",
         main="Survival & percentiles of development kernel", 
         ylab="Development", xlab="X")
    lines(1:nTempsP, P01)
    if ( 1 != dUnimodal(P99) * dUnimodal(P01) * dUnimodal(parasP[,3]) ) 
        stop("WARNING: Unimodality constraint violated in initial conditions")
    else
        print("Unimodality constraints satisfied by initial parameters.")
    if (FALSE) { ## TRUE
        lines(1:nTempsP, parasP[,3])
        mtext("Survival", 4)
    }
    rm(P01, P99)
}



###########################
## Lists for nimbleModel ##
###########################
Constants <- list(
    thresh  = thresh,     
    ## LARVAE/PUPAE
    nLP              = nLP, 
    nLP1             = nLP1,
    nTempsL          = nTempsL,
    SampSizeLP       = SampSizeLP,
    LPRowToTempIndex = LPRowToTempIndex, 
    ## PUPAE ONLY (IMPUTATION)
    nTempsP                = nTempsP,
    MullensRowToTempsIndex = MullensRowToTempsIndex,
    VaughanRowToTempsIndex = VaughanRowToTempsIndex,
    nTempsMullens          = nTempsMullens,
    nTempsVaughan          = nTempsVaughan, 
    VaughanTSSPre          = sum(VaughanPupaeData$ImpPreMortSS), 
    ImputationPrecision    = Precision4ImputingVaughan
)

Inits <- list(
    ## LARVAE
    logit_parasL   = logit(parasL),
    parasL         = parasL,
    resL           = resL,
    alphasL        = alphasL <- t(apply(parasL[,1:2], 1, function(x) nf_muVar2alpha(nf_muSc2muVar(x[1:2])))), 
    P01_L          = qbeta(0.01, alphasL[,1], alphasL[,2]),
    P99_L          = qbeta(0.99, alphasL[,1], alphasL[,2]),
    ## PUPAE
    logit_parasP   = logit(parasP),
    parasP         = parasP,
    resP           = resP,
    alphasP        = alphasP <- t(apply(parasP[,1:2], 1, function(x) nf_muVar2alpha(nf_muSc2muVar(x[1:2])))), 
    P01_P          = qbeta(0.01, alphasP[,1], alphasP[,2]),
    P99_P          = qbeta(0.99, alphasP[,1], alphasP[,2]),
    muSigSurvP     = muSigSurvP,
    MullensSSPost  = MullensPupaeData$SampleSizePost,
    ## IMPUTATION VECTORS
    pVaughanSSPre      = rep(1/nTempsVaughan, nTempsVaughan), 
    VaughanSSPre       = as.numeric(VaughanPupaeData$ImpPreMortSS), 
    VaughanSSPost      = as.numeric(VaughanPupaeData$ImpPostMortSS), 
    VaughanSurvivalObs = myRound(VaughanPupaeData$Survival, Precision4ImputingVaughan)
)

Data <- list(
    ## LARVAE-PUPAE
    mullensLP17 = mullensLP17, 
    mullensLP20 = mullensLP20,  
    mullensLP23 = mullensLP23,  
    mullensLP27 = mullensLP27,  
    mullensLP30 = mullensLP30,    
    ## PUPAE
    MullensPupaeMuSig = as.matrix(MullensPupaeData[,c("muOb","sigOb")]), 
    VaughanPupaeMuSig = as.matrix(VaughanPupaeData[,c("muOb","sigOb")]), 
    ## IMPUTATION CONSTRAINTS
    constraintDataL   = 1, 
    constraintDataP   = 1, 
    VaughanPseudoData = rep(1, 5)
) 

