######################
## Vaughan Egg Data ##
######################
VaughanEgg <- as.data.frame(
    cbind(temp          = c(20,23,27,30,35),
          meanOb        = c(2.65, 2.695833, 2.5583333, 2.120833, 2.379167),
          sigOb         = c(0.3416667, 0.325, 0.1708333, 0.5375, 0.325),
          ESampSizePre  = c(811, 811, 811, 811, 811),
          ESampSizePost = c(186.53, 514.174, 476.868, 551.48, 143.547),
          Surv          = c(0.23, 0.634, 0.588, 0.68, 0.177)))
temps          <- sort(unique(c(VaughanEgg$temp)))
nTempsVaughan  <- length(VaughanEgg[["temp"]])
VaughanMuSig   <- as.matrix(VaughanEgg[,c("meanOb","sigOb")])
VaughanSurvObs <- VaughanEgg[["Surv"]]
VaughanESSPre  <- VaughanEgg[["ESampSizePre"]]         ## ESS = Expected sample size
VaughanESSPost <- round(VaughanEgg[["ESampSizePost"]]) ## ESS = Expected sample size
VaughanTSSPre  <- sum(VaughanESSPre)                   ## TSS = Total sample size
VaughanSurvInd <- 1 + 0 * VaughanSurvObs               ## Ind = Indicator
pVaughanSSpre  <- rep(1/nTempsVaughan, nTempsVaughan)  ## Probabilities for imputing sample sizes

######################
## Mullens Egg Data ##
######################
MullensEgg <- as.data.frame(
    cbind(temp          = c(17, 20, 23, 27, 30), 
          meanOb        = c(5.633333, 3.233333, 2.895833, 2.141667, 1.833333), 
          sigOb         = c(0.26666667, 0.14166667, 0.14166667, 0.09166667, 0.30000000), 
          ESampSizePost = c(293, 293, 293, 293, 293)))
temps          <- sort(unique(c(temps, MullensEgg$temp)))
nTempsMullens  <- length(MullensEgg[["temp"]])
MullensMuSig   <- as.matrix(MullensEgg[,c("meanOb","sigOb")])
MullensESSPost <- MullensEgg[["ESampSizePost"]] ## ESS = Expected sample size
MullensTSSPost <- sum(MullensESSPost)           ## TSS = Total sample size
pSSPostMullens <- rep(1, nTempsMullens) / nTempsMullens

##########################
## Insectarium Egg Data ##
##########################
load("~/IPLM_paper/nimble/project_culicoides/egg/InsectariumData.Rdata") 
## InsectariumEggs$T15$eggs    ## SAMPLE SIZES PER EXPERIMENTAL BATCH
## InsectariumEggs$T15$rawData ## CUMULATIVE NUMBERS OF MATURATIONS OVER SUCCESSIVE DAYS. NA represents censoring.

(InsectariumEgg15 <- as.matrix(InsectariumEggs$T15$rawData[,-1])) ## Remove day 0
InsectariumEgg15[is.na(InsectariumEgg15)] <- -99                  ## Replace NAs with a -ve "missing data" indicator
InsectariumEgg15

(sampSizesInsectariumEgg15 <- InsectariumEggs$T15$eggs)

InsectariumEggTemps <- 15
(temps <- sort(unique(c(temps, InsectariumEggTemps))))

nTempsInsectarium <- length(InsectariumEggTemps) ## temperature not present in either Mullens' or Vaughan studies


##############################################
## Indexing for linking data to temperature ##
##############################################
temps
nTempsE <- length(temps)
## VaughanRowOfTemp ## Row index of Vaughan data at a given temperature
VaughanRowToTempsIndex  <- sapply(VaughanEgg$temp, function(x, temps) {which(temps == x)}, temps=temps)
VaughanRowToTempsIndex  
MullensRowToTempsIndex  <- sapply(MullensEgg$temp, function(x, temps) {which(temps == x)}, temps=temps)
MullensRowToTempsIndex  
InsectariumToTempsIndex <- sapply(InsectariumEggTemps, function(x, temps) {which(temps == x)}, temps=temps)
InsectariumToTempsIndex


#####################
## MODEL CONSTANTS ##
#####################
P      <- c(0.01, 0.99) ## Percentiles for assessing shape constraint on kernels as a function of T
thresh <- 1E-6          ## Threshold probability at which right censor is applied. Should be set to a negligible value. 

########################
## INITIAL PARAMETERS ##
########################
(paras  <- cbind(mu   = seq(0.45, 0.46, l=nTempsE),  ## Expected development
                 sc   = seq(0.3,  0.31, l=nTempsE),  ## Scale: variance = scale * Mu * (1 - Mu)
                 surv = seq(0.98, 0.99, l=nTempsE))) ## Daily survival
resE <- 32 ## 6 

###############################################
## Initial Values for muSigSurv and pDevDead ## 
###############################################
muSigSurv <- 0 * paras
colnames(muSigSurv) <- c("mu","sig","surv")
for (tt in 1:nTempsE) {
    muSigSurv[tt,1:3] <- nf_TW1_MuSigSurv(paras = paras[tt,1:3], res = resE, thresh = thresh)
}
for (tt in InsectariumToTempsIndex) {    
    (pDevDead <- travellingWave_paras2pDevDead(paras = paras[tt,1:3], res = resE, nIter = ncol(InsectariumEgg15)))
}
muSigSurv
pDevDead

#######################
## Check on pDevDead ##
#######################
if (min(pDevDead[,1]) > 0) {
    print("Check passed")
} else {
    stop("For initial values ensure all pDevDead[,1] > 0")
}

####################################
## Initial Value Checks : Vaughan ##
####################################
for (ii in 1:nTempsVaughan) {
    print(as.numeric(dObMuSig(VaughanMuSig[ii,1:2],
                              mu=muSigSurv[VaughanRowToTempsIndex[ii],1],
                              sig=muSigSurv[VaughanRowToTempsIndex[ii],2],
                              N = VaughanESSPost[ii],log=TRUE)))
}

####################################
## Initial Value Checks : Mullens ##
####################################
for (ii in 1:nTempsMullens) {
    print(as.numeric(dObMuSig(MullensMuSig[ii,1:2],
                              mu=muSigSurv[MullensRowToTempsIndex[ii],1],
                              sig=muSigSurv[MullensRowToTempsIndex[ii],2],
                              N = MullensESSPost[ii],log=TRUE)))
}

########################################
## Initial Value Checks : Insectarium ##
########################################
for (ii in 1:nTempsInsectarium) {    
    pDevDead <- travellingWave_paras2pDevDead (paras = paras[InsectariumToTempsIndex[ii], 1:3], res = resE, nIter = ncol(InsectariumEgg15))
    for (i in 1:nrow(InsectariumEgg15)) {
        print(dMaturation(InsectariumEgg15[i,], sSize = sampSizesInsectariumEgg15[i], pDev = pDevDead[, 1], log=TRUE))
    }
}

###################################
## Initial Value Checks : Alphas ##
###################################
for (ii in 1:nTempsE) { 
    print(alphas <- nf_muVar2alpha( nf_muSc2muVar ( paras[ii,1:2] ))) ## Initial alphas
    if (min(alphas) <= 1)
        stop("Warning: Ensure initial alphas are > 1")
}

##########################################
## TEST UNIMODALITY OF KERNEL QUANTILES ##
##########################################
if (FALSE) { ## TRUE
    (P01 <- nf_get_quantiles_vec(paras[,1:2], P=P[1]))
    (P99 <- nf_get_quantiles_vec(paras[,1:2], P=P[2]))
    plot(1:nTempsE, P99, ylim=0:1, typ="l",
         main="Survival & percentiles of development kernel", 
         ylab="Development", xlab="X")
    lines(1:nTempsE, P01)
    if ( 1 != dUnimodal(P99) * dUnimodal(P01) * dUnimodal(paras[,3]) ) 
        stop("WARNING: Unimodality constraint violated in initial conditions")
    else
        print("Unimodality constraints satisfied by initial parameters.")
    if (FALSE) { ## TRUE
        lines(1:nTempsE, paras[,3])
        mtext("Survival", 4)
    }
    rm(P01, P99)
}

###########################
## Lists for nimbleModel ##
###########################
Constants <- list(
    thresh                  = thresh,
    temps                   = temps, 
    nTempsE                 = nTempsE,
    ## Vaughan 
    nTempsVaughan           = nTempsVaughan,
    pVaughanSSpre           = pVaughanSSpre,
    VaughanRowToTempsIndex  = VaughanRowToTempsIndex,
    VaughanTSSPre           = VaughanTSSPre,
    ## Mullens 
    nTempsMullens           = nTempsMullens,
    pSSPostMullens          = pSSPostMullens,
    MullensRowToTempsIndex  = MullensRowToTempsIndex,
    MullensTSSPost          = MullensTSSPost, 
    ## Insectarium
    nTempsInsectarium       = nTempsInsectarium, 
    nGroupsInsectarium      = nrow(InsectariumEgg15),
    nStepsInsectarium       = ncol(InsectariumEgg15),
    InsectariumToTempsIndex = InsectariumToTempsIndex    
)

Inits <- list(    
    logit_paras   = logit(paras),
    paras         = paras,
    resE          = resE,
    muSigSurv     = muSigSurv,
    ## Alternative parameterisations    
    alphas = alphas <- t(apply(paras[,1:2], 1, function(x) nf_muVar2alpha(nf_muSc2muVar(x[1:2])))), 
    P01    = qbeta(0.01, alphas[,1], alphas[,2]),
    P99    = qbeta(0.99, alphas[,1], alphas[,2]),
    ## Vaughan 
    VaughanSSPre   = VaughanESSPre,
    VaughanSSPost  = VaughanESSPost,
    VaughanSurvObs = VaughanSurvObs,
    ## Mullens 
    MullensSSPost = MullensESSPost,
    ## Insectarium
    sampSizesInsectariumEgg15 = sampSizesInsectariumEgg15
)

Data <- list(
    constraintData = 1, 
    ## Vaughan 
    VaughanMuSig   = VaughanMuSig,    
    VaughanSurvInd = VaughanSurvInd,
    ## Mullens 
    MullensMuSig   = MullensMuSig, 
    ## Insectarium
    InsectariumEgg15 = InsectariumEgg15
) 


