########################################
## Insectarium Gonotrophic Cycle Data ##
########################################
dataDir <- "~/IPLM_paper/nimble/project_culicoides/gonotrophicCycle/Data"
gc15 <- read.csv(paste(dataDir, "Insectarium-gc-15C-uniqueEventPerRow.csv", sep="/"), header=TRUE)
gc20 <- read.csv(paste(dataDir, "Insectarium-gc-20C-uniqueEventPerRow.csv", sep="/"), header=TRUE)
gc25 <- read.csv(paste(dataDir, "Insectarium-gc-25C-uniqueEventPerRow.csv", sep="/"), header=TRUE)
## Nb. days into a cycle of a laying or death event
dayInCycle15 <- gc15$dayInCycle
dayInCycle20 <- gc20$dayInCycle
dayInCycle25 <- gc25$dayInCycle
## Nb. eggs layed. Zero indicates a death.
eggs15 <- gc15$eggs
eggs20 <- gc20$eggs
eggs25 <- gc25$eggs
## Round zeros (death on same day as laying) up to one time step (minimum given temporal resolution of discrete time model)
dayInCycle15[dayInCycle15==0] <- 1
dayInCycle20[dayInCycle20==0] <- 1
dayInCycle25[dayInCycle25==0] <- 1
## Maximum number of time steps required per temperature
(nSteps15 <- max(dayInCycle15)) ## 23
(nSteps20 <- max(dayInCycle20)) ## 12
(nSteps25 <- max(dayInCycle25)) ##  9
(nSteps   <- c(nSteps15, nSteps20, nSteps25))
(nSteps2  <- nSteps * 2)
(nSteps21 <- nSteps * 2 + 1)
##
gonoData15 <- cbind(dayInCycle15, eggs15)
gonoData20 <- cbind(dayInCycle20, eggs20)
gonoData25 <- cbind(dayInCycle25, eggs25)
## Fecundity model data
totalEggs   <- c(sum(gonoData15[,"eggs15"]), sum(gonoData20[,"eggs20"]), sum(gonoData25[,"eggs25"]))       ## nb. eggs per temperature
totalOvipos <- c(sum(0<gonoData15[,"eggs15"]), sum(0<gonoData20[,"eggs20"]), sum(0<gonoData25[,"eggs25"])) ## nb. ovipositions per temperature
##
temps      <- c(15, 20, 25)
nTempsGC   <- 3
nObs       <- c(nrow(gonoData15), nrow(gonoData20), nrow(gonoData25))
EFecundity <- totalEggs / totalOvipos ## Initial value

##################################################
## Check all individuals die once and only once ##
##################################################
all(sort(unique(gc15$id)) == sort(unique(gc15$id[gc15$eggs==0])))
all(sort(unique(gc20$id)) == sort(unique(gc20$id[gc20$eggs==0])))
all(sort(unique(gc25$id)) == sort(unique(gc25$id[gc25$eggs==0])))
0 == sum(duplicated(gc15$id[gc15$eggs==0]))
0 == sum(duplicated(gc20$id[gc20$eggs==0]))
0 == sum(duplicated(gc25$id[gc25$eggs==0]))

#####################
## MODEL CONSTANTS ##
#####################
thresh <- 1E-6  ## Threshold probability at which right censor is applied. Should be set to a negligeable value. 

########################
## INITIAL PARAMETERS ##
########################
(paras  <- cbind(mu   = seq(0.2,  0.21, l=nTempsGC),  ## Expected development
                 sc   = seq(0.1,  0.11, l=nTempsGC),  ## Scale: variance = scale * Mu * (1 - Mu)
                 surv = seq(0.98, 0.99, l=nTempsGC))) ## Daily survival
resGC <- 7

###################################
## Initial Value Checks : Alphas ##
###################################
for (ii in 1:nTempsGC) {
    print( alphas <- nf_muVar2alpha( nf_muSc2muVar ( paras[ii,1:2] ))) ## Initial alphas
    if (min(alphas) <= 1)
        stop("Warning: Ensure initial alphas are > 1")
}

###############################################
## Initial Values for muSigSurv and pLayDie ## 
###############################################
(pLayDie15 <- travellingWave_paras2pDevDead(paras=paras[1,1:3], res=resGC, nIter=nSteps[1]))
(pLayDie20 <- travellingWave_paras2pDevDead(paras=paras[2,1:3], res=resGC, nIter=nSteps[2]))
(pLayDie25 <- travellingWave_paras2pDevDead(paras=paras[3,1:3], res=resGC, nIter=nSteps[3]))

pVec15 <- twoColumnMatrix2Vector(pLayDie15)
pVec20 <- twoColumnMatrix2Vector(pLayDie20)
pVec25 <- twoColumnMatrix2Vector(pLayDie25)

#######################
## Check on pLayDie ##
#######################
if (min(c(pLayDie15, pLayDie20, pLayDie25)) > 0) {
    print("Check passed")
} else {
    stop("For initial values ensure all pLayDie[,1] > 0")
}

##########################################
## TEST UNIMODALITY OF KERNEL QUANTILES ##
##########################################
if (FALSE) { ## TRUE
    (P01 <- nf_get_quantiles_vec(paras[,1:2], P=0.01))
    (P99 <- nf_get_quantiles_vec(paras[,1:2], P=0.99))
    plot(1:nTempsGC, P99, ylim=0:1, typ="l",
         main="Survival & percentiles of development kernel", 
         ylab="Development", xlab="X")
    lines(1:nTempsGC, P01)
    if ( 1 != dUnimodal(P99) * dUnimodal(P01) * dUnimodal(paras[,3]) ) 
        stop("WARNING: Unimodality constraint violated in initial conditions")
    else
        print("Unimodality constraints satisfied by initial parameters.")
    if (FALSE) { ## TRUE
        lines(1:nTempsGC, paras[,3])
        mtext("Survival", 4)
    }
    rm(P01, P99)
}

########################################
## Initial Value Checks : Insectarium ##
########################################
gonoData15Ori <- gonoData15
gonoData20Ori <- gonoData20 
gonoData25Ori <- gonoData25

dGonotrophicCycle (x = gonoData15Ori, pLayDie = pLayDie15, EFecundity = EFecundity[1], log = 1)
dGonotrophicCycle (x = gonoData20Ori, pLayDie = pLayDie20, EFecundity = EFecundity[2], log = 1)
dGonotrophicCycle (x = gonoData25Ori, pLayDie = pLayDie25, EFecundity = EFecundity[3], log = 1)

####################################
## Adjust first column for deaths ##
####################################
gonoData15[gonoData15[,2]==0, 1] <- gonoData15[gonoData15[,2]==0, 1] + nSteps[1]
gonoData20[gonoData20[,2]==0, 1] <- gonoData20[gonoData20[,2]==0, 1] + nSteps[2]
gonoData25[gonoData25[,2]==0, 1] <- gonoData25[gonoData25[,2]==0, 1] + nSteps[3]

LL <- 0
for (ii in 1:nObs[1]) {
    (LL <- LL + dcat(gonoData15[ii,1], pVec15[1:nSteps21[1]], log=1) + dpois(gonoData15[ii,2], EFecundity[1] * (gonoData15[ii,1] <= nSteps[1]), log=1))
}
LL 
LL <- 0
for (ii in 1:nObs[2]) { 
    LL <- LL + dcat(gonoData20[ii,1], pVec20[1:nSteps21[2]], log=1) + dpois(gonoData20[ii,2], EFecundity[2] * (gonoData20[ii,1] <= nSteps[2]), log=1) 
}
LL 
LL <- 0
for (ii in 1:nObs[3]) { 
    LL <- LL + dcat(gonoData25[ii,1], pVec25[1:nSteps21[3]], log=1) + dpois(gonoData25[ii,2], EFecundity[3] * (gonoData25[ii,1] <= nSteps[3]), log=1) 
} 
LL

###########################
## Lists for nimbleModel ##
###########################
Constants <- list(
    thresh   = thresh,
    temps    = temps, 
    nTempsGC = nTempsGC,
    nSteps   = nSteps,
    nSteps2  = nSteps2,
    nSteps21 = nSteps21,
    nObs     = nObs
)

Inits <- list(    
    logit_paras    = logit(paras),
    paras          = paras,
    resGC          = resGC,
    EFecundity     = EFecundity,
    pLayDie15      = pLayDie15, 
    pLayDie20      = pLayDie20, 
    pLayDie25      = pLayDie25,
    pVec15         = pVec15,
    pVec20         = pVec20,
    pVec25         = pVec25,
    ## Alternative parameterisations    
    alphas = alphas <- t(apply(paras[,1:2], 1, function(x) nf_muVar2alpha(nf_muSc2muVar(x[1:2])))), 
    P01    = qbeta(0.01, alphas[,1], alphas[,2]),
    P99    = qbeta(0.99, alphas[,1], alphas[,2])
)

Data <- list(
    constraintData = 1, 
    gonoData15  = gonoData15[,1], 
    gonoData20  = gonoData20[,1], 
    gonoData25  = gonoData25[,1],
    ## FECUNDITY DATA 
    totalEggs   = totalEggs,
    totalOvipos = totalOvipos
) 

