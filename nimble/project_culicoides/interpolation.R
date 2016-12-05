## ####################################################
## GENERATE INTERPOLATION FOR TWO CASES: CLM  & IPLM ## 

## FOR CLM,  take samples from mcmc dirs corresponding to res=1 for any stage
## FOR IPLM, take samples from mcmcSamples_[Stg].csv
## FOR MAP,  take samples from mcmc dirs corresponding to res_MAP for any stage
rm(list=ls())

options(width=165)
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
setwd(baseDir)
library(nimble) 
library(latex2exp)

source("../FUNCTIONS_NIMBLE.R")
source("../FUNCTIONS_R.R")
source("FUNCTIONS_CULICOIDES.R")

nLines <- 1000 ## Nb. of lines used to generate interpolation

## ############################
## SET THE CASE TO WORK WITH ##
CASE <- "IPLM" ## "CLM" ## "MAP"

## #############################
## SET THE STAGE TO WORK WITH ##
Stg <- "larvae-pupae" ## "egg" ##"gonotrophicCycle"
STG <- Stg            ## for plotting purposes & data filename
## #############################
## IF larvae-pupae, CHOOSE TOO:
if (Stg=="larvae-pupae") 
    subStg <- "PUPA" ##  "LARVA"

##################################
## SET DIRECTORY & READ SAMPLES ##
source("setDirectories.R")
## CHECK outputs
print(getwd())
dim(samples); class(samples) ## must be matrix
head(samples,5)
summary(samples)

if (Stg=="larvae-pupae") {
    summary(samplesL); summary(samplesP)
    RESL <- samples[,ncol(samples)-1]
    RESP <- samples[,ncol(samples)]
    table(RESL); table(RESP)
} else {
    RES <- samples[,ncol(samples)]
    table(RES)
}

## ################################################################################
## CREATE SUBSET LISTS FROM TAIL(SAMPLES) THAT GROUP PARAMETERS AT EACH TEMPERATURE
if (Stg=="larvae-pupae") {
    if (subStg=="LARVA") {
        subsets <- tailSamples2subsetList(samples=samplesL, nLines=nLines, parNames=parNamesL)
        temps   <- tempsL
    } else {
        subsets <- tailSamples2subsetList(samples=samplesP, nLines=nLines, parNames=parNamesP)
        temps   <- tempsP
    }
} else {
    ## For Eggs or GC:
    subsets <- tailSamples2subsetList(samples=samples, nLines=nLines, parNames=parNames)
}

## Re-define subsets for all parameters
if (Stg=="egg" | Stg=="larvae-pupae") { ## dim of subsets elements: nLines x temps, except for resVec 
    P01Mat  <- subsets[[1]]; head(P01Mat,2) ## quantile 0.01 % at empirical temperatures
    P99Mat  <- subsets[[2]]; head(P99Mat,2) ## quantile 0.99 % at empirical temperatures
    muMat   <- subsets[[3]]; head(muMat,2)  ## mu at empirical temperatures
    scMat   <- subsets[[4]]; head(scMat,2)  ## sc at empirical temperatures
    surMat  <- subsets[[5]]; head(surMat,2) ## surv at empirical temperatures
    resVec  <- as.numeric(subsets[[6]]); head(resVec,10) ## res, length: nLines
    if(CASE=="MAP")
        res <- resVec[1]
}

if (Stg=="gonotrophicCycle") {
    EFecMat <- subsets[[1]]; head(EFecMat,2) ## EFec at empirical temperatures
    P01Mat  <- subsets[[2]]; head(P01Mat,2) 
    P99Mat  <- subsets[[3]]; head(P99Mat,2) 
    muMat   <- subsets[[4]]; head(muMat,2) 
    scMat   <- subsets[[5]]; head(scMat,2) 
    surMat <- subsets[[6]]; head(surMat,2) 
    resVec  <- as.numeric(subsets[[7]]); head(resVec,10) 
    if(CASE=="MAP")
        res <- resVec[1]
}


####################
## OBTAIN SPLINES ##
setwd(baseDir)
source("getSplines.R")

setwd(baseDir)
setwd(paste(getwd(), "mcmcOutput4plots", sep= "/"))

#################
## Save Output ##
if (Stg=="larvae-pupae") {
    save(Stg, nLines, resVec, subsets, subStg,
         X, a1X, a2X, surX,
         q01X, q99X, meanX, 
         q01.CI95, q99.CI95, mean.CI95, 
         file=paste0("interp_", subStg, "_", CASE, "_", nLines, "Lines", ".Rdata"))
    ## 
    save(Stg, nLines, resVec, subStg,
         X.15.25, a1X.15.25, a2X.15.25, surX.15.25,
         X.15.30, a1X.15.30, a2X.15.30, surX.15.30,
         file=paste0("interp_cosine_", subStg, "_", CASE, "_", nLines, "Lines",".Rdata"))
} else if (Stg=="gonotrophicCycle") {
    save(Stg, nLines, resVec, subsets,
         X, a1X, a2X, surX, fecX, 
         q01X, q99X, meanX, 
         q01.CI95, q99.CI95, mean.CI95, 
         file=paste0("interp_", Stg, "_", CASE, "_", nLines, "Lines", ".Rdata"))
    ## 
    save(Stg, nLines, resVec, 
         X.15.25, a1X.15.25, a2X.15.25, surX.15.25, fecX.15.25,
         X.15.30, a1X.15.30, a2X.15.30, surX.15.30, fecX.15.30,
         file=paste0("interp_cosine_", Stg, "_", CASE, "_", nLines, "Lines",".Rdata"))
} else { ## Stg=="egg"
    save(Stg, nLines, resVec, subsets,
         X, a1X, a2X, surX, 
         q01X, q99X, meanX, 
         q01.CI95, q99.CI95, mean.CI95, 
         file=paste0("interp_", Stg, "_", CASE, "_", nLines, "Lines", ".Rdata"))
    ## 
    save(Stg, nLines, resVec,
         X.15.25, a1X.15.25, a2X.15.25, surX.15.25,
         X.15.30, a1X.15.30, a2X.15.30, surX.15.30,
         file=paste0("interp_cosine_", Stg, "_", CASE, "_", nLines, "Lines",".Rdata"))
}



