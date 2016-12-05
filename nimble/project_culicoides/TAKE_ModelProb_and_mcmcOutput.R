#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ########################################################################################################################
## Script to obtain model probability  & MCMC output from IPLM model by integrating different (fixed) resolution outputs ##
## ########################################################################################################################

rm(list=ls())
## Set directories and basic parameters
options(width=160)
library(nimble)
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
setwd(baseDir)
source("../FUNCTIONS_R.R")
source("../FUNCTIONS_NIMBLE.R")
source("FUNCTIONS_CULICOIDES.R")


## #############################
## SET THE STAGE TO WORK WITH ##
Stg  <- "gonotrophicCycle"  ## "egg" ## "larvae-pupae" ##  
CASE <- Stg ## for plotting purposes
## ############################
## Read mcmc output directories
setwd(baseDir)
stgDir <- paste(getwd(), Stg, sep='/') 
setwd(stgDir) 
MCMC <- "mcmc"
setwd(MCMC)
## Find mcmc results
(mcmcDirs <- system("ls -1t|head -n 750", TRUE)) ## 750 is the maximum we could find (L-P case)
## Define lists to use
postLikByRes <- list() ## takes samples2 at every res
samplesList  <- list()  ## takes samples at every res

## Read 1st mcmc to take ncols of samples & samples2
tmpDir <- paste(getwd(),mcmcDirs[3],sep="/")
setwd(tmpDir)
(ncolLogLik  <- ncol(read.csv("mcmc_samples2.csv")))
samples      <- read.csv("mcmc_samples.csv")
(nParams     <- ncol(samples))
(colnSamples <- colnames(samples))
dim(samples)
rm(samples)

## ############################
## READ OUTPUT FOR EVERY MODEL. 
Res <- numeric()
if(Stg=="larvae-pupae" ) {
    resLP <- ResP <- ResL <- numeric()
    rm(Res)
}

(length(mcmcDirs))

## ######################################
## READ MCMC outpus - THIS CAN BE SLOW ##
imcmc <- 0
for (mcmcDir in mcmcDirs) {
    setwd(baseDir)
    setwd(stgDir)
    setwd(MCMC)
    newdir <- paste(getwd(),mcmcDir,sep="/")
    setwd(newdir)
    if(length(dir())>1) {
        imcmc <- imcmc + 1
        ## Read data
        postLikByRes[[imcmc]] <- as.matrix(read.csv("mcmc_samples2.csv"))
        samplesList[[imcmc]]  <- as.matrix(read.csv("mcmc_samples.csv"))
        if(Stg=="larvae-pupae") {
            (resLPchar    <- tail(strsplit(getwd(), "_")[[1]], 2))
            (resLPnum     <- as.numeric(resLPchar))
            (resLP[imcmc] <- paste0(resLPnum, collapse="_")) ## to append as colnames to modelProbs
            (ResL[imcmc]  <- resLPnum[1])
            (ResP[imcmc]  <- resLPnum[2])
        } else {
            Res[imcmc]    <- samplesList[[imcmc]][1, nParams] ## store mcmcDir resolution
        }
    }
}
    
imcmc


## ######################
## Take model probability
modelProbs <- mProb(L=postLikByRes) 
nlines     <- nrow(modelProbs)
## Check modelProbs sums to one by line
ans <- apply(modelProbs, 1, function(x) sum(x))
all(round(ans,13)==1) ## round, otherwise it gives FALSE because of rounding error

## Attach res to every column of modelProbs
if(Stg=="larvae-pupae" ) {
    colnames(modelProbs) <- resLP
} else {
    colnames(modelProbs) <- as.character(Res)
}

## #########################################################
## Store most likely posteriors based on model probabilities

## Build an index vector 
indxMcmcOrder <- apply(modelProbs,1, function(x) rcat(1,x))
table(indxMcmcOrder) ## note, the order doesn't correspond to res value but to mcmcmDirs order
post <- matrix(0, nlines, nParams)
for ( i in 1:nlines) 
    post[i,] <- samplesList[[indxMcmcOrder[i]]][i, ]
##
colnames(post) <- colnSamples ## this preserves original colnames 

## Put resolution output apart for plotting resolution histograms 
if(Stg=="larvae-pupae" ) {
    resSamples <- post[ ,c((ncol(post)-1), ncol(post))]
    dim(resSamples)
} else {
    resSamples <- post[,ncol(post)]
}



## ###########################################
## Set output directory & Save MCMC samples ##
setwd(baseDir)
setwd(paste(getwd(), "mcmcOutput4plots", sep='/') )

mcmcFile     <- paste0('mcmcSamples_', Stg ,'.csv')
modProbsFile <- paste0('modelProbs_', Stg ,'.csv')
resSampFile  <-  paste0('resSamples_', Stg, '_', nlines,'.csv')

print(paste0("Current working directory is: ", getwd()))

print(paste0("MCMC output going to file: ", mcmcFile))
print(paste0("Model probabilities going to file: ", modProbsFile))

    
print("START WRITING MCMC SAMPLES")
write.table(post, file=mcmcFile, col.names=TRUE, row.names = FALSE, sep=",")

print("START WRITING MODEL PROBABILITIES")
write.table(modelProbs, file=modProbsFile, col.names=TRUE, row.names = FALSE, sep=",")

print("START WRITING RESOLUTION SAMPLE TO USE FOR HISTOGRAMS")
write.table(resSamples, file=resSampFile, col.names=TRUE, row.names = FALSE, sep=",")

