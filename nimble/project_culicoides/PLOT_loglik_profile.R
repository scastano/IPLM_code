#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ####################################################################
## Plots case study II - loglik profile in terms of resolution values #
## ####################################################################
rm(list=ls())
## Set directories and basic parameters
options(width=160)
library(latex2exp)
library(nimble)
library(colorRamps)
library(lattice)
library(scatterplot3d)
library(akima)

Stg <- "larvae-pupae" ## "gonotrophicCycle"  ## egg
PDF <- FALSE ## TRUE

baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
MCMC    <- "mcmc"
setwd(baseDir)
source("../FUNCTIONS_R.R")
source("../FUNCTIONS_NIMBLE.R")
source("FUNCTIONS_CULICOIDES.R")


## ############################
## Read mcmc output directories
setwd(baseDir)
stgDir <- paste(getwd(), Stg, sep='/') 
setwd(stgDir)
##
if (Stg=="egg") 
    source("eggDataConstantsInitial.R")
if (Stg=="gonotrophicCycle") 
    source("gonoDataConstantsInitial.R")
if (Stg=="larvae-pupae") { 
    source("larvaePupaeDataConstantsInitial.R")
    setwd(baseDir)
    setwd(paste(getwd(),"mcmcOutput4plots",sep="/"))
    resLPmat <- as.matrix(read.csv("resSamples_larvae-pupae_10000.csv"))
}

setwd(baseDir)
setwd(paste(getwd(), Stg, MCMC, sep="/"))

## Find mcmc results
(mcmcDirs <- system("ls -1t|head -n 750", TRUE)) ## 750 is the maximum we could find (L-P case)
length(mcmcDirs)

## Check all directories contain mcmc ouput 
for (mcmcDir in mcmcDirs) { 
    setwd(baseDir)
    setwd(paste(getwd(), Stg, MCMC, mcmcDir, sep="/"))
    if(is.na(system("find |grep FINISHED", TRUE)[1]))
        print(paste("WARNING", mcmcDir,"doesn't have mcmc output"))
}

## Define lists to use
postLoglikList <- list() ## takes samples2 at every resolution
samplesList    <- list() ## takes samples at every resolutionres


## ##################################################
## Read 1st mcmc to take ncols of samples & samples2
setwd(baseDir)
setwd(paste(getwd(), Stg, MCMC, sep="/"))

tmpDir <- paste(getwd(), mcmcDirs[1], sep="/")
setwd(tmpDir)
samples2 <- read.csv("mcmc_samples2.csv")
(nLines  <- nrow(samples2))
(ncolLogLik  <- ncol(samples2))
samples      <- read.csv("mcmc_samples.csv")
(nParams     <- ncol(samples))
(nLines      <- nrow(samples2)) ## it should match to line belowk
nrow(samples)
rm(samples)

## ###############################
## READ OUTPUT FOR EVERY MODEL. ##
## ############################### This can be slow.

## ###########################
## EGG OR GONOTROPHIC CYCLE ##
## ###########################
if (Stg!="larvae-pupae") { 
    imcmc <- 0
    Res <- numeric()
    for (mcmcDir in mcmcDirs) {
        setwd(baseDir)
        setwd(stgDir)
        setwd(MCMC)
        newdir <- paste(getwd(),mcmcDir,sep="/")
        setwd(newdir)
        imcmc <- imcmc + 1
        ## Read data
        postLoglikList[[imcmc]] <- as.matrix(read.csv("mcmc_samples2.csv"))
        Res[imcmc] <- as.numeric(tail(strsplit(getwd(), "_")[[1]], 1))
    }
}

if (Stg=="larvae-pupae") { 
    resLPmat <- matrix(0,length(mcmcDirs),2)
    ## READ OUTPUT FOR EVERY MODEL
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
            postLoglikList[[imcmc]] <- as.matrix(read.csv("mcmc_samples2.csv"))
            (resLPchar <- tail(strsplit(getwd(), "_")[[1]], 2))
            (resLPmat[imcmc,] <- as.numeric(resLPchar))
        }
    }
    ## Remove zero lines if there are combinations of LP lacking
    resLPmat <- resLPmat[1:imcmc,]
}


## For every element of postLoglikList, sum partial loglikMat 
nbRes     <- length(postLoglikList) ## nb of models/(resolutions for E/GC, combinations for LP) to consider
loglikMat <- matrix(0,nLines,nbRes) 

for (j in 1:nbRes) 
    loglikMat[,j] <- rowSums(postLoglikList[[j]]) ## sum loglikMat of samples2 by line
  
## Attach resolution
if (Stg!="larvae-pupae")
    colnames(loglikMat) <- Res

## Order
if (Stg=="larvae-pupae") {
    ## take mean of logllikMat
    dim(loglikMat); head(loglikMat)
    LLM <- apply(loglikMat,2, function(x) mean(x))
    ## order first column(resL) and then second column (resP)
    ordered <- cbind(resLPmat[,2], resLPmat[,1], LLM)[order(resLPmat[,1],resLPmat[,2]), ]
} else {
    ## Order LLMat in terms of resolution
    LLMat <- loglikMat[, order(as.integer(colnames(loglikMat)))]
}


## ###################################################################
## Come back to the stage directory to put figure in figures directory
if (Stg=="egg") {
    setwd(baseDir)
    stgDir <- paste(getwd(), Stg, sep='/') 
    setwd(stgDir) 
}
if (Stg=="gonotrophicCycle") {
    setwd(baseDir)
    stgDir <- paste(getwd(), Stg, sep='/') 
    setwd(stgDir) 
}
if (Stg=="larvae-pupae") {
    setwd(baseDir)
    stgDir <- paste(getwd(), Stg, sep='/') 
    setwd(stgDir) 
}

print(getwd())

## Set directory for plotting
if (!is.element("figures", dir())) {
    system("mkdir figures")
}
setwd("figures")
print(getwd())


#############################################
## Plot loglik profile in terms of resolution
HEIGHT <- 10 ## 4 inches for pdf

if (Stg=="larvae-pupae") {
    if (PDF)
        pdf("LP-Mloglik-surface.pdf")
    par.set <- list(axis.line = list(col = "transparent"), clip = list(panel = "off"),
                    axis.xlab.padding = 5)
    MAR <- c(0, 0.5, 0.2, 0.5)
    CEXAXIS  <- 1.7
    CEX      <- 1.7
    par(mar=MAR)
    reggrid <- interp(ordered[,2],ordered[,1],ordered[,3],linear=T,extrap=F,nx=15,ny=50,
                      duplicate="mean")
    x.ticks <- c(5,10,15)
    y.ticks <- c(1,seq(10,max(reggrid$y),10))
    z.ticks <- c(-2500,-2000,-1500)
    wireframe(reggrid$z,
              xlim=c(1,max(reggrid$x)), ylim=c(1,max(reggrid$y)),
              xlab=list(label="Res' Larval Dev'",rot=25, cex=CEX),
              ylab=list(label="Res' Pupal Development",rot=-35, cex=CEX),
              zlab=list(label="MPLL",rot=90, cex=CEX),
              par.settings = par.set,
              scales=list(arrows=FALSE,
                          distance=c(0.8, 1.9, 1.2),
                          tck=c(2.4, 0.6, 1.3),
                          x=list(at=x.ticks, labels=x.ticks, cex=CEXAXIS),
                          y=list(at=y.ticks, labels=y.ticks, cex=CEXAXIS),
                          z=list(at=z.ticks, labels=z.ticks, cex=1.2),
                          colorkey=FALSE, aspect=c(3.2,1), panel.aspect=0.9))
    if (PDF)
        dev.off()
}


if (Stg=="gonotrophicCycle") {
    if (PDF)
        pdf("GC-loglik-post_vs_res.pdf", height=HEIGHT, width=2*HEIGHT)
    MAR <- c(3.75, 3.85, 1.2, 0.5)
    LINEXLAB <- 2.4
    LINEYLAB <- 2.3
    CEXAXIS  <- 1.8
    CEXLAB   <- 1.9
    CEX      <- 0.2
    LWD      <- 0.5
    par(mar=MAR)
    ## remove higher resolution values (50 res-values for most of stages)
    CLIP <- -c(31:50)
    boxplot(LLMat[, CLIP], ylab="", xlab="", axes = FALSE, las=TRUE, cex=CEX, lwd=LWD)
    box(lwd=0.7)
    axis(1, tcl=-0.25, cex.axis=CEXAXIS, cex.lab=CEXLAB, padj=-0.2, lwd=LWD)
    axis(2, tcl=-0.25, cex.axis=CEXAXIS, cex.lab=CEXLAB, padj= 0.5, lwd=LWD)
    mtext("Resolution of Gonotrophic Cycle", 1, cex=CEXLAB, line=LINEXLAB)
    mtext("Posterior Log Likelihood",        2, cex=CEXLAB, line=LINEYLAB)
    if (PDF)
        dev.off()
}


if (Stg=="egg") {
    if (PDF)
        pdf("Egg-loglik-post_vs_res.pdf", height=HEIGHT, width=2*HEIGHT)
    MAR <- c(3.75, 3.85, 1.2, 0.5)
    LINEXLAB <- 2.4
    LINEYLAB <- 2.3
    CEXAXIS  <- 1.8
    CEXLAB   <- 1.9
    CEX      <- 0.2
    LWD      <- 0.5
    par(mar=MAR)
    ## remove higher resolution values (50 res-values for most of stages)
    CLIP <- -c(31:50)
    boxplot(LLMat[, CLIP], ylab="", xlab="", axes = FALSE, las=TRUE, cex=CEX, lwd=LWD)
    box(lwd=0.7)
    axis(1, tcl=-0.25, cex.axis=CEXAXIS, cex.lab=CEXLAB, padj=-0.2, lwd=LWD)
    axis(2, tcl=-0.25, cex.axis=CEXAXIS, cex.lab=CEXLAB, padj= 0.5, lwd=LWD)
    mtext("Resolution of Egg Development", 1, cex=CEXLAB, line=LINEXLAB)
    mtext("Posterior Log Likelihood",      2, cex=CEXLAB, line=LINEYLAB)
    if (PDF)
        dev.off()
}

