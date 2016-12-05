#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## #####################################################################
## Plots case study II -  Interpolated developmental rates & survival ##
## #####################################################################
rm(list=ls())
options(width=165)
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
setwd(baseDir)
library(nimble)
library(coda)
library(colorRamps)
library(latex2exp)
source("../FUNCTIONS_NIMBLE.R")
source("../FUNCTIONS_R.R")
source("FUNCTIONS_CULICOIDES.R")
## Basic parameters
nLines <- 1000 ## nb lines taken for plots 
PDF    <- TRUE ## FALSE

## #########
## DEFINE ##
STAGE <- c("egg", "larvae", "pupae", "gonotrophicCycle")
CASE  <- "MAP" ## it is always the case here


## #######################################
## PLOT INTERPOLATED KERNEL & SURVIVAL  ##
## #######################################
for (Stg in STAGE) { 
    print(paste("Starting stage", Stg))
    if(Stg=="egg")
        STG <- "EGG"
    if(Stg=="gonotrophicCycle")
        STG <- "GC"
    ##
    if(Stg=="larvae") {
        Stg    <- "larvae-pupae"
        subStg <- "LARVA"
        STG    <- subStg ## for plotting purposes
    }
    if(Stg=="pupae") {
        Stg    <- "larvae-pupae"
        subStg <- "PUPA"
        STG    <- subStg ## for plotting purposes
    }   
    setwd(baseDir)
    source("setDirectories.R")
    ## ####################################################################
    ## BURN IN & SUBSET SAMPLES, DONE BY tailSamples2subsetList function ##
    if (Stg=="larvae-pupae") {
        if (subStg=="LARVA") {
            subsets <- tailSamples2subsetList(samples=samplesL, nLines=nLines, parNames=parNamesL)
            temps <- tempsL
        } 
        if (subStg=="PUPA") {
            subsets <- tailSamples2subsetList(samples=samplesP, nLines=nLines, parNames=parNamesP)
            temps <- tempsP
        }
    } else { ## for eggs or GC:
        subsets <- tailSamples2subsetList(samples=samples, nLines=nLines, parNames=parNames)
    }
    ## Re-define subsets for all parameters
    if (Stg=="egg" | Stg=="larvae-pupae") {## dim of subsets elements: nLines x temps, except for resVec 
        P01Mat  <- subsets[[1]] ## quantile 0.01 % at empirical temperatures
        P99Mat  <- subsets[[2]] ## quantile 0.99 % at empirical temperatures
        muMat   <- subsets[[3]] ## mu at empirical temperatures
        scMat   <- subsets[[4]] ## sc at empirical temperatures
        surMat  <- subsets[[5]] ## surv at empirical temperatures
        resVec  <- as.numeric(subsets[[6]]); head(resVec,10) ## res, length: nLines
    }
    ##
    if (Stg=="gonotrophicCycle") {
        EFecMat <- subsets[[1]] ## EFec at empirical temperatures
        P01Mat  <- subsets[[2]]
        P99Mat  <- subsets[[3]]
        muMat   <- subsets[[4]]
        scMat   <- subsets[[5]]
        surMat  <- subsets[[6]]
        resVec  <- as.numeric(subsets[[7]])
    }
    ## Define res
    res <- resVec[1]
    ## #################
    ## OBTAIN SPLINES ##
    setwd(baseDir)
    source("getSplines.R")
    ##
    if(Stg==STAGE[1]) {
        ## ##########################
        ## Set directory for plotting
        setwd(stgDir)
        if (!is.element("figures", dir())) {
            system("mkdir figures")
        }
        setwd("figures")
        print(getwd())
        ## ########
        ## PLOTS ##
        ## ########
        ##
        ## ########################
        ## Plot interpolated kernel
        if(PDF)
            pdf(paste0("Interp_MAP_nLines", nLines, ".pdf", sep=""),
                height=14, width=31)
        Mar=c(.5, 5, 0.5, 1.4)
        Omi=c(1.2, 0.9,.6, 0.1)
        cAxis <- 3.9
        cLab  <- 3
        cMain <- 3.5
        Padj  <- 1.2
        layout(matrix(1:8, 2, 4, byrow=FALSE), heights=c(3,2.5), widths=c(1,1)) 
        par(cex.main=3, cex.lab=cLab, cex.axis=cAxis, mar=Mar, omi=Omi) 
        iSub <- X >= 10 & X <= 40
    }
    plot(X[iSub], X[iSub], ylim=c(0,1.05), xlab="", ylab="", las=1, typ="n", xaxt="n")
    if(Stg==STAGE[1]) 
        mtext("Development", 2, cex=cLab, line=7.5)
    LWD <- 3
    Grey <- rgb(0,0,0,0.1)
    for (TT in temps) { lines(c(TT,TT),c(0,1.05), col=Grey) }
    lines(c(0,50),c(0,0),col=Grey) 
    ALPHA <- 10/nLines
    for (jj in 1:nLines) {
        lines(X[iSub], q01X[jj, iSub], col=rgb(0,0,1,ALPHA), lwd=LWD)
        lines(X[iSub], q99X[jj, iSub], col=rgb(0,0,1,ALPHA), lwd=LWD)
        lines(X[iSub], meanX[jj,iSub], col=rgb(1,0,0,ALPHA), lwd=LWD)
    }
    LWD <- 1
    ## Quantiles for q01
    lines(X[iSub], q01.CI95[1, iSub], lty=2, lwd=LWD)
    lines(X[iSub], q01.CI95[2, iSub], lty=1, lwd=LWD)
    lines(X[iSub], q01.CI95[3, iSub], lty=2, lwd=LWD)
    ## Quantiles for q99
    lines(X[iSub], q99.CI95[1,iSub], lty=2, lwd=LWD)
    lines(X[iSub], q99.CI95[2,iSub], lty=1, lwd=LWD)
    lines(X[iSub], q99.CI95[3,iSub], lty=2, lwd=LWD)
    ## Quantiles for mean
    lines(X[iSub], mean.CI95[1,iSub], lty=2, lwd=LWD)
    lines(X[iSub], mean.CI95[2,iSub], lty=1, lwd=LWD)
    lines(X[iSub], mean.CI95[3,iSub], lty=2, lwd=LWD)
    ## Legend
    legend("topleft", ## "topright"
           col=c(rgb(0,0,1),rgb(1,0,0), rgb(0,0,0), rgb(0,0,0)),
           legend=c(TeX("Sampled $P_{01}$ & $P_{99}$ "),
                    TeX("Mean | sampled $P_{01}$ & $P_{99}$ "),
                    "Mean of samples", 
                    "95% CI of samples"), bty="n",
           lty=c(1,1,1,2), bg="white", cex=2.4, lwd=4)
    ## Title
    if (Stg=="egg") { 
        mtext("Egg", line=1, cex=cMain)
    } else if (Stg=="gonotrophicCycle") {
        mtext("Adult", line=1, cex=cMain)
    } else if  (Stg=="larvae-pupae"){
        if (subStg=="LARVA")
            mtext("Larva", line=1, cex=cMain)
        if (subStg=="PUPA")
            mtext("Pupa", line=1, cex=cMain)
    }
    ## ##########################
    ## Plot Interpolated Survival 
    iSub <- X >= 10 & X <= 40
    plot(X[iSub], X[iSub], ylim=0:1, xlab="", ylab="", las=1, typ="n", xaxt="n")
    axis(1, at=seq(10,40,5), labels=seq(10,40,5), cex.axis=cAxis, padj=Padj)
    mtext("Temperature, °C", 1, cex=cLab, line=8)
    if(Stg==STAGE[1])
        mtext("Daily survival probability", 2, cex=cLab, line=7.5)
    LWD <- 3; Grey <- rgb(0,0,0,0.1)
    for (TT in temps) { lines(c(TT,TT),0:1, col=Grey) }
    lines(c(0,50),c(0,0),col=Grey) 
    ALPHA <- 10/nLines
    mycol <- rainbow(10, s=.9, alpha = ALPHA)
    for (jj in 1:nLines) {
        lines(X[iSub], surX[jj,iSub ], col=mycol, lwd=LWD)
    }
    LWD <- 1
    Quants <- matrix(0, ncol=3, nrow=length(X))
    ## Quantiles for q01, mean and q99 (95% CI)
    for (jj in 1:length(X))
        Quants[jj,] <- quantile(surX[,jj], probs=c(0.025, 0.5, 0.975))
    lines(X[iSub], Quants[iSub,1], lty=2, lwd=LWD)
    lines(X[iSub], Quants[iSub,2], lty=1, lwd=LWD)
    lines(X[iSub], Quants[iSub,3], lty=2, lwd=LWD)
    if(Stg==STAGE[length(STAGE)]) {
        if (PDF)
            dev.off()
    }
}    


## #######################################
## PLOT INTERPOLATED EXPECTED FECUNDITY ##
## #######################################

CASE  <- "IPLM" ## it is always the case here
setwd(baseDir)
Stg   <- "gonotrophicCycle"
STAGE <- Stg
STG   <- "GC"
source("setDirectories.R")
subsets <- tailSamples2subsetList(samples=samples, nLines=nLines, parNames=parNames)
EFecMat <- subsets[[1]]; head(EFecMat,2) ## EFec at empirical temperatures
P01Mat  <- subsets[[2]]; head(P01Mat,2) 
P99Mat  <- subsets[[3]]; head(P99Mat,2) 
muMat   <- subsets[[4]]; head(muMat,2) 
scMat   <- subsets[[5]]; head(scMat,2) 
surMat  <- subsets[[6]]; head(surMat,2) 
resVec  <- as.numeric(subsets[[7]]); head(resVec,10)
## Define res
res <- resVec[1]
## #################
## OBTAIN SPLINES ##
setwd(baseDir)
source("getSplines.R")
setwd(baseDir)
setwd(paste(getwd(), "gonotrophicCycle/figures", sep= "/"))
print(getwd())
## ######
## Plot #
if(PDF)
    pdf(paste0(STG,"_interp-ExpecFecundity_", CASE, "_nLines", nLines, ".pdf"))
par(mfrow=c(1,1), mar=c(5,6,1,1))
cAxis <- 2.3
cLab  <- 2.5
iSub  <- X >= 10 & X <= 35
plot(X[iSub], X[iSub], ylim=c(0, 130), typ="n", xlab="", 
     ylab="", cex.axis=cAxis, cex.lab=cLab, yaxt="n")
axis(2, at=c(0,60,120), las=2, cex.axis=cAxis)
mtext("Temperature", 1, line=3.3, cex=cLab)
mtext("E[Fecundity|Y]", 2, line=3.8, cex=cLab)
ALPHA <- 10/nLines
for (ii in 1:nLines)
    lines(X[iSub], fecX[ii,iSub], col=rgb(1,0,0, ALPHA)) 
Grey <- rgb(0,0,0,0.1)
for (TT in temps) { lines(c(TT,TT),c(0,130), col=Grey) }
Quants <- matrix(0, ncol=3, nrow=length(X))
## Quantiles for q01, mean and q99 (95% CI)
for (jj in 1:length(X)) 
    Quants[jj,] <- quantile(fecX[,jj], probs=c(0.025, 0.5, 0.975))
lines(X[iSub], Quants[iSub,1], lty=2)
lines(X[iSub], Quants[iSub,2], lty=1)
lines(X[iSub], Quants[iSub,3], lty=2)
if(PDF)
    dev.off()








