#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ####################################################
## Plots case study II - Stable population structure ##
## ####################################################
rm(list=ls())
options(width=165)
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
interpDir <- "mcmcOutput4plots"
figDir    <- "figures"
setwd(baseDir)
library(nimble)
library(latex2exp)
source("../FUNCTIONS_NIMBLE.R")
source("../FUNCTIONS_R.R")
source("FUNCTIONS_CULICOIDES.R")

## Basic parameters
nLinesPlot  <- 1000 ## nb lines taken for plots 
PROP_FEMALE <- 1/2  ## 50% of eggs are female
DIY         <- 365  ## Days In Year    
TRAJ_LENGTH <- DIY * 1.5
PLOT        <- "SteadyState" ## used within LOAD_interpData.R to load T correctly
CASE        <- "CLM" ## "IPLM" ## "MAP"
PDF         <- TRUE  ## FALSE
print(paste("Working with case", INI_POP, ",", CASE))

## ###########################
## LOAD INTERPOLATION Rdata ##
source("LOAD_interpData.R")

## Define
resE   <- resE[1]
resL   <- resL[1]
resP   <- resP[1]
resGC  <- resGC[1]
resTot <- resE + resL + resP + resGC

## #########################
## SET  FIGURES DIRECTORY ##
setwd(baseDir)
setwd(figDir)

## for plotting purposes
if(CASE=="IPLM")
    Case <- "MAP" 
if(CASE=="CLM")
    Case <- "CLM"


## #######
## PLOT ##
## #######
if (PDF)
    pdf(file = paste0("MAP_Stable-pop_", CASE, "_nLines", nLinesPlot, ".pdf"))
par(mar=c(5,5,1,1), mfrow=c(3,2), mai=c(0.4,0.4,0.3,0.2), omi=c(0.3,0.4,0,0.1))
for (tt in T) {
    (jj  <- which(T==tt))
    print(tt)
    EigenX      <- matrix(NA, nrow=nLinesPlot, lT)
    Proportions <- matrix(0, nrow=nLinesPlot, ncol=resTot)
    plot(0:resTot, 0*(0:resTot),
         main=paste(T[jj],"°C"), cex.main=2, 
         ylim=c(0,1), 
         typ="n", cex.axis=2, cex.lab=1.8, xlab="", ylab="")
    title(xlab = paste("Developmental sub-stage", Case),
          ylab = "",
          outer = TRUE, line = 0, cex.lab= 2.5)
    if (is.element (tt, c(10,20,30)))
        mtext("Proportion", 2, line= 2.6, cex=1.7)
    lines(rep(0, 2), c(0,2), col="grey")
    lines(rep(resE-1, 2), c(0,2), col="grey")
    lines(rep(resE+resL-2, 2), c(0,2), col="grey")
    lines(rep(resE+resL+resP-3, 2), c(0,2), col="grey")
    lines(rep(resE+resL+resP+resGC-4,2), c(0,2), col="grey")
    abline(0, 0, col="grey") ## Longterm Growth Threshold
    lines(c(0,0), c(-1,1)*pi, col="grey")        
    for (ii in 1:nLinesPlot) {
        jj    <- which(T==tt)
        Paras <- rbind(c(muE[ii,jj],  scE[ii,jj],  surE[ii,jj]),
                       c(muL[ii,jj],  scL[ii,jj],  surL[ii,jj]),
                       c(muP[ii,jj],  scP[ii,jj],  surP[ii,jj]),
                       c(muGC[ii,jj], scGC[ii,jj], surGC[ii,jj]))
        Femfec <- fec[ii,jj] * PROP_FEMALE
        Res <- c(resE, resL, resP, resGC)
        M <- nf_IPLM(paras=Paras, res=Res, femfec=Femfec, gCycle=1)
        ## Eigen values and vectors
        Eigen      <- eigen(M)
        (DomEigVal <- Mod(Eigen$values[1]))
        (w         <- Mod(Eigen$vectors[,1]))
        (v         <- Mod(eigen(t(M))$vectors[,1]))
        domEVecij  <- Eigen$vectors[,1]
        domEVecij  <- domEVecij / sum(domEVecij)
        Proportions[ii,] <- domEVecij
        for (ll in 1:resTot)
            lines((ll-1):ll, rep(domEVecij[ll],2), col=rgb(1,0,0,alpha=0.05))
    }
    ## Quantiles
    Quants <- matrix(0, ncol=3, nrow=resTot)
    for (ll in 1:resTot) {
        Quants[ll,] <- quantile(Proportions[,ll], c(0.025, 0.5, 0.975))
        lines((ll-1):ll, rep(Quants[ll,1],2))
        lines((ll-1):ll, rep(Quants[ll,2],2), lwd=1.5)
        lines((ll-1):ll, rep(Quants[ll,3],2))
    }
}
if (PDF)
    dev.off()


