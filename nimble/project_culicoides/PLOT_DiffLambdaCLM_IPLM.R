#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ####################################################################################
## Plots case study II - Difference in interpolated growth rate for CLM & IPLM models #
## ####################################################################################
rm(list=ls())
options(width=165)
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
interpDir <- "mcmcOutput4plots"
figDir    <- "figures"
setwd(baseDir)
library(nimble)
library(latex2exp)
library(colorRamps)
source("../FUNCTIONS_NIMBLE.R")
source("../FUNCTIONS_R.R")
source("FUNCTIONS_CULICOIDES.R")

## Basic parameters
nLinesPlot  <- 1000 ## nb lines taken for plots 
PROP_FEMALE <- 1/2  ## 50% of eggs are female
DIY         <- 365  ## Days In Year    
TRAJ_LENGTH <- 250 
empiricT    <- c(15,17,20,23,27,30,35) ## experimental temperatures
PDF         <- TRUE ## FALSE
PLOT        <- "PopGrowthRate" ## used within LOAD_interpData.R to load T correctly
print(paste("Taking differences between CLM-lambda  & IPLM-lambda, nlines=", nLinesPlot))

## #############
## Firs load CLM, rename some objects & then load for IPLM
CASE    <- "CLM" 
source("LOAD_interpData.R")
## CLM parameters
muE.clm   <- muE
scE.clm   <- scE
surE.clm  <- surE
resE.clm  <- resE
##
muL.clm   <- muL
scL.clm   <- scL
surL.clm  <- surL
resL.clm  <- resL
##
muP.clm   <- muP
scP.clm   <- scP
surP.clm  <- surP
resP.clm  <- resP
##
muGC.clm  <- muGC
scGC.clm  <- scGC
surGC.clm <- surGC
fec.clm   <- fec
resGC.clm <- resGC

setwd(baseDir)
CASE    <- "IPLM"
source("LOAD_interpData.R")
## IPLM parameters
muE.tdlm   <- muE
scE.tdlm   <- scE
surE.tdlm  <- surE
resE.tdlm  <- resE
##
muL.tdlm   <- muL
scL.tdlm   <- scL
surL.tdlm  <- surL
resL.tdlm  <- resL
##
muP.tdlm   <- muP
scP.tdlm   <- scP
surP.tdlm  <- surP
resP.tdlm  <- resP
##
muGC.tdlm  <- muGC
scGC.tdlm  <- scGC
surGC.tdlm <- surGC
fec.tdlm   <- fec
resGC.tdlm <- resGC

## #########################
## SET  FIGURES DIRECTORY ##
setwd(baseDir)
setwd(figDir)

## ########
## PLOT  ##
EigenT.clm     <- matrix(NA, nrow=nLinesPlot, lT)
ModDomEig.clm  <- matrix(0, nrow=nLinesPlot, ncol=lT)
EigenT.tdlm    <- matrix(NA, nrow=nLinesPlot, lT)
ModDomEig.tdlm <- matrix(0, nrow=nLinesPlot, ncol=lT)
print(nLinesPlot)
if (PDF)
    pdf(file = paste("Diff_absolute_CLM-IPLM_", "nLinesPlot", nLinesPlot,".pdf",sep=""))
par(mar=c(5,5,3,1))
yLim <- c(-0.8,0.8)
plot(T,T, ylim=yLim, main=TeX("$\\lambda_1^{CLM} - \\lambda_1^{IPLM}$"), cex.main=2.2, lwd=8,
     typ="n", xlab="", ylab="", cex.axis=1.8, cex.lab=2)
Grey <- rgb(0,0,0,0.1)
for (n in empiricT) { lines(c(n,n), yLim, col=Grey) }
mtext("Temperature, °C", 1, line=2.8, cex=2)
mtext("Difference", 2, line=2.8, cex=2)
abline(0, 0, lty=2) ## Long-term Growth Threshold
Alpha <- 0.04
for (ii in 1:nLinesPlot) { 
    for (tt in T) {
        Res.clm   <- c(resE.clm[ii],resL.clm[ii], resP.clm[ii],resGC.clm[ii])
        lres.clm  <- sum(Res.clm)
        Res.tdlm  <- c(resE.tdlm[ii],resL.tdlm[ii], resP.tdlm[ii],resGC.tdlm[ii])
        lres.tdlm <- sum(Res.tdlm)
        (jj  <- which(T==tt)) 
        Paras.clm <- rbind(c(muE.clm[ii,jj],  scE.clm[ii,jj],  surE.clm[ii,jj]),
                           c(muL.clm[ii,jj],  scL.clm[ii,jj],  surL.clm[ii,jj]),
                           c(muP.clm[ii,jj],  scP.clm[ii,jj],  surP.clm[ii,jj]),
                           c(muGC.clm[ii,jj], scGC.clm[ii,jj], surGC.clm[ii,jj]))
        Femfec.clm <- fec.clm[ii,jj] * PROP_FEMALE
        M.clm      <- nf_IPLM(paras=Paras.clm, res=Res.clm, femfec=Femfec.clm, gCycle=1)
        ## Eigen values and vectors
        EigenT.clm[ii,jj]  <- eigen(M.clm)$values[1]            
        ModDomEig.clm[ii,] <- Mod(EigenT.clm[ii,])
        Paras.tdlm <- rbind(c(muE.tdlm[ii,jj],  scE.tdlm[ii,jj],  surE.tdlm[ii,jj]),
                            c(muL.tdlm[ii,jj],  scL.tdlm[ii,jj],  surL.tdlm[ii,jj]),
                            c(muP.tdlm[ii,jj],  scP.tdlm[ii,jj],  surP.tdlm[ii,jj]),
                            c(muGC.tdlm[ii,jj], scGC.tdlm[ii,jj], surGC.tdlm[ii,jj]))
        Femfec.tdlm <- fec.tdlm[ii,jj] * PROP_FEMALE
        M.tdlm      <- nf_IPLM(paras=Paras.tdlm, res=Res.tdlm, femfec=Femfec.tdlm, gCycle=1)
        ## Eigen values and vectors
        EigenT.tdlm[ii,jj]  <- eigen(M.tdlm)$values[1]            
        ModDomEig.tdlm[ii,] <- Mod(EigenT.tdlm[ii,])
    }
    ## Plot differences
    lines(T, ModDomEig.clm[ii,] - ModDomEig.tdlm[ii,], col=rgb(1,0,0,Alpha))
}
Diff <- ModDomEig.clm - ModDomEig.tdlm
## Quantiles
quants <- matrix(0, nrow=lT, ncol=3)
for (jj in 1:lT)
    quants[jj,] <- quantile(Diff[,jj], c(0.025, 0.5, 0.975))
lines(T, quants[,1], lty=2)
lines(T, quants[,2])
lines(T, quants[,3], lty=2)
if (PDF)
    dev.off()
 
