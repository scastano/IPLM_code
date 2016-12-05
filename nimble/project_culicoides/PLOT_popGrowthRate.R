#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ############################################################################
## Plots case study II - Growth rate in a range of interpolated temperatures ## 
## ############################################################################
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
PROP_FEMALE <- 1/2  ##  50% of eggs are female
DIY         <- 365  ## Days In Year    
TRAJ_LENGTH <- 250 
PDF         <- TRUE ## FALSE
PLOT        <- "PopGrowthRate" ## used within LOAD_interpData.R to load T correctly
CASE        <- "IPLM" ## "CLM 
print(paste("Working with case", CASE))

## ###########################
## LOAD INTERPOLATION Rdata ##
source("LOAD_interpData.R")


## #########################
## SET  FIGURES DIRECTORY ##
setwd(baseDir)
setwd(figDir)

empiricT <- c(15,17,20,23,27,30,35)
## #######
## PLOT ##
EigenT    <- matrix(NA, nrow=nLinesPlot, lT)
ModDomEig <- matrix(0, nrow=nLinesPlot, ncol=lT)
print(nLinesPlot)
if (PDF)
    pdf(file = paste("Population_growth_rate_", CASE, "_", "nLinesPlot", nLinesPlot,".pdf",sep=""))
par(mar=c(5,5,3,1))
yLim <- c(0.7,1.3)
plot(T,T, ylim=yLim, main=CASE,cex.main=2,
     typ="n", xlab="", ylab="", cex.axis=1.8, cex.lab=2)
Grey <- rgb(0,0,0,0.1)
for (n in empiricT) { lines(c(n,n), yLim, col=Grey) }
mtext("Temperature, °C", 1, line=2.8, cex=2)
mtext(TeX("$\\lambda_1$"), 2, line=2.8, cex=2)
abline(1, 0, lty=2) ## Long-term Growth Threshold
for (ii in 1:nLinesPlot) { 
    for (tt in T) {
        Res   <- c(resE[ii],resL[ii], resP[ii],resGC[ii])
        lres  <- sum(Res)
        jj    <- which(T==tt)
        Paras <- rbind(c(muE[ii,jj],  scE[ii,jj],  surE[ii,jj]),
                       c(muL[ii,jj],  scL[ii,jj],  surL[ii,jj]),
                       c(muP[ii,jj],  scP[ii,jj],  surP[ii,jj]),
                       c(muGC[ii,jj], scGC[ii,jj], surGC[ii,jj]))
        Femfec <- fec[ii,jj] * PROP_FEMALE
        M      <- nf_IPLM(paras=Paras, res=Res, femfec=Femfec, gCycle=1)
        ## Eigen values and vectors
        EigenT[ii,jj]  <- eigen(M)$values[1]            
        ModDomEig[ii,] <- Mod(EigenT[ii,])
    }
    lines(T, ModDomEig[ii,], col=rgb(1,0,0,0.04))
}
## Quantiles
quants <- matrix(0, nrow=lT, ncol=3)
for (jj in 1:lT)
    quants[jj,] <- quantile(ModDomEig[,jj], c(0.025, 0.5, 0.975))
lines(T, quants[,1], lty=2)
lines(T, quants[,2])
lines(T, quants[,3], lty=2)
if (PDF)
    dev.off()
