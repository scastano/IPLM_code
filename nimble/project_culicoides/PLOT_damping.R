#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ###################################################################
## Plots case study II - Damping ratio at several fixed temperatures #
## ###################################################################
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
TRAJ_LENGTH <- 250 
PDF  <- TRUE ## FALSE
PLOT <- "Damping" ## used in LOAD_interpData.R to load T correctly
CASE <- "CLM" ## IPLM 
print(paste("Working with case", CASE))

## ###########################
## LOAD INTERPOLATION Rdata ##
source("LOAD_interpData.R")

## #########################
## SET  FIGURES DIRECTORY ##
setwd(baseDir)
setwd(figDir)


##########
## PLOT ##
if (PDF)
    pdf(file = paste0("damping_", CASE, "_nlinesPlot", nLinesPlot, ".pdf"), width = 11, height = 2.4)
(subT  <- seq(15, 35, by=5))
MAR <- c(2.4, 0, 1.2, 1.4)
CEX.AXIS <- 1.5
CEX.MAIN <- 1.6
YMAX <- 1
Padj <- 0.15
par(mar=MAR, mfrow=c(1,5))
if(CASE=="CLM") 
    par(mai=c(0.3, 0.5, 0.3, 0),  omi=c(0, 0.3, 0, 0.25)) 
if(CASE=="IPLM" | CASE=="MAP")
    par(mai=c(0.57, 0.5, 0.0, 0), omi=c(0, 0.3, 0.05, 0.25)) 
for (tt in subT) { 
    print(tt)
    if (tt>15) {
        TRAJ_LENGTH <- 200
    } else {
        TRAJ_LENGTH <- 800}
    if (tt==15) {
        if(CASE=="CLM") {
            plot(0:TRAJ_LENGTH, 0:TRAJ_LENGTH, las=1,
                 ylim=c(0, YMAX),
                 main=paste(tt, "°C", sep=""),
                 typ="n", xlab="", ylab="", xaxt="n",
                 cex.axis=CEX.AXIS, cex.main=CEX.MAIN)
        }
        if(CASE=="IPLM" | CASE=="MAP") {
            plot(0:TRAJ_LENGTH, 0:TRAJ_LENGTH, las=1,
                 ylim=c(0, YMAX),
                 typ="n", xlab="Days", ylab="", 
                 cex.axis=CEX.AXIS, cex.lab=2,
                 xaxt="n")
        }
    } else {
        if(CASE=="CLM") {
            plot(0:TRAJ_LENGTH, 0:TRAJ_LENGTH, las=1,
                 ylim=c(0, YMAX),
                 main=paste(tt, "°C", sep=""),
                 typ="n", xlab="", ylab="", xaxt="n",
                 cex.axis=CEX.AXIS, cex.main=CEX.MAIN)
        }
        if(CASE=="IPLM" | CASE=="MAP") {
            plot(0:TRAJ_LENGTH, 0:TRAJ_LENGTH, las=1,
                 ylim=c(0, YMAX),
                 typ="n", xlab="Days", ylab="",
                 cex.axis=CEX.AXIS, cex.main=CEX.MAIN, cex.lab=2,
                 xaxt="n")
        }
    }
    ##
    if(tt==15) {
        AT <- seq(0, TRAJ_LENGTH, by= 400)
        axis(1, at=AT, lab=AT, las=1, cex.axis=CEX.AXIS, padj=Padj)
    } else {
        AT <- seq(0, TRAJ_LENGTH, by= 100)
        axis(1, at=AT, lab=AT, las=1, cex.axis=CEX.AXIS, padj=Padj)
    }
    ## X and Y labels
    Cex.OuterL <- 1.8
    Cex.OuterR <- 1.45
    if (CASE=="CLM") {
            mtext(TeX("$\\rho^{-t}$"), outer=TRUE, side=2, line=-0.7, cex=Cex.OuterL,las=2)
        if(tt==subT[length(subT)])
            mtext(CASE, outer=TRUE, side=4, line=0.75, cex=Cex.OuterR)
    }
    if (CASE=="IPLM" | CASE=="MAP") {
        mtext(TeX("$\\rho^{-t}$"), outer=TRUE, side=2, line=-0.7, cex=Cex.OuterL,las=2, padj=-0.2)
        if(tt==subT[length(subT)])
            mtext(CASE, outer=TRUE, side=4, line=0.75, cex=Cex.OuterR, adj=0.65)
    }
    ## 
    abline(1,0,col="grey")
    invDamping <- matrix(0, nrow=nLinesPlot, ncol=TRAJ_LENGTH+1)
    for (ii in 1:nLinesPlot) {
        Res   <- c(resE[ii],resL[ii], resP[ii],resGC[ii])
        lres  <- sum(Res)
        (jj   <- which(T==tt)) 
        Paras <- rbind(c(muE[ii,jj],  scE[ii,jj],  surE[ii,jj]),
                       c(muL[ii,jj],  scL[ii,jj],  surL[ii,jj]),
                       c(muP[ii,jj],  scP[ii,jj],  surP[ii,jj]),
                       c(muGC[ii,jj], scGC[ii,jj], surGC[ii,jj]))
        Femfec <- fec[ii,jj] * PROP_FEMALE
        M <- nf_IPLM(paras=Paras, res=Res, femfec=Femfec, gCycle=1)
        ## Eigen values and vectors
        Eigen           <- eigen(M)
        (DomEigVal      <- Mod(Eigen$values[1]))
        (modEigVal2     <- Mod(Eigen$values[2]))
        (dampingRatio   <- DomEigVal / modEigVal2)
        invDamping[ii,] <- dampingRatio ^ -seq(0, TRAJ_LENGTH, by=1)            
        lines(0:TRAJ_LENGTH, invDamping[ii,], col=rgb(1,0,0,0.05)) ## Adults - black
    }
    ## Quantiles
    quants <- matrix(0, nrow=TRAJ_LENGTH+1, ncol=3)
    for (jj in 1:(TRAJ_LENGTH+1)) {
        quants[jj,] <- quantile(invDamping[,jj], c(0.025, 0.5, 0.975))
    }
    lines(0:TRAJ_LENGTH, quants[,1], lty=2)
    lines(0:TRAJ_LENGTH, quants[,2])
    lines(0:TRAJ_LENGTH, quants[,3], lty=2)
}
if (PDF)
    dev.off()

