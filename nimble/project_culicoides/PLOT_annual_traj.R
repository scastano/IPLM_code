#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## #########################################################################
## Plots case study II - Annual trajectories at several fixed temperatures #
## #########################################################################
options(width=165)
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
interpDir <- "mcmcOutput4plots"
figDir    <- "figures"

setwd(baseDir)
library(nimble)
library(latex2exp)
library(Matrix)

source("../FUNCTIONS_NIMBLE.R")
source("../FUNCTIONS_R.R")
source("FUNCTIONS_CULICOIDES.R")

## Basic parameters
nLinesPlot  <- 1000 ## nb lines taken for plots 
PROP_FEMALE <- 1/2  ## Currently assume 50% of eggs are female
TRAJ_LENGTH <- 150
PDF  <- TRUE ## FALSE
PLOT <- "AnnualTraj"
CASE <- "IPLM" ## "CLM
print(paste("Working with case", CASE))

## ###########################
## LOAD INTERPOLATION Rdata ##
source("LOAD_interpData.R")

## ########################
## SET FIGURES DIRECTORY ##
setwd(baseDir)
setwd(figDir)


W0_10 <- vector("list", nLinesPlot)
## obtain W0 at 10C
for (ii in 1:nLinesPlot) {  ## ii = 1
    Res   <- c(resE[ii],resL[ii], resP[ii],resGC[ii])
    lres  <- sum(Res)
    (jj   <- which(T==10)) 
    Paras <- rbind(c(muE[ii,jj],  scE[ii,jj], surE[ii,jj]),
                   c(muL[ii,jj],  scL[ii,jj], surL[ii,jj]),
                   c(muP[ii,jj],  scP[ii,jj], surP[ii,jj]),
                   c(muGC[ii,jj], scGC[ii,jj], surGC[ii,jj]))
    Femfec <- fec[ii,jj] * PROP_FEMALE
    M <- nf_IPLM(paras=Paras, res=Res, femfec=Femfec, gCycle=1)
    ## Eigen values and vectors
    Eigen      <- eigen(M)
    (DomEigVal <- Mod(Eigen$values[1])) ## Should be real, but small numerical errors occur when stable_adult_proportion ~= 0
    (w         <- Mod(Eigen$vectors[,1]))
    (v         <- Mod(eigen(t(M))$vectors[,1]))
    (sumW  <- sum(w))
    (sumVW <- as.numeric(v %*% w))
    W0_10[[ii]] <- w / sum(w)
}    


#########################################
## Relative Population Density - Total ##
if (PDF) 
    pdf(file = paste0("Transient_Traj_", CASE, "_", "_nlinesPlot", nLinesPlot, ".pdf"),
        width = 11, height = 2.4)
MAR <- c(2.4, 0, 1.2, 1.4)
CEX.AXIS <- 1.5 
CEX.MAIN <- 1.6 
YMAX <- 1
par(mar=MAR, mfrow=c(1,5))
if(CASE=="CLM") 
    par(mai=c(0.3, 0.4, 0.3, 0),  omi=c(0, 0.2, 0, 0.25))
if(CASE=="IPLM")
    par(mai=c(0.35, 0.4, 0.0, 0), omi=c(0, 0.2, 0.13, 0.25))
##
subT <- T[-1]
for (tt in subT) { 
    print(tt)
    if(tt==10) {}
    RelDensity        <- matrix(1, nrow=nLinesPlot, ncol=TRAJ_LENGTH+1)
    Amplification     <- matrix(1, nrow=nLinesPlot, ncol=TRAJ_LENGTH+1)
    Attenuation       <- matrix(1, nrow=nLinesPlot, ncol=TRAJ_LENGTH+1)
    AmplifiedInertia  <- rep(1, nLinesPlot)
    AttenuatedInertia <- rep(1, nLinesPlot)
    for (ii in 1:nLinesPlot) { 
        if (ii %% floor(nLinesPlot / 20) == 0)
            print(ii)
        Res  <- c(resE[ii],resL[ii], resP[ii],resGC[ii])
        lres <- sum(Res)
        (jj   <- which(subT==tt) + 1) 
        Paras <- rbind(c(muE[ii,jj],  scE[ii,jj], surE[ii,jj]),
                       c(muL[ii,jj],  scL[ii,jj], surL[ii,jj]),
                       c(muP[ii,jj],  scP[ii,jj], surP[ii,jj]),
                       c(muGC[ii,jj], scGC[ii,jj], surGC[ii,jj]))
        Femfec <- fec[ii,jj] * PROP_FEMALE
        M <- nf_IPLM(paras=Paras, res=Res, femfec=Femfec, gCycle=1)
        ## Define population density vector N
        N <- matrix(0, ncol=lres, nrow=TRAJ_LENGTH+1)
        ## Eigen values and vectors
        Eigen      <- eigen(M)
        (DomEigVal <- Mod(Eigen$values[1]))
        (w     <- Mod(Eigen$vectors[,1]))
        (v     <- Mod(eigen(t(M))$vectors[,1]))
        (sumW  <- sum(w))
        (sumVW <- as.numeric(v %*% w))
        ## Inertia
        (AmplifiedInertia[ii]  <- max(v) * sumW / sumVW)
        (AttenuatedInertia[ii] <- min(v) * sumW / sumVW)
        ## Standardise Population Matrix
        stM   <- M / DomEigVal
        stMtt <- diag(1, lres, lres)
        ## Set initial conditions
        N0    <- W0_10[[ii]]
        N[1,] <- N0
        ## Loop on time steps
        for (ll in 1:TRAJ_LENGTH) { 
            ## Generate standardised matrix for step ll
            (stMtt <- stMtt %*% stM)
            ## Project population 1 time step
            N[ll+1,] <- stMtt %*% N0
            ## Obtain bounds on transients
            CS <- colSums(stMtt)
            Amplification[ii,ll+1] <- max(CS)
            Attenuation[ii,ll+1]   <- min(CS)                
        }
        RelDensity[ii,] <- rowSums(N)
    } 
    ## #######
    ## Plot ##
    YLIM  <- c(-2,2)
    BY <- 0.1
    (AT <- sort(seq(-2, 2, by=1)))
    if (CASE=="CLM") {
        (AT <- seq(floor(YLIM[1]), ceiling(YLIM[2]), by=BY))
        while (length(AT) > 10) {        
            (AT <- sort(seq(ceiling(YLIM[2]), floor(YLIM[1]), by=-BY)))
            (BY <- BY * 2)
        }
        (AT <- sort(pretty(c(YLIM)))) 
    }
    (LAB  <- rep("", length(AT)))
    (iLAB <- c(1,length(LAB)))
    LAB[ iLAB ] <- sapply(AT[ iLAB ], function(i) as.expression(bquote(10^ .(i))))
    LAB[3] <- "0"
    (YLIM <- range(AT))
    length_plot <- TRAJ_LENGTH 
    if(CASE=="CLM") {
        plot(0:length_plot, 0:length_plot, las=1,
             ylim=YLIM, yaxt="n", xaxt="n",
             main=paste(tt, "°C", sep=""),
             typ="n", xlab="", ylab="", cex.axis=CEX.AXIS, cex.main=CEX.MAIN)
        ## Ylab
        title(xlab = "", ylab = TeX("Density$/ \\lambda_1^t$"),
              outer = TRUE, line = -1, cex.lab= 1.9)
    } else {
        plot(0:length_plot, 0:length_plot, las=1,
             ylim=YLIM, yaxt="n", xaxt="n",
             typ="n", xlab="", ylab="", cex.axis=CEX.AXIS, cex.main=CEX.MAIN)
        ## Ylab
        title(xlab = "", ylab = TeX("Density$/ \\lambda_1^t$"),
              outer = TRUE, line = -1, cex.lab= 1.9, adj=0.7)
    }
    ## ADD 'CASE' in right side
    if(tt==T[lT]) {
        if (CASE=="CLM")
            mtext(CASE, outer=TRUE, side=4, line=0.5, cex.lab=1.3)
        if (CASE=="IPLM")
            mtext(CASE, outer=TRUE, side=4, line=0.5, cex.lab=1.3, adj=0.65)
    }
    ## ADD X-AXIS TICKS
    (atX <- pretty(c(0, TRAJ_LENGTH), n=3))
    axis(1, at=atX, labels=rep("",length(atX)), lwd=0.8, adj=0.95)
    ## ADD X-AXIS LABELS
    Is.Odd <- !!(length(atX)%%2)
    atX <- range(atX)
    if (Is.Odd) { ## If length is odd add the midpoint
        (atX <- unique(sort(c(atX, median(atX)))))
    } 
    axis(1, at=atX, labels=atX, las=3, cex.axis=CEX.AXIS, las=1, padj=-0.15)
    ## ADD Y-AXIS TICKS & LABELS
    axis(2, at=AT, labels=LAB, cex.axis=CEX.AXIS, las=3, hadj=0.3, padj=0.3)   
    ## 
    if (CASE=="IPLM") mtext("Days", 1, line=1.5, cex=1.2, cex.lab=1.8)
    for (ii in 1:nLinesPlot) { 
        try(abline(h=log10(AmplifiedInertia[ii]),  col=rgb(0,1,0,0.01)))     ## Maximum Inertia - green
        try(abline(h=log10(AttenuatedInertia[ii]), col=rgb(0,1,0,0.01)))     ## Minimum Inertia - green
        lines(0:TRAJ_LENGTH, log10(Attenuation[ii,]), col=rgb(1,0,0,0.05))   ## Lower Bound on Attenuation   - red
        lines(0:TRAJ_LENGTH, log10(Amplification[ii,]), col=rgb(1,0,0,0.05)) ## Upper Bound on Amplification - red
        lines(0:TRAJ_LENGTH, log10(RelDensity[ii,]), col=rgb(0,0,1,0.05))    ## Trajectory given N0          - blue
    }
    abline(log10(1), 0 , col=rgb(0,0,0), lty=9, lwd=1) ## Stable state solution
}
if(PDF)
    dev.off()



