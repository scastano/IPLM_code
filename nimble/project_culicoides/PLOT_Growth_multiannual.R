#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## #######################################################################################################
## Plots case study II - Growht rate & Multiannual trajectories under cosine annual temperature profile ##
## #######################################################################################################
rm(list=ls())
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

## SET UP BEFORE RUNNING THE SCRIPT
CASE        <- "IPLM"  ## "CLM"
TempAmpl    <- "15-30" ## "15-25" 
nLinesPlot  <- 500
PDF         <- TRUE ## FALSE
PROP_FEMALE <- 1/2  ## 50% of eggs are female
DIY         <- 365  ## Days In Year    
TRAJ_LENGTH <- DIY * 2 
print(paste("Plotting", nLinesPlot, "lines, temperature amplitude" ,TempAmpl, "for case", CASE))

## ###############################################
## LOAD & FILTER nLinesPlot INTERPOLATION Rdata ##

source("LOAD_CosineInterpData.R") ## this can be slow
## First day of year with T > some thresholds
(plotDays <- c(which(T >= 10)[1], ## first day of the year at which temperature is >=10 degrees
               which(T >= 15)[1],
               which(T >= 30)[1]) [1])
(lengthPlotDays <- length(plotDays))

## #################################
## Set directory to store figures ##
setwd(baseDir)
setwd(figDir)

## ################################################################
## Relative Population PopDensity - Total - MULTI-ANNUAL CYCLES ##
Mlist <- vector("list", DIY)
Ylist <- vector("list", lengthPlotDays)
PopDensity    <- matrix(1, nrow=nLinesPlot, ncol=TRAJ_LENGTH) ## Relative population density - N0 = Multi-Annual steady state at 1st July
AmplifMax     <- matrix(1, nrow=nLinesPlot, ncol=TRAJ_LENGTH) ## Maximum amplification at each time step
AttenMin      <- matrix(1, nrow=nLinesPlot, ncol=TRAJ_LENGTH) ## Minimum attenuation at each time step
SeasonLength  <- numeric(length=nLinesPlot)
YearlyGrowth  <- numeric(length=nLinesPlot)
EigenM1_lines <- matrix(0, nrow=nLinesPlot, ncol=lT)
EigenT        <- matrix(NA, nrow=nLinesPlot, lT)
lambda1       <- matrix(0, nrow=nLinesPlot, ncol=lT)
## 
for (ii in 1:nLinesPlot) { 
    if (ii %% floor(nLinesPlot / 10) == 0)
        print(ii)
    ## DAILY PROJECTION MATRIX FOR EACH DAY OF YEAR
    EigenM1 <- numeric(length=lT)
    for (jj in 1:lT) {  ## jj=93 
        Res    <- c(resE[ii],resL[ii], resP[ii],resGC[ii])
        lres   <- sum(Res) ##   (resTot <- resE + resL + resP + resGC)
        (Tjj   <- T[jj])                        
        Paras  <- rbind(c(muE[ii,jj],  scE[ii,jj], surE[ii,jj]),
                        c(muL[ii,jj],  scL[ii,jj], surL[ii,jj]),
                        c(muP[ii,jj],  scP[ii,jj], surP[ii,jj]),
                        c(muGC[ii,jj], scGC[ii,jj], surGC[ii,jj]))
        Femfec <- fec[ii,jj] * PROP_FEMALE
        M      <- nf_IPLM(paras=Paras, res=Res, femfec=Femfec, gCycle=1)
        ## Stock in list : DAILY PROJECTION MATRIX 
        Mlist[[jj]] <- M
        ## Eigen values and vectors : DAILY PROJECTION MATRIX
        EigenM <- eigen(M)
        ## for plotting lambda 1 below:
        EigenT[ii,jj] <- EigenM$values[1]
        ## 
        (lambdaMod                        <- Mod(EigenM$values)) ## First should be real, but small numerical errors occur when stable_adult_proportion ~= 0 
        (lambdaArg                        <- Arg(EigenM$values)) ## First should be real, but small numerical errors occur when stable_adult_proportion ~= 0            
        (w                                <- Mod(EigenM$vectors[,1]))
        (v                                <- Mod(eigen(t(M))$vectors[,1]))
        attributes(Mlist[[jj]])$lambdaMod <- lambdaMod
        attributes(Mlist[[jj]])$lambdaArg <- lambdaArg
        attributes(Mlist[[jj]])$w         <- w
        attributes(Mlist[[jj]])$v         <- v
        ##
        EigenM1[jj] <- lambdaMod[1]
    }
    ## for plotting lambda 1 below:
    (lambda1[ii,] <- Mod(EigenT[ii,]))
    ## ###################################################
    ## YEARLY PROJECTION MATRIX FOR CERTAIN DAY OF YEAR ##
    for (plotDay in plotDays) { ## jj = 1
        (jj <- which(plotDay == plotDays))
        Ylist[[jj]] <- Mlist[[jj]] 
        for (JJ in (jj:(jj+DIY-2)) %% DIY + 1) { ## Loop over subsequent 364 days
            Ylist[[jj]] <- Mlist[[JJ]] %*% Ylist[[jj]] 
        }
        ## Eigen values and vectors 
        EigenY     <- eigen(Ylist[[jj]])
        (lambdaMod <- Mod(EigenY$values))          ## First should be real, but small numerical errors occur when stable_adult_proportion ~= 0 
        (lambdaArg <- Arg(EigenY$values))          ## First should be real, but small numerical errors occur when stable_adult_proportion ~= 0 
        (w         <- Mod(EigenY$vectors[,1])) 
        (v         <- Mod(eigen(t(Ylist[[jj]]))$vectors[,1])) 
        attributes(Ylist[[jj]])$lambdaMod <- lambdaMod 
        attributes(Ylist[[jj]])$lambdaArg <- lambdaArg 
        attributes(Ylist[[jj]])$w         <- w 
        attributes(Ylist[[jj]])$v         <- v
        YearlyGrowth[ii]                  <- lambdaMod[1]
    }
    (Elambda.Annual <- attributes(Ylist[[1]])$lambdaMod[1]) ## dominant Eigen value for 1 year       
    (Elambda.Daily  <- Elambda.Annual ^ (1/DIY)) ## lambda of an average day
    ## Initialise population vectors
    if (INI_POP=="ENDEMIC") {
        ## take stable pop of temperature corresponding to Ylis[[1]] 
        (N0 <- Mod(eigen(Ylist[[1]])$vectors[,1]))
        (N0 <- N0 / sum(N0))                                      
    }
    N     <- matrix(0, nrow=TRAJ_LENGTH, ncol=lres)
    N[1,] <- N0 ## first row: initial condition
    ## Initialise projection matrix
    stMtt <- diag(1, lres, lres)
    ## Loop on time steps
    for (tt in 2:TRAJ_LENGTH) { 
        ## Obtain day of year
        (jj <- (tt-1) %% DIY + 1)
        ## Generate standardised matrix for step tt
        (stM   <- Mlist[[jj]] / Elambda.Daily)
        (stMtt <- stM %*% stMtt )
        ## Project population 1 time step
        N[tt,] <- stMtt %*% N0
        ## Obtain bounds on transients
        CS <- colSums(stMtt)
        AmplifMax[ii,tt] <- max(CS)
        AttenMin[ii,tt]  <- min(CS)                
    }
    PopDensity[ii,]    <- rowSums(N)
    EigenM1_lines[ii,] <- EigenM1
}


## #################
## INITALISE PLOT ##
if (PDF)
    pdf(file = paste0("FIG_6_", TempAmpl, "_", CASE, "_",INI_POP, "_nLines",
                      nLinesPlot, ".pdf"), width = 5, height = 5.5)
par(mar=c(2,1,0,1), omi=c(0,0, 0, 0)) 
layout(matrix(c(1,2,3), 3, 1, byrow=TRUE), widths=c(2,2,2), heights=c(3,2,3/2))
ALPHA <- 0.025 
## ##############################
## PLOT BI-ANNUAL TRAJECTORIES ##
trajVec   <- c(0,1:TRAJ_LENGTH)
(TRAJ_L   <- trajVec[ trajVec %% 90 == 0])
TRAJ_L[1] <- 1
months   <- c("0","3","6","9","12","15","18","21","24")
if (TempAmpl=="15-25") {
    YMIN  <- -2.1
    YMAX  <- 2.6
}
if (TempAmpl=="15-30") {
    YMIN  <- -3.1
    YMAX  <- 3.8
}
YLIM <- c(YMIN, YMAX)
BY   <- 0.1
AT   <- seq(floor(YLIM[1]), ceiling(YLIM[2]), by=BY)
while (length(AT) > 10) {        
    (AT <- seq(floor(YLIM[1]), ceiling(YLIM[2]), by=BY))
    (BY <- BY * 2)
}
(AT  <- pretty(sort(unique(c(YMIN, YMAX, AT) ))))
(LAB <- rep("", length(AT)))
LAB[AT==round(AT)] <- sapply(AT[AT==round(AT)], function(i) as.expression(bquote(10^ .(i))))
CET.LAB  <- 1.7
CET.AXIS <- 1.6
par(mai=c(0.05, 0.7, 0.4, 0.1))
plot(1:TRAJ_LENGTH, 1:TRAJ_LENGTH, typ="n", ylim=YLIM , xaxt="n", yaxt="n", ylab="", xlab="",
     cex.axis=CET.AXIS, cex.lab=CET.LAB)## , cex.main=2)    
abline(0,0, col="grey")
title(CASE, line=0.3, cex.main=2.1)
mtext(TeX("Pop. Density$/ \\lambda_Y^t$"), 2, line=2.3, cex=1.5)
axis(2, at=AT, lab=LAB, cex.axis=CET.AXIS)
## Loop on lines
for (ii in 1:nLinesPlot) {         
    lines(1:TRAJ_LENGTH, log10(PopDensity[ii,]), col=rgb(0,0,1,ALPHA)) ## blue - trajectory
    lines(1:TRAJ_LENGTH, log10(AmplifMax[ii,]), col=rgb(1,0,0,ALPHA))  ## red  - Max Amplification
    lines(1:TRAJ_LENGTH, log10(AttenMin[ii,]), col=rgb(1,0,0,ALPHA))   ## red  - Min Attenuation
}
## ###########
## Lambda 1 ## 
par(mai=c(0.05, 0.7, 0, 0.1))
if (TempAmpl=="15-25") 
    yLim <- c(0.99,1.18)
if (TempAmpl=="15-30") 
    yLim <- c(0.95,1.25)
plot(1:(2*lT),1:(2*lT), ylim=yLim, 
     typ="n", xlab="", xaxt="n", yaxt="n", ylab="", cex.axis=CET.AXIS, cex.lab=1.8)
abline(h=1, col="gray")
if (TempAmpl=="15-25") 
    axis(2, at=c(1,1.05, 1.1, 1.15),lab=c("1","1.05","1.1","1.15"),  cex.axis=CET.AXIS)
if (TempAmpl=="15-30")
    axis(2, at=c(1,1.1,1.2),lab=c("1","1.1","1.2"),  cex.axis=CET.AXIS)
mtext(TeX("$\\lambda_1$"), 2, line=2.3, cex=CET.LAB)
for (ii in 1:nLinesPlot) 
    lines(1:(2*lT),c(lambda1[ii,],lambda1[ii,]), col=rgb(1,0,0,ALPHA))
##
## ########################################
## THIRD ROW: cosine temperature profile ##
par(mai=c(0.5, 0.7, 0, 0.1)) 
plot(1:(2*lT), c(T,T), ylim=c(15,32) , xaxt="n", yaxt="n", xlab="", ylab="",
     cex.axis=CET.AXIS, cex.lab=2.1, cex=0.1)
axis(1, at=TRAJ_L, lab=months,  cex.axis=CET.AXIS)
mtext("Months", 1, line=2.6, cex=CET.LAB)
if (TempAmpl=="15-25")
    axis(2, at=c(15,25),lab=c("15","25"),  cex.axis=CET.AXIS)
if (TempAmpl=="15-30")
    axis(2, at=c(15,30),lab=c("15","30"),  cex.axis=CET.AXIS)
mtext("Temperature", 2, line=2.55, cex=1.23, at=22)
if (PDF)
    dev.off()







