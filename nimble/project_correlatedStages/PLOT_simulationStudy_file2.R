#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ############################################
## Plots associated to case study I - file 2 ##
## ############################################
rm(list=ls())

## Set directories and basic parameters
options(width=160)
library(latex2exp)
library(nimble)
library(plotrix)

baseDir <- "~/IPLM_paper/nimble/project_correlatedStages/"
MCMC    <- "mcmc"
setwd(baseDir)
source("../FUNCTIONS_R.R")
source("../FUNCTIONS_NIMBLE.R")
source("../project_culicoides/FUNCTIONS_CULICOIDES.R")

## Compile nf_TW1_MuSigSurv
thresh         <- 1E-6
cTW1_MuSigSurv <- compileNimble(nf_TW1_MuSigSurv)

PDF <- FALSE ## TRUE

## ############################
## Read mcmc output directories
setwd(baseDir)
mcmcDir <- paste(getwd(), MCMC, sep='/') 
setwd(mcmcDir)
## Find mcmc results
(mcmcDirs <- system("ls -1t|head -n 1000", TRUE)) ## > 500 simulations
nsims <- length(mcmcDirs)
## Define matrix to store true qualities by sim (row)
Tqual <- matrix(0, nsims, 50)
## Set a minimum of fully developed individuals
mindev <- 35
## Define matrices to store (by row) true, median, 95%CI (lower & upper bounds) 
qTotDieStg1 <- qTotDieStg2 <- qTotCenStg1 <- qTotCenStg2 <- qTotFullDev <- matrix(0,4,1)

imcmc <- 0
######################
## THIS CAN BE SLOW ##
for (mcmcDir in mcmcDirs) { 
    setwd(baseDir)
    setwd(MCMC)
    newdir <- paste(getwd(),mcmcDir,sep="/")
    setwd(newdir)
    ## Read data
    if(length(dir())>1) {
        load(system("ls", TRUE)[2] )
        if (fulldev >= mindev) { 
            imcmc <- imcmc + 1
            ## Just checking 
            if (FALSE) { ## TRUE  
                ## these 3 below must match
                fulldev
                length(qual_history$qual_FullDev)
                length(qualOrder_history$qualOrder_FullDev)
                ## these 2 below must match
                Inits[["qual"]][qualOrder_history$qualOrder_DieStg1]
                qual_history$qual_DieStg1
            }
            ## EXTRACT parameters & data
            N <- Constants[["N"]]
            Tqual[imcmc,] <- Inits[["qual"]]
            ## Burn-in Samples
            samples <- tail(samples, nrow(samples) - round(nrow(samples)/2))
            (nlines <- nrow(samples))
            colnames(samples)
            ## TAKE ESTIMATED QUALITIES
            Q.imcmc <- samples[,paste0("qual" , 1:N)]
            ## ##################################################
            ## TAKE DATA OF Q-INDIVIDUALS THAT DIED IN STAGE 1 ##
            ##
            ## Take estimates & true values of 'q' associated to individual/s that died in Stg 1
            if (length(qualOrder_history$qualOrder_DieStg1)==0) {
                qDieStg1  <- matrix(99,nrow(samples),1)
                TqdieStg1 <- 99
                storeqDieStg1 <- matrix(rep(99,4),ncol=1) ## by row: TqdieStg1, qdieStg1.M, qdieStg1.L, qdieStg1.U
                ##
            } else if (length(qualOrder_history$qualOrder_DieStg1)==1) { ## vector
                qDieStg1   <- Q.imcmc[,qualOrder_history$qualOrder_DieStg1]
                (TqdieStg1 <- qual_history$qual_DieStg1) ## true value
                ## Take median and 95%CI of estimated qualities
                (qdieStg1.M <- median(qDieStg1))
                (qdieStg1.L <- quantile(qDieStg1, 0.025))
                (qdieStg1.U <- quantile(qDieStg1, 0.975))
                storeqDieStg1 <- matrix(c(TqdieStg1,qdieStg1.M,qdieStg1.L,qdieStg1.U),ncol=1) ## by row: TqdieStg1, qdieStg1.M, qdieStg1.L, qdieStg1.U
            } else { ## matrix
                qDieStg1  <- Q.imcmc[,qualOrder_history$qualOrder_DieStg1]
                TqdieStg1 <- qual_history$qual_DieStg1
                ## Take median and 95%CI of estimated qualities
                qdieStg1.M <- apply(qDieStg1, 2, function(x) median(x))
                qdieStg1.L <- apply(qDieStg1, 2, function(x) quantile(x, 0.025))
                qdieStg1.U <- apply(qDieStg1, 2, function(x) quantile(x, 0.975))
                storeqDieStg1 <- rbind(TqdieStg1,qdieStg1.M,qdieStg1.L,qdieStg1.U)
            }
            ## store
            qTotDieStg1 <- cbind(qTotDieStg1, storeqDieStg1)
            ## ##################################################
            ## TAKE DATA OF Q-INDIVIDUALS THAT DIED IN STAGE 2 ##
            ##
            ## Take estimates & true values of 'q' associated to individual/s that died in Stg 2
            if (length(qualOrder_history$qualOrder_DieStg2)==0) {
                qDieStg2  <- matrix(99,nrow(samples),1)
                TqdieStg2 <- 99
                storeqDieStg2 <- matrix(rep(99,4),ncol=1) ## by row: TqdieStg1, qdieStg1.M, qdieStg1.L, qdieStg1.U
                ##
            } else if (length(qualOrder_history$qualOrder_DieStg2)==1) { ## vector
                qDieStg2   <- Q.imcmc[,qualOrder_history$qualOrder_DieStg2]
                (TqdieStg2 <- qual_history$qual_DieStg2) ## true value
                ## Take median and 95%CI of estimated qualities
                (qdieStg2.M <- median(qDieStg2))
                (qdieStg2.L <- quantile(qDieStg2, 0.025))
                (qdieStg2.U <- quantile(qDieStg2, 0.975))
                storeqDieStg2 <- matrix(c(TqdieStg2,qdieStg2.M,qdieStg2.L,qdieStg2.U),ncol=1) ## by row: TqdieStg2, qdieStg2.M, qdieStg2.L, qdieStg2.U
            } else { ## matrix
                qDieStg2  <- Q.imcmc[,qualOrder_history$qualOrder_DieStg2]
                TqdieStg2 <- qual_history$qual_DieStg2
                ## Take median and 95%CI of estimated qualities
                qdieStg2.M <- apply(qDieStg2, 2, function(x) median(x))
                qdieStg2.L <- apply(qDieStg2, 2, function(x) quantile(x, 0.025))
                qdieStg2.U <- apply(qDieStg2, 2, function(x) quantile(x, 0.975))
                storeqDieStg2 <- rbind(TqdieStg2,qdieStg2.M,qdieStg2.L,qdieStg2.U)
            }
            ## store
            qTotDieStg2 <- cbind(qTotDieStg2, storeqDieStg2)
            ## ##########################################################
            ## TAKE DATA OF Q-INDIVIDUALS THAT GOT CENSORED IN STAGE 1 ##
            ## 
            ## Take estimates & true values of 'q' associated to individual/s that died in Stg 1
            if (length(qualOrder_history$qualOrder_CensStg1)==0) {
                qCenStg1  <- matrix(99,nrow(samples),1)
                TqdieStg1 <- 99
                storeqCenStg1 <- matrix(rep(99,4),ncol=1) ## by row: TqdieStg1, qdieStg1.M, qdieStg1.L, qdieStg1.U
                ##
            } else if (length(qualOrder_history$qualOrder_CensStg1)==1) { ## vector
                qCenStg1   <- Q.imcmc[,qualOrder_history$qualOrder_CensStg1]
                (TqdieStg1 <- qual_history$qual_CensStg1) ## true value
                ## Take median and 95%CI of estimated qualities
                (qdieStg1.M <- median(qCenStg1))
                (qdieStg1.L <- quantile(qCenStg1, 0.025))
                (qdieStg1.U <- quantile(qCenStg1, 0.975))
                storeqCenStg1 <- matrix(c(TqdieStg1,qdieStg1.M,qdieStg1.L,qdieStg1.U),ncol=1) ## by row: TqdieStg1, qdieStg1.M, qdieStg1.L, qdieStg1.U
            } else { ## matrix
                qCenStg1  <- Q.imcmc[,qualOrder_history$qualOrder_CensStg1]
                TqdieStg1 <- qual_history$qual_CensStg1
                ## Take median and 95%CI of estimated qualities
                qdieStg1.M <- apply(qCenStg1, 2, function(x) median(x))
                qdieStg1.L <- apply(qCenStg1, 2, function(x) quantile(x, 0.025))
                qdieStg1.U <- apply(qCenStg1, 2, function(x) quantile(x, 0.975))
                storeqCenStg1 <- rbind(TqdieStg1,qdieStg1.M,qdieStg1.L,qdieStg1.U)
            }
            ## store
            qTotCenStg1 <- cbind(qTotCenStg1, storeqCenStg1)
            ## ##########################################################
            ## TAKE DATA OF Q-INDIVIDUALS THAT GOT CENSORED IN STAGE 2 ##
            ## 
            ## Take estimates & true values of 'q' associated to individual/s that died in Stg 1
            if (length(qualOrder_history$qualOrder_CensStg2)==0) {
                qCenStg2  <- matrix(99,nrow(samples),1)
                TqdieStg2 <- 99
                storeqCenStg2 <- matrix(rep(99,4),ncol=1) ## by row: TqdieStg2, qdieStg2.M, qdieStg2.L, qdieStg2.U
                ##
            } else if (length(qualOrder_history$qualOrder_CensStg2)==1) { ## vector
                qCenStg2   <- Q.imcmc[,qualOrder_history$qualOrder_CensStg2]
                (TqdieStg2 <- qual_history$qual_CensStg2) ## true value
                ## Take median and 95%CI of estimated qualities
                (qdieStg2.M <- median(qCenStg2))
                (qdieStg2.L <- quantile(qCenStg2, 0.025))
                (qdieStg2.U <- quantile(qCenStg2, 0.975))
                storeqCenStg2 <- matrix(c(TqdieStg2,qdieStg2.M,qdieStg2.L,qdieStg2.U),ncol=1) ## by row: TqdieStg2, qdieStg2.M, qdieStg2.L, qdieStg2.U
            } else { ## matrix
                qCenStg2  <- Q.imcmc[,qualOrder_history$qualOrder_CensStg2]
                TqdieStg2 <- qual_history$qual_CensStg2
                ## Take median and 95%CI of estimated qualities
                qdieStg2.M <- apply(qCenStg2, 2, function(x) median(x))
                qdieStg2.L <- apply(qCenStg2, 2, function(x) quantile(x, 0.025))
                qdieStg2.U <- apply(qCenStg2, 2, function(x) quantile(x, 0.975))
                storeqCenStg2 <- rbind(TqdieStg2,qdieStg2.M,qdieStg2.L,qdieStg2.U)
            }
            ## store
            qTotCenStg2 <- cbind(qTotCenStg2, storeqCenStg2)
            ## #######################################################
            ## TAKE DATA OF Q-INDIVIDUALS THAT COMPLETE BOTH STAGES ##
            ## 
            ## Take quality true values AND estimates of the quality associated to individual/s that fully developed
            (TfullDev <- qual_history$qual_FullDev)
            qFullDev  <- Q.imcmc[,qualOrder_history$qualOrder_FullDev]
            dim(qFullDev); length(TfullDev)
            ## Take median and 95%CI of estimated qualities here
            (qfullDevM <- apply(qFullDev, 2, function(x) median(x)))
            (qfullDevL <- apply(qFullDev, 2, function(x) quantile(x, 0.025)))
            (qfullDevU <- apply(qFullDev, 2, function(x) quantile(x, 0.975)))
            ## store
            qTotFullDev <- cbind(qTotFullDev, rbind(TfullDev,qfullDevM,qfullDevL,qfullDevU))
        }
    }
}

## Remove first column of zeros of qTotDieStg1 & qTotFullDev
qTotDieStg1[,1:5]
qTotDieStg1 <-qTotDieStg1[,-1]
qTotDieStg2 <-qTotDieStg2[,-1]
qTotCenStg1 <-qTotCenStg1[,-1]
qTotCenStg2 <-qTotCenStg2[,-1]
qTotFullDev <- qTotFullDev[,-1]
## Check
qTotDieStg1[,1:5]
qTotFullDev[,1:5]

## Check
dim(qTotDieStg1) ## ncols: total indiv that die in Stg1 over the imcmc simulations
dim(qTotDieStg2) ## ncols: total indiv that die in Stg2 over the imcmc simulations
dim(qTotFullDev) ## ncols: total indiv that fully develop over the imcmc simulations
dim(Tqual)       ## ncol : total number of individuals, equals imcmc X N


## ###########################################
## COUNT THE NUMBER OF OUTLIERS IN EVERY CASE

## Total number of individuals that Die in  Stg 1 & nb. of outliers
ncol(qTotDieStg1) ## total
outl1 <- sum(apply(qTotDieStg1, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1) 
outl2 <- sum(apply(qTotDieStg1, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1 + outl2)/ncol(qTotDieStg1)

## Total number of individuals that Die in  Stg 2 & nb. of outliers
ncol(qTotDieStg2) ## total
outl1 <- sum(apply(qTotDieStg2, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1) 
outl2 <- sum(apply(qTotDieStg2, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1 + outl2)/ncol(qTotDieStg2)

## Total number of individuals that are censored in Stg1  & nb. of outliers
ncol(qTotCenStg1) ## total
outl1 <- sum(apply(qTotCenStg1, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1) 
outl2 <- sum(apply(qTotCenStg1, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1 + outl2)/ncol(qTotCenStg1)

## Total number of individuals that are censored in Stg2  & nb. of outliers
ncol(qTotCenStg2) ## total
outl1 <- sum(apply(qTotCenStg2, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1) 
outl2 <- sum(apply(qTotCenStg2, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1 + outl2)/ncol(qTotCenStg2)

## Total number of individuals that fully develop & nb. of outliers
ncol(qTotFullDev) ## total
outl1 <- sum(apply(qTotFullDev, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1) 
outl2 <- sum(apply(qTotFullDev, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1)
outl1; outl2
## proportion:
(outl1 + outl2)/ncol(qTotFullDev)


##########################
## PREPARE FOR PLOTTING ##
##########################

## #############
## Q-die Stg1 ##
## check reordering is OK: 
trueQdieStg1 <- sort(qTotDieStg1[1,])
as.matrix(qTotDieStg1[,which(qTotDieStg1[1,]==trueQdieStg1[1])])
qTotDieStg1[,order(qTotDieStg1[1,])][,1:5] ## first column must match the previous line
##
## ORDER the matrix & REMOVE the latest columns full of 99 (i.e. when no indiv dies in a given simulation)
qTotDieStg1 <- qTotDieStg1[,order(trueQdieStg1)]
tail(trueQdieStg1) ## next, remove columns with 99s
dim(qTotDieStg1)
qTotDieStg1 <- qTotDieStg1[,-which(qTotDieStg1[1,]==99)]
dim(qTotDieStg1)

## #############
## Q-die Stg2 ##
## check reordering is OK: 
trueQdieStg2 <- sort(qTotDieStg2[1,])
as.matrix(qTotDieStg2[,which(qTotDieStg2[1,]==trueQdieStg2[1])])
qTotDieStg2[,order(qTotDieStg2[1,])][,1:5] ## first column must match the previous line
##
## ORDER the matrix & REMOVE the latest columns full of 99 (i.e. when no indiv dies in a given simulation)
qTotDieStg2 <- qTotDieStg2[,order(trueQdieStg2)]
tail(trueQdieStg2) ## next, remove columns with 99s
dim(qTotDieStg2)
qTotDieStg2 <- qTotDieStg2[,-which(qTotDieStg2[1,]==99)]
dim(qTotDieStg2)


## #############
## Q-cen Stg1 ##
## check reordering is OK: 
trueQcenStg1 <- sort(qTotCenStg1[1,])
as.matrix(qTotCenStg1[,which(qTotCenStg1[1,]==trueQcenStg1[1])])
qTotCenStg1[,order(qTotCenStg1[1,])][,1:5] ## first column must match the previous line
##
## ORDER the matrix & REMOVE the latest columns full of 99 (i.e. when no indiv cens in a given simulation)
qTotCenStg1 <- qTotCenStg1[,order(trueQcenStg1)]
tail(trueQcenStg1) ## next, remove columns with 99s
dim(qTotCenStg1)
qTotCenStg1 <- qTotCenStg1[,-which(qTotCenStg1[1,]==99)]
dim(qTotCenStg1)
## #############
## Q-cen Stg2 ##
## check reordering is OK: 
trueQcenStg2 <- sort(qTotCenStg2[1,])
as.matrix(qTotCenStg2[,which(qTotCenStg2[1,]==trueQcenStg2[1])])
qTotCenStg2[,order(qTotCenStg2[1,])][,1:5] ## first column must match the previous line
##
## ORDER the matrix & REMOVE the latest columns full of 99 (i.e. when no indiv cens in a given simulation)
qTotCenStg2 <- qTotCenStg2[,order(trueQcenStg2)]
tail(trueQcenStg2) ## next, remove columns with 99s
dim(qTotCenStg2)
qTotCenStg2 <- qTotCenStg2[,-which(qTotCenStg2[1,]==99)]
dim(qTotCenStg2)



## #############
## Q-Full Dev ##
## check reordering is OK: 
trueQfullDev <- sort(qTotFullDev[1,])
as.matrix(qTotFullDev[,which(qTotFullDev[1,]==trueQfullDev[1])])
qTotFullDev[,order(qTotFullDev[1,])][,1:5] ## first column must match the previous line

## ########
## PLOTS ##
setwd(baseDir)
setwd(paste(getwd(), "figures",sep="/"))


## ######################################
## ## PLOT estimated vs true parameters #
if (PDF)
    pdf(paste0("qEstimates_byCases.pdf"), width = 25, height = 7)
par(mfrow=c(1,5), mai=c(1, 0.8, 1, 0.1), omi=c(0.1, 0.1, 0.1, 0.15))
cexLab  <- 2.7
cexMain <- 2.6
Cex    <- 0.6
LWD    <- 1.2
LX     <- 3.9
LY     <- 3.3
LZ     <- 1.5
LYnu   <- 3.3 
Caxis  <- 3
PadjX  <- 0.45
PadjY  <- 0.05
Hadj   <- 0.5
## #################
## INITIALISE PLOT 1 
plot(qTotDieStg1[1,], qTotDieStg1[1,], pch=19, cex=Cex, type="n", xlab="",
     ylab="", cex.axis=Caxis, col="darkgreen", xaxt="n", yaxt="n")
AT <- seq(0, 1, 0.2)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis,  padj=PadjX)
axis(2, at=AT, lab=AT, las=3, cex.axis=Caxis,  padj=PadjY)
## plot median & 95%CI (lower & upper bounds)
points(qTotDieStg1[1,], qTotDieStg1[2,], pch=19, cex=0.2, col="darkgreen")
points(qTotDieStg1[1,], qTotDieStg1[3,], pch=19, cex=0.2, col="blue")
points(qTotDieStg1[1,], qTotDieStg1[4,], pch=19, cex=0.2, col="red")
lines(seq(0,1,0.01), seq(0,1,0.01), type="l", lty=2, lwd=3)
# nb. of individuals plotted here: ncol(qTotDieStg1)
mtext("Mortality Stage 1", 3, line=LZ, cex=cexMain)
mtext("q", 1, line=LX, cex=cexLab) 
mtext(TeX("$ \\hat{q} $"),  2, line=LY, cex=cexLab, las=2)
## #################
## INITIALISE PLOT 2 
plot(qTotDieStg2[1,], qTotDieStg2[1,], pch=19, cex=Cex, type="n", xlab="",
     ylab="", cex.axis=Caxis, col="darkgreen", xaxt="n", yaxt="n")
AT <- seq(0, 1, 0.2)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis,  padj=PadjX)
axis(2, at=AT, lab=AT, las=3, cex.axis=Caxis,  padj=PadjY)
## plot median & 95%CI (lower & upper bounds)
points(qTotDieStg2[1,], qTotDieStg2[2,], pch=19, cex=0.2, col="darkgreen")
points(qTotDieStg2[1,], qTotDieStg2[3,], pch=19, cex=0.2, col="blue")
points(qTotDieStg2[1,], qTotDieStg2[4,], pch=19, cex=0.2, col="red")
mtext("Mortality Stage 2", 3, line=LZ, cex=cexMain)
lines( seq(0,1,0.01), seq(0,1,0.01), type="l", lty=2, lwd=3)
mtext("q", 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{q} $"),  2, line=LY, cex=cexLab, las=2)
## ## #################
## ## INITIALISE PLOT 3 
plot(qTotCenStg1[1,], qTotCenStg1[1,], pch=19, cex=Cex, type="n", xlim=c(0,1), ylim=c(0,1), xlab="",
     ylab="", cex.axis=Caxis, col="darkgreen", xaxt="n", yaxt="n")
AT <- seq(0, 1, 0.2)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis,  padj=PadjX)
axis(2, at=AT, lab=AT, las=3, cex.axis=Caxis,  padj=PadjY)
## plot median & 95%CI (lower & upper bounds)
points(qTotCenStg1[1,], qTotCenStg1[2,], pch=19, cex=0.5, col="darkgreen")
points(qTotCenStg1[1,], qTotCenStg1[3,], pch=19, cex=0.5, col="blue")
points(qTotCenStg1[1,], qTotCenStg1[4,], pch=19, cex=0.5, col="red")
lines( seq(0,1,0.01), seq(0,1,0.01), type="l", lty=2, lwd=3)
mtext("Censored Stage 1", 3, line=LZ, cex=cexMain)
mtext("q", 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{q} $"),  2, line=LY, cex=cexLab, las=2)
## #################
## INITIALISE PLOT 4 
plot(qTotCenStg2[1,], qTotCenStg2[1,], pch=19, cex=Cex, type="n", xlim=c(0,1), ylim=c(0,1), xlab="",
     ylab="", cex.axis=Caxis, col="darkgreen", xaxt="n", yaxt="n")
AT <- seq(0, 1, 0.2)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis,  padj=PadjX)
axis(2, at=AT, lab=AT, las=3, cex.axis=Caxis,  padj=PadjY)
## plot median & 95%CI (lower & upper bounds)
points(qTotCenStg2[1,], qTotCenStg2[2,], pch=19, cex=0.5, col="darkgreen")
points(qTotCenStg2[1,], qTotCenStg2[3,], pch=19, cex=0.5, col="blue")
points(qTotCenStg2[1,], qTotCenStg2[4,], pch=19, cex=0.5, col="red")
lines( seq(0,1,0.01), seq(0,1,0.01), type="l", lty=2, lwd=3)
mtext("Censored Stage 2", 3, line=LZ, cex=cexMain)
mtext("q", 1, line=LX, cex=cexLab) #TeX("$ \\q $")
mtext(TeX("$ \\hat{q} $"),  2, line=LY, cex=cexLab, las=2)
## #################
## INITIALISE PLOT 5 
####################
plot(qTotFullDev[1,], qTotFullDev[1,], pch=19, cex=Cex, type="n", xlab="",
     ylab="", cex.axis=Caxis, col="darkgreen", xaxt="n", yaxt="n")
AT <- seq(0, 1, 0.2)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis,  padj=PadjX)
axis(2, at=AT, lab=AT, las=3, cex.axis=Caxis,  padj=PadjY)
## plot median & 95%CI (lower & upper bounds)
points(qTotFullDev[1,], qTotFullDev[3,], pch=19, cex=0.2, col="blue")
points(qTotFullDev[1,], qTotFullDev[4,], pch=19, cex=0.2, col="red")
points(qTotFullDev[1,], qTotFullDev[2,], pch=19, cex=0.2, col="darkgreen")
lines( seq(0,1,0.01), seq(0,1,0.01), type="l", lty=2, lwd=3)
mtext("Fully developed", 3, line=LZ, cex=cexMain)
mtext("q", 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{q} $"),  2, line=LY, cex=cexLab, las=2)
if (PDF)
    dev.off()




















