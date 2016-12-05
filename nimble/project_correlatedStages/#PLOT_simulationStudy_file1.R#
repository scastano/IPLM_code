#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ############################################
## Plots associated to case study I - file 1 ##
## ############################################
rm(list=ls())

## Set directories and basic parameters
options(width=160)
library(latex2exp)
library(nimble)
library(plotrix)

baseDir <- "~/IPLM_paper/nimble/project_correlatedStages/"
MCMC    <- "mcmc"

PDF <- TRUE ## FALSE
PLOT_KernelEstimates <- TRUE ## FALSE

setwd(baseDir)
source("../FUNCTIONS_R.R")
source("../FUNCTIONS_NIMBLE.R")
source("../project_culicoides/FUNCTIONS_CULICOIDES.R")

## Compile nf_TW1_MuSigSurv
thresh         <- 1E-6
cTW1_MuSigSurv <- compileNimble(nf_TW1_MuSigSurv)

## ############################
## Read mcmc output directories
setwd(baseDir)
mcmcDir <- paste(getwd(), MCMC, sep='/') 
setwd(mcmcDir)
## Find mcmc results
(mcmcDirs <- system("ls -1t|head -n 1000", TRUE)) ## >500 simulations
nsims     <- length(mcmcDirs)
    
## Define matrices to store median and 95%CI of estimated mu,sig and surv (by stage)
medianStg1 <- qLow_Stg1 <- qUpp_Stg1 <- matrix(0, nsims, 3)
medianStg2 <- qLow_Stg2 <- qUpp_Stg2 <- matrix(0, nsims, 3)
colnames(medianStg1) <- colnames(qLow_Stg1) <- colnames(qUpp_Stg1) <- c("mu1","sig1","surv1")
colnames(medianStg2) <- colnames(qLow_Stg2) <- colnames(qUpp_Stg2) <- c("mu2","sig2","surv2")
## Define matrices to store median and 95%CI of estimated res1 and res2 (two stages)
Tres <- medianRes <- qLowRes <- qUppRes <- matrix(0, nsims, 2)
colnames(Tres) <- colnames(medianRes) <- colnames(qLowRes) <- colnames(qUppRes) <- c("res1","res2")
## Define matrix to store TRUE mu,sig and surv of both stages
trueMuSigSurv <- matrix(0, nsims, 6)
colnames(trueMuSigSurv) <- c("mu1","sig1","surv1","mu2","sig2","surv2")
## Define vectors to store median and 95%CI of TRUE  rho
TrueRho <- medianRho <- qLowRho <- qUppRho <- numeric()
## Define vector to store number of individuals fully developed in every simulated data
Fulldev <- numeric()
## Define matrix to store true qualities by simulation (row)
Tqual  <- matrix(0, nsims, 50)
## Define matrices to store median and 95%CI of estimated res1 and res2 (two stages)
Mq <- Uq <- Lq <- matrix(0, nsims, 50)
mindev <- 35

######################
## THIS CAN BE SLOW ##
imcmc <- 0

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
            print(imcmc)
            ## Get parameters & data
            N    <- Constants[["N"]]
            S    <- Constants[["S"]]
            TrueRho[imcmc] <- Inits[["rho"]]
            Tres[imcmc,]   <- Inits[["res"]]
            Tqual[imcmc,]  <- Inits[["qual"]]
            ##
            Fulldev[imcmc] <- fulldev
            ## Burn-in Samples
            samples <- tail(samples, nrow(samples) - round(nrow(samples)/2))  ## Ad-hoc burn-in for working with running chains
            (nlines <- nrow(samples))
            ##
            colnames(samples)
            ## Take estimated qualities
            Q.imcmc <- samples[,paste0("qual" , 1:N)]
            ## Take median and 95%CI of estimated qualities
            Mq[imcmc,] <- apply(Q.imcmc, 2, function(x) median(x))
            Lq[imcmc,] <- apply(Q.imcmc, 2, function(x) quantile(x, 0.025))
            Uq[imcmc,] <- apply(Q.imcmc, 2, function(x) quantile(x, 0.975))
            ## #######################
            ## Take samples by stage #
            ##Define set of params by stage
            parStg1 <- c("mu1", "sc1", "sur1", "res1")
            parStg2 <- c("mu2", "sc2", "sur2", "res2")
            lpar <- length(parStg1)
            parMatStg1 <- samples[,parStg1]
            parMatStg2 <- samples[,parStg2]
            ## ###############################################
            ## Store resolution of every stage by simulation #
            resmat <- samples[,c("res1","res2")]
            ## the same for rho
            rhovec <- samples[,"rho"]
            ## ###################################################################
            ## Define matrices to store moments of travelling wave by simulation #
            muSigSurvStg1 <- matrix(0,nrow(samples), 3)
            muSigSurvStg2 <- matrix(0,nrow(samples), 3)
            colnames(muSigSurvStg1) <- c("mu","sig","surv")
            colnames(muSigSurvStg2) <- c("mu","sig","surv")
            for (i in 1:nrow(samples)) {
                muSigSurvStg1[i,] <- cTW1_MuSigSurv(paras = parMatStg1[i,-lpar], res = parMatStg1[i,lpar], thresh = thresh)
                muSigSurvStg2[i,] <- cTW1_MuSigSurv(paras = parMatStg2[i,-lpar], res = parMatStg2[i,lpar], thresh = thresh)
            }
            ## #################################
            ## Calculate moments for true values 
            truePar <- Inits$paras
            trueRes <- Inits$res
            trueMuSigSurv[imcmc,] <- c(cTW1_MuSigSurv(paras = truePar[1,], res = trueRes[1], thresh = thresh),
                                       cTW1_MuSigSurv(paras = truePar[2,], res = trueRes[2], thresh = thresh))
            ## Get mean and 95%CI of moments and survival & res
            medianStg1[imcmc,] <- t(apply(muSigSurvStg1, 2, function(x) median(x)))
            medianStg2[imcmc,] <- t(apply(muSigSurvStg2, 2, function(x) median(x)))
            qLow_Stg1[imcmc,]  <- t(apply(muSigSurvStg1, 2, function(x) quantile(x, 0.025)))
            qLow_Stg2[imcmc,]  <- t(apply(muSigSurvStg2, 2, function(x) quantile(x, 0.025)))
            qUpp_Stg1[imcmc,]  <- t(apply(muSigSurvStg1, 2, function(x) quantile(x, 0.975)))
            qUpp_Stg2[imcmc,]  <- t(apply(muSigSurvStg2, 2, function(x) quantile(x, 0.975)))
            ## the same for resolution
            medianRes[imcmc,] <- c(t(apply(resmat, 2, function(x) median(x))))
            qLowRes[imcmc,]  <- t(apply(resmat, 2, function(x) quantile(x, 0.025)))
            qUppRes[imcmc,]  <- t(apply(resmat, 2, function(x) quantile(x, 0.975)))
            ## the same for rho
            medianRho[imcmc] <- median(rhovec)
            qLowRho[imcmc]   <- quantile(rhovec, 0.025)
            qUppRho[imcmc]   <- quantile(rhovec, 0.975)
        }
    }
}


## ########
## PLOTS ##
setwd(baseDir)
setwd(paste(getwd(), "figures",sep="/"))

## ################################
## REORDER FOR PLOTTING PURPOSES ##
## ################################

## #########
## STAGE 1 #
## mu
Tmu1 <- sort(trueMuSigSurv[,"mu1"])
Mmu1 <- medianStg1[, "mu1"][order(trueMuSigSurv[,"mu1"])]
Lmu1 <- qLow_Stg1[ , "mu1"][order(trueMuSigSurv[,"mu1"])]
Umu1 <- qUpp_Stg1[ , "mu1"][order(trueMuSigSurv[,"mu1"])]
## sigma
Tsig1 <- sort(trueMuSigSurv[,"sig1"])
Msig1 <- medianStg1[, "sig1"][order(trueMuSigSurv[,"sig1"])]
Lsig1 <- qLow_Stg1[ , "sig1"][order(trueMuSigSurv[,"sig1"])]
Usig1 <- qUpp_Stg1[ , "sig1"][order(trueMuSigSurv[,"sig1"])]
## survival
Tsurv1 <- sort(trueMuSigSurv[,"surv1"])
Msurv1 <- medianStg1[, "surv1"][order(trueMuSigSurv[,"surv1"])]
Lsurv1 <- qLow_Stg1[ , "surv1"][order(trueMuSigSurv[,"surv1"])]
Usurv1 <- qUpp_Stg1[ , "surv1"][order(trueMuSigSurv[,"surv1"])]
## res
Tres1 <- sort(Tres[,"res1"])
Mres1 <- medianRes[,"res1"][order(Tres[,"res1"])]
Lres1 <- qLowRes[ , "res1"][order(Tres[,"res1"])]
Ures1 <- qUppRes[ , "res1"][order(Tres[,"res1"])]

## #########
## STAGE 2 #
## mu
Tmu2 <- sort(trueMuSigSurv[,"mu2"])
Mmu2 <- medianStg2[, "mu2"][order(trueMuSigSurv[,"mu2"])]
Lmu2 <- qLow_Stg2[ , "mu2"][order(trueMuSigSurv[,"mu2"])]
Umu2 <- qUpp_Stg2[ , "mu2"][order(trueMuSigSurv[,"mu2"])]
## sigma
Tsig2 <- sort(trueMuSigSurv[,"sig2"])
Msig2 <- medianStg2[, "sig2"][order(trueMuSigSurv[,"sig2"])]
Lsig2 <- qLow_Stg2[ , "sig2"][order(trueMuSigSurv[,"sig2"])]
Usig2 <- qUpp_Stg2[ , "sig2"][order(trueMuSigSurv[,"sig2"])]
## survival
Tsurv2 <- sort(trueMuSigSurv[,"surv2"])
Msurv2 <- medianStg2[, "surv2"][order(trueMuSigSurv[,"surv2"])]
Lsurv2 <- qLow_Stg2[ , "surv2"][order(trueMuSigSurv[,"surv2"])]
Usurv2 <- qUpp_Stg2[ , "surv2"][order(trueMuSigSurv[,"surv2"])]
## res
Tres2 <- sort(Tres[,"res2"])
Mres2 <- medianRes[,"res2"][order(Tres[,"res2"])]
Lres2 <- qLowRes[ , "res2"][order(Tres[,"res2"])]
Ures2 <- qUppRes[ , "res2"][order(Tres[,"res2"])]

## #####
## RHO #
Trho <- sort(TrueRho)
Mrho <- medianRho[order(TrueRho)]
Lrho <- qLowRho[order(TrueRho)]
Urho <- qUppRho[order(TrueRho)]

## ###########
## QUALITIES #
Tsortedqual <- sort(Tqual[1:imcmc,])
Mq <- Mq[order(Tqual)]
Lq <- Lq[order(Tqual)]
Uq <- Uq[order(Tqual)]

## Count for reporting
median(Mres1)
median(Lres1)
median(Ures1)
median(Mres2)
median(Lres2)
median(Ures2)


## #################
## COUNT OUTLIERS ##

## STAGE 1
Mu1 <- rbind(Tmu1, Mmu1, Lmu1, Umu1)
## Total number of individuals that Die in Stg 1 & nb. of outliers
outl1 <- sum(apply(Mu1, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1)
outl2 <- sum(apply(Mu1, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1+outl2)/imcmc

Sig1 <- rbind(Tsig1, Msig1, Lsig1, Usig1)
## Total number of individuals that Die in  Stg 1 & nb. of outliers
outl1 <- sum(apply(Sig1, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1)
outl2 <- sum(apply(Sig1, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1+outl2)/imcmc

Surv1 <- rbind(Tsurv1, Msurv1, Lsurv1, Usurv1)
## Total number of individuals that Die in  Stg 1 & nb. of outliers
outl1 <- sum(apply(Surv1, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1)
outl2 <- sum(apply(Surv1, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1+outl2)/imcmc

## STAGE 2
Mu2 <- rbind(Tmu2, Mmu2, Lmu2, Umu2)
## Total number of individuals that Die in  Stg 1 & nb. of outliers
outl1 <- sum(apply(Mu2, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1)
outl2 <- sum(apply(Mu2, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1+outl2)/imcmc

Sig2 <- rbind(Tsig2, Msig2, Lsig2, Usig2)
## Total number of individuals that Die in  Stg 1 & nb. of outliers
outl1 <- sum(apply(Sig2, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1)
outl2 <- sum(apply(Sig2, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1+outl2)/imcmc

Surv2 <- rbind(Tsurv2, Msurv2, Lsurv2, Usurv2)
## Total number of individuals that Die in  Stg 1 & nb. of outliers
outl1 <- sum(apply(Surv2, 2, function(x) x[1] > x[4])) ## upper bound (elem4) lower than true value (elem1)
outl2 <- sum(apply(Surv2, 2, function(x) x[1] < x[3])) ## lower bound (elem3) greater than true value (elem1) 
outl1; outl2
## proportion:
(outl1+outl2)/imcmc



## ######################################################################################
## ## PLOT estimated WTD parameters vs true  WTD parameters (mean, sigma, total survival) 
if (PDF)
    pdf("MuSigSurv2_against_True.pdf", width = 12, height = 7.5)
MAR <- c(3, 20, 3, 3)
par(mar=MAR, mfrow=c(2,3))
par(mai=c(0.6, 0.9, 0.2, 0.1),  omi=c(0.1, 0.3, 0, 0.1)) 
cexLab <- 2.2
cexM   <- 0.4
cexUL  <- 0.4
LWD    <- 1.2
LX     <- 3.9
LY     <- 2.9
LYnu   <- 3.3 
Caxis  <- 2.2
## plot MU1
plot(  Tmu1, Mmu1, pch=19, cex=cexUL, type="n", xlab="", ylab="", cex.axis=Caxis, col="darkgreen")
points(Tmu1, Umu1, pch=19, cex=cexUL, col="red")
points(Tmu1, Lmu1, pch=19, cex=cexUL, col="blue")
points(Tmu1, Mmu1, pch=19, cex=cexUL, col="darkgreen")
lines(Tmu1,   Tmu1,  type="l", lty=2, lwd=LWD)
mtext(TeX("$ \\tilde{\\mu}_1 $"), 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\tilde{\\mu}}_1 $"),  2, line=LY, cex=cexLab, las=2)
## plot SIG1
plot(  Tsig1, Msig1, pch=19, cex=cexUL, type="n", xlab="", ylab="", cex.axis=Caxis,
     xaxt="n", yaxt="n", col="darkgreen")
AT <- seq(0, 2, 1)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Tsig1, Usig1, pch=19, cex=cexUL, col="red")
points(Tsig1, Lsig1, pch=19, cex=cexUL, col="blue")
points(Tsig1, Msig1, pch=19, cex=cexUL, col="darkgreen") 
lines(Tsig1,  Tsig1, type="l", lty=2, lwd=LWD)
mtext(TeX("$ \\tilde{\\sigma}_1 $"), 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\tilde{\\sigma}}_1 $"), 2, line=LY, cex=cexLab, las=2)
## plot SURV1
plot(Tsurv1, Msurv1, pch=19, cex=cexUL, type="n", xlab="", ylab="", cex.axis=Caxis,
     col="darkgreen", yaxt="n", ylim=c(0.55, 1)) ## ylim=c(0.65,1.05))
AT <- seq(0.5, 1.1, 0.1)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)#, cex.axis=CET.AXIS)
points(Tsurv1, Usurv1, pch=19, cex=cexUL, col="red")
points(Tsurv1, Lsurv1, pch=19, cex=cexUL, col="blue")
points(Tsurv1, Msurv1, pch=19, cex=cexUL, col="darkgreen")
lines(Tsurv1,   Tsurv1,  type="l", lty=2, lwd=LWD)
## for (i in 1:N)
##     abline(h=i/N, lty=2)
mtext(TeX("$ \\tilde{\\nu}_1 $"),     1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\tilde{\\nu}}_1 $"), 2, line=LYnu, cex=cexLab, las=2)
## plot MU2
plot(  Tmu2, Mmu2, pch=19, cex=cexUL, type="n", xlab="", ylab="", cex.axis=Caxis,
     col="darkgreen", xaxt="n", yaxt="n")#, xlim=c(1,1)) 
AT <- seq(2, 12, 2)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Tmu2, Umu2, pch=19, cex=cexUL, col="red")
points(Tmu2, Lmu2, pch=19, cex=cexUL, col="blue")
points(Tmu2, Mmu2, pch=19, cex=cexUL, col="darkgreen") 
lines(Tmu2,   Tmu2,  type="l", lty=2, lwd=LWD)
mtext(TeX("$ \\tilde{\\mu}_{2} $"),     1, line=LX, cex=cexLab)#, ylab=seq(2,14,2))
mtext(TeX("$ \\hat{\\tilde{\\mu}}_2 $"), 2, line=LY, cex=cexLab, las=2)
## plot SIG2 
plot(  Tsig2, Msig2, pch=19, cex=cexUL, type="n", xlab="", ylab="", cex.axis=Caxis, 
     xaxt="n", yaxt="n", col="darkgreen")# , xlim=c(0,7), ylim=c(0,7))
AT <- seq(0, 6, 1)
axis(1, at=AT, lab=AT, las=1, cex.axis=Caxis)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Tsig2, Usig2, pch=19, cex=cexUL, col="red")
points(Tsig2, Lsig2, pch=19, cex=cexUL, col="blue")
points(Tsig2, Msig2, pch=19, cex=cexUL, col="darkgreen")
lines(Tsig2,  Tsig2, type="l", lty=2, lwd=LWD)
mtext(TeX("$ \\tilde{\\sigma}_2 $"), 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\tilde{\\sigma}}_2 $"), 2, line=LY, cex=cexLab, las=2)
## plot SURV2
plot(Tsurv2, Msurv2, pch=19, cex=cexUL, type="n", xlab="", ylab="", cex.axis=Caxis,
     col="darkgreen", yaxt="n", ylim=c(0.6, 1))
AT <- seq(0.4, 1.1, 0.1)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Tsurv2, Usurv2, pch=19, cex=cexUL, col="red")
points(Tsurv2, Lsurv2, pch=19, cex=cexUL, col="blue")
points(Tsurv2, Msurv2, pch=19, cex=cexUL, col="darkgreen")
## for (i in 1:N)
##     abline(h=i/N, lty=2)
lines(Tsurv2, Tsurv2,  type="l", lty=2, lwd=LWD)
mtext(TeX("$ \\tilde{\\nu}_2 $"),     1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\tilde{\\nu}}_2 $"), 2, line=LYnu, cex=cexLab, las=2)
##
if (PDF)
    dev.off()



## PLOT  for supplementary material
if (PDF)
 pdf("Res_Rho_FullDevHistog1.pdf", width = 12, height = 10)
MAR <- c(3, 20, 3, 3)
par(mar=MAR, mfrow=c(2,2))
par(mai=c(0.9, 1.1, 0.2, 0.1),  omi=c(0.2, 0.2, 0.1, 0.15)) 
cexLab <- 2.8
cexM   <- 0.4
cexUL  <- 0.6
LWD    <- 1.2
LX     <- 3.5
LY     <- 3.4
LYnu   <- 3.3 
Caxis  <- 2.2
cexHist <- 2.5
## plot RES 
plot(Tres1,   Mres1, pch=19, cex=cexUL, xlab="", ylab="", ylim=c(0,100),
     col="darkgreen", cex.axis=Caxis, yaxt="n")
AT <- seq(0, 100, 20)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Tres1, Ures1, pch=19, cex=cexUL, col="red")
points(Tres1, Lres1, pch=19, cex=cexUL, col="blue")
mtext(TeX("$ r_{1} $"), 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\r}_1 $"), 2, line=LY, cex=cexLab, las=2)
lines(Tres1,   Tres1,  type="l", lty=2, lwd=LWD)
##
plot(Tres2,   Mres2, pch=19, cex=cexUL, xlab="", ylab="", ylim=c(0,100),
     col="darkgreen", cex.axis=Caxis, yaxt="n")
AT <- seq(0, 100, 20)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Tres2, Ures2, pch=19, cex=cexUL, col="red")
points(Tres2, Lres2, pch=19, cex=cexUL, col="blue")
mtext(TeX("$ r_{2} $"), 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\r}_2 $"), 2, line=LY, cex=cexLab, las=2)
lines(Tres2,   Tres2,  type="l", lty=2, lwd=LWD)
##
## plot RHO
plot(Trho,   Mrho, pch=19, cex=cexUL,  xlab="", ylab="", ylim=c(0,1),
     col="darkgreen", cex.axis=Caxis, yaxt="n")
AT <- seq(0, 1, 0.2)
axis(2, at=AT, lab=AT, las=1, cex.axis=Caxis)
points(Trho, Urho, pch=19, cex=cexUL, col="red")
points(Trho, Lrho, pch=19, cex=cexUL, col="blue")
lines(Trho,Trho,  type="l", lty=2, lwd=1.4)
mtext(TeX("$ \\rho $"), 1, line=LX, cex=cexLab)
mtext(TeX("$ \\hat{\\rho} $"), 2, line=3.7, cex=cexLab, las=2)
## plot HISTOGRAM FULL DEV
hist(Fulldev, breaks=(mindev-1):N, main="", cex=cexLab,   xlab="", ylab="", cex.axis=Caxis)
mtext("Nb. individuals fully developed", 1, line=3.5, cex=2)
mtext("Frequency",      2, line=3.5, cex=cexHist)
if (PDF)
    dev.off()














