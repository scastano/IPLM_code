#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

## ##############################################
## Plots case study II - Resolution histograms ##
## ##############################################
rm(list=ls())
options(width=160)
PDF <- TRUE ## FALSE
baseDir <- "~/IPLM_paper/nimble/project_culicoides/"
setwd(baseDir)

## ############################
## Read mcmc output directories
setwd(baseDir)
readDir <- paste(getwd(), "mcmcOutput4plots", sep='/') 
setwd(readDir)
##
resE  <- as.matrix(read.csv("resSamples_egg_10000.csv"))
resA  <- as.matrix(read.csv("resSamples_gonotrophicCycle_10000.csv"))
resLP <- as.matrix(read.csv("resSamples_larvae-pupae_10000.csv"))

resE <- as.numeric(resE)
resA <- as.numeric(resA)
resL <- resLP[,1]
resP <- resLP[,2]

## ############################################################
## Take the maximum a posteriori  (MAP) values of resolutions #
(resE_MAP  <- sort(table(resE), decreasing=T)[1])
(resL_MAP  <- sort(table(resL), decreasing=T)[1])
(resP_MAP  <- sort(table(resP), decreasing=T)[1])
(resGC_MAP <- sort(table(resA), decreasing=T)[1])

## Check same number of lines
(mcmcLines <- length(resE))
(mcmcLines <- length(resA))
(mcmcLines <- length(resL))


## ############################
## Set directory for plotting #
setwd(baseDir)
if (!is.element("figures", dir())) {
    system("mkdir figures")
}
setwd("figures")
print(getwd())

## #######
## PLOT ##
## #######

if (PDF)
    pdf("RES_histograms.pdf", height=7, width=25)
MAR <- c(5, 10, 1.2, 1)
par(mar=MAR, mfrow=c(1,4))
par(mai=c(0.1, 0.8, 0.2, 0.1),  omi=c(1, 0.7, 0.3, 0.2)) 
CEX.AXIS <- 3.8
CEX.LAB  <- 3
LineX    <- 6
LineY    <- 5
## Egg
hist(resE, include.lowest=TRUE, breaks=1:50,
     axes=FALSE,
     ylim=c(0,mcmcLines), xlab="", ylab="", main="", xaxt="n",
     cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
mtext("Resolution Egg", 1, line=LineX, cex=CEX.LAB)
mtext("Frequency", 2, line=LineY, cex=CEX.LAB)
axis(1, at=seq(0,50,10), cex.axis=CEX.AXIS, padj=.6)
axis(2, at=seq(0,10000,by=1000), cex.axis=CEX.AXIS)
## Larva 
hist(resL, breaks=1:15, axes=FALSE,
     ylim=c(0,mcmcLines), xlab="", ylab="", main="",
     cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
mtext("Resolution Larva", 1, line=LineX, cex=CEX.LAB)
axis(1, at=seq(0,15,5), cex.axis=CEX.AXIS, padj=.6)
axis(2, at=seq(0,10000,by=1000), cex.axis=CEX.AXIS)
## Pupa
hist(resP, breaks=1:max(resP), axes=FALSE,
     ylim=c(0,400), xlab="", ylab="", main="",
     cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
mtext("Resolution Pupa", 1, line=LineX, cex=CEX.LAB)
axis(1, at=seq(0,50,10), cex.axis=CEX.AXIS, padj=.6)
axis(2, at=seq(0,400,by=100), cex.axis=CEX.AXIS)
## Adult
hist(resA, breaks=1:max(resA), axes=FALSE,
     ylim=c(0,800), xlab="", ylab="", main="",
     cex.lab=CEX.LAB, cex.axis=CEX.AXIS, xaxt="n")
mtext("Resolution Adult", 1, line=LineX, cex=CEX.LAB)
axis(1, at=seq(0,50,10), cex.axis=CEX.AXIS, padj=.6)
axis(2, at=seq(0,800,by=100), cex.axis=CEX.AXIS)
if (PDF)
    dev.off()



## ## Good values for y-axes if I'm using burnin'
## if (PDF)
##     pdf("RES_histograms.pdf", height=7, width=25)
## MAR <- c(5, 10, 1.2, 1)
## par(mar=MAR, mfrow=c(1,4))
## par(mai=c(0.1, 0.8, 0.2, 0.1),  omi=c(1, 0.7, 0.3, 0.2)) 
## CEX.AXIS <- 3.8
## CEX.LAB  <- 3
## LineX    <- 6
## LineY    <- 5
## ## Egg
## hist(resE, include.lowest=TRUE, breaks=1:50,
##      axes=FALSE,
##      ylim=c(0,mcmcLines), xlab="", ylab="", main="", xaxt="n",
##      cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
## mtext("Resolution Egg", 1, line=LineX, cex=CEX.LAB)
## mtext("Frequency", 2, line=LineY, cex=CEX.LAB)
## axis(1, at=seq(0,50,10), cex.axis=CEX.AXIS, padj=.6)
## axis(2, at=seq(0,5000,by=1000), cex.axis=CEX.AXIS)#, tick=FALSE)
## ## Larva 
## hist(resL, breaks=1:15, axes=FALSE,
##      ylim=c(0,mcmcLines), xlab="", ylab="", main="",
##      cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
## mtext("Resolution Larva", 1, line=LineX, cex=CEX.LAB)
## axis(1, at=seq(0,15,5), cex.axis=CEX.AXIS, padj=.6)
## axis(2, at=seq(0,5000,by=1000), cex.axis=CEX.AXIS)
## ## Pupa
## hist(resP, breaks=1:max(resP), axes=FALSE,
##      ylim=c(0,200), xlab="", ylab="", main="",
##      cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
## mtext("Resolution Pupa", 1, line=LineX, cex=CEX.LAB)
## axis(1, at=seq(0,50,10), cex.axis=CEX.AXIS, padj=.6)
## axis(2, at=seq(0,200,by=50), cex.axis=CEX.AXIS)
## ## Adult
## hist(resA, breaks=1:max(resA), axes=FALSE,
##      ylim=c(0,400), xlab="", ylab="", main="",
##      cex.lab=CEX.LAB, cex.axis=CEX.AXIS, xaxt="n")
## mtext("Resolution Adult", 1, line=LineX, cex=CEX.LAB)
## axis(1, at=seq(0,50,10), cex.axis=CEX.AXIS, padj=.6)
## axis(2, at=seq(0,400,by=50), cex.axis=CEX.AXIS)
## ##
## if (PDF)
##     dev.off()


