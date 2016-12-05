## ################
## OBTAIN SPLINES #
## This script is called from interpolation.R
require(nimble)

#############################################################################
## Set temperatures (X variable) at which interpolation is to be performed ##
Delta <- 1/5
minX  <-  0
maxX  <- 50
X     <- seq(minX, maxX, by=Delta)
lX    <- length(X)
extrX <- c(0, 50) ## Extreme limits used by rejection samplers to reject crossing q01 ~ q99 extrapolations

## STARTING THE COLDEST DAY
## For cosine profile at different ranges
X.15.25  <- 20 -  5 * cos (2*pi*(0:364)/364)
X.15.30  <- 22.5 -  7.5 * cos (2*pi*(0:364)/364)
lX.15.25 <- length(X.15.25)
lX.15.30 <- length(X.15.30)
if (FALSE) {
    plot(X.15.25)
    plot(X.15.30)
}

#####################################################
## Initialise matrices for output of interpolation ##

## sur : daily survival
surX       <- matrix(0, nLines, lX)
surX.15.25 <- matrix(0, nLines, lX.15.25)
surX.15.30 <- matrix(0, nLines, lX.15.30)

## fecundity
if (Stg=="gonotrophicCycle"){
    fecX       <- matrix(0, nLines, lX)
    fecX.15.25 <- matrix(0, nLines, lX.15.25)
    fecX.15.30 <- matrix(0, nLines, lX.15.30)
}

## (a1,a2): beta distribution parameters
a1X <- a2X <- meanX <- surX  
## For cosine profile at different ranges
a1X.15.25 <- a2X.15.25 <- surX.15.25
a1X.15.30 <- a2X.15.30 <- surX.15.30

## For storing quantiles (P01 & P99 in the paper)
q01X <- surX
q99X <- surX
q01X.15.30 <- surX.15.30
q99X.15.30 <- surX.15.30
q01X.15.25 <- surX.15.25
q99X.15.25 <- surX.15.25

########################################
## Initialise List for SUCH functions ##
splineList <- vector("list", nLines)

###############################
## Loop on MCMC output lines ##
for (jj in 1:nLines) { ## jj=1
    if(jj%%50==0)
        print(jj)
    if (Stg=="gonotrophicCycle") 
        (EFec <- EFecMat[jj,])
    (q01 <- P01Mat[jj,]) 
    (q99 <- P99Mat[jj,]) 
    (sur <- surMat[jj,])
    (res <- resVec[jj])
    ## SURVIVAL; perform interpolations on logit scale
    surSpline        <- SUCH(temps, logit(sur)) 
    (surX[jj,]       <- ilogit(surSpline(X)))
    (surX.15.30[jj,] <- ilogit(surSpline(X.15.30)))
    (surX.15.25[jj,] <- ilogit(surSpline(X.15.25)))
    ## FECUNDITY; perform interpolations on log scale
    if (Stg=="gonotrophicCycle")  {
        fecSpline <- SUCH(temps, log(EFec))
        fecX[jj,] <- exp(fecSpline(X))
        fecX.15.30[jj,] <- exp(fecSpline(X.15.30))
        fecX.15.25[jj,] <- exp(fecSpline(X.15.25))
    }
    ## KERNEL; perform interpolations on logit scale
    fgSpline    <- DoubleSUCH(fOri=logit(q01), gOri=logit(q99), TempsOri=temps, RangeTempsNew=extrX)
    (q01X[jj,] <- ilogit(fgSpline[[1]](X))) 
    (q99X[jj,] <- ilogit(fgSpline[[2]](X))) 
    ## Apply splines to cosine temperature profile
    (q01X.15.25[jj,] <- ilogit(fgSpline[[1]](X.15.25))) 
    (q99X.15.25[jj,] <- ilogit(fgSpline[[2]](X.15.25))) 
    (q01X.15.30[jj,] <- ilogit(fgSpline[[1]](X.15.30))) 
    (q99X.15.30[jj,] <- ilogit(fgSpline[[2]](X.15.30))) 
    ## Set lower & upper limits on extrapolated values of q
    (LowLim   <- 1E-4)
    (UpperLim <- 1 - 1E-10)
    ## ###########################################
    ## Estimate alpha via Nelder-Mead algorithm ##
    for (tt in 1:length(X)) { 
        (qn <- c(q01X[jj,tt], q99X[jj,tt])) 
        ## Enforce lower limit on q 
        qn[qn < LowLim] <- LowLim  ## Much smaller than this and optim crashes, but this is well below 1/res
        ## Run Nelder-Mead 
        (aIni        <- IniValAlpha(p=c(0.01, 0.99), qn=qn))
        aIni[aIni<1] <- 1
        (op          <- optim(par=as.numeric(aIni), fn=sumSqError, target=qn, p=c(0.01, 0.99)))
        (aHat        <- op$par)
        count <- 0
        ## Run optim several times since one alone can fail to converge 
        while (any(abs(aIni - aHat) > 1E-5)) {
            aIni  <- aHat
            (op   <- optim(par=as.numeric(aIni), fn=sumSqError, target=qn, p=c(0.01, 0.99)))        
            (aHat <- op$par)
            count <- count + 1
        }
        (a1X[jj,tt] <- aHat[1]) 
        (a2X[jj,tt] <- aHat[2])
        (meanX[jj,tt]   <- aHat[1] / (aHat[1] + aHat[2]))
    }
    ## T in (15,25) 
    for (tt in 1:length(X.15.25)) {
        (qn <- c(q01X.15.25[jj,tt], q99X.15.25[jj,tt]))
        ## Enforce lower limit
        qn[qn < 1E-4] <- 1E-4 ## Much smaller than this and optim crashes, but this is well bellow 1/res
        ## Run Nelder-Mead 
        (aIni        <- IniValAlpha(p=c(0.01, 0.99), qn=qn))
        aIni[aIni<1] <- 1
        (op          <- optim(par=as.numeric(aIni), fn=sumSqError, target=qn, p=c(0.01, 0.99)))
        (aHat        <- op$par)
        count <- 0
        ## Run optim several times since one alone can fail to converge 
        while (any(abs(aIni - aHat) > 1E-5)) {
            aIni  <- aHat
            (op   <- optim(par=as.numeric(aIni), fn=sumSqError, target=qn, p=c(0.01, 0.99)))        
            (aHat <- op$par)
            count <- count + 1
            if(count %% 10 ==0)
                print(paste("optim problems at count", count, "with  qn =", qn, sep=" "))
        }
        (a1X.15.25[jj,tt] <- aHat[1]) 
        (a2X.15.25[jj,tt] <- aHat[2])
    } 
    ## T in (15,30) 
    for (tt in 1:length(X.15.30)) { 
        (qn <- c(q01X.15.30[jj,tt], q99X.15.30[jj,tt]))
        ## Enforce lower limit
        qn[qn < 1E-4] <- 1E-4
        ## Run Nelder-Mead 
        (aIni        <- IniValAlpha(p=c(0.01, 0.99), qn=qn))
        aIni[aIni<1] <- 1
        (op          <- optim(par=as.numeric(aIni), fn=sumSqError, target=qn, p=c(0.01, 0.99)))
        (aHat        <- op$par)
        count <- 0
        ## Run optim several times since one alone can fail to converge 
        while (any(abs(aIni - aHat) > 1E-5)) {
            aIni  <- aHat
            (op   <- optim(par=as.numeric(aIni), fn=sumSqError, target=qn, p=c(0.01, 0.99)))        
            (aHat <- op$par)
            count <- count + 1
            if(count %% 10 ==0)
                print(paste("optim problems at count", count, "with  qn =", qn, sep=" "))
        }
        (a1X.15.30[jj,tt] <- aHat[1]) 
        (a2X.15.30[jj,tt] <- aHat[2])
    }
}
    
    
## ####################################
## Obtain 95% CI on random quantiles ##
q01.CI95  <- matrix(0,nrow=3,ncol=lX)
mean.CI95 <- matrix(0,nrow=3,ncol=lX)
q99.CI95  <- matrix(0,nrow=3,ncol=lX)
Quantiles <- c(0.025,0.5,0.975)
for (tt in 1:length(X)) {
    q01.CI95[,tt]  <- quantile(q01X[,tt],  Quantiles)
    mean.CI95[,tt] <- quantile(meanX[,tt], Quantiles)
    q99.CI95[,tt]  <- quantile(q99X[,tt],  Quantiles)
}


